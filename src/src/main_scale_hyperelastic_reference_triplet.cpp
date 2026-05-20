#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "AppRuntimeSupport.h"
#include "FEMModel.h"
#include "HyperelasticScenarioSupport.h"
#include "ResultFileExporter.h"
#include "constants.h"
#include "loadFunction.h"
#include "meshgenerator.h"

using namespace Constants;

namespace {

constexpr double kYoungsModulus = 11.84;
constexpr double kPoissonsRatio = 0.48;
constexpr double kThickness = 1.0;

const Eigen::Vector2d kTireCenter(0.0, 0.0);
constexpr double kTireInnerRadius = 250.0;
constexpr double kTireOuterRadius = 300.0;
constexpr double kTireStartAngle = PI + 60.0 * DEG_TO_RAD;
constexpr double kTireEndAngle = kTireStartAngle + 60.0 * DEG_TO_RAD;
constexpr double kContactCenterAngle = 1.5 * PI;
constexpr double kExpectedContactHalfAngle = 5.0 * DEG_TO_RAD;
constexpr double kPlaneInitialGap = 2.0;
constexpr double kPrescribedInnerBoundaryDy = -10.0;

struct MeshVariant {
    std::string label;
    int radialLayers = 0;
    int circumferentialNodes = 0;
    double circumferentialRefinementStrength = 0.0;
    double radialRefinementStrength = 0.0;
    double circumferentialLocalizationPower = 2.0;
};

struct CaseOutcome {
    bool success = false;
    FEModel::PerformanceMetrics metrics;
    HyperelasticScenarioSupport::FiniteStrainDiagnostics diagnostics;
    double minimumSignedDistance = std::numeric_limits<double>::quiet_NaN();
    int nodeCount = 0;
    int elementCount = 0;
    int candidateFacetCount = 0;
    std::filesystem::path outputDirectory;
};

struct SurrogateLoadParameters {
    double maxPressure = 0.0;
    double contactHalfAngle = 0.0;
    double contactHalfWidth = 0.0;
    double contactCenterX = 0.0;
    double totalNormalForce = 0.0;
    int activeFacetCount = 0;
};

struct ContactDriveSettings {
    bool useMovingPlane = false;
    double movingPlaneNormalDisplacement = -kPrescribedInnerBoundaryDy;
    bool prescribeInnerBoundaryX = false;
    double prescribedInnerBoundaryDx = 0.0;
};

std::string environmentString(std::string_view variableName) {
    if (const char* value = std::getenv(std::string(variableName).c_str());
        value != nullptr) {
        return value;
    }
    return {};
}

double environmentDouble(std::string_view variableName, double defaultValue) {
    const std::string value = environmentString(variableName);
    if (value.empty()) {
        return defaultValue;
    }
    return std::stod(value);
}

int environmentInt(std::string_view variableName, int defaultValue) {
    const std::string value = environmentString(variableName);
    if (value.empty()) {
        return defaultValue;
    }
    return std::stoi(value);
}

std::string toLower(std::string value) {
    std::transform(value.begin(),
        value.end(),
        value.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::tolower(character));
        });
    return value;
}

bool environmentBool(std::string_view variableName, bool defaultValue) {
    const std::string value = toLower(environmentString(variableName));
    if (value.empty()) {
        return defaultValue;
    }
    return value == "1" || value == "true" || value == "yes" || value == "on";
}

ContactDriveSettings readContactDriveSettings() {
    ContactDriveSettings settings;
    const std::string contactDrive =
        toLower(environmentString("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_CONTACT_DRIVE"));
    settings.useMovingPlane =
        environmentBool("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_MOVE_PLANE", false) ||
        contactDrive == "plane" ||
        contactDrive == "moving_plane" ||
        contactDrive == "rigid_plane";
    settings.movingPlaneNormalDisplacement = environmentDouble(
        "FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_PLANE_DY",
        -kPrescribedInnerBoundaryDy);
    settings.prescribeInnerBoundaryX = environmentBool(
        "FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_PRESCRIBE_INNER_X", false);
    settings.prescribedInnerBoundaryDx = environmentDouble(
        "FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_INNER_DX", 0.0);
    return settings;
}

std::set<std::string> parseRequestedCases(std::string_view casesSpec) {
    std::set<std::string> requestedCases;
    std::stringstream stream{std::string(casesSpec)};
    std::string token;
    while (std::getline(stream, token, ',')) {
        token.erase(std::remove_if(token.begin(),
                          token.end(),
                          [](unsigned char character) {
                              return std::isspace(character) != 0;
                          }),
            token.end());
        if (!token.empty()) {
            requestedCases.insert(toLower(token));
        }
    }
    return requestedCases;
}

std::vector<MeshVariant> buildMeshVariants() {
    return {
        {"contact_focused_very_coarse", 61, 81, 6.0, 1.4, 2.0},
        {"contact_focused_coarse", 91, 121, 7.5, 1.8, 2.15},
        {"contact_focused_medium", 181, 241, 10.0, 2.5, 2.5},
        {"contact_focused_dense", 241, 321, 12.0, 3.0, 2.75},
        {"contact_focused_heavy", 301, 401, 14.0, 3.5, 3.0},
    };
}

MeshVariant resolveMeshVariant(std::string_view requestedLabel) {
    const auto variants = buildMeshVariants();
    const auto variantIt = std::find_if(variants.begin(),
        variants.end(),
        [requestedLabel](const MeshVariant& variant) {
            return variant.label == requestedLabel;
        });
    if (variantIt == variants.end()) {
        throw std::invalid_argument(
            "Unknown hyperelastic triplet mesh label: " + std::string(requestedLabel));
    }
    return *variantIt;
}

MeshGenerator::TireContactAnalysisControl buildContactControl(
    const MeshVariant& meshVariant,
    const ContactDriveSettings& driveSettings) {
    MeshGenerator::TireContactAnalysisControl control;
    control.mesh.center = kTireCenter;
    control.mesh.innerRadius = kTireInnerRadius;
    control.mesh.outerRadius = kTireOuterRadius;
    control.mesh.startAngle = kTireStartAngle;
    control.mesh.endAngle = kTireEndAngle;
    control.mesh.radialLayers = meshVariant.radialLayers;
    control.mesh.circumferentialNodes = meshVariant.circumferentialNodes;
    control.mesh.materialId = 1;
    control.mesh.refineCircumferentiallyNearContact = true;
    control.mesh.refineRadiallyToOuterSurface = true;
    control.mesh.localizeCircumferentialRefinementToOuterSurface = true;
    control.mesh.expectedContactCenterAngle = kContactCenterAngle;
    control.mesh.expectedContactHalfAngle = kExpectedContactHalfAngle;
    control.mesh.circumferentialRefinementStrength =
        meshVariant.circumferentialRefinementStrength;
    control.mesh.radialRefinementStrength = meshVariant.radialRefinementStrength;
    control.mesh.circumferentialLocalizationPower =
        meshVariant.circumferentialLocalizationPower;
    control.mesh.candidateFacetWindowScale = 2.5;
    control.mesh.outerRadiusTolerance = 1.0e-5;

    control.rigidPlane = RigidPlane2D{
        Eigen::Vector2d(0.0, 1.0),
        -(kTireOuterRadius + kPlaneInitialGap)
    };
    control.innerRadiusTolerance = 1.0e-5;
    if (driveSettings.useMovingPlane) {
        control.prescribeInnerBoundaryX = false;
        control.prescribeInnerBoundaryY = false;
        control.addInnerBoundaryAnchor = false;
    }
    else {
        control.prescribeInnerBoundaryX = driveSettings.prescribeInnerBoundaryX;
        control.prescribeInnerBoundaryY = true;
        control.prescribedInnerBoundaryDx = driveSettings.prescribedInnerBoundaryDx;
        control.prescribedInnerBoundaryDy = kPrescribedInnerBoundaryDy;
        control.addInnerBoundaryAnchor = !driveSettings.prescribeInnerBoundaryX;
    }
    control.anchorFixX = true;
    control.anchorFixY = false;
    control.anchorSelectMinimumX = true;
    control.anchorSelectContactCenterAngle = true;
    return control;
}

double computeMinimumSignedDistance(
    const FEModel& model,
    const std::shared_ptr<Assembly>& assembly,
    const RigidPlane2D& plane) {
    const auto& nodes = assembly->getNodes();
    const auto nodalDisplacements = model.getNodalDisplacements();

    double minimumSignedDistance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < nodes.size() && i < nodalDisplacements.size(); ++i) {
        const Eigen::Vector2d deformedPosition =
            nodes[i]->getCoordinates() + nodalDisplacements[i];
        minimumSignedDistance =
            std::min(minimumSignedDistance, plane.signedDistance(deformedPosition));
    }
    return minimumSignedDistance;
}

RigidPlane2D movePlaneAlongNormal(
    const RigidPlane2D& plane,
    double normalDisplacement) {
    RigidPlane2D movedPlane = plane;
    movedPlane.offset += normalDisplacement;
    return movedPlane;
}

void fixInnerBoundaryForMovingPlaneDrive(
    const HyperelasticScenarioSupport::TireContactCase& tireCase) {
    for (int nodeId : tireCase.innerBoundaryNodeIds) {
        tireCase.assembly->addFixedNode(nodeId, true, true);
    }
}

std::string formatNumber(double value) {
    std::ostringstream stream;
    stream << std::setprecision(10) << value;
    return stream.str();
}

std::string wrapCaseName(
    std::string_view prefix,
    const MeshVariant& meshVariant,
    int loadSteps) {
    std::ostringstream stream;
    stream << prefix << "_" << meshVariant.label << "_steps_" << loadSteps;
    return stream.str();
}

void addSharedExportMetrics(ResultFileExportOptions& exportOptions,
    const MeshVariant& meshVariant,
    double loadSteps,
    const MeshGenerator::RingMeshDiagnostics& diagnostics) {
    exportOptions.extraNumericMetrics.insert(exportOptions.extraNumericMetrics.end(), {
        {"hyperelastic_load_steps", loadSteps},
        {"youngs_modulus", kYoungsModulus},
        {"poissons_ratio", kPoissonsRatio},
        {"thickness", kThickness},
        {"ring_center_x", kTireCenter.x()},
        {"ring_center_y", kTireCenter.y()},
        {"ring_inner_radius", kTireInnerRadius},
        {"ring_outer_radius", kTireOuterRadius},
        {"ring_start_angle_rad", kTireStartAngle},
        {"ring_end_angle_rad", kTireEndAngle},
        {"ring_contact_center_angle_rad", kContactCenterAngle},
        {"ring_contact_half_angle_rad", kExpectedContactHalfAngle},
        {"ring_radial_layers", static_cast<double>(meshVariant.radialLayers)},
        {"ring_circumferential_nodes", static_cast<double>(meshVariant.circumferentialNodes)},
        {"mesh_circumferential_refinement_strength",
            meshVariant.circumferentialRefinementStrength},
        {"mesh_radial_refinement_strength", meshVariant.radialRefinementStrength},
        {"mesh_circumferential_localization_power",
            meshVariant.circumferentialLocalizationPower},
        {"mesh_min_radial_step", diagnostics.minRadialStep},
        {"mesh_max_radial_step", diagnostics.maxRadialStep},
        {"mesh_min_angular_step", diagnostics.minAngularStep},
        {"mesh_max_angular_step", diagnostics.maxAngularStep},
        {"mesh_min_outer_arc_step", diagnostics.minOuterArcStep},
        {"mesh_max_outer_arc_step", diagnostics.maxOuterArcStep},
        {"mesh_min_aspect_ratio", diagnostics.minAspectRatio},
        {"mesh_max_aspect_ratio", diagnostics.maxAspectRatio},
    });
}

SurrogateLoadParameters extractSurrogateLoadParameters(
    const std::vector<ContactFacetResult>& facetResults) {
    SurrogateLoadParameters parameters;
    if (facetResults.empty()) {
        throw std::runtime_error("Cannot extract surrogate load: no contact facet results");
    }

    double weightedXSum = 0.0;
    for (const auto& facetResult : facetResults) {
        if (!facetResult.active || facetResult.activeLength <= 0.0 ||
            facetResult.averagePressure <= 0.0) {
            continue;
        }

        parameters.activeFacetCount += 1;
        parameters.totalNormalForce += facetResult.integratedNormalForce;
        parameters.maxPressure = std::max(parameters.maxPressure, facetResult.averagePressure);
        weightedXSum += facetResult.integratedNormalForce * facetResult.referenceMidpoint.x();
    }

    if (parameters.activeFacetCount <= 0 || parameters.totalNormalForce <= 0.0 ||
        parameters.maxPressure <= 0.0) {
        throw std::runtime_error(
            "Cannot extract surrogate load: active contact patch was not detected");
    }

    parameters.contactCenterX = weightedXSum / parameters.totalNormalForce;

    for (const auto& facetResult : facetResults) {
        if (!facetResult.active || facetResult.activeLength <= 0.0 ||
            facetResult.averagePressure <= 0.0) {
            continue;
        }

        const Eigen::Vector2d relativePoint = facetResult.referenceMidpoint - kTireCenter;
        const double angle =
            std::atan2(relativePoint.y(), relativePoint.x()) - kContactCenterAngle;
        const double wrappedAngle = std::atan2(std::sin(angle), std::cos(angle));
        parameters.contactHalfAngle =
            std::max(parameters.contactHalfAngle, std::abs(wrappedAngle));
        parameters.contactHalfWidth =
            std::max(parameters.contactHalfWidth,
                std::abs(facetResult.referenceMidpoint.x() - parameters.contactCenterX));
    }

    if (parameters.contactHalfAngle <= 0.0) {
        throw std::runtime_error(
            "Cannot extract surrogate load: contact half-angle is not positive");
    }
    if (parameters.contactHalfWidth <= 0.0) {
        parameters.contactHalfWidth = kTireOuterRadius * std::sin(parameters.contactHalfAngle);
    }

    return parameters;
}

void exportCase(
    FEModel& model,
    const MeshVariant& meshVariant,
    const std::filesystem::path& outputDirectory,
    std::string_view caseName,
    std::string_view studyFamily,
    std::string_view contactMethod,
    std::string_view contactParameterName,
    double contactParameterValue,
    const MeshGenerator::RingMeshDiagnostics& diagnostics,
    const HyperelasticScenarioSupport::FiniteStrainDiagnostics& finiteStrainDiagnostics,
    double loadSteps,
    double minimumSignedDistance,
    const std::vector<std::pair<std::string, double>>& extraNumericMetrics = {}) {
    ResultFileExportOptions exportOptions;
    exportOptions.outputDirectory = outputDirectory;
    exportOptions.baseName = "solution";
    exportOptions.extraStringMetrics = {
        {"case_name", std::string(caseName)},
        {"mesh_label", meshVariant.label},
        {"geometry_family", "tire_ring"},
        {"study_family", std::string(studyFamily)},
        {"contact_method", std::string(contactMethod)},
        {"contact_parameter_name", std::string(contactParameterName)}
    };
    exportOptions.extraNumericMetrics = {
        {"contact_parameter_value", contactParameterValue},
        {"minimum_jacobian_determinant", finiteStrainDiagnostics.minJ},
        {"maximum_jacobian_determinant", finiteStrainDiagnostics.maxJ},
        {"maximum_strain_energy_density", finiteStrainDiagnostics.maxStrainEnergyDensity},
        {"gauss_point_count", static_cast<double>(finiteStrainDiagnostics.gaussPointCount)}
    };

    if (std::isfinite(minimumSignedDistance)) {
        exportOptions.extraNumericMetrics.push_back(
            {"minimum_signed_distance", minimumSignedDistance});
    }

    addSharedExportMetrics(exportOptions, meshVariant, loadSteps, diagnostics);
    for (const auto& [name, value] : extraNumericMetrics) {
        exportOptions.extraNumericMetrics.push_back({name, value});
    }

    ResultFileExporter::exportSolution(model, exportOptions);
    AppRuntimeSupport::runPostprocessorIfEnabled(outputDirectory);
}

CaseOutcome runPenaltyContactCase(
    const MeshVariant& meshVariant,
    double penalty,
    int loadSteps,
    int maxAdaptiveLoadSteps,
    int maxIterations,
    const std::filesystem::path& outputRoot) {
    const ContactDriveSettings driveSettings = readContactDriveSettings();
    const auto control = buildContactControl(meshVariant, driveSettings);
    const auto tireCase = HyperelasticScenarioSupport::buildFiniteStrainTireContactCase(
        control, kYoungsModulus, kPoissonsRatio, kThickness);
    if (driveSettings.useMovingPlane) {
        fixInnerBoundaryForMovingPlaneDrive(tireCase);
    }
    const RigidPlane2D finalContactPlane = driveSettings.useMovingPlane
        ? movePlaneAlongNormal(
              tireCase.rigidPlane,
              driveSettings.movingPlaneNormalDisplacement)
        : tireCase.rigidPlane;

    FEModel model;
    model.setAssembly(tireCase.assembly);
    model.setSolverTolerance(1.0e-8);
    model.setMaxIterations(maxIterations);
    model.setHyperelasticLoadSteps(loadSteps);
    model.setMaxAdaptiveHyperelasticLoadSteps(maxAdaptiveLoadSteps);
    model.configureRigidPlanePenaltyContact(
        tireCase.rigidPlane, tireCase.mesh.candidateContactFacets, penalty);
    if (driveSettings.useMovingPlane) {
        model.setContactPlaneLoadRamp(tireCase.rigidPlane, finalContactPlane);
    }

    CaseOutcome outcome;
    outcome.success = model.solveHyperelastic();
    outcome.metrics = model.getPerformanceMetrics();
    outcome.diagnostics =
        HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);
    outcome.minimumSignedDistance =
        computeMinimumSignedDistance(model, tireCase.assembly, model.getContactPlane());
    outcome.nodeCount = tireCase.assembly->getNodeCount();
    outcome.elementCount =
        static_cast<int>(tireCase.assembly->getFiniteStrainElements().size());
    outcome.candidateFacetCount =
        static_cast<int>(tireCase.mesh.candidateContactFacets.size());
    outcome.outputDirectory = outputRoot / "penalty_contact";

    exportCase(model,
        meshVariant,
        outcome.outputDirectory,
        wrapCaseName("penalty_contact", meshVariant, loadSteps),
        "main_scale_hyperelastic_reference_triplet",
        "penalty",
        "penalty_parameter",
        penalty,
        tireCase.mesh.diagnostics,
        outcome.diagnostics,
        static_cast<double>(loadSteps),
        outcome.minimumSignedDistance,
        {
            {"penalty_parameter", penalty},
            {"plane_initial_gap", kPlaneInitialGap},
            {"contact_drive_moving_plane", driveSettings.useMovingPlane ? 1.0 : 0.0},
            {"plane_prescribed_normal_displacement",
                driveSettings.useMovingPlane
                    ? driveSettings.movingPlaneNormalDisplacement
                    : 0.0},
            {"plane_initial_offset", tireCase.rigidPlane.offset},
            {"plane_final_offset", model.getContactPlane().offset},
            {"inner_boundary_prescribe_x", control.prescribeInnerBoundaryX ? 1.0 : 0.0},
            {"inner_boundary_prescribed_dx", control.prescribedInnerBoundaryDx},
            {"inner_boundary_prescribe_y", control.prescribeInnerBoundaryY ? 1.0 : 0.0},
            {"inner_boundary_prescribed_dy",
                control.prescribeInnerBoundaryY ? control.prescribedInnerBoundaryDy : 0.0},
            {"inner_boundary_fixed_x", driveSettings.useMovingPlane ? 1.0 : 0.0},
            {"inner_boundary_fixed_y", driveSettings.useMovingPlane ? 1.0 : 0.0},
            {"inner_boundary_anchor_enabled", control.addInnerBoundaryAnchor ? 1.0 : 0.0},
            {"rigid_plane_normal_x", model.getContactPlane().normal.x()},
            {"rigid_plane_normal_y", model.getContactPlane().normal.y()},
            {"rigid_plane_offset", model.getContactPlane().offset},
        });

    return outcome;
}

std::pair<CaseOutcome, SurrogateLoadParameters> runAugmentedLagrangianContactCase(
    const MeshVariant& meshVariant,
    double augmentationScalingFactor,
    int loadSteps,
    int maxAdaptiveLoadSteps,
    int maxIterations,
    const std::filesystem::path& outputRoot) {
    const ContactDriveSettings driveSettings = readContactDriveSettings();
    const auto control = buildContactControl(meshVariant, driveSettings);
    const auto tireCase = HyperelasticScenarioSupport::buildFiniteStrainTireContactCase(
        control, kYoungsModulus, kPoissonsRatio, kThickness);
    if (driveSettings.useMovingPlane) {
        fixInnerBoundaryForMovingPlaneDrive(tireCase);
    }
    const RigidPlane2D finalContactPlane = driveSettings.useMovingPlane
        ? movePlaneAlongNormal(
              tireCase.rigidPlane,
              driveSettings.movingPlaneNormalDisplacement)
        : tireCase.rigidPlane;

    FEModel model;
    model.setAssembly(tireCase.assembly);
    model.setSolverTolerance(1.0e-8);
    model.setMaxIterations(maxIterations);
    model.setHyperelasticLoadSteps(loadSteps);
    model.setMaxAdaptiveHyperelasticLoadSteps(maxAdaptiveLoadSteps);

    AugmentedLagrangianSettings settings;
    settings.augmentationParameter = 1.0;
    settings.useAutomaticAugmentationScaling = true;
    settings.automaticScalingFactor = augmentationScalingFactor;
    model.configureRigidPlaneAugmentedLagrangianContact(
        tireCase.rigidPlane, tireCase.mesh.candidateContactFacets, settings);
    if (driveSettings.useMovingPlane) {
        model.setContactPlaneLoadRamp(tireCase.rigidPlane, finalContactPlane);
    }

    CaseOutcome outcome;
    outcome.success = model.solveHyperelastic();
    outcome.metrics = model.getPerformanceMetrics();
    outcome.diagnostics =
        HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);
    outcome.minimumSignedDistance =
        computeMinimumSignedDistance(model, tireCase.assembly, model.getContactPlane());
    outcome.nodeCount = tireCase.assembly->getNodeCount();
    outcome.elementCount =
        static_cast<int>(tireCase.assembly->getFiniteStrainElements().size());
    outcome.candidateFacetCount =
        static_cast<int>(tireCase.mesh.candidateContactFacets.size());
    outcome.outputDirectory = outputRoot / "augmented_lagrangian_contact";

    const SurrogateLoadParameters surrogateParameters =
        extractSurrogateLoadParameters(model.getContactFacetResults());

    exportCase(model,
        meshVariant,
        outcome.outputDirectory,
        wrapCaseName("augmented_lagrangian_contact", meshVariant, loadSteps),
        "main_scale_hyperelastic_reference_triplet",
        "augmented_lagrangian",
        "augmentation_scaling_factor",
        augmentationScalingFactor,
        tireCase.mesh.diagnostics,
        outcome.diagnostics,
        static_cast<double>(loadSteps),
        outcome.minimumSignedDistance,
        {
            {"augmentation_scaling_factor", augmentationScalingFactor},
            {"plane_initial_gap", kPlaneInitialGap},
            {"contact_drive_moving_plane", driveSettings.useMovingPlane ? 1.0 : 0.0},
            {"plane_prescribed_normal_displacement",
                driveSettings.useMovingPlane
                    ? driveSettings.movingPlaneNormalDisplacement
                    : 0.0},
            {"plane_initial_offset", tireCase.rigidPlane.offset},
            {"plane_final_offset", model.getContactPlane().offset},
            {"inner_boundary_prescribe_x", control.prescribeInnerBoundaryX ? 1.0 : 0.0},
            {"inner_boundary_prescribed_dx", control.prescribedInnerBoundaryDx},
            {"inner_boundary_prescribe_y", control.prescribeInnerBoundaryY ? 1.0 : 0.0},
            {"inner_boundary_prescribed_dy",
                control.prescribeInnerBoundaryY ? control.prescribedInnerBoundaryDy : 0.0},
            {"inner_boundary_fixed_x", driveSettings.useMovingPlane ? 1.0 : 0.0},
            {"inner_boundary_fixed_y", driveSettings.useMovingPlane ? 1.0 : 0.0},
            {"inner_boundary_anchor_enabled", control.addInnerBoundaryAnchor ? 1.0 : 0.0},
            {"rigid_plane_normal_x", model.getContactPlane().normal.x()},
            {"rigid_plane_normal_y", model.getContactPlane().normal.y()},
            {"rigid_plane_offset", model.getContactPlane().offset},
            {"surrogate_max_pressure", surrogateParameters.maxPressure},
            {"surrogate_contact_half_angle_rad", surrogateParameters.contactHalfAngle},
            {"surrogate_contact_half_width", surrogateParameters.contactHalfWidth},
            {"surrogate_contact_center_x", surrogateParameters.contactCenterX},
            {"surrogate_total_normal_force", surrogateParameters.totalNormalForce},
            {"surrogate_active_facets", static_cast<double>(surrogateParameters.activeFacetCount)},
        });

    return {outcome, surrogateParameters};
}

CaseOutcome runNoContactSurrogateCase(
    const MeshVariant& meshVariant,
    const SurrogateLoadParameters& surrogateParameters,
    int loadSteps,
    int maxAdaptiveLoadSteps,
    int maxIterations,
    const std::filesystem::path& outputRoot) {
    auto control = buildContactControl(meshVariant, readContactDriveSettings());
    control.mesh.expectedContactHalfAngle =
        std::max(surrogateParameters.contactHalfAngle, 0.5 * DEG_TO_RAD);
    control.mesh.candidateFacetWindowScale = 1.5;
    control.prescribeInnerBoundaryX = false;
    control.prescribeInnerBoundaryY = false;
    control.addInnerBoundaryAnchor = false;

    auto tireCase = HyperelasticScenarioSupport::buildFiniteStrainTireContactCase(
        control, kYoungsModulus, kPoissonsRatio, kThickness);
    for (int nodeId : tireCase.innerBoundaryNodeIds) {
        tireCase.assembly->addFixedNode(nodeId, true, true);
    }

    auto surrogateLoad = std::make_shared<LoadFunction>(
        LoadFunction::parabolicPressure(
            surrogateParameters.maxPressure,
            surrogateParameters.contactHalfWidth,
            surrogateParameters.contactCenterX));
    for (const auto& facet : tireCase.mesh.candidateContactFacets) {
        tireCase.assembly->addSurfaceLoad(
            facet.elementId, facet.surfaceIndex, surrogateLoad);
    }

    FEModel model;
    model.setAssembly(tireCase.assembly);
    model.setSolverTolerance(1.0e-8);
    model.setMaxIterations(maxIterations);
    model.setHyperelasticLoadSteps(loadSteps);
    model.setMaxAdaptiveHyperelasticLoadSteps(maxAdaptiveLoadSteps);

    CaseOutcome outcome;
    outcome.success = model.solveHyperelastic();
    outcome.metrics = model.getPerformanceMetrics();
    outcome.diagnostics =
        HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);
    outcome.nodeCount = tireCase.assembly->getNodeCount();
    outcome.elementCount =
        static_cast<int>(tireCase.assembly->getFiniteStrainElements().size());
    outcome.candidateFacetCount =
        static_cast<int>(tireCase.mesh.candidateContactFacets.size());
    outcome.outputDirectory = outputRoot / "no_contact_surrogate";

    exportCase(model,
        meshVariant,
        outcome.outputDirectory,
        wrapCaseName("no_contact_surrogate", meshVariant, loadSteps),
        "main_scale_hyperelastic_reference_triplet",
        "none",
        "surrogate_parabolic_pressure",
        surrogateParameters.maxPressure,
        tireCase.mesh.diagnostics,
        outcome.diagnostics,
        static_cast<double>(loadSteps),
        std::numeric_limits<double>::quiet_NaN(),
        {
            {"surrogate_max_pressure", surrogateParameters.maxPressure},
            {"surrogate_contact_half_angle_rad", surrogateParameters.contactHalfAngle},
            {"surrogate_contact_half_width", surrogateParameters.contactHalfWidth},
            {"surrogate_contact_center_x", surrogateParameters.contactCenterX},
            {"surrogate_total_normal_force", surrogateParameters.totalNormalForce},
            {"surrogate_active_facets", static_cast<double>(surrogateParameters.activeFacetCount)},
            {"surrogate_loaded_facets",
                static_cast<double>(tireCase.mesh.candidateContactFacets.size())},
        });

    return outcome;
}

void writeTripletSummary(
    const std::filesystem::path& outputRoot,
    const std::optional<CaseOutcome>& penaltyOutcome,
    const std::optional<CaseOutcome>& alOutcome,
    const std::optional<CaseOutcome>& surrogateOutcome,
    const std::optional<SurrogateLoadParameters>& surrogateParameters) {
    std::ofstream stream(outputRoot / "triplet_summary.csv");
    if (!stream.is_open()) {
        throw std::runtime_error("Failed to open triplet summary CSV output");
    }

    stream << "case,success,load_steps,converged_load_steps,nonlinear_iterations,"
              "max_penetration,total_time_seconds,contact_force_norm,minimum_signed_distance,"
              "minimum_jacobian_determinant,maximum_jacobian_determinant,linear_solver_backend\n";

    const auto writeRow = [&stream](std::string_view caseName, const CaseOutcome& outcome) {
        stream << caseName << ','
               << std::boolalpha << outcome.success << ','
               << outcome.metrics.loadSteps << ','
               << outcome.metrics.convergedLoadSteps << ','
               << outcome.metrics.nonlinearIterations << ','
               << outcome.metrics.maxPenetration << ','
               << outcome.metrics.totalTimeSeconds << ','
               << outcome.metrics.contactForceNorm << ','
               << formatNumber(outcome.minimumSignedDistance) << ','
               << outcome.diagnostics.minJ << ','
               << outcome.diagnostics.maxJ << ','
               << outcome.metrics.linearSolverBackend << '\n';
    };

    if (penaltyOutcome.has_value()) {
        writeRow("penalty_contact", *penaltyOutcome);
    }
    if (alOutcome.has_value()) {
        writeRow("augmented_lagrangian_contact", *alOutcome);
    }
    if (surrogateOutcome.has_value()) {
        writeRow("no_contact_surrogate", *surrogateOutcome);
    }

    if (surrogateParameters.has_value()) {
        stream << "\nparameter,value\n"
               << "surrogate_max_pressure," << surrogateParameters->maxPressure << '\n'
               << "surrogate_contact_half_angle_rad,"
               << surrogateParameters->contactHalfAngle << '\n'
               << "surrogate_contact_half_angle_deg,"
               << surrogateParameters->contactHalfAngle * RAD_TO_DEG << '\n'
               << "surrogate_contact_half_width," << surrogateParameters->contactHalfWidth
               << '\n'
               << "surrogate_contact_center_x," << surrogateParameters->contactCenterX << '\n'
               << "surrogate_total_normal_force," << surrogateParameters->totalNormalForce
               << '\n'
               << "surrogate_active_facets," << surrogateParameters->activeFacetCount << '\n';
    }
}

} // namespace

int main() {
    try {
        const std::string meshLabel =
            environmentString("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_MESH").empty()
            ? "contact_focused_very_coarse"
            : environmentString("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_MESH");
        const double penalty =
            environmentDouble("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_PENALTY", 5.0e2);
        const double alScalingFactor =
            environmentDouble("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_AL_SCALING", 20.0);
        const int contactLoadSteps =
            std::max(2, environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_CONTACT_STEPS", 12));
        const int maxAdaptiveContactLoadSteps =
            std::max(contactLoadSteps,
                environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_MAX_CONTACT_STEPS", 120));
        const int surrogateLoadSteps =
            std::max(2, environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_SURROGATE_STEPS", 12));
        const int maxAdaptiveSurrogateLoadSteps =
            std::max(surrogateLoadSteps,
                environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_MAX_SURROGATE_STEPS", 120));
        const int maxIterations =
            std::max(8, environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_MAX_ITERATIONS", 30));
        const int alLoadSteps = std::max(2,
            environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_AL_STEPS", contactLoadSteps));
        const int maxAdaptiveAlLoadSteps = std::max(alLoadSteps,
            environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_AL_MAX_STEPS",
                maxAdaptiveContactLoadSteps));
        const int alMaxIterations = std::max(8,
            environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_AL_MAX_ITERATIONS",
                maxIterations));
        const std::string outputSubdirectory =
            environmentString("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_OUTPUT_SUBDIRECTORY");
        const std::string requestedCasesSpec =
            environmentString("FEM_MAIN_SCALE_HYPERELASTIC_TRIPLET_CASES");
        const std::set<std::string> requestedCases =
            requestedCasesSpec.empty()
            ? std::set<std::string>{"penalty", "al", "surrogate"}
            : parseRequestedCases(requestedCasesSpec);

        if (requestedCases.empty()) {
            throw std::invalid_argument("No triplet cases requested");
        }
        if (requestedCases.contains("surrogate") && !requestedCases.contains("al")) {
            throw std::invalid_argument(
                "The no-contact surrogate case requires augmented Lagrangian contact "
                "to run in the same invocation");
        }

        const MeshVariant meshVariant = resolveMeshVariant(meshLabel);
        const std::filesystem::path outputRoot =
            AppRuntimeSupport::caseOutputDirectory(
                !outputSubdirectory.empty()
                ? outputSubdirectory
                : "main_scale_hyperelastic_reference_triplet");
        std::filesystem::create_directories(outputRoot);

        std::optional<CaseOutcome> penaltyOutcome;
        std::optional<CaseOutcome> alOutcome;
        std::optional<CaseOutcome> surrogateOutcome;
        std::optional<SurrogateLoadParameters> surrogateParameters;

        if (requestedCases.contains("penalty")) {
            penaltyOutcome = runPenaltyContactCase(
                meshVariant,
                penalty,
                contactLoadSteps,
                maxAdaptiveContactLoadSteps,
                maxIterations,
                outputRoot);
        }

        if (requestedCases.contains("al")) {
            auto [alCaseOutcome, extractedSurrogateParameters] =
                runAugmentedLagrangianContactCase(
                    meshVariant,
                    alScalingFactor,
                    alLoadSteps,
                    maxAdaptiveAlLoadSteps,
                    alMaxIterations,
                    outputRoot);
            alOutcome = std::move(alCaseOutcome);
            surrogateParameters = std::move(extractedSurrogateParameters);
        }

        if (requestedCases.contains("surrogate")) {
            surrogateOutcome = runNoContactSurrogateCase(
                meshVariant,
                *surrogateParameters,
                surrogateLoadSteps,
                maxAdaptiveSurrogateLoadSteps,
                maxIterations,
                outputRoot);
        }

        writeTripletSummary(
            outputRoot,
            penaltyOutcome,
            alOutcome,
            surrogateOutcome,
            surrogateParameters);

        std::cout << std::boolalpha
                  << "triplet_root=" << outputRoot << '\n'
                  << "mesh=" << meshVariant.label << '\n'
                  << "penalty_success="
                  << (penaltyOutcome.has_value() ? penaltyOutcome->success : true) << '\n'
                  << "al_success=" << (alOutcome.has_value() ? alOutcome->success : true) << '\n'
                  << "surrogate_success="
                  << (surrogateOutcome.has_value() ? surrogateOutcome->success : true) << '\n';

        if (surrogateParameters.has_value()) {
            std::cout << "surrogate_max_pressure=" << surrogateParameters->maxPressure << '\n'
                      << "surrogate_contact_half_angle_deg="
                      << surrogateParameters->contactHalfAngle * RAD_TO_DEG << '\n'
                      << "surrogate_contact_half_width="
                      << surrogateParameters->contactHalfWidth << '\n'
                      << "surrogate_contact_center_x="
                      << surrogateParameters->contactCenterX << '\n';
        }
        std::cout << std::endl;

        const bool overallSuccess =
            (!penaltyOutcome.has_value() || penaltyOutcome->success) &&
            (!alOutcome.has_value() || alOutcome->success) &&
            (!surrogateOutcome.has_value() || surrogateOutcome->success);
        return overallSuccess ? 0 : 1;
    }
    catch (const std::exception& exception) {
        std::cerr << "main_scale_hyperelastic_reference_triplet_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
