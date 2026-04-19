#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "AppRuntimeSupport.h"
#include "FEMModel.h"
#include "ResultFileExporter.h"
#include "assembly.h"
#include "constants.h"
#include "material.h"
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

struct StudyResult {
    std::string meshLabel;
    double penalty = 0.0;
    bool success = false;
    int candidateFacetCount = 0;
    int nodeCount = 0;
    int elementCount = 0;
    double minimumSignedDistance = std::numeric_limits<double>::infinity();
    MeshGenerator::RingMeshDiagnostics meshDiagnostics;
    FEModel::PerformanceMetrics metrics;
};

bool quickModeEnabled() {
    if (const char* quickModeValue = std::getenv("FEM_MAIN_SCALE_CONTACT_QUICK");
        quickModeValue != nullptr && std::string_view(quickModeValue) == "1") {
        return true;
    }
    return false;
}

bool longModeEnabled() {
    if (const char* longModeValue = std::getenv("FEM_MAIN_SCALE_CONTACT_LONG");
        longModeValue != nullptr && std::string_view(longModeValue) == "1") {
        return true;
    }
    return false;
}

std::string environmentString(std::string_view variableName) {
    if (const char* value = std::getenv(std::string(variableName).c_str());
        value != nullptr) {
        return value;
    }
    return {};
}

bool matchesOptionalFilter(std::string_view filterValue, std::string_view actualValue) {
    return filterValue.empty() || filterValue == actualValue;
}

bool matchesOptionalNumericFilter(std::string_view filterValue, double actualValue) {
    if (filterValue.empty()) {
        return true;
    }

    try {
        const double requestedValue = std::stod(std::string(filterValue));
        return std::abs(requestedValue - actualValue) <=
            1.0e-12 * std::max({1.0, std::abs(requestedValue), std::abs(actualValue)});
    }
    catch (const std::exception&) {
        return false;
    }
}

std::vector<double> buildPenaltyValues(bool quickMode, bool longMode) {
    if (quickMode) {
        return {5.0e2};
    }
    if (longMode) {
        return {1.0e1, 2.5e1, 5.0e1, 1.0e2, 2.5e2, 5.0e2, 1.0e3, 2.0e3, 5.0e3, 1.0e4};
    }
    return {1.0e2, 5.0e2, 2.0e3};
}

std::vector<MeshVariant> buildMeshVariants(bool quickMode, bool longMode) {
    if (quickMode) {
        return {
            {"contact_focused_medium", 141, 181, 10.0, 2.5, 2.5}
        };
    }
    if (longMode) {
        return {
            {"contact_focused_very_coarse", 61, 81, 6.0, 1.4, 2.0},
            {"contact_focused_coarse", 91, 121, 7.5, 1.8, 2.15},
            {"contact_focused_medium", 181, 241, 10.0, 2.5, 2.5},
            {"contact_focused_dense", 241, 321, 12.0, 3.0, 2.75},
            {"contact_focused_heavy", 301, 401, 14.0, 3.5, 3.0},
        };
    }

    return {
        {"contact_focused_medium", 181, 241, 10.0, 2.5, 2.5},
        {"contact_focused_dense", 241, 321, 12.0, 3.0, 2.75},
        {"contact_focused_heavy", 301, 401, 14.0, 3.5, 3.0},
    };
}

std::string makeCaseToken(const std::string& meshLabel, double penalty) {
    std::ostringstream stream;
    stream << meshLabel << "_penalty_" << std::llround(penalty);
    return stream.str();
}

std::string buildStudyHeader() {
    return "contact_method,contact_parameter_name,contact_parameter_value,mesh,success,"
           "nodes,elements,candidate_facets,active_facets,active_gauss_points,"
           "nonlinear_iterations,linear_iterations,max_penetration,minimum_signed_distance,"
           "contact_force_norm_n,contact_force_norm_kn,max_normal_multiplier,"
           "mean_normal_multiplier,contact_state_relative_update_norm,"
           "assembly_time_seconds,solve_time_seconds,total_time_seconds,"
           "matrix_nonzeros,min_radial_step,max_radial_step,min_angular_step,max_angular_step,"
           "min_outer_arc_step,max_outer_arc_step,min_aspect_ratio,max_aspect_ratio,linear_solver_backend";
}

std::string buildStudyRow(const StudyResult& result) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(8)
           << "penalty" << ','
           << "penalty_parameter" << ','
           << result.penalty << ','
           << result.meshLabel << ','
           << std::boolalpha << result.success << ','
           << result.nodeCount << ','
           << result.elementCount << ','
           << result.candidateFacetCount << ','
           << result.metrics.activeSetSize << ','
           << result.metrics.activeContactGaussPoints << ','
           << result.metrics.nonlinearIterations << ','
           << result.metrics.linearIterations << ','
           << result.metrics.maxPenetration << ','
           << result.minimumSignedDistance << ','
           << result.metrics.contactForceNorm << ','
           << (1.0e-3 * result.metrics.contactForceNorm) << ','
           << result.metrics.maxNormalContactMultiplier << ','
           << result.metrics.meanNormalContactMultiplier << ','
           << result.metrics.contactStateRelativeUpdateNorm << ','
           << result.metrics.assemblyTimeSeconds << ','
           << result.metrics.solveTimeSeconds << ','
           << result.metrics.totalTimeSeconds << ','
           << result.metrics.matrixNonZeros << ','
           << result.meshDiagnostics.minRadialStep << ','
           << result.meshDiagnostics.maxRadialStep << ','
           << result.meshDiagnostics.minAngularStep << ','
           << result.meshDiagnostics.maxAngularStep << ','
           << result.meshDiagnostics.minOuterArcStep << ','
           << result.meshDiagnostics.maxOuterArcStep << ','
           << result.meshDiagnostics.minAspectRatio << ','
           << result.meshDiagnostics.maxAspectRatio << ','
           << result.metrics.linearSolverBackend;
    return stream.str();
}

double computeMinimumSignedDistance(const FEModel& model,
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

StudyResult runStudyCase(const std::filesystem::path& outputRoot,
    const MeshVariant& meshVariant,
    double penalty) {
    StudyResult result;
    result.meshLabel = meshVariant.label;
    result.penalty = penalty;

    auto assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<Material>(1, kYoungsModulus, kPoissonsRatio, kThickness);
    assembly->addMaterial(material);

    MeshGenerator meshGenerator(assembly);

    MeshGenerator::TireContactAnalysisControl control;
    control.mesh.center = kTireCenter;
    control.mesh.innerRadius = kTireInnerRadius;
    control.mesh.outerRadius = kTireOuterRadius;
    control.mesh.startAngle = kTireStartAngle;
    control.mesh.endAngle = kTireEndAngle;
    control.mesh.radialLayers = meshVariant.radialLayers;
    control.mesh.circumferentialNodes = meshVariant.circumferentialNodes;
    control.mesh.materialId = material->getId();
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
    control.prescribeInnerBoundaryX = false;
    control.prescribeInnerBoundaryY = true;
    control.prescribedInnerBoundaryDy = kPrescribedInnerBoundaryDy;
    control.addInnerBoundaryAnchor = true;
    control.anchorFixX = true;
    control.anchorFixY = false;
    control.anchorSelectMinimumX = true;

    const auto setup = meshGenerator.generateTireContactAnalysisSetup(control);

    FEModel model;
    model.setAssembly(assembly);
    model.setMaxIterations(60);
    model.setSolverTolerance(1.0e-9);
    model.configureRigidPlaneContact(setup.rigidPlane, setup.mesh.candidateContactFacets, penalty);

    result.success = model.solveContact();
    result.metrics = model.getPerformanceMetrics();
    result.meshDiagnostics = setup.mesh.diagnostics;
    result.candidateFacetCount = static_cast<int>(setup.mesh.candidateContactFacets.size());
    result.nodeCount = assembly->getNodeCount();
    result.elementCount = assembly->getElementCount();

    if (result.success) {
        result.minimumSignedDistance =
            computeMinimumSignedDistance(model, assembly, setup.rigidPlane);
    }

    ResultFileExportOptions exportOptions;
    exportOptions.outputDirectory = outputRoot / makeCaseToken(meshVariant.label, penalty);
    exportOptions.baseName = "solution";
    exportOptions.extraStringMetrics = {
        {"case_name", makeCaseToken(meshVariant.label, penalty)},
        {"mesh_label", meshVariant.label},
        {"geometry_family", "tire_ring"},
        {"study_family", "main_scale_contact"},
        {"contact_method", "penalty"},
        {"contact_parameter_name", "penalty_parameter"}
    };
    exportOptions.extraNumericMetrics = {
        {"contact_parameter_value", penalty},
        {"penalty_parameter", penalty},
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
        {"rigid_plane_normal_x", setup.rigidPlane.normal.x()},
        {"rigid_plane_normal_y", setup.rigidPlane.normal.y()},
        {"rigid_plane_offset", setup.rigidPlane.offset},
        {"inner_boundary_prescribed_dy", kPrescribedInnerBoundaryDy},
        {"plane_initial_gap", kPlaneInitialGap},
        {"mesh_circumferential_refinement_strength", meshVariant.circumferentialRefinementStrength},
        {"mesh_radial_refinement_strength", meshVariant.radialRefinementStrength},
        {"mesh_circumferential_localization_power", meshVariant.circumferentialLocalizationPower},
        {"mesh_min_radial_step", result.meshDiagnostics.minRadialStep},
        {"mesh_max_radial_step", result.meshDiagnostics.maxRadialStep},
        {"mesh_min_angular_step", result.meshDiagnostics.minAngularStep},
        {"mesh_max_angular_step", result.meshDiagnostics.maxAngularStep},
        {"mesh_min_outer_arc_step", result.meshDiagnostics.minOuterArcStep},
        {"mesh_max_outer_arc_step", result.meshDiagnostics.maxOuterArcStep},
        {"mesh_min_aspect_ratio", result.meshDiagnostics.minAspectRatio},
        {"mesh_max_aspect_ratio", result.meshDiagnostics.maxAspectRatio},
    };
    if (std::isfinite(result.minimumSignedDistance)) {
        exportOptions.extraNumericMetrics.push_back(
            {"minimum_signed_distance", result.minimumSignedDistance});
    }

    ResultFileExporter::exportSolution(model, exportOptions);
    return result;
}

} // namespace

int main() {
    try {
        const bool quickMode = quickModeEnabled();
        const bool longMode = longModeEnabled();
        const std::string meshFilter = environmentString("FEM_MAIN_SCALE_CONTACT_MESH_FILTER");
        const std::string penaltyFilter =
            environmentString("FEM_MAIN_SCALE_CONTACT_PENALTY_FILTER");
        const std::string outputSubdirectory =
            environmentString("FEM_MAIN_SCALE_CONTACT_OUTPUT_SUBDIRECTORY");
        const auto penalties = buildPenaltyValues(quickMode, longMode);
        const auto meshVariants = buildMeshVariants(quickMode, longMode);

        const std::filesystem::path outputRoot =
            AppRuntimeSupport::caseOutputDirectory(
                !outputSubdirectory.empty()
                ? outputSubdirectory
                : (longMode
                    ? "main_scale_contact_penalty_mesh_sensitivity"
                    : "main_scale_contact_study"));
        std::filesystem::create_directories(outputRoot);

        std::ofstream summaryStream(outputRoot / "study_summary.csv");
        if (!summaryStream.is_open()) {
            throw std::runtime_error("Failed to open study summary CSV file");
        }

        const std::string header = buildStudyHeader();
        std::cout << header << std::endl;
        summaryStream << header << '\n';

        for (const auto& meshVariant : meshVariants) {
            if (!matchesOptionalFilter(meshFilter, meshVariant.label)) {
                continue;
            }
            for (double penalty : penalties) {
                if (!matchesOptionalNumericFilter(penaltyFilter, penalty)) {
                    continue;
                }
                const StudyResult result = runStudyCase(outputRoot, meshVariant, penalty);
                const std::string row = buildStudyRow(result);
                std::cout << row << std::endl;
                summaryStream << row << '\n';
            }
        }

        AppRuntimeSupport::runPostprocessorIfEnabled(outputRoot);
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "main_scale_contact_study_failed: " << exception.what() << std::endl;
        return 1;
    }
}
