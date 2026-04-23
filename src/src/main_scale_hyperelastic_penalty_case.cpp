#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
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
    const auto variantIt = std::find_if(variants.begin(), variants.end(),
        [requestedLabel](const MeshVariant& variant) {
            return variant.label == requestedLabel;
        });
    if (variantIt == variants.end()) {
        throw std::invalid_argument(
            "Unknown hyperelastic main-scale mesh label: " + std::string(requestedLabel));
    }
    return *variantIt;
}

MeshGenerator::TireContactAnalysisControl buildControl(const MeshVariant& meshVariant) {
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
    control.prescribeInnerBoundaryX = false;
    control.prescribeInnerBoundaryY = true;
    control.prescribedInnerBoundaryDy = kPrescribedInnerBoundaryDy;
    control.addInnerBoundaryAnchor = true;
    control.anchorFixX = true;
    control.anchorFixY = false;
    control.anchorSelectMinimumX = true;
    return control;
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

std::string makeCaseName(std::string_view meshLabel, double penalty, int loadSteps) {
    std::ostringstream stream;
    stream << meshLabel
           << "_penalty_" << std::llround(penalty)
           << "_steps_" << loadSteps;
    return stream.str();
}

} // namespace

int main() {
    try {
        const std::string meshLabel =
            environmentString("FEM_MAIN_SCALE_HYPERELASTIC_PENALTY_MESH").empty()
            ? "contact_focused_very_coarse"
            : environmentString("FEM_MAIN_SCALE_HYPERELASTIC_PENALTY_MESH");
        const double penalty =
            environmentDouble("FEM_MAIN_SCALE_HYPERELASTIC_PENALTY_VALUE", 5.0e2);
        const int loadSteps =
            std::max(2, environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_LOAD_STEPS", 12));
        const int maxAdaptiveLoadSteps =
            std::max(loadSteps,
                environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_MAX_ADAPTIVE_LOAD_STEPS", 120));
        const int maxIterations =
            std::max(8, environmentInt("FEM_MAIN_SCALE_HYPERELASTIC_MAX_ITERATIONS", 30));
        const std::string outputSubdirectory =
            environmentString("FEM_MAIN_SCALE_HYPERELASTIC_OUTPUT_SUBDIRECTORY");

        const MeshVariant meshVariant = resolveMeshVariant(meshLabel);
        const auto control = buildControl(meshVariant);
        const auto tireCase = HyperelasticScenarioSupport::buildFiniteStrainTireContactCase(
            control, kYoungsModulus, kPoissonsRatio, kThickness);

        FEModel model;
        model.setAssembly(tireCase.assembly);
        model.setSolverTolerance(1.0e-8);
        model.setMaxIterations(maxIterations);
        model.setHyperelasticLoadSteps(loadSteps);
        model.setMaxAdaptiveHyperelasticLoadSteps(maxAdaptiveLoadSteps);
        model.configureRigidPlanePenaltyContact(
            tireCase.rigidPlane, tireCase.mesh.candidateContactFacets, penalty);

        const bool success = model.solveHyperelastic();
        const auto metrics = model.getPerformanceMetrics();
        const auto diagnostics =
            HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);
        const double minimumSignedDistance =
            computeMinimumSignedDistance(model, tireCase.assembly, tireCase.rigidPlane);

        const std::string caseName = makeCaseName(meshVariant.label, penalty, loadSteps);
        const std::filesystem::path outputDirectory =
            AppRuntimeSupport::caseOutputDirectory(
                !outputSubdirectory.empty()
                ? outputSubdirectory
                : "main_scale_hyperelastic_penalty_case");
        std::filesystem::create_directories(outputDirectory);

        ResultFileExportOptions exportOptions;
        exportOptions.outputDirectory = outputDirectory;
        exportOptions.baseName = "solution";
        exportOptions.extraStringMetrics = {
            {"case_name", caseName},
            {"mesh_label", meshVariant.label},
            {"geometry_family", "tire_ring"},
            {"study_family", "main_scale_hyperelastic_penalty_case"},
            {"contact_method", "penalty"},
            {"contact_parameter_name", "penalty_parameter"}
        };
        exportOptions.extraNumericMetrics = {
            {"contact_parameter_value", penalty},
            {"penalty_parameter", penalty},
            {"hyperelastic_load_steps", static_cast<double>(loadSteps)},
            {"max_adaptive_hyperelastic_load_steps", static_cast<double>(maxAdaptiveLoadSteps)},
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
            {"rigid_plane_normal_x", tireCase.rigidPlane.normal.x()},
            {"rigid_plane_normal_y", tireCase.rigidPlane.normal.y()},
            {"rigid_plane_offset", tireCase.rigidPlane.offset},
            {"inner_boundary_prescribed_dy", kPrescribedInnerBoundaryDy},
            {"plane_initial_gap", kPlaneInitialGap},
            {"mesh_circumferential_refinement_strength",
                meshVariant.circumferentialRefinementStrength},
            {"mesh_radial_refinement_strength", meshVariant.radialRefinementStrength},
            {"mesh_circumferential_localization_power",
                meshVariant.circumferentialLocalizationPower},
            {"mesh_min_radial_step", tireCase.mesh.diagnostics.minRadialStep},
            {"mesh_max_radial_step", tireCase.mesh.diagnostics.maxRadialStep},
            {"mesh_min_angular_step", tireCase.mesh.diagnostics.minAngularStep},
            {"mesh_max_angular_step", tireCase.mesh.diagnostics.maxAngularStep},
            {"mesh_min_outer_arc_step", tireCase.mesh.diagnostics.minOuterArcStep},
            {"mesh_max_outer_arc_step", tireCase.mesh.diagnostics.maxOuterArcStep},
            {"mesh_min_aspect_ratio", tireCase.mesh.diagnostics.minAspectRatio},
            {"mesh_max_aspect_ratio", tireCase.mesh.diagnostics.maxAspectRatio},
            {"minimum_signed_distance", minimumSignedDistance},
            {"minimum_jacobian_determinant", diagnostics.minJ},
            {"maximum_jacobian_determinant", diagnostics.maxJ},
            {"maximum_strain_energy_density", diagnostics.maxStrainEnergyDensity},
            {"gauss_point_count", static_cast<double>(diagnostics.gaussPointCount)}
        };

        ResultFileExporter::exportSolution(model, exportOptions);

        std::cout << std::boolalpha
                  << "case_name=" << caseName
                  << ", success=" << success
                  << ", nodes=" << tireCase.assembly->getNodeCount()
                  << ", elements="
                  << tireCase.assembly->getFiniteStrainElements().size()
                  << ", candidate_facets="
                  << tireCase.mesh.candidateContactFacets.size()
                  << ", load_steps=" << metrics.loadSteps
                  << ", converged_load_steps=" << metrics.convergedLoadSteps
                  << ", nonlinear_iterations=" << metrics.nonlinearIterations
                  << ", max_penetration=" << metrics.maxPenetration
                  << ", minimum_signed_distance=" << minimumSignedDistance
                  << ", min_J=" << diagnostics.minJ
                  << ", max_J=" << diagnostics.maxJ
                  << ", strain_energy=" << metrics.strainEnergy
                  << ", total_time_seconds=" << metrics.totalTimeSeconds
                  << ", linear_solver_backend=" << metrics.linearSolverBackend
                  << std::endl;

        AppRuntimeSupport::runPostprocessorIfEnabled(outputDirectory);
        return success ? 0 : 1;
    }
    catch (const std::exception& exception) {
        std::cerr << "main_scale_hyperelastic_penalty_case_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
