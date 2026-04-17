#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>

#include "FEMModel.h"
#include "ResultFileExporter.h"
#include "assembly.h"
#include "constants.h"
#include "material.h"
#include "meshgenerator.h"

using namespace Constants;

namespace {

const Eigen::Vector2d kTireCenter(0.0, 0.58);
constexpr double kTireInnerRadius = 0.25;
constexpr double kTireOuterRadius = 0.50;
constexpr double kTireStartAngle = -120.0 * DEG_TO_RAD;
constexpr double kTireEndAngle = -60.0 * DEG_TO_RAD;
constexpr int kTireRadialLayers = 9;
constexpr int kTireCircumferentialNodes = 61;
constexpr double kContactCenterAngle = -90.0 * DEG_TO_RAD;
constexpr double kContactHalfAngle = 12.0 * DEG_TO_RAD;
constexpr double kPenaltyParameter = 1.0e7;

struct TireContactCase {
    std::shared_ptr<Assembly> assembly;
    MeshGenerator::RingMeshDiagnostics diagnostics;
    std::vector<ContactFacet> contactFacets;
    RigidPlane2D rigidPlane;
};

TireContactCase buildTireContactCase() {
    TireContactCase tireCase;
    tireCase.assembly = std::make_shared<Assembly>();

    auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
    tireCase.assembly->addMaterial(material);

    MeshGenerator meshGenerator(tireCase.assembly);

    MeshGenerator::TireContactAnalysisControl control;
    control.mesh.center = kTireCenter;
    control.mesh.innerRadius = kTireInnerRadius;
    control.mesh.outerRadius = kTireOuterRadius;
    control.mesh.startAngle = kTireStartAngle;
    control.mesh.endAngle = kTireEndAngle;
    control.mesh.radialLayers = kTireRadialLayers;
    control.mesh.circumferentialNodes = kTireCircumferentialNodes;
    control.mesh.materialId = material->getId();
    control.mesh.refineCircumferentiallyNearContact = true;
    control.mesh.refineRadiallyToOuterSurface = true;
    control.mesh.expectedContactCenterAngle = kContactCenterAngle;
    control.mesh.expectedContactHalfAngle = kContactHalfAngle;
    control.mesh.circumferentialRefinementStrength = 6.0;
    control.mesh.radialRefinementStrength = 2.5;
    control.mesh.candidateFacetWindowScale = 3.0;
    control.rigidPlane = RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0};
    control.prescribedInnerBoundaryDy = -0.12;

    const auto setup = meshGenerator.generateTireContactAnalysisSetup(control);
    tireCase.diagnostics = setup.mesh.diagnostics;
    tireCase.contactFacets = setup.mesh.candidateContactFacets;
    tireCase.rigidPlane = setup.rigidPlane;

    return tireCase;
}

double computeMinimumSignedDistance(const FEModel& model) {
    const auto signedDistances = model.getNodalContactSignedDistances();
    if (signedDistances.empty()) {
        return 0.0;
    }
    return *std::min_element(signedDistances.begin(), signedDistances.end());
}

} // namespace

int main() {
    try {
        auto tireCase = buildTireContactCase();

        FEModel model;
        model.setAssembly(tireCase.assembly);
        model.setSolverTolerance(1.0e-8);
        model.setMaxIterations(40);
        model.configureRigidPlaneContact(tireCase.rigidPlane, tireCase.contactFacets, kPenaltyParameter);

        const bool success = model.solveContact();
        if (!success) {
            std::cerr << "Tire contact single-case solve failed" << std::endl;
            return 1;
        }

        ResultFileExportOptions exportOptions;
        exportOptions.outputDirectory = std::filesystem::path("results") / "tire_contact_single_case";
        exportOptions.baseName = "solution";
        exportOptions.extraStringMetrics = {
            {"case_name", "tire_contact_single_case"},
            {"mesh_label", "graded"},
            {"geometry_family", "tire_ring"}
        };
        exportOptions.extraNumericMetrics = {
            {"penalty_parameter", kPenaltyParameter},
            {"candidate_contact_facets", static_cast<double>(tireCase.contactFacets.size())},
            {"ring_center_x", kTireCenter.x()},
            {"ring_center_y", kTireCenter.y()},
            {"ring_inner_radius", kTireInnerRadius},
            {"ring_outer_radius", kTireOuterRadius},
            {"ring_start_angle_rad", kTireStartAngle},
            {"ring_end_angle_rad", kTireEndAngle},
            {"ring_contact_center_angle_rad", kContactCenterAngle},
            {"ring_contact_half_angle_rad", kContactHalfAngle},
            {"ring_radial_layers", static_cast<double>(kTireRadialLayers)},
            {"ring_circumferential_nodes", static_cast<double>(kTireCircumferentialNodes)},
            {"rigid_plane_normal_x", tireCase.rigidPlane.normal.x()},
            {"rigid_plane_normal_y", tireCase.rigidPlane.normal.y()},
            {"rigid_plane_offset", tireCase.rigidPlane.offset},
            {"mesh_min_radial_step", tireCase.diagnostics.minRadialStep},
            {"mesh_max_radial_step", tireCase.diagnostics.maxRadialStep},
            {"mesh_min_angular_step", tireCase.diagnostics.minAngularStep},
            {"mesh_max_angular_step", tireCase.diagnostics.maxAngularStep},
            {"mesh_min_aspect_ratio", tireCase.diagnostics.minAspectRatio},
            {"mesh_max_aspect_ratio", tireCase.diagnostics.maxAspectRatio}
        };

        const auto exportArtifacts = ResultFileExporter::exportSolution(model, exportOptions);
        const auto& metrics = model.getPerformanceMetrics();

        std::cout << "Tire contact single-case example" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "nodes=" << tireCase.assembly->getNodeCount() << std::endl;
        std::cout << "elements=" << tireCase.assembly->getElementCount() << std::endl;
        std::cout << "candidate_contact_facets=" << tireCase.contactFacets.size() << std::endl;
        std::cout << "active_contact_facets=" << metrics.activeSetSize << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        std::cout << "max_penetration=" << metrics.maxPenetration << std::endl;
        std::cout << "minimum_signed_distance=" << computeMinimumSignedDistance(model) << std::endl;
        std::cout << "assembly_time_seconds=" << metrics.assemblyTimeSeconds << std::endl;
        std::cout << "solve_time_seconds=" << metrics.solveTimeSeconds << std::endl;
        std::cout << "total_time_seconds=" << metrics.totalTimeSeconds << std::endl;
        std::cout << "matrix_nonzeros=" << metrics.matrixNonZeros << std::endl;
        std::cout << "linear_solver_backend=" << metrics.linearSolverBackend << std::endl;
        std::cout << "results_vtu=" << exportArtifacts.vtuPath.string() << std::endl;
        std::cout << "results_metrics_json=" << exportArtifacts.metricsJsonPath.string() << std::endl;

        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "Tire contact single-case example failed: " << exception.what() << std::endl;
        return 1;
    }
}
