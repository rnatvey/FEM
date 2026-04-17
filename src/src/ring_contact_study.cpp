#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

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

struct StudyResult {
    std::string meshLabel;
    double penalty = 0.0;
    bool success = false;
    int contactFacetCount = 0;
    double minimumSignedDistance = std::numeric_limits<double>::infinity();
    MeshGenerator::RingMeshDiagnostics meshDiagnostics;
    FEModel::PerformanceMetrics metrics;
};

std::string makeCaseToken(const std::string& meshLabel, double penalty) {
    std::ostringstream stream;
    stream << meshLabel << "_penalty_" << std::llround(penalty);
    return stream.str();
}

std::string buildStudyHeader() {
    return "mesh,penalty,success,candidate_facets,active_facets,nonlinear_iterations,"
           "max_penetration,minimum_signed_distance,contact_force_norm,"
           "assembly_time_seconds,solve_time_seconds,equilibrium_residual_norm,"
           "linear_residual_norm,min_radial_step,max_radial_step,min_angular_step,"
           "max_angular_step,min_aspect_ratio,max_aspect_ratio,linear_solver_backend";
}

std::string buildStudyRow(const StudyResult& result) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(8)
           << result.meshLabel << ','
           << result.penalty << ','
           << std::boolalpha << result.success << ','
           << result.contactFacetCount << ','
           << result.metrics.activeSetSize << ','
           << result.metrics.nonlinearIterations << ','
           << result.metrics.maxPenetration << ','
           << result.minimumSignedDistance << ','
           << result.metrics.contactForceNorm << ','
           << result.metrics.assemblyTimeSeconds << ','
           << result.metrics.solveTimeSeconds << ','
           << result.metrics.equilibriumResidualNorm << ','
           << result.metrics.linearResidualNorm << ','
           << result.meshDiagnostics.minRadialStep << ','
           << result.meshDiagnostics.maxRadialStep << ','
           << result.meshDiagnostics.minAngularStep << ','
           << result.meshDiagnostics.maxAngularStep << ','
           << result.meshDiagnostics.minAspectRatio << ','
           << result.meshDiagnostics.maxAspectRatio << ','
           << result.metrics.linearSolverBackend;
    return stream.str();
}

std::shared_ptr<Assembly> buildRingAssembly(bool graded,
    MeshGenerator::RingMeshDiagnostics& diagnostics,
    int& contactFacetCount,
    std::vector<ContactFacet>& contactFacets,
    RigidPlane2D& rigidPlane) {
    auto assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
    assembly->addMaterial(material);

    MeshGenerator meshGenerator(assembly);

    MeshGenerator::TireContactAnalysisControl control;
    control.mesh.center = kTireCenter;
    control.mesh.innerRadius = kTireInnerRadius;
    control.mesh.outerRadius = kTireOuterRadius;
    control.mesh.startAngle = kTireStartAngle;
    control.mesh.endAngle = kTireEndAngle;
    control.mesh.radialLayers = kTireRadialLayers;
    control.mesh.circumferentialNodes = kTireCircumferentialNodes;
    control.mesh.materialId = material->getId();
    control.mesh.refineCircumferentiallyNearContact = graded;
    control.mesh.refineRadiallyToOuterSurface = graded;
    control.mesh.expectedContactCenterAngle = kContactCenterAngle;
    control.mesh.expectedContactHalfAngle = kContactHalfAngle;
    control.mesh.circumferentialRefinementStrength = 6.0;
    control.mesh.radialRefinementStrength = 2.5;
    control.mesh.candidateFacetWindowScale = 3.0;
    control.rigidPlane = RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0};
    control.prescribedInnerBoundaryDy = -0.12;

    const auto analysisSetup = meshGenerator.generateTireContactAnalysisSetup(control);
    diagnostics = analysisSetup.mesh.diagnostics;
    contactFacets = analysisSetup.mesh.candidateContactFacets;
    contactFacetCount = static_cast<int>(contactFacets.size());
    rigidPlane = analysisSetup.rigidPlane;

    return assembly;
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
    const std::string& meshLabel,
    bool graded,
    double penalty) {
    StudyResult result;
    result.meshLabel = meshLabel;
    result.penalty = penalty;

    std::vector<ContactFacet> contactFacets;
    RigidPlane2D rigidPlane;
    auto assembly = buildRingAssembly(
        graded,
        result.meshDiagnostics,
        result.contactFacetCount,
        contactFacets,
        rigidPlane);

    FEModel model;
    model.setAssembly(assembly);
    model.setMaxIterations(40);
    model.setSolverTolerance(1.0e-8);
    model.configureRigidPlaneContact(rigidPlane, contactFacets, penalty);

    result.success = model.solveContact();
    result.metrics = model.getPerformanceMetrics();
    if (result.success) {
        result.minimumSignedDistance = computeMinimumSignedDistance(model, assembly, rigidPlane);
    }

    const std::filesystem::path caseOutputDirectory =
        outputRoot / makeCaseToken(meshLabel, penalty);
    ResultFileExportOptions exportOptions;
    exportOptions.outputDirectory = caseOutputDirectory;
    exportOptions.baseName = "solution";
    exportOptions.extraStringMetrics = {
        {"mesh_label", meshLabel},
        {"geometry_family", "tire_ring"}
    };
    exportOptions.extraNumericMetrics = {
        {"penalty_parameter", penalty},
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
        {"rigid_plane_normal_x", rigidPlane.normal.x()},
        {"rigid_plane_normal_y", rigidPlane.normal.y()},
        {"rigid_plane_offset", rigidPlane.offset},
        {"mesh_min_radial_step", result.meshDiagnostics.minRadialStep},
        {"mesh_max_radial_step", result.meshDiagnostics.maxRadialStep},
        {"mesh_min_angular_step", result.meshDiagnostics.minAngularStep},
        {"mesh_max_angular_step", result.meshDiagnostics.maxAngularStep},
        {"mesh_min_aspect_ratio", result.meshDiagnostics.minAspectRatio},
        {"mesh_max_aspect_ratio", result.meshDiagnostics.maxAspectRatio}
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
        const std::vector<double> penalties = {1.0e6, 1.0e7, 1.0e8};
        const std::filesystem::path outputRoot =
            std::filesystem::path("results") / "ring_contact_study";
        std::filesystem::create_directories(outputRoot);
        std::ofstream summaryStream(outputRoot / "study_summary.csv");
        if (!summaryStream.is_open()) {
            throw std::runtime_error("Failed to open study summary CSV file");
        }

        const std::string header = buildStudyHeader();
        std::cout << header << std::endl;
        summaryStream << header << '\n';

        for (double penalty : penalties) {
            const auto row = buildStudyRow(runStudyCase(outputRoot, "uniform", false, penalty));
            std::cout << row << std::endl;
            summaryStream << row << '\n';
        }
        for (double penalty : penalties) {
            const auto row = buildStudyRow(runStudyCase(outputRoot, "graded", true, penalty));
            std::cout << row << std::endl;
            summaryStream << row << '\n';
        }

        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "ring_contact_study_failed: " << exception.what() << std::endl;
        return 1;
    }
}
