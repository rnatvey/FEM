#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "FEMModel.h"
#include "assembly.h"
#include "constants.h"
#include "material.h"
#include "meshgenerator.h"

using namespace Constants;

namespace {

struct StudyResult {
    std::string meshLabel;
    double penalty = 0.0;
    bool success = false;
    int contactFacetCount = 0;
    double minimumSignedDistance = std::numeric_limits<double>::infinity();
    MeshGenerator::RingMeshDiagnostics meshDiagnostics;
    FEModel::PerformanceMetrics metrics;
};

std::shared_ptr<Assembly> buildRingAssembly(bool graded,
    MeshGenerator::RingMeshDiagnostics& diagnostics,
    int& contactFacetCount,
    std::vector<ContactFacet>& contactFacets,
    RigidPlane2D& rigidPlane) {
    auto assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
    assembly->addMaterial(material);

    MeshGenerator meshGenerator(assembly);
    const Eigen::Vector2d center(0.0, 0.58);
    const double innerRadius = 0.25;
    const double outerRadius = 0.50;
    const double startAngle = -120.0 * DEG_TO_RAD;
    const double endAngle = -60.0 * DEG_TO_RAD;

    MeshGenerator::TireContactAnalysisControl control;
    control.mesh.center = center;
    control.mesh.innerRadius = innerRadius;
    control.mesh.outerRadius = outerRadius;
    control.mesh.startAngle = startAngle;
    control.mesh.endAngle = endAngle;
    control.mesh.radialLayers = 9;
    control.mesh.circumferentialNodes = 61;
    control.mesh.materialId = material->getId();
    control.mesh.refineCircumferentiallyNearContact = graded;
    control.mesh.refineRadiallyToOuterSurface = graded;
    control.mesh.expectedContactCenterAngle = -90.0 * DEG_TO_RAD;
    control.mesh.expectedContactHalfAngle = 12.0 * DEG_TO_RAD;
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

StudyResult runStudyCase(const std::string& meshLabel, bool graded, double penalty) {
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

    return result;
}

void printStudyRow(const StudyResult& result) {
    std::cout << std::fixed << std::setprecision(8)
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
              << result.metrics.linearSolverBackend
              << std::endl;
}

} // namespace

int main() {
    try {
        const std::vector<double> penalties = {1.0e6, 1.0e7, 1.0e8};

        std::cout << "mesh,penalty,success,candidate_facets,active_facets,nonlinear_iterations,"
                     "max_penetration,minimum_signed_distance,contact_force_norm,"
                     "assembly_time_seconds,solve_time_seconds,equilibrium_residual_norm,"
                     "linear_residual_norm,min_radial_step,max_radial_step,min_angular_step,"
                     "max_angular_step,min_aspect_ratio,max_aspect_ratio,linear_solver_backend"
                  << std::endl;

        for (double penalty : penalties) {
            printStudyRow(runStudyCase("uniform", false, penalty));
        }
        for (double penalty : penalties) {
            printStudyRow(runStudyCase("graded", true, penalty));
        }

        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "ring_contact_study_failed: " << exception.what() << std::endl;
        return 1;
    }
}
