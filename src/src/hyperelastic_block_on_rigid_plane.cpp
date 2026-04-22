#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "AppRuntimeSupport.h"
#include "HyperelasticScenarioSupport.h"
#include "ResultFileExporter.h"

namespace {

constexpr int kNodesX = 7;
constexpr int kNodesY = 7;
constexpr double kYoungsModulus = 1.0e3;
constexpr double kPoissonsRatio = 0.30;
constexpr double kThickness = 1.0;
constexpr double kPrescribedTopEdgeDy = -0.08;
constexpr double kPenaltyParameter = 1.0e6;
constexpr double kCoordinateTolerance = 1.0e-9;
constexpr double kMaximumAllowedPenetration = 2.0e-2;
constexpr double kMinimumAllowedJacobianDeterminant = 0.90;

struct HyperelasticContactBlockCase {
    std::shared_ptr<Assembly> assembly;
    int anchorNodeId = -1;
    int topEdgeNodeCount = 0;
};

HyperelasticContactBlockCase buildHyperelasticContactBlockCase() {
    HyperelasticContactBlockCase blockCase;
    blockCase.assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<NeoHookeanMaterial>(
        NeoHookeanMaterial::fromYoungsModulusAndPoissonsRatio(
            1, kYoungsModulus, kPoissonsRatio, kThickness));
    blockCase.assembly->addFiniteStrainMaterial(material);

    const Eigen::Vector2d lowerCorner(-0.5, 0.05);
    const Eigen::Vector2d upperCorner(0.5, 1.05);

    std::vector<std::vector<int>> nodeGrid(
        static_cast<size_t>(kNodesY),
        std::vector<int>(static_cast<size_t>(kNodesX), 0));

    int nextNodeId = 1;
    for (int j = 0; j < kNodesY; ++j) {
        const double eta = static_cast<double>(j) / static_cast<double>(kNodesY - 1);
        const double y = lowerCorner.y() + eta * (upperCorner.y() - lowerCorner.y());
        for (int i = 0; i < kNodesX; ++i) {
            const double xi = static_cast<double>(i) / static_cast<double>(kNodesX - 1);
            const double x = lowerCorner.x() + xi * (upperCorner.x() - lowerCorner.x());
            blockCase.assembly->addNode(std::make_shared<Node>(nextNodeId, x, y));
            nodeGrid[static_cast<size_t>(j)][static_cast<size_t>(i)] = nextNodeId++;
        }
    }

    int nextElementId = 1;
    for (int j = 0; j < kNodesY - 1; ++j) {
        for (int i = 0; i < kNodesX - 1; ++i) {
            blockCase.assembly->addFiniteStrainElement(
                std::make_shared<FiniteStrainQ4Element>(
                    nextElementId++,
                    std::vector<int>{
                        nodeGrid[static_cast<size_t>(j + 1)][static_cast<size_t>(i + 1)],
                        nodeGrid[static_cast<size_t>(j + 1)][static_cast<size_t>(i)],
                        nodeGrid[static_cast<size_t>(j)][static_cast<size_t>(i)],
                        nodeGrid[static_cast<size_t>(j)][static_cast<size_t>(i + 1)]},
                    material->getId()));
        }
    }

    const auto& nodes = blockCase.assembly->getNodes();
    double maxY = -std::numeric_limits<double>::infinity();
    double smallestAnchorX = std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    for (const auto& node : nodes) {
        const Eigen::Vector2d coordinates = node->getCoordinates();
        if (std::abs(coordinates.y() - maxY) > kCoordinateTolerance) {
            continue;
        }

        blockCase.assembly->addPrescribedDisplacementY(node->getId(), kPrescribedTopEdgeDy);
        ++blockCase.topEdgeNodeCount;
        if (coordinates.x() < smallestAnchorX) {
            smallestAnchorX = coordinates.x();
            blockCase.anchorNodeId = node->getId();
        }
    }

    if (blockCase.anchorNodeId < 0) {
        throw std::runtime_error("Failed to select anchor node in hyperelastic contact block case");
    }

    blockCase.assembly->addFixedNode(blockCase.anchorNodeId, true, false);
    return blockCase;
}

double computeMinimumSignedDistanceToPlane(const FEModel& model,
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

} // namespace

int main() {
    try {
        const auto blockCase = buildHyperelasticContactBlockCase();
        const auto contactFacets =
            HyperelasticScenarioSupport::collectBoundaryFacetsByCoordinate(
                blockCase.assembly, 1, true, 1.0e-8);
        if (contactFacets.empty()) {
            std::cerr << "No candidate contact facets found for hyperelastic block case"
                      << std::endl;
            return 1;
        }

        FEModel model;
        model.setAssembly(blockCase.assembly);
        model.setSolverTolerance(1.0e-8);
        model.setMaxIterations(25);

        const RigidPlane2D plane{Eigen::Vector2d(0.0, 1.0), 0.0};
        model.configureRigidPlanePenaltyContact(plane, contactFacets, kPenaltyParameter);
        if (!model.solveHyperelastic()) {
            std::cerr << "Hyperelastic block-on-rigid-plane solve failed" << std::endl;
            return 1;
        }

        const auto& metrics = model.getPerformanceMetrics();
        const auto finiteStrainDiagnostics =
            HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);
        const double minimumSignedDistance =
            computeMinimumSignedDistanceToPlane(model, blockCase.assembly, plane);

        if (metrics.activeSetSize <= 0) {
            std::cerr << "Hyperelastic contact case did not activate any contact facets"
                      << std::endl;
            return 1;
        }
        if (!(metrics.maxPenetration > 0.0) ||
            metrics.maxPenetration >= kMaximumAllowedPenetration) {
            std::cerr << "Hyperelastic contact penetration is outside accepted range: "
                      << metrics.maxPenetration << std::endl;
            return 1;
        }
        if (!(metrics.contactForceNorm > 0.0)) {
            std::cerr << "Hyperelastic contact case produced zero contact force norm"
                      << std::endl;
            return 1;
        }
        if (finiteStrainDiagnostics.minJ <= kMinimumAllowedJacobianDeterminant) {
            std::cerr << "Hyperelastic contact case has unacceptable minimum det(F): "
                      << finiteStrainDiagnostics.minJ << std::endl;
            return 1;
        }

        ResultFileExportOptions exportOptions;
        exportOptions.outputDirectory =
            AppRuntimeSupport::caseOutputDirectory("hyperelastic_block_on_rigid_plane");
        exportOptions.baseName = "solution";
        exportOptions.extraStringMetrics = {
            {"case_name", "hyperelastic_block_on_rigid_plane"},
            {"verification_category", "hyperelastic_contact_penalty"}
        };
        exportOptions.extraNumericMetrics = {
            {"penalty_parameter", kPenaltyParameter},
            {"prescribed_top_edge_dy", kPrescribedTopEdgeDy},
            {"candidate_contact_facets", static_cast<double>(contactFacets.size())},
            {"active_contact_facets", static_cast<double>(metrics.activeSetSize)},
            {"minimum_signed_distance", minimumSignedDistance},
            {"minimum_jacobian_determinant", finiteStrainDiagnostics.minJ},
            {"maximum_jacobian_determinant", finiteStrainDiagnostics.maxJ}
        };

        const auto exportArtifacts = ResultFileExporter::exportSolution(model, exportOptions);
        AppRuntimeSupport::runPostprocessorIfEnabled(exportOptions.outputDirectory);

        std::cout << "Hyperelastic block on rigid plane" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "candidate_contact_facets=" << contactFacets.size() << std::endl;
        std::cout << "active_contact_facets=" << metrics.activeSetSize << std::endl;
        std::cout << "max_penetration=" << metrics.maxPenetration << std::endl;
        std::cout << "contact_force_norm=" << metrics.contactForceNorm << std::endl;
        std::cout << "minimum_signed_distance=" << minimumSignedDistance << std::endl;
        std::cout << "minimum_jacobian_determinant=" << finiteStrainDiagnostics.minJ
                  << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        std::cout << "results_vtu=" << exportArtifacts.vtuPath.string() << std::endl;
        std::cout << "results_metrics_json=" << exportArtifacts.metricsJsonPath.string()
                  << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "hyperelastic_block_on_rigid_plane_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
