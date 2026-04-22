#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "HyperelasticScenarioSupport.h"
#include "NeoHookeanMaterial.h"
#include "assembly.h"
#include "node.h"
#include "planeisometric/FiniteStrainQ4Element.h"

namespace {

constexpr int kNodesX = 7;
constexpr int kNodesY = 7;
constexpr double kYoungsModulus = 1.0e3;
constexpr double kPoissonsRatio = 0.30;
constexpr double kThickness = 1.0;
constexpr double kPrescribedTopEdgeDy = -0.08;
constexpr double kAugmentationParameter = 1.0e6;
constexpr double kMaximumAllowedPenetration = 1.0e-3;
constexpr double kMinimumAllowedJacobianDeterminant = 0.95;
constexpr double kCoordinateTolerance = 1.0e-9;

struct HyperelasticContactBlockCase {
    std::shared_ptr<Assembly> assembly;
    int anchorNodeId = -1;
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
        if (coordinates.x() < smallestAnchorX) {
            smallestAnchorX = coordinates.x();
            blockCase.anchorNodeId = node->getId();
        }
    }

    if (blockCase.anchorNodeId < 0) {
        throw std::runtime_error("Failed to select anchor node in hyperelastic AL smoke test");
    }

    blockCase.assembly->addFixedNode(blockCase.anchorNodeId, true, false);
    return blockCase;
}

} // namespace

int main() {
    try {
        const auto hyperelasticCase = buildHyperelasticContactBlockCase();
        const auto contactFacets =
            HyperelasticScenarioSupport::collectBoundaryFacetsByCoordinate(
                hyperelasticCase.assembly, 1, true, 1.0e-8);

        if (contactFacets.empty()) {
            std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                      << "no contact facets collected" << std::endl;
            return 1;
        }

        FEModel model;
        model.setAssembly(hyperelasticCase.assembly);
        model.setSolverTolerance(1.0e-8);
        model.setMaxIterations(25);
        AugmentedLagrangianSettings settings;
        settings.augmentationParameter = kAugmentationParameter;
        model.configureRigidPlaneAugmentedLagrangianContact(
            RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0},
            contactFacets,
            settings);

        if (!model.solveHyperelastic()) {
            std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                      << "solveHyperelastic() returned false" << std::endl;
            return 1;
        }

        const auto& metrics = model.getPerformanceMetrics();
        const auto diagnostics =
            HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);

        if (metrics.activeSetSize <= 0) {
            std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                      << "no active contact facets"
                      << " max_penetration=" << metrics.maxPenetration
                      << " contact_force_norm=" << metrics.contactForceNorm
                      << " nonlinear_iterations=" << metrics.nonlinearIterations
                      << " lambda_relative_update=" << metrics.contactStateRelativeUpdateNorm
                      << std::endl;
            return 1;
        }
        if (!(metrics.maxPenetration > 0.0) ||
            metrics.maxPenetration >= kMaximumAllowedPenetration) {
            std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                      << "penetration out of range: " << metrics.maxPenetration << std::endl;
            return 1;
        }
        if (!(metrics.maxNormalContactMultiplier > 0.0)) {
            std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                      << "non-positive normal multiplier" << std::endl;
            return 1;
        }
        if (diagnostics.minJ <= kMinimumAllowedJacobianDeterminant) {
            std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                      << "minimum det(F) too low: " << diagnostics.minJ << std::endl;
            return 1;
        }

        std::cout << "Hyperelastic augmented-Lagrangian smoke" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "candidate_contact_facets=" << contactFacets.size() << std::endl;
        std::cout << "active_contact_facets=" << metrics.activeSetSize << std::endl;
        std::cout << "active_gauss_points=" << metrics.activeContactGaussPoints << std::endl;
        std::cout << "max_penetration=" << metrics.maxPenetration << std::endl;
        std::cout << "max_normal_multiplier=" << metrics.maxNormalContactMultiplier << std::endl;
        std::cout << "lambda_relative_update="
                  << metrics.contactStateRelativeUpdateNorm << std::endl;
        std::cout << "minimum_jacobian_determinant=" << diagnostics.minJ << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "hyperelastic_augmented_lagrangian_smoke_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
