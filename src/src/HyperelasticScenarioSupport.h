#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "FEMModel.h"
#include "NeoHookeanMaterial.h"
#include "assembly.h"
#include "material.h"
#include "node.h"
#include "ContactTypes.h"
#include "planeisometric/FiniteStrainQ4Element.h"
#include "planeisometric/Planeisoparametric.h"
#include "meshgenerator.h"

namespace HyperelasticScenarioSupport {

constexpr double kCoordinateTolerance = 1.0e-9;

struct BlockCase {
    std::shared_ptr<Assembly> assembly;
    int anchorNodeId = -1;
    int topEdgeNodeCount = 0;
};

struct TireContactCase {
    std::shared_ptr<Assembly> assembly;
    MeshGenerator::TireContactMeshResult mesh;
    RigidPlane2D rigidPlane = RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0};
    std::vector<int> innerBoundaryNodeIds;
    int anchorNodeId = -1;
};

inline BlockCase buildFiniteStrainBlockCase(int nodesX,
    int nodesY,
    const Eigen::Vector2d& lowerCorner,
    const Eigen::Vector2d& upperCorner,
    double youngsModulus,
    double poissonsRatio,
    double thickness,
    double prescribedTopEdgeDy) {
    if (nodesX < 2 || nodesY < 2) {
        throw std::invalid_argument("Finite-strain block case requires at least a 2x2 node grid");
    }

    BlockCase blockCase;
    blockCase.assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<NeoHookeanMaterial>(
        NeoHookeanMaterial::fromYoungsModulusAndPoissonsRatio(
            1, youngsModulus, poissonsRatio, thickness));
    blockCase.assembly->addFiniteStrainMaterial(material);

    std::vector<std::vector<int>> nodeGrid(
        static_cast<size_t>(nodesY),
        std::vector<int>(static_cast<size_t>(nodesX), 0));

    int nextNodeId = 1;
    for (int j = 0; j < nodesY; ++j) {
        const double eta = static_cast<double>(j) / static_cast<double>(nodesY - 1);
        const double y = lowerCorner.y() + eta * (upperCorner.y() - lowerCorner.y());
        for (int i = 0; i < nodesX; ++i) {
            const double xi = static_cast<double>(i) / static_cast<double>(nodesX - 1);
            const double x = lowerCorner.x() + xi * (upperCorner.x() - lowerCorner.x());
            blockCase.assembly->addNode(std::make_shared<Node>(nextNodeId, x, y));
            nodeGrid[static_cast<size_t>(j)][static_cast<size_t>(i)] = nextNodeId++;
        }
    }

    int nextElementId = 1;
    for (int j = 0; j < nodesY - 1; ++j) {
        for (int i = 0; i < nodesX - 1; ++i) {
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
    double minY = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        minY = std::min(minY, node->getCoordinates().y());
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    double smallestAnchorX = std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        const Eigen::Vector2d coordinates = node->getCoordinates();
        if (std::abs(coordinates.y() - maxY) > kCoordinateTolerance) {
            if (std::abs(coordinates.y() - minY) <= kCoordinateTolerance) {
                blockCase.assembly->addFixedNode(node->getId(), false, true);
                if (coordinates.x() < smallestAnchorX) {
                    smallestAnchorX = coordinates.x();
                    blockCase.anchorNodeId = node->getId();
                }
            }
            continue;
        }

        blockCase.assembly->addPrescribedDisplacementY(node->getId(), prescribedTopEdgeDy);
        ++blockCase.topEdgeNodeCount;
    }

    if (blockCase.anchorNodeId < 0) {
        throw std::runtime_error("Failed to select an anchor node for the finite-strain block");
    }

    blockCase.assembly->addFixedNode(blockCase.anchorNodeId, true, false);
    return blockCase;
}

inline BlockCase buildLinearBlockCase(int nodesX,
    int nodesY,
    const Eigen::Vector2d& lowerCorner,
    const Eigen::Vector2d& upperCorner,
    double youngsModulus,
    double poissonsRatio,
    double thickness,
    double prescribedTopEdgeDy) {
    if (nodesX < 2 || nodesY < 2) {
        throw std::invalid_argument("Linear block case requires at least a 2x2 node grid");
    }

    BlockCase blockCase;
    blockCase.assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<Material>(1, youngsModulus, poissonsRatio, thickness);
    blockCase.assembly->addMaterial(material);

    std::vector<std::vector<int>> nodeGrid(
        static_cast<size_t>(nodesY),
        std::vector<int>(static_cast<size_t>(nodesX), 0));

    int nextNodeId = 1;
    for (int j = 0; j < nodesY; ++j) {
        const double eta = static_cast<double>(j) / static_cast<double>(nodesY - 1);
        const double y = lowerCorner.y() + eta * (upperCorner.y() - lowerCorner.y());
        for (int i = 0; i < nodesX; ++i) {
            const double xi = static_cast<double>(i) / static_cast<double>(nodesX - 1);
            const double x = lowerCorner.x() + xi * (upperCorner.x() - lowerCorner.x());
            blockCase.assembly->addNode(std::make_shared<Node>(nextNodeId, x, y));
            nodeGrid[static_cast<size_t>(j)][static_cast<size_t>(i)] = nextNodeId++;
        }
    }

    int nextElementId = 1;
    for (int j = 0; j < nodesY - 1; ++j) {
        for (int i = 0; i < nodesX - 1; ++i) {
            blockCase.assembly->addElement(
                std::make_shared<PlaneIsoparametricElement>(
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
    double minY = std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        minY = std::min(minY, node->getCoordinates().y());
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    double smallestAnchorX = std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        const Eigen::Vector2d coordinates = node->getCoordinates();
        if (std::abs(coordinates.y() - maxY) > kCoordinateTolerance) {
            if (std::abs(coordinates.y() - minY) <= kCoordinateTolerance) {
                blockCase.assembly->addFixedNode(node->getId(), false, true);
                if (coordinates.x() < smallestAnchorX) {
                    smallestAnchorX = coordinates.x();
                    blockCase.anchorNodeId = node->getId();
                }
            }
            continue;
        }

        blockCase.assembly->addPrescribedDisplacementY(node->getId(), prescribedTopEdgeDy);
        ++blockCase.topEdgeNodeCount;
    }

    if (blockCase.anchorNodeId < 0) {
        throw std::runtime_error("Failed to select an anchor node for the linear block");
    }

    blockCase.assembly->addFixedNode(blockCase.anchorNodeId, true, false);
    return blockCase;
}

inline double computeMaximumTopEdgeDisplacementError(const FEModel& model,
    const std::shared_ptr<Assembly>& assembly,
    double prescribedTopEdgeDy) {
    const auto& nodes = assembly->getNodes();
    const auto nodalDisplacements = model.getNodalDisplacements();

    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    double maxError = 0.0;
    for (size_t i = 0; i < nodes.size() && i < nodalDisplacements.size(); ++i) {
        if (std::abs(nodes[i]->getCoordinates().y() - maxY) > kCoordinateTolerance) {
            continue;
        }
        maxError = std::max(maxError, std::abs(nodalDisplacements[i].y() - prescribedTopEdgeDy));
    }

    return maxError;
}

inline Eigen::Vector2d computeReactionForceSum(const FEModel& model) {
    Eigen::Vector2d reactionSum = Eigen::Vector2d::Zero();
    for (const auto& nodalReaction : model.getNodalReactionForces()) {
        reactionSum += nodalReaction;
    }
    return reactionSum;
}

struct FiniteStrainDiagnostics {
    double minJ = 1.0;
    double maxJ = 1.0;
    double maxStrainEnergyDensity = 0.0;
    double maxSecondPiolaStressNorm = 0.0;
    int gaussPointCount = 0;
};

inline FiniteStrainDiagnostics evaluateFiniteStrainDiagnostics(const FEModel& model) {
    FiniteStrainDiagnostics diagnostics;
    const auto responses = model.evaluateCurrentFiniteStrainElementResponses();
    bool firstGaussPoint = true;
    for (const auto& response : responses) {
        for (const auto& gaussPoint : response.gaussPointResults) {
            if (firstGaussPoint) {
                diagnostics.minJ = gaussPoint.jacobianDeterminant;
                diagnostics.maxJ = gaussPoint.jacobianDeterminant;
                firstGaussPoint = false;
            }
            diagnostics.minJ = std::min(diagnostics.minJ, gaussPoint.jacobianDeterminant);
            diagnostics.maxJ = std::max(diagnostics.maxJ, gaussPoint.jacobianDeterminant);
            diagnostics.maxStrainEnergyDensity =
                std::max(diagnostics.maxStrainEnergyDensity, gaussPoint.strainEnergyDensity);
            diagnostics.maxSecondPiolaStressNorm =
                std::max(diagnostics.maxSecondPiolaStressNorm,
                    gaussPoint.secondPiolaStressVoigt.norm());
            ++diagnostics.gaussPointCount;
        }
    }

    return diagnostics;
}

inline std::vector<ContactFacet> collectBoundaryFacetsByCoordinate(
    const std::shared_ptr<Assembly>& assembly,
    int axis,
    bool selectMinimum,
    double tolerance = 1.0e-8) {
    if (!assembly || assembly->getNodes().empty()) {
        return {};
    }

    double boundaryCoordinate = assembly->getNodes().front()->getCoordinates()(axis);
    for (const auto& node : assembly->getNodes()) {
        boundaryCoordinate = selectMinimum
            ? std::min(boundaryCoordinate, node->getCoordinates()(axis))
            : std::max(boundaryCoordinate, node->getCoordinates()(axis));
    }

    static constexpr int surfaceNodePairs[4][2] = {
        {2, 3},
        {0, 3},
        {0, 1},
        {1, 2}
    };

    std::vector<ContactFacet> facets;
    auto facetBelongsToBoundary = [&](const std::vector<int>& nodeIds, int surfaceIndex) {
        const int localNodeA = surfaceNodePairs[surfaceIndex][0];
        const int localNodeB = surfaceNodePairs[surfaceIndex][1];
        const auto nodeA = assembly->getNode(nodeIds[static_cast<size_t>(localNodeA)]);
        const auto nodeB = assembly->getNode(nodeIds[static_cast<size_t>(localNodeB)]);
        if (!nodeA || !nodeB) {
            return false;
        }

        const double coordinateA = nodeA->getCoordinates()(axis);
        const double coordinateB = nodeB->getCoordinates()(axis);
        return std::abs(coordinateA - boundaryCoordinate) <= tolerance &&
            std::abs(coordinateB - boundaryCoordinate) <= tolerance;
    };

    const auto& finiteStrainElements = assembly->getFiniteStrainElements();
    if (!finiteStrainElements.empty()) {
        for (const auto& element : finiteStrainElements) {
            if (!element) {
                continue;
            }
            for (int surfaceIndex = 0; surfaceIndex < 4; ++surfaceIndex) {
                if (facetBelongsToBoundary(element->getNodeIds(), surfaceIndex)) {
                    facets.push_back(ContactFacet{element->getId(), surfaceIndex});
                }
            }
        }
        return facets;
    }

    for (const auto& element : assembly->getElements()) {
        if (!element) {
            continue;
        }
        for (int surfaceIndex = 0; surfaceIndex < 4; ++surfaceIndex) {
            if (facetBelongsToBoundary(element->getNodeIds(), surfaceIndex)) {
                facets.push_back(ContactFacet{element->getId(), surfaceIndex});
            }
        }
    }

    return facets;
}

inline TireContactCase buildFiniteStrainTireContactCase(
    const MeshGenerator::TireContactAnalysisControl& control,
    double youngsModulus,
    double poissonsRatio,
    double thickness) {
    auto temporaryAssembly = std::make_shared<Assembly>();
    auto temporaryMaterial = std::make_shared<Material>(
        1, youngsModulus, poissonsRatio, thickness);
    temporaryAssembly->addMaterial(temporaryMaterial);

    MeshGenerator temporaryMeshGenerator(temporaryAssembly);
    const auto setup = temporaryMeshGenerator.generateTireContactAnalysisSetup(control);

    TireContactCase tireCase;
    tireCase.assembly = std::make_shared<Assembly>();
    auto finiteStrainMaterial = std::make_shared<NeoHookeanMaterial>(
        NeoHookeanMaterial::fromYoungsModulusAndPoissonsRatio(
            1, youngsModulus, poissonsRatio, thickness));
    tireCase.assembly->addFiniteStrainMaterial(finiteStrainMaterial);

    std::vector<std::shared_ptr<Node>> nodes;
    nodes.reserve(temporaryAssembly->getNodes().size());
    for (const auto& node : temporaryAssembly->getNodes()) {
        nodes.push_back(std::make_shared<Node>(
            node->getId(),
            node->getCoordinates().x(),
            node->getCoordinates().y()));
    }
    tireCase.assembly->addNodes(nodes);

    std::vector<std::shared_ptr<FiniteStrainQ4Element>> finiteStrainElements;
    finiteStrainElements.reserve(temporaryAssembly->getElements().size());
    for (const auto& element : temporaryAssembly->getElements()) {
        finiteStrainElements.push_back(std::make_shared<FiniteStrainQ4Element>(
            element->getId(), element->getNodeIds(), finiteStrainMaterial->getId()));
    }
    tireCase.assembly->addFiniteStrainElements(finiteStrainElements);

    for (int nodeId : setup.innerBoundaryNodeIds) {
        if (control.prescribeInnerBoundaryX && control.prescribeInnerBoundaryY) {
            tireCase.assembly->addPrescribedDisplacement(
                nodeId,
                control.prescribedInnerBoundaryDx,
                control.prescribedInnerBoundaryDy);
        }
        else if (control.prescribeInnerBoundaryX) {
            tireCase.assembly->addPrescribedDisplacementX(
                nodeId, control.prescribedInnerBoundaryDx);
        }
        else if (control.prescribeInnerBoundaryY) {
            tireCase.assembly->addPrescribedDisplacementY(
                nodeId, control.prescribedInnerBoundaryDy);
        }
    }

    if (control.addInnerBoundaryAnchor) {
        if (setup.anchorNodeId < 0) {
            throw std::runtime_error(
                "Failed to select anchor node for finite-strain tire-contact case");
        }
        tireCase.assembly->addFixedNode(
            setup.anchorNodeId, control.anchorFixX, control.anchorFixY);
    }

    tireCase.mesh = setup.mesh;
    tireCase.rigidPlane = setup.rigidPlane;
    tireCase.innerBoundaryNodeIds = setup.innerBoundaryNodeIds;
    tireCase.anchorNodeId = setup.anchorNodeId;
    return tireCase;
}

} // namespace HyperelasticScenarioSupport
