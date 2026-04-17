#include "RigidPlaneAugmentedLagrangianContactSolver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "assembly.h"
#include "material.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kGaussPoint = 0.577350269189626;
const std::vector<double> kGaussPoints = {-kGaussPoint, kGaussPoint};
const std::vector<double> kGaussWeights = {1.0, 1.0};
constexpr int kGaussPointCountPerFacet = 2;

} // namespace

RigidPlaneAugmentedLagrangianContactSolver::RigidPlaneAugmentedLagrangianContactSolver(
    std::shared_ptr<Assembly> assembly,
    const RigidPlane2D& plane,
    const AugmentedLagrangianSettings& settings)
    : assembly_(std::move(assembly)), plane_(plane), settings_(settings) {
    if (!assembly_) {
        throw std::invalid_argument("Assembly cannot be null for contact solver");
    }
    if (plane_.normal.norm() == 0.0) {
        throw std::invalid_argument("Rigid plane normal cannot be zero");
    }
    if (settings_.augmentationParameter <= 0.0) {
        throw std::invalid_argument("Augmented Lagrangian parameter must be positive");
    }
    if (settings_.multiplierTolerance < 0.0) {
        throw std::invalid_argument("Multiplier tolerance must be non-negative");
    }

    plane_.normal.normalize();
}

void RigidPlaneAugmentedLagrangianContactSolver::setContactFacets(std::vector<ContactFacet> facets) {
    facets_ = std::move(facets);
    initializeGaussPointStates();
}

void RigidPlaneAugmentedLagrangianContactSolver::setPlane(const RigidPlane2D& plane) {
    plane_ = plane;
    if (plane_.normal.norm() == 0.0) {
        throw std::invalid_argument("Rigid plane normal cannot be zero");
    }
    plane_.normal.normalize();
}

void RigidPlaneAugmentedLagrangianContactSolver::setSettings(
    const AugmentedLagrangianSettings& settings) {
    if (settings.augmentationParameter <= 0.0) {
        throw std::invalid_argument("Augmented Lagrangian parameter must be positive");
    }
    if (settings.multiplierTolerance < 0.0) {
        throw std::invalid_argument("Multiplier tolerance must be non-negative");
    }

    settings_ = settings;
}

void RigidPlaneAugmentedLagrangianContactSolver::assembleContact(
    const Eigen::VectorXd& fullDisplacements,
    Eigen::SparseMatrix<double>& contactK,
    Eigen::VectorXd& contactF,
    ContactState& state) const {
    if (gaussPointStates_.size() != facets_.size() * kGaussPointCountPerFacet) {
        throw std::runtime_error("Augmented Lagrangian GP state is not initialized");
    }

    const int totalDof = assembly_->getTotalDofCount();
    contactF = Eigen::VectorXd::Zero(totalDof);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(facets_.size() * 128);

    state = {};
    state.method = getMethod();
    state.activeFacetIds.reserve(facets_.size());
    state.facetResults.reserve(facets_.size());
    state.gaussPointStates.reserve(facets_.size() * kGaussPointCountPerFacet);

    for (size_t facetIndex = 0; facetIndex < facets_.size(); ++facetIndex) {
        const auto& facet = facets_[facetIndex];
        auto element = assembly_->getElement(facet.elementId);
        auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
        if (!planeElement) {
            continue;
        }

        auto material = assembly_->getMaterial(planeElement->getMaterialId());
        if (!material) {
            continue;
        }

        std::vector<std::shared_ptr<Node>> elementNodes = getElementNodes(facet.elementId);
        std::vector<int> dofIndices = assembly_->getElementFullDofIndices(facet.elementId);

        ContactFacetResult facetResult;
        facetResult.facetId = static_cast<int>(facetIndex);
        facetResult.elementId = facet.elementId;
        facetResult.surfaceIndex = facet.surfaceIndex;
        facetResult.normal = plane_.normal;
        facetResult.thickness = material->getThickness();
        bool facetActive = false;

        auto [midXi, midEta] = mapSurfaceCoordinates(facet.surfaceIndex, 0.0);
        Eigen::Vector4d midpointShapeFunctions =
            evaluateScalarShapeFunctions(*planeElement, midXi, midEta);
        for (int localNode = 0; localNode < 4; ++localNode) {
            facetResult.referenceMidpoint +=
                midpointShapeFunctions(localNode) * elementNodes[localNode]->getCoordinates();
            const int dofX = dofIndices[2 * localNode];
            const int dofY = dofIndices[2 * localNode + 1];
            if (dofX >= 0 && dofY >= 0 &&
                dofX < fullDisplacements.size() && dofY < fullDisplacements.size()) {
                facetResult.deformedMidpoint.x() +=
                    midpointShapeFunctions(localNode) * fullDisplacements(dofX);
                facetResult.deformedMidpoint.y() +=
                    midpointShapeFunctions(localNode) * fullDisplacements(dofY);
            }
        }
        facetResult.deformedMidpoint += facetResult.referenceMidpoint;

        for (int gpIndex = 0; gpIndex < kGaussPointCountPerFacet; ++gpIndex) {
            const int stateIndex =
                static_cast<int>(facetIndex) * kGaussPointCountPerFacet + gpIndex;
            const ContactGaussPointState& storedState = gaussPointStates_[stateIndex];

            auto [xi, eta] =
                mapSurfaceCoordinates(facet.surfaceIndex, kGaussPoints[gpIndex]);
            Eigen::Vector4d N = evaluateScalarShapeFunctions(*planeElement, xi, eta);

            Eigen::Vector2d xGp = Eigen::Vector2d::Zero();
            Eigen::Vector2d uGp = Eigen::Vector2d::Zero();
            for (int localNode = 0; localNode < 4; ++localNode) {
                xGp += N(localNode) * elementNodes[localNode]->getCoordinates();
                const int dofX = dofIndices[2 * localNode];
                const int dofY = dofIndices[2 * localNode + 1];
                if (dofX >= 0 && dofY >= 0 &&
                    dofX < fullDisplacements.size() && dofY < fullDisplacements.size()) {
                    uGp.x() += N(localNode) * fullDisplacements(dofX);
                    uGp.y() += N(localNode) * fullDisplacements(dofY);
                }
            }

            const double thickness = std::max(material->getThickness(), 1.0e-15);
            const double ds = computeSurfaceJacobian(
                facet.surfaceIndex, *planeElement, elementNodes, xi, eta) *
                kGaussWeights[gpIndex] *
                thickness;
            const double lengthContribution = ds / thickness;
            const double gap = plane_.signedDistance(xGp + uGp);
            const double penetration = std::max(0.0, -gap);
            const double augmentedPressure = std::max(
                0.0,
                storedState.lambdaN - settings_.augmentationParameter * gap);

            ContactGaussPointState gaussPointState = storedState;
            gaussPointState.gapN = gap;
            gaussPointState.penetrationN = penetration;
            gaussPointState.pressureN = augmentedPressure;
            gaussPointState.active = augmentedPressure > 0.0;
            state.gaussPointStates.push_back(gaussPointState);

            facetResult.integratedArea += ds;
            facetResult.facetLength += lengthContribution;
            facetResult.averageGap += gap * ds;
            facetResult.averagePenetration += penetration * ds;
            facetResult.maximumPenetration = std::max(facetResult.maximumPenetration, penetration);

            if (augmentedPressure <= 0.0) {
                continue;
            }

            facetActive = true;
            facetResult.active = true;
            facetResult.activeGaussPointCount += 1;
            state.activeGaussPointCount += 1;
            facetResult.activeArea += ds;
            facetResult.activeLength += lengthContribution;
            state.maxPenetration = std::max(state.maxPenetration, penetration);

            const Eigen::Vector2d traction = augmentedPressure * plane_.normal;
            facetResult.integratedNormalForce += augmentedPressure * ds;

            for (int i = 0; i < 4; ++i) {
                contactF(dofIndices[2 * i]) += N(i) * traction.x() * ds;
                contactF(dofIndices[2 * i + 1]) += N(i) * traction.y() * ds;

                for (int j = 0; j < 4; ++j) {
                    Eigen::Matrix2d kij = settings_.augmentationParameter * N(i) * N(j) *
                        (plane_.normal * plane_.normal.transpose()) * ds;
                    triplets.emplace_back(dofIndices[2 * i], dofIndices[2 * j], kij(0, 0));
                    triplets.emplace_back(dofIndices[2 * i], dofIndices[2 * j + 1], kij(0, 1));
                    triplets.emplace_back(dofIndices[2 * i + 1], dofIndices[2 * j], kij(1, 0));
                    triplets.emplace_back(dofIndices[2 * i + 1], dofIndices[2 * j + 1], kij(1, 1));
                }
            }

            state.contactForceNorm += augmentedPressure * ds;
        }

        if (facetActive) {
            state.activeFacetIds.push_back(static_cast<int>(facetIndex));
        }

        if (facetResult.integratedArea > 0.0) {
            facetResult.averageGap /= facetResult.integratedArea;
            facetResult.averagePenetration /= facetResult.integratedArea;
            facetResult.averagePressure =
                facetResult.integratedNormalForce / facetResult.integratedArea;
        }

        state.facetResults.push_back(facetResult);
    }

    contactK.resize(totalDof, totalDof);
    contactK.setFromTriplets(triplets.begin(), triplets.end());
    contactK.makeCompressed();
}

ContactSolverUpdateInfo RigidPlaneAugmentedLagrangianContactSolver::updateState(
    const Eigen::VectorXd& fullDisplacements) {
    if (gaussPointStates_.size() != facets_.size() * kGaussPointCountPerFacet) {
        initializeGaussPointStates();
    }

    ContactSolverUpdateInfo updateInfo;
    updateInfo.converged = true;

    double deltaLambdaSquaredNorm = 0.0;
    double lambdaSquaredNorm = 0.0;
    double lambdaSum = 0.0;

    for (size_t facetIndex = 0; facetIndex < facets_.size(); ++facetIndex) {
        const auto& facet = facets_[facetIndex];
        auto element = assembly_->getElement(facet.elementId);
        auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
        if (!planeElement) {
            continue;
        }

        std::vector<std::shared_ptr<Node>> elementNodes = getElementNodes(facet.elementId);
        std::vector<int> dofIndices = assembly_->getElementFullDofIndices(facet.elementId);

        for (int gpIndex = 0; gpIndex < kGaussPointCountPerFacet; ++gpIndex) {
            const int stateIndex =
                static_cast<int>(facetIndex) * kGaussPointCountPerFacet + gpIndex;
            ContactGaussPointState& gaussPointState = gaussPointStates_[stateIndex];

            auto [xi, eta] =
                mapSurfaceCoordinates(facet.surfaceIndex, kGaussPoints[gpIndex]);
            Eigen::Vector4d N = evaluateScalarShapeFunctions(*planeElement, xi, eta);

            Eigen::Vector2d xGp = Eigen::Vector2d::Zero();
            Eigen::Vector2d uGp = Eigen::Vector2d::Zero();
            for (int localNode = 0; localNode < 4; ++localNode) {
                xGp += N(localNode) * elementNodes[localNode]->getCoordinates();
                const int dofX = dofIndices[2 * localNode];
                const int dofY = dofIndices[2 * localNode + 1];
                if (dofX >= 0 && dofY >= 0 &&
                    dofX < fullDisplacements.size() && dofY < fullDisplacements.size()) {
                    uGp.x() += N(localNode) * fullDisplacements(dofX);
                    uGp.y() += N(localNode) * fullDisplacements(dofY);
                }
            }

            const double gap = plane_.signedDistance(xGp + uGp);
            const double penetration = std::max(0.0, -gap);
            const double previousLambda = gaussPointState.lambdaN;
            const double updatedLambda = std::max(
                0.0,
                previousLambda - settings_.augmentationParameter * gap);

            gaussPointState.gapN = gap;
            gaussPointState.penetrationN = penetration;
            gaussPointState.lambdaN = updatedLambda;
            gaussPointState.pressureN = updatedLambda;
            gaussPointState.active = updatedLambda > 0.0;

            const double deltaLambda = updatedLambda - previousLambda;
            deltaLambdaSquaredNorm += deltaLambda * deltaLambda;
            lambdaSquaredNorm += updatedLambda * updatedLambda;

            if (gaussPointState.active) {
                updateInfo.activeGaussPointCount += 1;
                updateInfo.maxNormalMultiplier =
                    std::max(updateInfo.maxNormalMultiplier, updatedLambda);
                lambdaSum += updatedLambda;
            }
        }
    }

    updateInfo.stateUpdateNorm = std::sqrt(deltaLambdaSquaredNorm);
    updateInfo.relativeStateUpdateNorm =
        updateInfo.stateUpdateNorm / (std::sqrt(lambdaSquaredNorm) + 1.0e-15);
    updateInfo.meanNormalMultiplier =
        updateInfo.activeGaussPointCount > 0
        ? lambdaSum / static_cast<double>(updateInfo.activeGaussPointCount)
        : 0.0;
    updateInfo.converged =
        updateInfo.relativeStateUpdateNorm <= settings_.multiplierTolerance;
    return updateInfo;
}

void RigidPlaneAugmentedLagrangianContactSolver::resetState() {
    initializeGaussPointStates();
}

std::pair<double, double> RigidPlaneAugmentedLagrangianContactSolver::mapSurfaceCoordinates(
    int surfaceIndex,
    double surfaceParameter) {
    switch (surfaceIndex) {
    case 0:
        return {surfaceParameter, -1.0};
    case 1:
        return {1.0, surfaceParameter};
    case 2:
        return {surfaceParameter, 1.0};
    case 3:
        return {-1.0, surfaceParameter};
    default:
        throw std::invalid_argument("Surface index must be in [0, 3]");
    }
}

double RigidPlaneAugmentedLagrangianContactSolver::computeSurfaceJacobian(
    int surfaceIndex,
    const PlaneIsoparametricElement& element,
    const std::vector<std::shared_ptr<Node>>& elementNodes,
    double xi,
    double eta) {
    Eigen::Matrix2d jacobian = element.jacobian(xi, eta, elementNodes);

    switch (surfaceIndex) {
    case 0:
    case 2:
        return std::sqrt(jacobian(0, 0) * jacobian(0, 0) + jacobian(0, 1) * jacobian(0, 1));
    case 1:
    case 3:
        return std::sqrt(jacobian(1, 0) * jacobian(1, 0) + jacobian(1, 1) * jacobian(1, 1));
    default:
        throw std::invalid_argument("Surface index must be in [0, 3]");
    }
}

void RigidPlaneAugmentedLagrangianContactSolver::initializeGaussPointStates() {
    gaussPointStates_.clear();
    gaussPointStates_.reserve(facets_.size() * kGaussPointCountPerFacet);

    for (size_t facetIndex = 0; facetIndex < facets_.size(); ++facetIndex) {
        const auto& facet = facets_[facetIndex];
        for (int gpIndex = 0; gpIndex < kGaussPointCountPerFacet; ++gpIndex) {
            auto [xi, eta] =
                mapSurfaceCoordinates(facet.surfaceIndex, kGaussPoints[gpIndex]);

            ContactGaussPointState gaussPointState;
            gaussPointState.facetId = static_cast<int>(facetIndex);
            gaussPointState.gaussPointIndex = gpIndex;
            gaussPointState.surfaceParameter = kGaussPoints[gpIndex];
            gaussPointState.xi = xi;
            gaussPointState.eta = eta;
            gaussPointStates_.push_back(gaussPointState);
        }
    }
}

std::vector<std::shared_ptr<Node>> RigidPlaneAugmentedLagrangianContactSolver::getElementNodes(
    int elementId) const {
    auto element = assembly_->getElement(elementId);
    if (!element) {
        throw std::invalid_argument("Element not found: " + std::to_string(elementId));
    }

    std::vector<std::shared_ptr<Node>> elementNodes;
    elementNodes.reserve(element->getNodeIds().size());
    for (int nodeId : element->getNodeIds()) {
        auto node = assembly_->getNode(nodeId);
        if (!node) {
            throw std::invalid_argument("Node not found: " + std::to_string(nodeId));
        }
        elementNodes.push_back(node);
    }

    return elementNodes;
}

Eigen::Vector4d RigidPlaneAugmentedLagrangianContactSolver::evaluateScalarShapeFunctions(
    const PlaneIsoparametricElement& element,
    double xi,
    double eta) const {
    Eigen::MatrixXd shapeMatrix = element.shapeFunctions(xi, eta);
    Eigen::Vector4d N = Eigen::Vector4d::Zero();
    for (int i = 0; i < 4; ++i) {
        N(i) = shapeMatrix(0, 2 * i);
    }
    return N;
}
