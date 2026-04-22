#include "RigidPlaneAugmentedLagrangianContactSolver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "FiniteStrainMaterial.h"
#include "NeoHookeanMaterial.h"
#include "assembly.h"
#include "node.h"
#include "planeisometric/FiniteStrainQ4Element.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kGaussPoint = 0.577350269189626;
const std::vector<double> kGaussPoints = {-kGaussPoint, kGaussPoint};
const std::vector<double> kGaussWeights = {1.0, 1.0};
constexpr int kGaussPointCountPerFacet = 2;

std::pair<double, double> mapSurfaceCoordinatesLocal(int surfaceIndex, double surfaceParameter) {
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

struct ContactElementView {
    int elementId = -1;
    int materialId = -1;
    double thickness = 0.0;
    double characteristicYoungsModulus = 0.0;
    std::vector<int> nodeIds;
    std::shared_ptr<PlaneIsoparametricElement> linearElement;
    std::shared_ptr<FiniteStrainQ4Element> finiteStrainElement;
};

ContactElementView resolveContactElementView(
    const std::shared_ptr<Assembly>& assembly,
    int elementId) {
    ContactElementView elementView;
    elementView.elementId = elementId;

    auto linearElement = assembly->getElement(elementId);
    auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(linearElement);
    if (planeElement) {
        const auto material = assembly->getMaterial(planeElement->getMaterialId());
        if (!material) {
            throw std::invalid_argument(
                "Linear material not found for contact element " + std::to_string(elementId));
        }

        elementView.materialId = planeElement->getMaterialId();
        elementView.thickness = material->getThickness();
        elementView.characteristicYoungsModulus = material->getYoungsModulus();
        elementView.nodeIds = planeElement->getNodeIds();
        elementView.linearElement = std::move(planeElement);
        return elementView;
    }

    auto finiteStrainElement = assembly->getFiniteStrainElement(elementId);
    if (finiteStrainElement) {
        const auto material =
            assembly->getFiniteStrainMaterial(finiteStrainElement->getMaterialId());
        if (!material) {
            throw std::invalid_argument(
                "Finite-strain material not found for contact element " +
                std::to_string(elementId));
        }
        const auto neoHookeanMaterial =
            std::dynamic_pointer_cast<NeoHookeanMaterial>(material);
        if (!neoHookeanMaterial) {
            throw std::invalid_argument(
                "Finite-strain augmented Lagrangian contact currently expects "
                "Neo-Hookean material for characteristic stiffness scaling on element " +
                std::to_string(elementId));
        }

        elementView.materialId = finiteStrainElement->getMaterialId();
        elementView.thickness = material->getThickness();
        elementView.characteristicYoungsModulus =
            neoHookeanMaterial->getEquivalentYoungsModulus();
        elementView.nodeIds = finiteStrainElement->getNodeIds();
        elementView.finiteStrainElement = std::move(finiteStrainElement);
        return elementView;
    }

    throw std::invalid_argument(
        "Contact element not found in either linear or finite-strain assembly branch: " +
        std::to_string(elementId));
}

std::vector<std::shared_ptr<Node>> getElementNodes(
    const std::shared_ptr<Assembly>& assembly,
    const std::vector<int>& nodeIds) {
    std::vector<std::shared_ptr<Node>> elementNodes;
    elementNodes.reserve(nodeIds.size());
    for (int nodeId : nodeIds) {
        auto node = assembly->getNode(nodeId);
        if (!node) {
            throw std::invalid_argument("Node not found: " + std::to_string(nodeId));
        }
        elementNodes.push_back(node);
    }
    return elementNodes;
}

std::vector<int> getElementFullDofIndices(
    const std::shared_ptr<Assembly>& assembly,
    const std::vector<int>& nodeIds) {
    std::vector<int> dofIndices;
    dofIndices.reserve(nodeIds.size() * 2);
    for (int nodeId : nodeIds) {
        dofIndices.push_back(assembly->getGlobalDofIndex(nodeId, 0));
        dofIndices.push_back(assembly->getGlobalDofIndex(nodeId, 1));
    }
    return dofIndices;
}

Eigen::Vector4d evaluateScalarShapeFunctions(
    const ContactElementView& elementView,
    double xi,
    double eta) {
    if (elementView.linearElement) {
        Eigen::MatrixXd shapeMatrix = elementView.linearElement->shapeFunctions(xi, eta);
        Eigen::Vector4d N = Eigen::Vector4d::Zero();
        for (int i = 0; i < 4; ++i) {
            N(i) = shapeMatrix(0, 2 * i);
        }
        return N;
    }

    if (elementView.finiteStrainElement) {
        return elementView.finiteStrainElement->shapeFunctionsLocal(xi, eta);
    }

    throw std::runtime_error("Unsupported contact element view without element backend");
}

double computeSurfaceJacobian(int surfaceIndex,
    const ContactElementView& elementView,
    const std::vector<std::shared_ptr<Node>>& elementNodes,
    double xi,
    double eta) {
    Eigen::Matrix2d jacobian = Eigen::Matrix2d::Zero();
    if (elementView.linearElement) {
        jacobian = elementView.linearElement->jacobian(xi, eta, elementNodes);
    }
    else if (elementView.finiteStrainElement) {
        jacobian = elementView.finiteStrainElement->referenceJacobian(xi, eta, elementNodes);
    }
    else {
        throw std::runtime_error("Unsupported contact element view without element backend");
    }

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

double computeReferenceFacetLength(int surfaceIndex,
    const ContactElementView& elementView,
    const std::vector<std::shared_ptr<Node>>& elementNodes) {
    double facetLength = 0.0;
    for (size_t gpIndex = 0; gpIndex < kGaussPoints.size(); ++gpIndex) {
        const auto [xi, eta] = mapSurfaceCoordinatesLocal(surfaceIndex, kGaussPoints[gpIndex]);
        facetLength += computeSurfaceJacobian(
            surfaceIndex, elementView, elementNodes, xi, eta) * kGaussWeights[gpIndex];
    }
    return facetLength;
}

double effectiveAugmentationParameter(const AugmentedLagrangianSettings& settings,
    const ContactElementView& elementView,
    double characteristicFacetLength) {
    if (!settings.useAutomaticAugmentationScaling) {
        return settings.augmentationParameter;
    }

    const double safeFacetLength = std::max(characteristicFacetLength, 1.0e-12);
    const double safeYoungsModulus = std::max(elementView.characteristicYoungsModulus, 1.0e-12);
    return settings.automaticScalingFactor * safeYoungsModulus / safeFacetLength;
}

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
    if (settings_.absoluteMultiplierTolerance < 0.0) {
        throw std::invalid_argument("Absolute multiplier tolerance must be non-negative");
    }
    if (settings_.penetrationTolerance < 0.0) {
        throw std::invalid_argument("Penetration tolerance must be non-negative");
    }
    if (settings_.relativePenetrationToleranceFactor < 0.0) {
        throw std::invalid_argument(
            "Relative penetration tolerance factor must be non-negative");
    }
    if (!(settings_.multiplierRelaxation > 0.0) || settings_.multiplierRelaxation > 1.0) {
        throw std::invalid_argument("Multiplier relaxation must lie in (0, 1]");
    }
    if (settings_.useAutomaticAugmentationScaling &&
        !(settings_.automaticScalingFactor > 0.0)) {
        throw std::invalid_argument("Automatic augmentation scaling factor must be positive");
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
    if (settings.absoluteMultiplierTolerance < 0.0) {
        throw std::invalid_argument("Absolute multiplier tolerance must be non-negative");
    }
    if (settings.penetrationTolerance < 0.0) {
        throw std::invalid_argument("Penetration tolerance must be non-negative");
    }
    if (settings.relativePenetrationToleranceFactor < 0.0) {
        throw std::invalid_argument(
            "Relative penetration tolerance factor must be non-negative");
    }
    if (!(settings.multiplierRelaxation > 0.0) || settings.multiplierRelaxation > 1.0) {
        throw std::invalid_argument("Multiplier relaxation must lie in (0, 1]");
    }
    if (settings.useAutomaticAugmentationScaling &&
        !(settings.automaticScalingFactor > 0.0)) {
        throw std::invalid_argument("Automatic augmentation scaling factor must be positive");
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
        const ContactElementView elementView =
            resolveContactElementView(assembly_, facet.elementId);
        std::vector<std::shared_ptr<Node>> elementNodes =
            getElementNodes(assembly_, elementView.nodeIds);
        std::vector<int> dofIndices =
            getElementFullDofIndices(assembly_, elementView.nodeIds);
        const double characteristicFacetLength =
            computeReferenceFacetLength(facet.surfaceIndex, elementView, elementNodes);
        const double augmentationParameter = effectiveAugmentationParameter(
            settings_, elementView, characteristicFacetLength);

        ContactFacetResult facetResult;
        facetResult.facetId = static_cast<int>(facetIndex);
        facetResult.elementId = facet.elementId;
        facetResult.surfaceIndex = facet.surfaceIndex;
        facetResult.normal = plane_.normal;
        facetResult.thickness = elementView.thickness;
        bool facetActive = false;

        auto [midXi, midEta] = mapSurfaceCoordinates(facet.surfaceIndex, 0.0);
        Eigen::Vector4d midpointShapeFunctions =
            evaluateScalarShapeFunctions(elementView, midXi, midEta);
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
            Eigen::Vector4d N = evaluateScalarShapeFunctions(elementView, xi, eta);

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

            const double thickness = std::max(elementView.thickness, 1.0e-15);
            const double ds = computeSurfaceJacobian(
                facet.surfaceIndex, elementView, elementNodes, xi, eta) *
                kGaussWeights[gpIndex] *
                thickness;
            const double lengthContribution = ds / thickness;
            const double gap = plane_.signedDistance(xGp + uGp);
            const double penetration = std::max(0.0, -gap);
            const double augmentedPressure = std::max(
                0.0,
                storedState.lambdaN - augmentationParameter * gap);

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
                    Eigen::Matrix2d kij = augmentationParameter * N(i) * N(j) *
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
    double characteristicFacetLengthSum = 0.0;
    int characteristicFacetCount = 0;

    for (size_t facetIndex = 0; facetIndex < facets_.size(); ++facetIndex) {
        const auto& facet = facets_[facetIndex];
        const ContactElementView elementView =
            resolveContactElementView(assembly_, facet.elementId);
        std::vector<std::shared_ptr<Node>> elementNodes =
            getElementNodes(assembly_, elementView.nodeIds);
        std::vector<int> dofIndices =
            getElementFullDofIndices(assembly_, elementView.nodeIds);
        const double characteristicFacetLength =
            computeReferenceFacetLength(facet.surfaceIndex, elementView, elementNodes);
        characteristicFacetLengthSum += characteristicFacetLength;
        characteristicFacetCount += 1;
        const double augmentationParameter = effectiveAugmentationParameter(
            settings_, elementView, characteristicFacetLength);

        for (int gpIndex = 0; gpIndex < kGaussPointCountPerFacet; ++gpIndex) {
            const int stateIndex =
                static_cast<int>(facetIndex) * kGaussPointCountPerFacet + gpIndex;
            ContactGaussPointState& gaussPointState = gaussPointStates_[stateIndex];

            auto [xi, eta] =
                mapSurfaceCoordinates(facet.surfaceIndex, kGaussPoints[gpIndex]);
            Eigen::Vector4d N = evaluateScalarShapeFunctions(elementView, xi, eta);

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
            const double trialLambda = std::max(
                0.0,
                previousLambda - augmentationParameter * gap);
            const double updatedLambda =
                previousLambda +
                settings_.multiplierRelaxation * (trialLambda - previousLambda);

            gaussPointState.gapN = gap;
            gaussPointState.penetrationN = penetration;
            gaussPointState.lambdaN = updatedLambda;
            gaussPointState.pressureN = updatedLambda;
            gaussPointState.active = updatedLambda > 0.0;

            const double deltaLambda = updatedLambda - previousLambda;
            deltaLambdaSquaredNorm += deltaLambda * deltaLambda;
            lambdaSquaredNorm += updatedLambda * updatedLambda;
            updateInfo.maxStateUpdateMagnitude =
                std::max(updateInfo.maxStateUpdateMagnitude, std::abs(deltaLambda));
            updateInfo.maxPenetration =
                std::max(updateInfo.maxPenetration, penetration);

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
    updateInfo.relativeMaxStateUpdate =
        updateInfo.maxStateUpdateMagnitude / (updateInfo.maxNormalMultiplier + 1.0e-15);
    updateInfo.meanNormalMultiplier =
        updateInfo.activeGaussPointCount > 0
        ? lambdaSum / static_cast<double>(updateInfo.activeGaussPointCount)
        : 0.0;
    const double meanCharacteristicFacetLength =
        characteristicFacetCount > 0
        ? characteristicFacetLengthSum / static_cast<double>(characteristicFacetCount)
        : 0.0;
    const double effectivePenetrationTolerance = std::max(
        settings_.penetrationTolerance,
        settings_.relativePenetrationToleranceFactor * meanCharacteristicFacetLength);
    updateInfo.converged =
        updateInfo.maxPenetration <= effectivePenetrationTolerance &&
        (updateInfo.relativeStateUpdateNorm <= settings_.multiplierTolerance ||
            updateInfo.relativeMaxStateUpdate <= settings_.multiplierTolerance ||
            updateInfo.stateUpdateNorm <= settings_.absoluteMultiplierTolerance);
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
