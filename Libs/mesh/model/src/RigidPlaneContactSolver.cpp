#include "RigidPlaneContactSolver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "assembly.h"
#include "material.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

RigidPlaneContactSolver::RigidPlaneContactSolver(std::shared_ptr<Assembly> assembly,
    const RigidPlane2D& plane,
    double penaltyParameter)
    : assembly_(std::move(assembly)), plane_(plane), penaltyParameter_(penaltyParameter) {
    if (!assembly_) {
        throw std::invalid_argument("Assembly cannot be null for contact solver");
    }
    if (plane_.normal.norm() == 0.0) {
        throw std::invalid_argument("Rigid plane normal cannot be zero");
    }
    plane_.normal.normalize();
}

void RigidPlaneContactSolver::setContactFacets(std::vector<ContactFacet> facets) {
    facets_ = std::move(facets);
}

void RigidPlaneContactSolver::setPlane(const RigidPlane2D& plane) {
    plane_ = plane;
    if (plane_.normal.norm() == 0.0) {
        throw std::invalid_argument("Rigid plane normal cannot be zero");
    }
    plane_.normal.normalize();
}

void RigidPlaneContactSolver::setPenaltyParameter(double penaltyParameter) {
    penaltyParameter_ = penaltyParameter;
}

void RigidPlaneContactSolver::assembleContact(const Eigen::VectorXd& fullDisplacements,
    Eigen::SparseMatrix<double>& contactK,
    Eigen::VectorXd& contactF,
    ContactState& state) const {
    const int totalDof = assembly_->getTotalDofCount();
    contactF = Eigen::VectorXd::Zero(totalDof);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(facets_.size() * 32);

    static constexpr double gaussPoint = 0.577350269189626;
    static const std::vector<double> gaussPoints = {-gaussPoint, gaussPoint};
    static const std::vector<double> gaussWeights = {1.0, 1.0};

    state = {};

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
        bool facetActive = false;

        for (size_t gpIndex = 0; gpIndex < gaussPoints.size(); ++gpIndex) {
            auto [xi, eta] = mapSurfaceCoordinates(facet.surfaceIndex, gaussPoints[gpIndex]);
            Eigen::Vector4d N = evaluateScalarShapeFunctions(*planeElement, xi, eta);

            Eigen::Vector2d xGp = Eigen::Vector2d::Zero();
            Eigen::Vector2d uGp = Eigen::Vector2d::Zero();
            for (int localNode = 0; localNode < 4; ++localNode) {
                xGp += N(localNode) * elementNodes[localNode]->getCoordinates();
                int dofX = dofIndices[2 * localNode];
                int dofY = dofIndices[2 * localNode + 1];
                if (dofX >= 0 && dofY >= 0 &&
                    dofX < fullDisplacements.size() && dofY < fullDisplacements.size()) {
                    uGp.x() += N(localNode) * fullDisplacements(dofX);
                    uGp.y() += N(localNode) * fullDisplacements(dofY);
                }
            }

            double gap = plane_.signedDistance(xGp + uGp);
            if (gap >= 0.0) {
                continue;
            }

            facetActive = true;
            state.maxPenetration = std::max(state.maxPenetration, -gap);

            double ds = computeSurfaceJacobian(
                facet.surfaceIndex, *planeElement, elementNodes, xi, eta) *
                gaussWeights[gpIndex] *
                material->getThickness();

            Eigen::Vector2d traction = (-penaltyParameter_ * gap) * plane_.normal;

            for (int i = 0; i < 4; ++i) {
                contactF(dofIndices[2 * i]) += N(i) * traction.x() * ds;
                contactF(dofIndices[2 * i + 1]) += N(i) * traction.y() * ds;

                for (int j = 0; j < 4; ++j) {
                    Eigen::Matrix2d kij = penaltyParameter_ * N(i) * N(j) *
                        (plane_.normal * plane_.normal.transpose()) * ds;
                    triplets.emplace_back(dofIndices[2 * i], dofIndices[2 * j], kij(0, 0));
                    triplets.emplace_back(dofIndices[2 * i], dofIndices[2 * j + 1], kij(0, 1));
                    triplets.emplace_back(dofIndices[2 * i + 1], dofIndices[2 * j], kij(1, 0));
                    triplets.emplace_back(dofIndices[2 * i + 1], dofIndices[2 * j + 1], kij(1, 1));
                }
            }

            state.contactForceNorm += traction.norm() * ds;
        }

        if (facetActive) {
            state.activeFacetIds.push_back(static_cast<int>(facetIndex));
        }
    }

    contactK.resize(totalDof, totalDof);
    contactK.setFromTriplets(triplets.begin(), triplets.end());
    contactK.makeCompressed();
}

std::pair<double, double> RigidPlaneContactSolver::mapSurfaceCoordinates(int surfaceIndex,
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

double RigidPlaneContactSolver::computeSurfaceJacobian(int surfaceIndex,
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

std::vector<std::shared_ptr<Node>> RigidPlaneContactSolver::getElementNodes(int elementId) const {
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

Eigen::Vector4d RigidPlaneContactSolver::evaluateScalarShapeFunctions(
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
