#include "loadFunction.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "assembly.h"
#include "planeisometric/FiniteStrainQ4Element.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

std::vector<std::shared_ptr<Node>> loadElementNodes(
    const std::vector<int>& nodeIds,
    const std::shared_ptr<const Assembly>& assembly,
    int elementId,
    const char* elementKind) {
    std::vector<std::shared_ptr<Node>> nodes;
    nodes.reserve(nodeIds.size());
    for (int nodeId : nodeIds) {
        auto node = assembly->getNode(nodeId);
        if (node) {
            nodes.push_back(node);
        }
        else {
            std::cerr << "Assembly validation failed: " << elementKind << " element "
                      << elementId << " references non-existent node " << nodeId
                      << std::endl;
            break;
        }
    }
    return nodes;
}

Eigen::Vector2d legacyRigidPlaneNormal() {
    // Preserve the semantics used in the historical no-contact surrogate:
    // pressure acts along the rigid-plane normal, not the local outer-surface normal.
    return Eigen::Vector2d(0.0, 1.0);
}

double surfaceJacobianFromReferenceJacobian(
    const Eigen::Matrix2d& jacobian,
    int surfaceIndex) {
    switch (surfaceIndex) {
    case 0:
    case 2:
        return std::sqrt(jacobian(0, 0) * jacobian(0, 0) +
                         jacobian(0, 1) * jacobian(0, 1));
    case 1:
    case 3:
        return std::sqrt(jacobian(1, 0) * jacobian(1, 0) +
                         jacobian(1, 1) * jacobian(1, 1));
    default:
        throw std::invalid_argument("Invalid surface index");
    }
}

void surfaceLocalCoordinates(int surfaceIndex, double param, double& xi, double& eta) {
    switch (surfaceIndex) {
    case 0: xi = param; eta = -1.0; break;
    case 1: xi = 1.0; eta = param; break;
    case 2: xi = param; eta = 1.0; break;
    case 3: xi = -1.0; eta = param; break;
    default:
        throw std::invalid_argument("Invalid surface index");
    }
}

} // namespace

Eigen::VectorXd LoadFunction::applyToElementSurface(
    int elementId,
    int surfaceIndex,
    const std::shared_ptr<const Assembly>& assembly) {
    const double gp = 0.577350269189626;
    const std::vector<double> gaussPoints = {-gp, gp};
    const std::vector<double> weights = {1.0, 1.0};

    if (auto element = assembly->getElement(elementId)) {
        auto material = assembly->getMaterial(element->getMaterialId());
        const auto nodes =
            loadElementNodes(element->getNodeIds(), assembly, element->getId(), "linear");
        if (!material || nodes.size() != element->getNodeIds().size()) {
            throw std::invalid_argument("Invalid linear element or material");
        }

        auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
        if (!planeElement) {
            throw std::runtime_error("Element is not a PlaneIsoparametricElement");
        }

        Eigen::VectorXd surfaceForces = Eigen::VectorXd::Zero(element->getDofCount());
        const auto nodeCoords = planeElement->getNodalCoordinates(nodes);

        for (int gaussIndex = 0; gaussIndex < static_cast<int>(gaussPoints.size()); ++gaussIndex) {
            double xi = 0.0;
            double eta = 0.0;
            const double param = gaussPoints[static_cast<size_t>(gaussIndex)];
            const double weight = weights[static_cast<size_t>(gaussIndex)];
            surfaceLocalCoordinates(surfaceIndex, param, xi, eta);

            const Eigen::MatrixXd N = planeElement->shapeFunctions(xi, eta);
            Eigen::Vector2d point = Eigen::Vector2d::Zero();
            for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx) {
                const double shapeFunction = N(0, nodeIdx * 2);
                point += shapeFunction * nodeCoords.row(nodeIdx).transpose();
            }

            const Eigen::Vector2d load =
                distribution_(point.x(), point.y(), legacyRigidPlaneNormal());
            const Eigen::Matrix2d jacobian = planeElement->jacobian(xi, eta, nodes);
            const double surfaceJacobian =
                surfaceJacobianFromReferenceJacobian(jacobian, surfaceIndex);

            for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx) {
                const double shapeFunction = N(0, nodeIdx * 2);
                const int dofX = nodeIdx * 2;
                const int dofY = nodeIdx * 2 + 1;
                surfaceForces(dofX) += shapeFunction * load.x() * weight *
                    surfaceJacobian * material->getThickness();
                surfaceForces(dofY) += shapeFunction * load.y() * weight *
                    surfaceJacobian * material->getThickness();
            }
        }

        return surfaceForces;
    }

    if (auto element = assembly->getFiniteStrainElement(elementId)) {
        auto material = assembly->getFiniteStrainMaterial(element->getMaterialId());
        const auto nodes =
            loadElementNodes(element->getNodeIds(), assembly, element->getId(), "finite-strain");
        if (!material || nodes.size() != element->getNodeIds().size()) {
            throw std::invalid_argument("Invalid finite-strain element or material");
        }

        Eigen::VectorXd surfaceForces =
            Eigen::VectorXd::Zero(FiniteStrainQ4Element::kDofCount);

        for (int gaussIndex = 0; gaussIndex < static_cast<int>(gaussPoints.size()); ++gaussIndex) {
            double xi = 0.0;
            double eta = 0.0;
            const double param = gaussPoints[static_cast<size_t>(gaussIndex)];
            const double weight = weights[static_cast<size_t>(gaussIndex)];
            surfaceLocalCoordinates(surfaceIndex, param, xi, eta);

            const Eigen::Vector4d N = element->shapeFunctionsLocal(xi, eta);
            Eigen::Vector2d point = Eigen::Vector2d::Zero();
            for (int nodeIdx = 0; nodeIdx < FiniteStrainQ4Element::kNodeCount; ++nodeIdx) {
                point += N(nodeIdx) * nodes[static_cast<size_t>(nodeIdx)]->getCoordinates();
            }

            const Eigen::Vector2d load =
                distribution_(point.x(), point.y(), legacyRigidPlaneNormal());
            const Eigen::Matrix2d jacobian = element->referenceJacobian(xi, eta, nodes);
            const double surfaceJacobian =
                surfaceJacobianFromReferenceJacobian(jacobian, surfaceIndex);

            for (int nodeIdx = 0; nodeIdx < FiniteStrainQ4Element::kNodeCount; ++nodeIdx) {
                const int dofX = nodeIdx * 2;
                const int dofY = nodeIdx * 2 + 1;
                surfaceForces(dofX) += N(nodeIdx) * load.x() * weight *
                    surfaceJacobian * material->getThickness();
                surfaceForces(dofY) += N(nodeIdx) * load.y() * weight *
                    surfaceJacobian * material->getThickness();
            }
        }

        return surfaceForces;
    }

    throw std::invalid_argument(
        "Element with ID " + std::to_string(elementId) + " not found");
}

LoadFunction LoadFunction::constantLoad(double fx, double fy) {
    return LoadFunction([fx, fy](double x, double y, const Eigen::Vector2d& normal) {
        (void)x;
        (void)y;
        (void)normal;
        return Eigen::Vector2d(fx, fy);
    });
}

LoadFunction LoadFunction::linearPressure(double p0, double gradient) {
    return LoadFunction([p0, gradient](double x, double y, const Eigen::Vector2d& normal) {
        (void)y;
        const double pressure = p0 + gradient * x;
        return Eigen::Vector2d(pressure * normal);
    });
}

LoadFunction LoadFunction::hertzianPressure(double maxPressure, double contactWidth, double x0) {
    return LoadFunction(
        [maxPressure, contactWidth, x0](double x, double y, const Eigen::Vector2d& normal) {
            (void)y;
            const double dx = x - x0;
            if (std::abs(dx) > contactWidth) {
                return Eigen::Vector2d(0.0, 0.0);
            }

            const double relativeX = dx / contactWidth;
            const double pressure = maxPressure * std::sqrt(1.0 - relativeX * relativeX);
            return Eigen::Vector2d(pressure * normal);
        });
}

LoadFunction LoadFunction::sinusoidalLoad(double amplitude, double wavelength) {
    return LoadFunction(
        [amplitude, wavelength](double x, double y, const Eigen::Vector2d& normal) {
            (void)y;
            const double pressure = amplitude * std::sin(2 * M_PI * x / wavelength);
            return Eigen::Vector2d(pressure * normal);
        });
}

LoadFunction LoadFunction::parabolicPressure(double maxPressure,
    double contactHalfWidth,
    double x0) {
    return LoadFunction(
        [maxPressure, contactHalfWidth, x0](double x, double y, const Eigen::Vector2d& normal) {
            (void)y;
            const double dx = x - x0;
            if (std::abs(dx) > contactHalfWidth) {
                return Eigen::Vector2d(0.0, 0.0);
            }

            const double relativeX = dx / contactHalfWidth;
            const double pressure = maxPressure * (1.0 - relativeX * relativeX);
            return Eigen::Vector2d(pressure * normal);
        });
}
