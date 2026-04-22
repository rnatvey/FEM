#include "planeisometric/FiniteStrainQ4Element.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

constexpr double kGaussPoint = 0.5773502691896257;
constexpr std::array<double, 2> kGaussPoints = {-kGaussPoint, kGaussPoint};
constexpr std::array<double, 2> kGaussWeights = {1.0, 1.0};

Eigen::Matrix2d planeStrainStressVoigtToTensor(const Eigen::Vector3d& stressVoigt) {
    Eigen::Matrix2d stress = Eigen::Matrix2d::Zero();
    stress(0, 0) = stressVoigt(0);
    stress(1, 1) = stressVoigt(1);
    stress(0, 1) = stressVoigt(2);
    stress(1, 0) = stressVoigt(2);
    return stress;
}

Eigen::Vector3d greenLagrangeTensorToVoigt(const Eigen::Matrix2d& strainTensor) {
    return Eigen::Vector3d(strainTensor(0, 0), strainTensor(1, 1), strainTensor(0, 1));
}

double polygonArea(const Eigen::Matrix<double, 2, 4>& coordinates) {
    double twiceArea = 0.0;
    for (int i = 0; i < 4; ++i) {
        const int j = (i + 1) % 4;
        twiceArea +=
            coordinates(0, i) * coordinates(1, j) - coordinates(0, j) * coordinates(1, i);
    }
    return 0.5 * std::abs(twiceArea);
}

} // namespace

FiniteStrainQ4Element::FiniteStrainQ4Element(int id,
    const std::vector<int>& nodeIds,
    int materialId)
    : id_(id), nodeIds_(nodeIds), materialId_(materialId) {
    if (nodeIds_.size() != kNodeCount) {
        throw std::invalid_argument("FiniteStrainQ4Element requires exactly 4 nodes");
    }
}

FiniteStrainElementResponse FiniteStrainQ4Element::evaluateResponse(
    const std::vector<std::shared_ptr<Node>>& nodes,
    const Eigen::VectorXd& elementDisplacements,
    const FiniteStrainMaterial& material,
    bool computeTangent) const {
    FiniteStrainElementResponse response =
        evaluateResponseWithoutTangent(nodes, elementDisplacements, material);

    if (!computeTangent) {
        response.tangent = Eigen::MatrixXd::Zero(kDofCount, kDofCount);
        return response;
    }

    response.tangent = Eigen::MatrixXd::Zero(kDofCount, kDofCount);
    const Eigen::Matrix<double, 2, 4> coordinates = nodalCoordinates(nodes);
    const double characteristicLength = std::max(std::sqrt(std::max(polygonArea(coordinates), 0.0)), 1.0);

    for (int dofIndex = 0; dofIndex < kDofCount; ++dofIndex) {
        const double perturbation = finiteDifferenceStep(
            elementDisplacements(dofIndex), characteristicLength);

        Eigen::VectorXd displacedPlus = elementDisplacements;
        displacedPlus(dofIndex) += perturbation;
        const Eigen::VectorXd internalForcePlus =
            evaluateResponseWithoutTangent(nodes, displacedPlus, material).internalForce;

        Eigen::VectorXd displacedMinus = elementDisplacements;
        displacedMinus(dofIndex) -= perturbation;
        const Eigen::VectorXd internalForceMinus =
            evaluateResponseWithoutTangent(nodes, displacedMinus, material).internalForce;

        response.tangent.col(dofIndex) =
            (internalForcePlus - internalForceMinus) / (2.0 * perturbation);
    }

    return response;
}

Eigen::Vector4d FiniteStrainQ4Element::shapeFunctionsLocal(double xi, double eta) const {
    Eigen::Vector4d shapeFunctions;
    shapeFunctions(0) = 0.25 * (1.0 + xi) * (1.0 + eta);
    shapeFunctions(1) = 0.25 * (1.0 - xi) * (1.0 + eta);
    shapeFunctions(2) = 0.25 * (1.0 - xi) * (1.0 - eta);
    shapeFunctions(3) = 0.25 * (1.0 + xi) * (1.0 - eta);
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 2> FiniteStrainQ4Element::shapeFunctionsDerivativesLocal(
    double xi,
    double eta) const {
    Eigen::Matrix<double, 4, 2> derivatives = Eigen::Matrix<double, 4, 2>::Zero();
    derivatives(0, 0) = 0.25 * (1.0 + eta);
    derivatives(0, 1) = 0.25 * (1.0 + xi);
    derivatives(1, 0) = -0.25 * (1.0 + eta);
    derivatives(1, 1) = 0.25 * (1.0 - xi);
    derivatives(2, 0) = -0.25 * (1.0 - eta);
    derivatives(2, 1) = -0.25 * (1.0 - xi);
    derivatives(3, 0) = 0.25 * (1.0 - eta);
    derivatives(3, 1) = -0.25 * (1.0 + xi);
    return derivatives;
}

Eigen::Matrix2d FiniteStrainQ4Element::referenceJacobian(
    double xi,
    double eta,
    const std::vector<std::shared_ptr<Node>>& nodes) const {
    const Eigen::Matrix<double, 2, 4> coordinates = nodalCoordinates(nodes);
    const Eigen::Matrix<double, 4, 2> localDerivatives =
        shapeFunctionsDerivativesLocal(xi, eta);

    Eigen::Matrix2d jacobian = Eigen::Matrix2d::Zero();
    for (int nodeIndex = 0; nodeIndex < kNodeCount; ++nodeIndex) {
        jacobian(0, 0) += localDerivatives(nodeIndex, 0) * coordinates(0, nodeIndex);
        jacobian(0, 1) += localDerivatives(nodeIndex, 0) * coordinates(1, nodeIndex);
        jacobian(1, 0) += localDerivatives(nodeIndex, 1) * coordinates(0, nodeIndex);
        jacobian(1, 1) += localDerivatives(nodeIndex, 1) * coordinates(1, nodeIndex);
    }
    return jacobian;
}

FiniteStrainElementResponse FiniteStrainQ4Element::evaluateResponseWithoutTangent(
    const std::vector<std::shared_ptr<Node>>& nodes,
    const Eigen::VectorXd& elementDisplacements,
    const FiniteStrainMaterial& material) const {
    if (nodes.size() != kNodeCount) {
        throw std::invalid_argument("FiniteStrainQ4Element requires exactly 4 element nodes");
    }
    if (elementDisplacements.size() != kDofCount) {
        throw std::invalid_argument("FiniteStrainQ4Element expects 8 displacement DOFs");
    }

    FiniteStrainElementResponse response;
    response.internalForce = Eigen::VectorXd::Zero(kDofCount);
    response.tangent = Eigen::MatrixXd::Zero(kDofCount, kDofCount);
    response.gaussPointResults.reserve(kGaussPoints.size() * kGaussPoints.size());

    Eigen::Matrix<double, 2, 4> displacementMatrix = Eigen::Matrix<double, 2, 4>::Zero();
    for (int nodeIndex = 0; nodeIndex < kNodeCount; ++nodeIndex) {
        displacementMatrix(0, nodeIndex) = elementDisplacements(2 * nodeIndex);
        displacementMatrix(1, nodeIndex) = elementDisplacements(2 * nodeIndex + 1);
    }

    const double thickness = material.getThickness();

    for (double xi : kGaussPoints) {
        for (double eta : kGaussPoints) {
            double detJ0 = 0.0;
            const Eigen::Matrix<double, 4, 2> gradientsReference =
                shapeFunctionGradientsReference(xi, eta, nodes, detJ0);
            if (!(detJ0 > 0.0)) {
                throw std::runtime_error(
                    "FiniteStrainQ4Element encountered a non-positive reference Jacobian");
            }

            Eigen::Matrix2d deformationGradient = Eigen::Matrix2d::Identity();
            deformationGradient += displacementMatrix * gradientsReference;

            const MaterialPointResponse materialResponse =
                material.evaluatePlaneStrain(deformationGradient);
            const Eigen::Matrix2d secondPiolaStress =
                planeStrainStressVoigtToTensor(materialResponse.secondPiolaStressVoigt);
            const Eigen::Matrix2d firstPiolaStress = deformationGradient * secondPiolaStress;
            const Eigen::Matrix2d rightCauchyGreen =
                deformationGradient.transpose() * deformationGradient;
            const Eigen::Matrix2d greenLagrangeTensor =
                0.5 * (rightCauchyGreen - Eigen::Matrix2d::Identity());

            const double integrationWeight =
                kGaussWeights[0] * kGaussWeights[1] * detJ0 * thickness;

            for (int nodeIndex = 0; nodeIndex < kNodeCount; ++nodeIndex) {
                const Eigen::Vector2d gradient = gradientsReference.row(nodeIndex).transpose();
                const Eigen::Vector2d nodalForce = firstPiolaStress * gradient * integrationWeight;
                response.internalForce.segment<2>(2 * nodeIndex) += nodalForce;
            }

            FiniteStrainGaussPointResult gaussPointResult;
            gaussPointResult.xi = xi;
            gaussPointResult.eta = eta;
            gaussPointResult.integrationWeight = integrationWeight;
            gaussPointResult.referenceJacobianDeterminant = detJ0;
            gaussPointResult.jacobianDeterminant = materialResponse.jacobianDeterminant;
            gaussPointResult.strainEnergyDensity = materialResponse.strainEnergyDensity;
            gaussPointResult.deformationGradient = deformationGradient;
            gaussPointResult.greenLagrangeStrainVoigt =
                greenLagrangeTensorToVoigt(greenLagrangeTensor);
            gaussPointResult.secondPiolaStressVoigt =
                materialResponse.secondPiolaStressVoigt;
            response.gaussPointResults.push_back(gaussPointResult);
            response.totalStrainEnergy +=
                materialResponse.strainEnergyDensity * integrationWeight;
        }
    }

    return response;
}

Eigen::Matrix<double, 4, 2> FiniteStrainQ4Element::shapeFunctionGradientsReference(
    double xi,
    double eta,
    const std::vector<std::shared_ptr<Node>>& nodes,
    double& detJ0) const {
    const Eigen::Matrix<double, 4, 2> localDerivatives =
        shapeFunctionsDerivativesLocal(xi, eta);
    const Eigen::Matrix2d jacobian0 = referenceJacobian(xi, eta, nodes);
    detJ0 = jacobian0.determinant();
    if (!(detJ0 > 0.0)) {
        throw std::runtime_error(
            "FiniteStrainQ4Element encountered a non-positive reference Jacobian");
    }

    const Eigen::Matrix2d inverseJacobian0 = jacobian0.inverse();
    Eigen::Matrix<double, 4, 2> gradientsReference = Eigen::Matrix<double, 4, 2>::Zero();
    for (int nodeIndex = 0; nodeIndex < kNodeCount; ++nodeIndex) {
        const Eigen::Vector2d localGradient =
            localDerivatives.row(nodeIndex).transpose();
        gradientsReference.row(nodeIndex) =
            (inverseJacobian0 * localGradient).transpose();
    }

    return gradientsReference;
}

Eigen::Matrix<double, 2, 4> FiniteStrainQ4Element::nodalCoordinates(
    const std::vector<std::shared_ptr<Node>>& nodes) const {
    Eigen::Matrix<double, 2, 4> coordinates = Eigen::Matrix<double, 2, 4>::Zero();
    for (int nodeIndex = 0; nodeIndex < kNodeCount; ++nodeIndex) {
        if (!nodes[nodeIndex]) {
            throw std::invalid_argument("FiniteStrainQ4Element received a null node");
        }
        coordinates(0, nodeIndex) = nodes[nodeIndex]->getCoordinates().x();
        coordinates(1, nodeIndex) = nodes[nodeIndex]->getCoordinates().y();
    }
    return coordinates;
}

double FiniteStrainQ4Element::finiteDifferenceStep(double dofValue, double characteristicLength) const {
    const double scale = std::max({1.0, std::abs(dofValue), characteristicLength});
    return 1.0e-7 * scale;
}
