#include "NeoHookeanMaterial.h"

#include <cmath>
#include <stdexcept>

namespace {

constexpr double kNearZero = 1.0e-15;

double tangentComponent(const Eigen::Matrix3d& inverseRightCauchyGreen,
    double lameFirstParameter,
    double shearModulus,
    double logJacobian,
    int i,
    int j,
    int k,
    int l) {
    return lameFirstParameter * inverseRightCauchyGreen(i, j) * inverseRightCauchyGreen(k, l) +
        (shearModulus - lameFirstParameter * logJacobian) *
        (inverseRightCauchyGreen(i, k) * inverseRightCauchyGreen(j, l) +
            inverseRightCauchyGreen(i, l) * inverseRightCauchyGreen(j, k));
}

} // namespace

NeoHookeanMaterial::NeoHookeanMaterial(int id,
    double shearModulus,
    double bulkModulus,
    double thickness)
    : id_(id),
      shearModulus_(shearModulus),
      bulkModulus_(bulkModulus),
      thickness_(thickness) {
    if (!(shearModulus_ > 0.0)) {
        throw std::invalid_argument("Neo-Hookean shear modulus must be positive");
    }
    if (!(bulkModulus_ > 0.0)) {
        throw std::invalid_argument("Neo-Hookean bulk modulus must be positive");
    }
    if (!(thickness_ > 0.0)) {
        throw std::invalid_argument("Neo-Hookean thickness must be positive");
    }
}

NeoHookeanMaterial NeoHookeanMaterial::fromYoungsModulusAndPoissonsRatio(int id,
    double youngsModulus,
    double poissonsRatio,
    double thickness) {
    return NeoHookeanMaterial(id,
        computeShearModulus(youngsModulus, poissonsRatio),
        computeBulkModulus(youngsModulus, poissonsRatio),
        thickness);
}

MaterialPointResponse NeoHookeanMaterial::evaluatePlaneStrain(
    const PlaneStrainDeformationGradient& deformationGradient) const {
    const double jacobianDeterminant = deformationGradient.determinant();
    if (!(jacobianDeterminant > 0.0)) {
        throw std::runtime_error(
            "Neo-Hookean plane-strain deformation gradient must have positive determinant");
    }

    Eigen::Matrix3d deformationGradient3D = Eigen::Matrix3d::Identity();
    deformationGradient3D.topLeftCorner<2, 2>() = deformationGradient;

    const Eigen::Matrix3d rightCauchyGreen =
        deformationGradient3D.transpose() * deformationGradient3D;
    return evaluateFromRightCauchyGreen(rightCauchyGreen, jacobianDeterminant);
}

double NeoHookeanMaterial::getLameFirstParameter() const {
    return computeLameFirstParameter(shearModulus_, bulkModulus_);
}

double NeoHookeanMaterial::getEffectivePoissonsRatio() const {
    const double denominator = 2.0 * (3.0 * bulkModulus_ + shearModulus_);
    if (std::abs(denominator) <= kNearZero) {
        throw std::runtime_error(
            "Cannot recover effective Poisson ratio from Neo-Hookean moduli");
    }

    return (3.0 * bulkModulus_ - 2.0 * shearModulus_) / denominator;
}

double NeoHookeanMaterial::getBulkToShearRatio() const {
    if (std::abs(shearModulus_) <= kNearZero) {
        throw std::runtime_error("Cannot compute bulk-to-shear ratio for zero shear modulus");
    }

    return bulkModulus_ / shearModulus_;
}

bool NeoHookeanMaterial::isNearlyIncompressible(double threshold) const {
    return getEffectivePoissonsRatio() >= threshold;
}

double NeoHookeanMaterial::computeShearModulus(double youngsModulus, double poissonsRatio) {
    if (!(youngsModulus > 0.0)) {
        throw std::invalid_argument("Young's modulus must be positive");
    }

    const double denominator = 2.0 * (1.0 + poissonsRatio);
    if (std::abs(denominator) <= kNearZero) {
        throw std::invalid_argument("Invalid Poisson ratio for shear modulus conversion");
    }

    return youngsModulus / denominator;
}

double NeoHookeanMaterial::computeBulkModulus(double youngsModulus, double poissonsRatio) {
    if (!(youngsModulus > 0.0)) {
        throw std::invalid_argument("Young's modulus must be positive");
    }

    const double denominator = 3.0 * (1.0 - 2.0 * poissonsRatio);
    if (std::abs(denominator) <= kNearZero) {
        throw std::invalid_argument("Invalid Poisson ratio for bulk modulus conversion");
    }

    return youngsModulus / denominator;
}

double NeoHookeanMaterial::computeLameFirstParameter(double shearModulus, double bulkModulus) {
    return bulkModulus - (2.0 / 3.0) * shearModulus;
}

MaterialPointResponse NeoHookeanMaterial::evaluateFromRightCauchyGreen(
    const Eigen::Matrix3d& rightCauchyGreen,
    double jacobianDeterminant) const {
    const double lameFirstParameter =
        computeLameFirstParameter(shearModulus_, bulkModulus_);
    const double logJacobian = std::log(jacobianDeterminant);
    const double firstInvariant = rightCauchyGreen.trace();
    const Eigen::Matrix3d inverseRightCauchyGreen = rightCauchyGreen.inverse();

    MaterialPointResponse response;
    response.jacobianDeterminant = jacobianDeterminant;
    response.strainEnergyDensity =
        0.5 * shearModulus_ * (firstInvariant - 3.0 - 2.0 * logJacobian) +
        0.5 * lameFirstParameter * logJacobian * logJacobian;

    const Eigen::Matrix3d secondPiolaStress =
        shearModulus_ * (Eigen::Matrix3d::Identity() - inverseRightCauchyGreen) +
        lameFirstParameter * logJacobian * inverseRightCauchyGreen;

    response.secondPiolaStressVoigt <<
        secondPiolaStress(0, 0),
        secondPiolaStress(1, 1),
        secondPiolaStress(0, 1);

    Eigen::Matrix3d reducedTangent = Eigen::Matrix3d::Zero();
    reducedTangent(0, 0) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 0, 0, 0, 0);
    reducedTangent(0, 1) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 0, 0, 1, 1);
    reducedTangent(0, 2) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 0, 0, 0, 1);
    reducedTangent(1, 0) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 1, 1, 0, 0);
    reducedTangent(1, 1) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 1, 1, 1, 1);
    reducedTangent(1, 2) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 1, 1, 0, 1);
    reducedTangent(2, 0) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 0, 1, 0, 0);
    reducedTangent(2, 1) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 0, 1, 1, 1);
    reducedTangent(2, 2) = tangentComponent(
        inverseRightCauchyGreen, lameFirstParameter, shearModulus_, logJacobian, 0, 1, 0, 1);
    response.materialTangentVoigt = reducedTangent;

    return response;
}
