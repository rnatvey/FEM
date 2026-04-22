#pragma once

#include <string_view>

#include "FiniteStrainMaterial.h"

// Compressible Neo-Hookean baseline used for the first finite-strain branch.
//
// Stored parameters:
//   mu - shear modulus
//   K  - bulk modulus
//
// Internally the implementation uses
//   lambda = K - 2/3 mu
//
// and the 3D compressible Neo-Hookean energy
//   W = mu/2 * (I1 - 3 - 2 ln J) + lambda/2 * (ln J)^2
//
// with plane strain enforced by embedding the in-plane deformation gradient
// into a 3D state with F33 = 1.
//
// Important limitation of the current baseline:
//   when this material is used with the current displacement-only fully
//   integrated Q4 finite-strain element, nearly incompressible regimes may
//   suffer from volumetric locking. SRI and mixed u/p are intentionally left
//   for later steps instead of being introduced prematurely here.
class NeoHookeanMaterial final : public FiniteStrainMaterial {
public:
    NeoHookeanMaterial(int id,
        double shearModulus,
        double bulkModulus,
        double thickness = 0.1);

    static NeoHookeanMaterial fromYoungsModulusAndPoissonsRatio(int id,
        double youngsModulus,
        double poissonsRatio,
        double thickness = 0.1);

    int getId() const override { return id_; }
    std::string_view name() const override { return "neo_hookean"; }
    double getThickness() const override { return thickness_; }

    MaterialPointResponse evaluatePlaneStrain(
        const PlaneStrainDeformationGradient& deformationGradient) const override;
    double getShearModulus() const { return shearModulus_; }
    double getBulkModulus() const { return bulkModulus_; }
    double getLameFirstParameter() const;
    double getEquivalentYoungsModulus() const;
    double getEffectivePoissonsRatio() const;
    double getBulkToShearRatio() const;
    bool isNearlyIncompressible(double threshold = 0.475) const;

private:
    static double computeShearModulus(double youngsModulus, double poissonsRatio);
    static double computeBulkModulus(double youngsModulus, double poissonsRatio);
    static double computeLameFirstParameter(double shearModulus, double bulkModulus);

    MaterialPointResponse evaluateFromRightCauchyGreen(const Eigen::Matrix3d& rightCauchyGreen,
        double jacobianDeterminant) const;

    int id_ = -1;
    double shearModulus_ = 0.0;
    double bulkModulus_ = 0.0;
    double thickness_ = 0.1;
};
