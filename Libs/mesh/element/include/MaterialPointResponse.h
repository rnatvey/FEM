#pragma once

#include <Eigen/Dense>

// Finite-strain baseline for the current 2D plane-strain formulation:
// - stress measure: 2nd Piola-Kirchhoff stress
// - tangent measure: material tangent dS/dE
// - reduced Voigt order: [11, 22, 12] with tensorial shear component E12
//
// The constitutive model is evaluated in plane strain. Internally a material
// may reconstruct the full 3D state (for example with F33 = 1), but the
// element-facing response is kept in reduced 2D form to match the existing Q4
// workflow and keep the first hyperelastic branch focused and compact.
struct MaterialPointResponse {
    using PlaneStrainStress = Eigen::Vector3d;
    using PlaneStrainTangent = Eigen::Matrix3d;

    double strainEnergyDensity = 0.0;
    double jacobianDeterminant = 1.0;
    PlaneStrainStress secondPiolaStressVoigt = PlaneStrainStress::Zero();
    PlaneStrainTangent materialTangentVoigt = PlaneStrainTangent::Zero();
};
