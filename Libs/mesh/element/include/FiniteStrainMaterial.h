#pragma once

#include <Eigen/Dense>

#include <string_view>

#include "MaterialPointResponse.h"

class FiniteStrainMaterial {
public:
    using PlaneStrainDeformationGradient = Eigen::Matrix2d;

    virtual ~FiniteStrainMaterial() = default;

    // Evaluates a 2D plane-strain material point response for the given
    // in-plane deformation gradient. The returned response must be consistent
    // with the total Lagrangian formulation:
    // - stress: 2nd Piola-Kirchhoff, reduced Voigt [11, 22, 12]
    // - tangent: material tangent dS/dE in the same reduced basis
    virtual MaterialPointResponse evaluatePlaneStrain(
        const PlaneStrainDeformationGradient& deformationGradient) const = 0;

    virtual std::string_view name() const = 0;
};
