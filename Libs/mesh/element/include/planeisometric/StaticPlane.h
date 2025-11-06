#pragma once

#include "Planeisoparametric.h"

class StaticPlaneElement : public PlaneIsoparametricElement {
public:
    StaticPlaneElement(const Eigen::MatrixXd& nodeCoords,
        double youngsModulus = 2.0e11,
        double poissonsRatio = 0.3,
        double thickness = 0.1)
        : PlaneIsoparametricElement(nodeCoords, youngsModulus, poissonsRatio, thickness) {
    }

    std::string type() const override { return "StaticPlaneElement"; }

    // ƒополнительные методы дл€ статических задач
    Eigen::VectorXd computeEquivalentNodalForces(const Eigen::VectorXd& bodyForces) const;
};