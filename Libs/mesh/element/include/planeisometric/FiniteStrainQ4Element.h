#pragma once

#include <Eigen/Dense>

#include <memory>
#include <vector>

#include "FiniteStrainElementResponse.h"
#include "FiniteStrainMaterial.h"
#include "node.h"

class FiniteStrainQ4Element {
public:
    static constexpr int kNodeCount = 4;
    static constexpr int kDofCount = 8;

    // Current baseline element:
    // - total Lagrangian
    // - displacement-only Q4
    // - full 2x2 integration for both deviatoric and volumetric response
    //
    // This is intentionally kept simple for the first hyperelastic branch.
    // For nearly incompressible materials it may exhibit volumetric locking.
    // Selective reduced integration and mixed u/p are left as future work.
    FiniteStrainQ4Element(int id, const std::vector<int>& nodeIds, int materialId);

    int getId() const { return id_; }
    const std::vector<int>& getNodeIds() const { return nodeIds_; }
    int getMaterialId() const { return materialId_; }

    FiniteStrainElementResponse evaluateResponse(
        const std::vector<std::shared_ptr<Node>>& nodes,
        const Eigen::VectorXd& elementDisplacements,
        const FiniteStrainMaterial& material,
        bool computeTangent = true) const;

    Eigen::Vector4d shapeFunctionsLocal(double xi, double eta) const;
    Eigen::Matrix<double, 4, 2> shapeFunctionsDerivativesLocal(double xi, double eta) const;
    Eigen::Matrix2d referenceJacobian(double xi,
        double eta,
        const std::vector<std::shared_ptr<Node>>& nodes) const;

private:
    FiniteStrainElementResponse evaluateResponseWithoutTangent(
        const std::vector<std::shared_ptr<Node>>& nodes,
        const Eigen::VectorXd& elementDisplacements,
        const FiniteStrainMaterial& material) const;

    Eigen::Matrix<double, 4, 2> shapeFunctionGradientsReference(
        double xi,
        double eta,
        const std::vector<std::shared_ptr<Node>>& nodes,
        double& detJ0) const;

    Eigen::Matrix<double, 2, 4> nodalCoordinates(
        const std::vector<std::shared_ptr<Node>>& nodes) const;

    double finiteDifferenceStep(double dofValue, double characteristicLength) const;

    int id_ = -1;
    std::vector<int> nodeIds_;
    int materialId_ = -1;
};
