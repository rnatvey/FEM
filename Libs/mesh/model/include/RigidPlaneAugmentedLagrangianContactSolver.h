#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

#include "ContactTypes.h"
#include "IRigidPlaneContactSolver.h"

class Assembly;
class Node;
class PlaneIsoparametricElement;

class RigidPlaneAugmentedLagrangianContactSolver : public IRigidPlaneContactSolver {
public:
    RigidPlaneAugmentedLagrangianContactSolver(std::shared_ptr<Assembly> assembly,
        const RigidPlane2D& plane,
        const AugmentedLagrangianSettings& settings);

    RigidPlaneContactMethod getMethod() const override {
        return RigidPlaneContactMethod::AugmentedLagrangian;
    }
    std::string_view getMethodName() const override {
        return rigidPlaneContactMethodName(getMethod());
    }

    void setContactFacets(std::vector<ContactFacet> facets) override;
    void setPlane(const RigidPlane2D& plane) override;
    void setSettings(const AugmentedLagrangianSettings& settings);

    const RigidPlane2D& getPlane() const override { return plane_; }
    const std::vector<ContactFacet>& getContactFacets() const override { return facets_; }
    const AugmentedLagrangianSettings& getSettings() const { return settings_; }

    void assembleContact(const Eigen::VectorXd& fullDisplacements,
        Eigen::SparseMatrix<double>& contactK,
        Eigen::VectorXd& contactF,
        ContactState& state) const override;

    ContactSolverUpdateInfo updateState(const Eigen::VectorXd& fullDisplacements) override;
    void resetState() override;

private:
    static std::pair<double, double> mapSurfaceCoordinates(int surfaceIndex, double surfaceParameter);
    static double computeSurfaceJacobian(int surfaceIndex,
        const PlaneIsoparametricElement& element,
        const std::vector<std::shared_ptr<Node>>& elementNodes,
        double xi,
        double eta);

    void initializeGaussPointStates();
    std::vector<std::shared_ptr<Node>> getElementNodes(int elementId) const;
    Eigen::Vector4d evaluateScalarShapeFunctions(const PlaneIsoparametricElement& element,
        double xi,
        double eta) const;

    std::shared_ptr<Assembly> assembly_;
    RigidPlane2D plane_;
    AugmentedLagrangianSettings settings_;
    std::vector<ContactFacet> facets_;
    std::vector<ContactGaussPointState> gaussPointStates_;
};
