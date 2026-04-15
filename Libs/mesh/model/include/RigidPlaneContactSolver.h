#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

#include "ContactTypes.h"

class Assembly;
class Material;
class Node;
class PlaneIsoparametricElement;

class RigidPlaneContactSolver {
public:
    RigidPlaneContactSolver(std::shared_ptr<Assembly> assembly,
        const RigidPlane2D& plane,
        double penaltyParameter);

    void setContactFacets(std::vector<ContactFacet> facets);
    void setPlane(const RigidPlane2D& plane);
    void setPenaltyParameter(double penaltyParameter);

    const RigidPlane2D& getPlane() const { return plane_; }
    double getPenaltyParameter() const { return penaltyParameter_; }
    const std::vector<ContactFacet>& getContactFacets() const { return facets_; }

    void assembleContact(const Eigen::VectorXd& fullDisplacements,
        Eigen::SparseMatrix<double>& contactK,
        Eigen::VectorXd& contactF,
        ContactState& state) const;

private:
    static std::pair<double, double> mapSurfaceCoordinates(int surfaceIndex, double surfaceParameter);
    static double computeSurfaceJacobian(int surfaceIndex,
        const PlaneIsoparametricElement& element,
        const std::vector<std::shared_ptr<Node>>& elementNodes,
        double xi,
        double eta);

    std::vector<std::shared_ptr<Node>> getElementNodes(int elementId) const;
    Eigen::Vector4d evaluateScalarShapeFunctions(const PlaneIsoparametricElement& element,
        double xi,
        double eta) const;

    std::shared_ptr<Assembly> assembly_;
    RigidPlane2D plane_;
    double penaltyParameter_ = 1.0e6;
    std::vector<ContactFacet> facets_;
};
