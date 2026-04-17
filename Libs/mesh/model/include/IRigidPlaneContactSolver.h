#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string_view>
#include <vector>

#include "ContactTypes.h"

class IRigidPlaneContactSolver {
public:
    virtual ~IRigidPlaneContactSolver() = default;

    virtual RigidPlaneContactMethod getMethod() const = 0;
    virtual std::string_view getMethodName() const = 0;

    virtual void setContactFacets(std::vector<ContactFacet> facets) = 0;
    virtual void setPlane(const RigidPlane2D& plane) = 0;

    virtual const RigidPlane2D& getPlane() const = 0;
    virtual const std::vector<ContactFacet>& getContactFacets() const = 0;

    virtual void assembleContact(const Eigen::VectorXd& fullDisplacements,
        Eigen::SparseMatrix<double>& contactK,
        Eigen::VectorXd& contactF,
        ContactState& state) const = 0;

    virtual ContactSolverUpdateInfo updateState(const Eigen::VectorXd& fullDisplacements) = 0;
    virtual void resetState() = 0;
};
