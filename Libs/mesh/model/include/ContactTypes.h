#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

struct RigidPlane2D {
    Eigen::Vector2d normal = Eigen::Vector2d(0.0, 1.0);
    double offset = 0.0;

    double signedDistance(const Eigen::Vector2d& point) const {
        return normal.dot(point) - offset;
    }
};

struct ContactFacet {
    int elementId = -1;
    int surfaceIndex = -1;
};

struct ContactState {
    std::vector<int> activeFacetIds;
    double maxPenetration = 0.0;
    double contactForceNorm = 0.0;
};

struct ContactIterationInfo {
    ContactState state;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::VectorXd force;
};
