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

struct ContactFacetResult {
    int facetId = -1;
    int elementId = -1;
    int surfaceIndex = -1;
    bool active = false;
    int activeGaussPointCount = 0;
    Eigen::Vector2d referenceMidpoint = Eigen::Vector2d::Zero();
    Eigen::Vector2d deformedMidpoint = Eigen::Vector2d::Zero();
    Eigen::Vector2d normal = Eigen::Vector2d::Zero();
    double thickness = 0.0;
    double facetLength = 0.0;
    double activeLength = 0.0;
    double integratedArea = 0.0;
    double activeArea = 0.0;
    double averageGap = 0.0;
    double averagePenetration = 0.0;
    double maximumPenetration = 0.0;
    double integratedNormalForce = 0.0;
    double averagePressure = 0.0;
};

struct ContactState {
    std::vector<int> activeFacetIds;
    std::vector<ContactFacetResult> facetResults;
    double maxPenetration = 0.0;
    double contactForceNorm = 0.0;
};

struct ContactIterationInfo {
    ContactState state;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::VectorXd force;
};
