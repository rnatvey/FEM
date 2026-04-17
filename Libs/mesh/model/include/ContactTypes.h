#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string_view>
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

enum class RigidPlaneContactMethod {
    Penalty,
    AugmentedLagrangian
};

inline std::string_view rigidPlaneContactMethodName(const RigidPlaneContactMethod method) {
    switch (method) {
    case RigidPlaneContactMethod::Penalty:
        return "penalty";
    case RigidPlaneContactMethod::AugmentedLagrangian:
        return "augmented_lagrangian";
    }
    return "unknown";
}

struct AugmentedLagrangianSettings {
    double augmentationParameter = 1.0e6;
    double multiplierTolerance = 1.0e-8;
};

struct ContactGaussPointState {
    int facetId = -1;
    int gaussPointIndex = -1;
    double surfaceParameter = 0.0;
    double xi = 0.0;
    double eta = 0.0;
    double lambdaN = 0.0;
    double gapN = 0.0;
    double penetrationN = 0.0;
    double pressureN = 0.0;
    bool active = false;
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
    RigidPlaneContactMethod method = RigidPlaneContactMethod::Penalty;
    std::vector<int> activeFacetIds;
    std::vector<ContactFacetResult> facetResults;
    std::vector<ContactGaussPointState> gaussPointStates;
    double maxPenetration = 0.0;
    double contactForceNorm = 0.0;
    int activeGaussPointCount = 0;
};

struct ContactSolverUpdateInfo {
    bool converged = true;
    double stateUpdateNorm = 0.0;
    double relativeStateUpdateNorm = 0.0;
    int activeGaussPointCount = 0;
    double maxNormalMultiplier = 0.0;
    double meanNormalMultiplier = 0.0;
};

struct ContactIterationInfo {
    ContactState state;
    ContactSolverUpdateInfo updateInfo;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::VectorXd force;
};
