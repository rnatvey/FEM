#pragma once

#include <Eigen/Dense>

#include <vector>

struct FiniteStrainGaussPointResult {
    double xi = 0.0;
    double eta = 0.0;
    double integrationWeight = 0.0;
    double referenceJacobianDeterminant = 0.0;
    double jacobianDeterminant = 1.0;
    double strainEnergyDensity = 0.0;
    Eigen::Matrix2d deformationGradient = Eigen::Matrix2d::Identity();
    Eigen::Vector3d greenLagrangeStrainVoigt = Eigen::Vector3d::Zero();
    Eigen::Vector3d secondPiolaStressVoigt = Eigen::Vector3d::Zero();
};

struct FiniteStrainElementResponse {
    Eigen::VectorXd internalForce;
    Eigen::MatrixXd tangent;
    std::vector<FiniteStrainGaussPointResult> gaussPointResults;
    double totalStrainEnergy = 0.0;
};
