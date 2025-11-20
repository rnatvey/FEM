#pragma once
#include <Eigen/Dense>

class ConcentratedForce {
public:
    ConcentratedForce(int nodeId, double fx, double fy);

    int getNodeId() const { return nodeId_; }
    double getForceX() const { return fx_; }
    double getForceY() const { return fy_; }
    const Eigen::Vector2d& getForce() const { return force_; }

private:
    int nodeId_;
    Eigen::Vector2d force_;
    double fx_, fy_;
};