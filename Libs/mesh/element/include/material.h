#pragma once
#include <Eigen/Dense>

class Material {
public:
    Material(int id, double youngsModulus, double poissonsRatio, double thickness = 0.1)
        : id_(id), youngsModulus_(youngsModulus), poissonsRatio_(poissonsRatio), thickness_(thickness) {
    }

    int getId() const { return id_; }
    double getYoungsModulus() const { return youngsModulus_; }
    double getPoissonsRatio() const { return poissonsRatio_; }
    double getThickness() const { return thickness_; }

    Eigen::Matrix3d getElasticityMatrix() const;

private:
    int id_;
    double youngsModulus_;
    double poissonsRatio_;
    double thickness_;
};