#include "material.h"
#include <stdexcept>
#include <iostream>

Eigen::Matrix3d Material::getElasticityMatrix() const {
    double E = youngsModulus_;
    double nu = poissonsRatio_;

 /*   double factor = E / (1 - nu * nu);

    Eigen::Matrix3d D;
    D << 1, nu, 0,
        nu, 1, 0,
        0, 0, (1 - nu) / 2;

    D = factor * D;*/
    double factor = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
    Eigen::Matrix3d D;
    D << 1, nu/(1.0-nu), 0,
        nu / (1.0 - nu), 1, 0,
        0, 0, (1.0 - 2.0*nu) / 2/(1.0-nu);

    return D*factor;
}