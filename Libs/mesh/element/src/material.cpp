#include "material.h"
#include <stdexcept>

Eigen::Matrix3d Material::getElasticityMatrix() const {
    double E = youngsModulus_;
    double nu = poissonsRatio_;
    double factor = E / (1.0 - nu * nu);
    
    Eigen::Matrix3d D;
    D << 1.0,  nu, 0.0,
         nu, 1.0,  0.0,
         0.0,  0.0,  (1.0 - nu) / 2.0;
    
    return factor * D;
}