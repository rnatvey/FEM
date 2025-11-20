#include "material.h"
#include <stdexcept>
#include <iostream>

Eigen::Matrix3d Material::getElasticityMatrix() const {
    double E = youngsModulus_;
    double nu = poissonsRatio_;

    //std::cout << "=== Material Debug ===" << std::endl;
    //std::cout << "E=" << E << ", nu=" << nu << std::endl;

    // Äëÿ ïëîñêîãî ÍÀÏÐßÆÅÍÍÎÃÎ ñîñòîÿíèÿ
    double factor = E / (1 - nu * nu);
  /*  std::cout << "Factor = E/(1-nu^2) = " << factor << std::endl;*/

    Eigen::Matrix3d D;
    D << 1, nu, 0,
        nu, 1, 0,
        0, 0, (1 - nu) / 2;

    D = factor * D;

    //std::cout << "D-matrix:\n" << D << std::endl;
    //std::cout << "D-matrix norm: " << D.norm() << std::endl;

    return D;
}