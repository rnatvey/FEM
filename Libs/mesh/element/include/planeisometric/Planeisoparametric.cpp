
#include "Planeisoparametric.h"
#include "solver.h"
#include "node.h"
#include "material.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

LinearSolver solver;

PlaneIsoparametricElement::PlaneIsoparametricElement(int id, const std::vector<int>& nodeIds, int materialId)
    : BaseElement(id, nodeIds, materialId)
{
    if (nodeIds.size() != 4) {
        throw std::invalid_argument("PlaneIsoparametricElement requires exactly 4 nodes");
    }
}

Eigen::MatrixXd PlaneIsoparametricElement::computeStiffnessMatrix(const std::vector<std::shared_ptr<Node>>& nodes,
    const std::shared_ptr<Material>& material) const
{
    LinearSolver solver;

    auto integrand = [this, &nodes, &material](double xi, double eta) -> Eigen::MatrixXd {
        Eigen::MatrixXd B = strainDisplacementMatrix(xi, eta, nodes, material);
        Eigen::Matrix2d J = jacobian(xi, eta, nodes);
        double detJ = J.determinant();

        if (detJ <= 0) {
            throw std::runtime_error("Negative Jacobian determinant in PlaneIsoparametricElement");
        }

        return B.transpose() * material->getElasticityMatrix() * B * detJ * material->getThickness();
        };

    return solver.computeGaussIntegral(integrand, 2);
}

Eigen::MatrixXd PlaneIsoparametricElement::shapeFunctions(double xi, double eta) const {
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(2, 8);
    Eigen::Vector4d shapes = shapeFunctionsLocal(xi, eta);

    for (int i = 0; i < 4; ++i) {
        N(0, 2 * i) = shapes(i);
        N(1, 2 * i + 1) = shapes(i);
    }

    return N;
}

Eigen::Vector4d PlaneIsoparametricElement::shapeFunctionsLocal(double xi, double eta) const {
    Eigen::Vector4d N;
    N(0) = 0.25 * (1 - xi) * (1 - eta);
    N(1) = 0.25 * (1 + xi) * (1 - eta);
    N(2) = 0.25 * (1 + xi) * (1 + eta);
    N(3) = 0.25 * (1 - xi) * (1 + eta);
    return N;
}

Eigen::MatrixXd PlaneIsoparametricElement::shapeFunctionsDerivatives(double xi, double eta) const {
    return shapeFunctionsDerivativesLocal(xi, eta);
}

Eigen::MatrixXd PlaneIsoparametricElement::shapeFunctionsDerivativesLocal(double xi, double eta) const {
    Eigen::MatrixXd dN = Eigen::MatrixXd::Zero(4, 2);

    // dN/dxi
    dN(0, 0) = -0.25 * (1 - eta);
    dN(1, 0) = 0.25 * (1 - eta);
    dN(2, 0) = 0.25 * (1 + eta);
    dN(3, 0) = -0.25 * (1 + eta);

    // dN/deta
    dN(0, 1) = -0.25 * (1 - xi);
    dN(1, 1) = -0.25 * (1 + xi);
    dN(2, 1) = 0.25 * (1 + xi);
    dN(3, 1) = 0.25 * (1 - xi);

    return dN;
}

Eigen::Matrix2d PlaneIsoparametricElement::jacobian(double xi, double eta, const std::vector<std::shared_ptr<Node>>& nodes) const {
    Eigen::MatrixXd dN = shapeFunctionsDerivativesLocal(xi, eta);
    Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
    Eigen::MatrixXd nodeCoords = getNodalCoordinates(nodes);

    for (int i = 0; i < 4; ++i) {
        J(0, 0) += dN(i, 0) * nodeCoords(i, 0); // dx/dxi
        J(0, 1) += dN(i, 0) * nodeCoords(i, 1); // dy/dxi
        J(1, 0) += dN(i, 1) * nodeCoords(i, 0); // dx/deta
        J(1, 1) += dN(i, 1) * nodeCoords(i, 1); // dy/deta
    }

    return J;
}

Eigen::MatrixXd PlaneIsoparametricElement::strainDisplacementMatrix(double xi, double eta,
    const std::vector<std::shared_ptr<Node>>& nodes,
    const std::shared_ptr<Material>& material) const
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
    Eigen::MatrixXd dN = shapeFunctionsDerivativesLocal(xi, eta);
    Eigen::Matrix2d J = jacobian(xi, eta, nodes);

  /*  std::cout << "=== B-matrix Debug ===" << std::endl;
    std::cout << "xi=" << xi << ", eta=" << eta << std::endl;
    std::cout << "dN_local:\n" << dN << std::endl;
    std::cout << "Jacobian:\n" << J << std::endl;
    std::cout << "det(J)=" << J.determinant() << std::endl;*/

    if (J.determinant() <= 0) {
        throw std::runtime_error("Negative Jacobian determinant");
    }

    Eigen::Matrix2d invJ = J.inverse();
    /*std::cout << "inv(J):\n" << invJ << std::endl;*/

    // Преобразование производных в глобальные координаты
    for (int i = 0; i < 4; ++i) {
        double dNdx = invJ(0, 0) * dN(i, 0) + invJ(0, 1) * dN(i, 1);
        double dNdy = invJ(1, 0) * dN(i, 0) + invJ(1, 1) * dN(i, 1);

       /* std::cout << "Node " << i << ": dNdx=" << dNdx << ", dNdy=" << dNdy << std::endl;*/

        int col = 2 * i;
        B(0, col) = dNdx;        // epsilon_xx
        B(1, col + 1) = dNdy;        // epsilon_yy
        B(2, col) = dNdy;        // gamma_xy
        B(2, col + 1) = dNdx;
    }

    //std::cout << "B-matrix:\n" << B << std::endl;
    //std::cout << "B-matrix norm: " << B.norm() << std::endl;

    return B;
}

Eigen::Vector3d PlaneIsoparametricElement::computeStress(double xi, double eta,
    const Eigen::VectorXd& displacements,
    const std::vector<std::shared_ptr<Node>>& nodes,
    const std::shared_ptr<Material>& material) const
{
    Eigen::Vector3d strain = computeStrain(xi, eta, displacements, nodes, material);
    return material->getElasticityMatrix() * strain;
}

Eigen::Vector3d PlaneIsoparametricElement::computeStrain(double xi, double eta,
    const Eigen::VectorXd& displacements,
    const std::vector<std::shared_ptr<Node>>& nodes,
    const std::shared_ptr<Material>& material) const
{
    Eigen::MatrixXd B = strainDisplacementMatrix(xi, eta, nodes, material);
    return B * displacements;
}

Eigen::VectorXd PlaneIsoparametricElement::computeEquivalentNodalForces(const Eigen::VectorXd& bodyForces,
    const std::vector<std::shared_ptr<Node>>& nodes,
    const std::shared_ptr<Material>& material) const
{
    if (bodyForces.size() != 2) {
        throw std::invalid_argument("Body forces must be 2D vector [fx, fy]");
    }

    auto forceIntegrand = [this, &nodes, &material, &bodyForces](double xi, double eta) -> Eigen::MatrixXd {
        Eigen::MatrixXd N = shapeFunctions(xi, eta);
        Eigen::Matrix2d J = jacobian(xi, eta, nodes);
        double detJ = J.determinant();

        if (detJ <= 0) {
            throw std::runtime_error("Negative Jacobian determinant");
        }

        Eigen::Vector2d b(bodyForces[0], bodyForces[1]);
        Eigen::VectorXd fe = N.transpose() * b * detJ * material->getThickness();

        return fe;
        };

    Eigen::MatrixXd forceMatrix = solver.computeGaussIntegral(forceIntegrand, 2);
    return forceMatrix.col(0);
}

bool PlaneIsoparametricElement::isValid(const std::vector<std::shared_ptr<Node>>& nodes) const {
    // nodes должны быть узлами ЭТОГО элемента (4 узла), а не все узлы сборки

    // Проверяем, что передано ровно 4 узла (для Quad4 элемента)
    if (nodes.size() != 4) {
        std::cout << "WARNING: isValid called with " << nodes.size()
            << " nodes, expected 4 for Quad4 element" << std::endl;
        return false;
    }

    try {
        // Проверяем якобиан в центре элемента
        Eigen::Matrix2d J = jacobian(0, 0, nodes);
        double detJ = J.determinant();

        //std::cout << "Element " << id_ << " Jacobian det at center: " << detJ << std::endl;

        if (detJ <= 0) {
            std::cout << "ERROR: Negative or zero Jacobian determinant: " << detJ << std::endl;
            return false;
        }

        // Дополнительная проверка: якобиан в углах
        std::vector<std::pair<double, double>> corners = {
            {-1, -1}, {1, -1}, {1, 1}, {-1, 1}
        };

        for (const auto& corner : corners) {
            double xi = corner.first;
            double eta = corner.second;

            Eigen::Matrix2d Jcorner = jacobian(xi, eta, nodes);
            double detJcorner = Jcorner.determinant();

            if (detJcorner <= 0) {
                std::cout << "ERROR: Negative Jacobian at corner ("
                    << xi << ", " << eta << "): " << detJcorner << std::endl;
                return false;
            }
        }

        return true;

    }
    catch (const std::exception& e) {
        std::cout << "ERROR in element validation: " << e.what() << std::endl;
        return false;
    }
}

Eigen::Vector2d PlaneIsoparametricElement::getSurfaceNormal(int surfaceIndex, const std::vector<std::shared_ptr<Node>>& nodes) const {
    Eigen::MatrixXd nodeCoords = getNodalCoordinates(nodes);
    Eigen::Vector2d normal = Eigen::Vector2d::Zero();

    switch (surfaceIndex) {
    case 0: // нижняя (узлы 0-1)
        normal = (nodeCoords.row(1) - nodeCoords.row(0)).transpose();
        break;
    case 1: // правая (узлы 1-2)  
        normal = (nodeCoords.row(2) - nodeCoords.row(1)).transpose();
        break;
    case 2: // верхняя (узлы 2-3)
        normal = (nodeCoords.row(3) - nodeCoords.row(2)).transpose();
        break;
    case 3: // левая (узлы 3-0)
        normal = (nodeCoords.row(0) - nodeCoords.row(3)).transpose();
        break;
    default:
        throw std::invalid_argument("Invalid surface index");
    }

    normal = Eigen::Vector2d(-normal.y(), normal.x()).normalized();
    return normal;
}

std::vector<Eigen::Vector2d> PlaneIsoparametricElement::getSurfacePoints(int surfaceIndex, const std::vector<std::shared_ptr<Node>>& nodes) const {
    Eigen::MatrixXd nodeCoords = getNodalCoordinates(nodes);
    std::vector<Eigen::Vector2d> points;

    switch (surfaceIndex) {
    case 0: // нижняя
        points.push_back(nodeCoords.row(0).transpose());
        points.push_back(nodeCoords.row(1).transpose());
        break;
    case 1: // правая
        points.push_back(nodeCoords.row(1).transpose());
        points.push_back(nodeCoords.row(2).transpose());
        break;
    case 2: // верхняя
        points.push_back(nodeCoords.row(2).transpose());
        points.push_back(nodeCoords.row(3).transpose());
        break;
    case 3: // левая
        points.push_back(nodeCoords.row(3).transpose());
        points.push_back(nodeCoords.row(0).transpose());
        break;
    default:
        throw std::invalid_argument("Invalid surface index");
    }

    return points;
}

Eigen::MatrixXd PlaneIsoparametricElement::getNodalCoordinates(
    const std::vector<std::shared_ptr<Node>>& nodes) const {

    // nodes должны быть 4 узлами этого элемента
    if (nodes.size() != 4) {
        throw std::invalid_argument("PlaneIsoparametricElement requires exactly 4 nodes");
    }

    Eigen::MatrixXd coords(4, 2);

    // Просто берем координаты из переданных узлов
    // Предполагаем, что nodes переданы в правильном порядке
    for (int i = 0; i < 4; ++i) {
        coords.row(i) = nodes[i]->getCoordinates().transpose();
    }

    return coords;
}