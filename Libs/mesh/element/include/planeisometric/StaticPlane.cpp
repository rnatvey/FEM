
#include "StaticPlane.h"
#include "solver.h"
#include <stdexcept>

Eigen::VectorXd StaticPlaneElement::computeEquivalentNodalForces(const Eigen::VectorXd& bodyForces) const {
    if (bodyForces.size() != 2) {
        throw std::invalid_argument("Body forces must be 2D vector [fx, fy]");
    }

    LinearSolver solver;

    // Лямбда-функция для эквивалентных узловых сил от объемных сил
    auto forceIntegrand = [this, &bodyForces](double xi, double eta) -> Eigen::MatrixXd {
        Eigen::MatrixXd N = shapeFunctions(xi, eta);  // Матрица формы 2x8
        Eigen::Matrix2d J = jacobian(xi, eta);
        double detJ = J.determinant();

        if (detJ <= 0) {
            throw std::runtime_error("Negative Jacobian determinant in StaticPlaneElement");
        }

        // fe = N^T * b * t * |J|, где b - вектор объемных сил
        Eigen::Vector2d b(bodyForces[0], bodyForces[1]);
        Eigen::VectorXd fe = N.transpose() * b * detJ * thickness_;

        return fe; // Возвращаем как матрицу для совместимости с computeGaussIntegral
        };

    // Интегрируем с помощью LinearSolver
    Eigen::MatrixXd forceMatrix = solver.computeGaussIntegral(forceIntegrand, 2);

    // Преобразуем результат в вектор
    return forceMatrix.col(0); // forceMatrix будет 8x1 в данном случае
}