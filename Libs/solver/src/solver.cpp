#include "solver.h"
#include <stdexcept>

Eigen::MatrixXd LinearSolver::computeGaussIntegral(
    const std::function<Eigen::MatrixXd(double, double)>& matFunc,
    int numGaussPoints)
{
    auto gaussPoints = generateGaussPoints(numGaussPoints);
    if (gaussPoints.empty()) {
        throw std::runtime_error("No Gauss points generated");
    }

    Eigen::MatrixXd testMat = matFunc(0, 0);
    Eigen::MatrixXd integral = Eigen::MatrixXd::Zero(testMat.rows(), testMat.cols());

    for (const auto& gp : gaussPoints) {
        integral += matFunc(gp.xi, gp.eta) * gp.weight;
    }

    return integral;
}

Eigen::VectorXd LinearSolver::solveSystem(
    const Eigen::SparseMatrix<double>& systemMatrix,
    const Eigen::VectorXd& rightHandSide)
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(systemMatrix);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }
    return solver.solve(rightHandSide);
}

std::vector<LinearSolver::GaussPoint> LinearSolver::generateGaussPoints(int order) {
    std::vector<GaussPoint> points;

    if (order == 2) {
        const double pos = 0.577350269189626;
        const double w = 1.0;

        points = {
            {-pos, -pos, w * w},
            {-pos,  pos, w * w},
            { pos, -pos, w * w},
            { pos,  pos, w * w}
        };
    }

    return points;
}

void LinearSolver::applyBoundaryConditions(Eigen::SparseMatrix<double>& systemMatrix,
    Eigen::VectorXd& rightHandSide,
    const std::vector<int>& fixedDofs) const {
    for (int dof : fixedDofs) {
        // Обнуляем строку и столбец
        for (int k = 0; k < systemMatrix.outerSize(); ++k) {
            systemMatrix.coeffRef(dof, k) = 0.0;
            systemMatrix.coeffRef(k, dof) = 0.0;
        }
        // Ставим 1 на диагонали
        systemMatrix.coeffRef(dof, dof) = 1.0;
        rightHandSide[dof] = 0.0;
    }
}