#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <vector>

class LinearSolver {
public:
    LinearSolver() = default;
    ~LinearSolver() = default;

    // Вычисление интеграла методом Гаусса
    Eigen::MatrixXd computeGaussIntegral(
        const std::function<Eigen::MatrixXd(double, double)>& matFunc,
        int numGaussPoints);

    // Решение СЛАУ методом LDLT разложения
    Eigen::VectorXd solveSystem(
        const Eigen::SparseMatrix<double>& systemMatrix,
        const Eigen::VectorXd& rightHandSide);

    void applyBoundaryConditions(Eigen::SparseMatrix<double>& systemMatrix,
        Eigen::VectorXd& rightHandSide,
        const std::vector<int>& fixedDofs) const;

private:
    struct GaussPoint {
        double xi;
        double eta;
        double weight;
    };

    std::vector<GaussPoint> generateGaussPoints(int order);
};