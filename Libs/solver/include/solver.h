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

    void reduceSystem(const Eigen::SparseMatrix<double>& fullK,
        const Eigen::VectorXd& fullF,
        const std::vector<int>& fixedDofs,
        Eigen::SparseMatrix<double>& reducedK,
        Eigen::VectorXd& reducedF,
        std::vector<int>& activeDofs) const;

    // Восстановление полного вектора перемещений
    void expandSolution(const Eigen::VectorXd& reducedU,
        const std::vector<int>& fixedDofs,
        const std::vector<int>& activeDofs,
        Eigen::VectorXd& fullU) const;

    void applyPrescribedDisplacements(Eigen::SparseMatrix<double>& K,
        Eigen::VectorXd& F,
        const std::vector<int>& prescribedDofs,
        const std::vector<double>& prescribedValues,
        Eigen::VectorXd& reactions) const;

private:
    struct GaussPoint {
        double xi;
        double eta;
        double weight;
    };

    std::vector<GaussPoint> generateGaussPoints(int order);




};