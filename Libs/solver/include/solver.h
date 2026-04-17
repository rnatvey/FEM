#pragma once

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#ifndef EIGEN_USE_MKL
#define EIGEN_USE_MKL
#endif

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/PardisoSupport>
#include <Eigen/Sparse>
#include <functional>
#include <string>
#include <vector>

class LinearSolver {
public:
    LinearSolver() = default;
    ~LinearSolver() = default;

    struct SolveStats {
        int iterations = 0;
        int maxIterations = 0;
        double tolerance = 0.0;
        double estimatedError = 0.0;
        double residualNorm = 0.0;
        double relativeResidualNorm = 0.0;
        double solveTimeSeconds = 0.0;
        bool converged = false;
        bool usedDirectSolver = false;
        std::string backendName = "uninitialized";
    };

    Eigen::MatrixXd computeGaussIntegral(
        const std::function<Eigen::MatrixXd(double, double)>& matFunc,
        int numGaussPoints);

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

    void reduceSystem(const Eigen::SparseMatrix<double>& fullK,
        const Eigen::VectorXd& fullF,
        const std::vector<int>& constrainedDofs,
        const std::vector<int>& prescribedDofs,
        const std::vector<double>& prescribedValues,
        Eigen::SparseMatrix<double>& reducedK,
        Eigen::VectorXd& reducedF,
        std::vector<int>& activeDofs) const;

    void expandSolution(const Eigen::VectorXd& reducedU,
        const std::vector<int>& fixedDofs,
        const std::vector<int>& activeDofs,
        Eigen::VectorXd& fullU) const;

    void applyPrescribedDisplacements(Eigen::SparseMatrix<double>& K,
        Eigen::VectorXd& F,
        const std::vector<int>& prescribedDofs,
        const std::vector<double>& prescribedValues,
        Eigen::VectorXd& reactions) const;

    const SolveStats& getLastSolveStats() const { return lastSolveStats_; }

private:
    struct GaussPoint {
        double xi;
        double eta;
        double weight;
    };

    std::vector<GaussPoint> generateGaussPoints(int order);

    SolveStats lastSolveStats_;
};
