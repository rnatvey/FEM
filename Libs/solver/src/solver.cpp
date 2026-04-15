#include "solver.h"

#include <chrono>
#include <stdexcept>

Eigen::MatrixXd LinearSolver::computeGaussIntegral(
    const std::function<Eigen::MatrixXd(double, double)>& matFunc,
    int numGaussPoints)
{
    auto gaussPoints = generateGaussPoints(numGaussPoints);
    if (gaussPoints.empty()) {
        throw std::runtime_error("No Gauss points generated");
    }

    Eigen::MatrixXd testMat = matFunc(0.0, 0.0);
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
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
        Eigen::Lower,
        Eigen::IncompleteCholesky<double>> solver;

    constexpr double tolerance = 1.0e-10;
    constexpr int maxIterations = 2000;

    solver.setTolerance(tolerance);
    solver.setMaxIterations(maxIterations);

    const auto startTime = std::chrono::high_resolution_clock::now();
    solver.compute(systemMatrix);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }

    Eigen::VectorXd solution = solver.solve(rightHandSide);
    const auto endTime = std::chrono::high_resolution_clock::now();

    lastSolveStats_.iterations = solver.iterations();
    lastSolveStats_.maxIterations = maxIterations;
    lastSolveStats_.tolerance = tolerance;
    lastSolveStats_.estimatedError = solver.error();
    lastSolveStats_.solveTimeSeconds =
        std::chrono::duration<double>(endTime - startTime).count();
    lastSolveStats_.converged = solver.info() == Eigen::Success;

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
    }

    return solution;
}

std::vector<LinearSolver::GaussPoint> LinearSolver::generateGaussPoints(int order) {
    std::vector<GaussPoint> points;

    if (order == 2) {
        constexpr double pos = 0.577350269189626;
        constexpr double w = 1.0;

        points = {
            {-pos, -pos, w * w},
            {pos, -pos, w * w},
            {pos, pos, w * w},
            {-pos, pos, w * w}
        };
    }

    if (order == 3) {
        constexpr double p1 = -0.774596669241483;
        constexpr double p2 = 0.0;
        constexpr double p3 = 0.774596669241483;

        constexpr double w1 = 5.0 / 9.0;
        constexpr double w2 = 8.0 / 9.0;
        constexpr double w3 = 5.0 / 9.0;

        points = {
            {p1, p1, w1 * w1},
            {p1, p2, w1 * w2},
            {p1, p3, w1 * w3},
            {p2, p1, w2 * w1},
            {p2, p2, w2 * w2},
            {p2, p3, w2 * w3},
            {p3, p1, w3 * w1},
            {p3, p2, w3 * w2},
            {p3, p3, w3 * w3}
        };
    }

    return points;
}

void LinearSolver::applyBoundaryConditions(Eigen::SparseMatrix<double>& systemMatrix,
    Eigen::VectorXd& rightHandSide,
    const std::vector<int>& fixedDofs) const {
    for (int dof : fixedDofs) {
        for (int k = 0; k < systemMatrix.outerSize(); ++k) {
            systemMatrix.coeffRef(dof, k) = 0.0;
            systemMatrix.coeffRef(k, dof) = 0.0;
        }

        systemMatrix.coeffRef(dof, dof) = 0.0;
        rightHandSide[dof] = 0.0;
    }
}

void LinearSolver::reduceSystem(const Eigen::SparseMatrix<double>& fullK,
    const Eigen::VectorXd& fullF,
    const std::vector<int>& fixedDofs,
    Eigen::SparseMatrix<double>& reducedK,
    Eigen::VectorXd& reducedF,
    std::vector<int>& activeDofs) const {
    int totalDof = fullK.rows();

    if (totalDof == 0) {
        throw std::runtime_error("Empty matrix in reduceSystem");
    }

    if (totalDof != fullF.size()) {
        throw std::runtime_error("Matrix and vector size mismatch in reduceSystem");
    }

    std::vector<char> isFixed(totalDof, 0);
    for (int dof : fixedDofs) {
        if (dof < 0 || dof >= totalDof) {
            throw std::runtime_error("Invalid fixed DOF index: " + std::to_string(dof));
        }
        isFixed[dof] = 1;
    }

    activeDofs.clear();
    activeDofs.reserve(totalDof - static_cast<int>(fixedDofs.size()));
    for (int i = 0; i < totalDof; ++i) {
        if (!isFixed[i]) {
            activeDofs.push_back(i);
        }
    }

    int activeDofCount = static_cast<int>(activeDofs.size());
    if (activeDofCount == 0) {
        throw std::runtime_error("No active DOFs after applying boundary conditions");
    }

    std::vector<int> globalToActive(totalDof, -1);
    for (int i = 0; i < activeDofCount; ++i) {
        globalToActive[activeDofs[i]] = i;
    }

    reducedF.resize(activeDofCount);
    for (int i = 0; i < activeDofCount; ++i) {
        reducedF(i) = fullF[activeDofs[i]];
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(fullK.nonZeros());

    for (int k = 0; k < fullK.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(fullK, k); it; ++it) {
            int globalI = it.row();
            int globalJ = it.col();

            if (!isFixed[globalI] && !isFixed[globalJ]) {
                triplets.emplace_back(globalToActive[globalI], globalToActive[globalJ], it.value());
            }
        }
    }

    reducedK.resize(activeDofCount, activeDofCount);
    reducedK.setFromTriplets(triplets.begin(), triplets.end());
    reducedK.makeCompressed();
}

void LinearSolver::expandSolution(const Eigen::VectorXd& reducedU,
    const std::vector<int>& fixedDofs,
    const std::vector<int>& activeDofs,
    Eigen::VectorXd& fullU) const {
    (void)fixedDofs;

    fullU.setZero();
    for (size_t i = 0; i < activeDofs.size(); ++i) {
        fullU(activeDofs[i]) = reducedU(static_cast<Eigen::Index>(i));
    }
}

void LinearSolver::applyPrescribedDisplacements(Eigen::SparseMatrix<double>& K,
    Eigen::VectorXd& F,
    const std::vector<int>& prescribedDofs,
    const std::vector<double>& prescribedValues,
    Eigen::VectorXd& reactions) const {
    if (prescribedDofs.size() != prescribedValues.size()) {
        throw std::invalid_argument("Prescribed DOFs and values size mismatch");
    }

    int totalDof = K.rows();
    if (totalDof == 0) {
        throw std::runtime_error("Empty matrix in applyPrescribedDisplacements");
    }

    if (totalDof != F.size()) {
        throw std::runtime_error("Matrix and vector size mismatch in applyPrescribedDisplacements");
    }

    reactions.resize(totalDof);
    reactions.setZero();

    std::vector<char> isPrescribed(totalDof, 0);
    for (int prescribedDof : prescribedDofs) {
        if (prescribedDof < 0 || prescribedDof >= totalDof) {
            throw std::runtime_error("Invalid prescribed DOF index: " + std::to_string(prescribedDof));
        }
        isPrescribed[prescribedDof] = 1;
    }

    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        int prescribedDof = prescribedDofs[i];
        double prescribedValue = prescribedValues[i];

        for (int j = 0; j < totalDof; ++j) {
            if (!isPrescribed[j]) {
                F(j) -= K.coeff(j, prescribedDof) * prescribedValue;
            }
        }
    }

    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        int prescribedDof = prescribedDofs[i];
        double prescribedValue = prescribedValues[i];

        for (int j = 0; j < totalDof; ++j) {
            K.coeffRef(prescribedDof, j) = 0.0;
            K.coeffRef(j, prescribedDof) = 0.0;
        }

        K.coeffRef(prescribedDof, prescribedDof) = 1.0;
        F(prescribedDof) = prescribedValue;
    }
}
