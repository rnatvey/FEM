#include "solver.h"

#include <algorithm>
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
    constexpr double tolerance = 1.0e-10;
    constexpr int maxIterations = 2000;
    const auto startTime = std::chrono::high_resolution_clock::now();

    auto computeResidualNorm = [&systemMatrix, &rightHandSide](const Eigen::VectorXd& solution) {
        return (systemMatrix * solution - rightHandSide).norm();
    };

    auto finalizeStats = [this, &startTime, &rightHandSide](int iterations,
        int requestedMaxIterations,
        double requestedTolerance,
        double estimatedError,
        double residualNorm,
        bool converged,
        bool usedDirectSolver,
        const char* backendName) {
            const auto endTime = std::chrono::high_resolution_clock::now();
            lastSolveStats_.iterations = iterations;
            lastSolveStats_.maxIterations = requestedMaxIterations;
            lastSolveStats_.tolerance = requestedTolerance;
            lastSolveStats_.estimatedError = estimatedError;
            lastSolveStats_.residualNorm = residualNorm;
            const double rhsNorm = rightHandSide.norm();
            lastSolveStats_.relativeResidualNorm =
                residualNorm / std::max(rhsNorm, 1.0e-30);
            lastSolveStats_.solveTimeSeconds =
                std::chrono::duration<double>(endTime - startTime).count();
            lastSolveStats_.converged = converged;
            lastSolveStats_.usedDirectSolver = usedDirectSolver;
            lastSolveStats_.backendName = backendName;
        };

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
        Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<double>> iterativeSolver;

    iterativeSolver.setTolerance(tolerance);
    iterativeSolver.setMaxIterations(maxIterations);
    iterativeSolver.compute(systemMatrix);
    if (iterativeSolver.info() == Eigen::Success) {
        Eigen::VectorXd solution = iterativeSolver.solve(rightHandSide);
        if (iterativeSolver.info() == Eigen::Success) {
            const double residualNorm = computeResidualNorm(solution);
            finalizeStats(iterativeSolver.iterations(),
                maxIterations,
                tolerance,
                iterativeSolver.error(),
                residualNorm,
                true,
                false,
                "ConjugateGradient");
            return solution;
        }
    }

    // Fallback to a direct sparse solver for poorly conditioned systems
    // such as near-incompressible materials and penalty contact steps.
#ifdef SOLVER_USE_MKL
    Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> directSolver;
#else
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> directSolver;
#endif

    directSolver.compute(systemMatrix);
    if (directSolver.info() != Eigen::Success) {
        finalizeStats(iterativeSolver.iterations(),
            maxIterations,
            tolerance,
            iterativeSolver.error(),
            0.0,
            false,
            true,
#ifdef SOLVER_USE_MKL
            "PardisoLDLT");
#else
            "SimplicialLDLT");
#endif
        throw std::runtime_error("Matrix decomposition failed");
    }

    Eigen::VectorXd solution = directSolver.solve(rightHandSide);
    const bool converged = directSolver.info() == Eigen::Success;
    const double residualNorm = converged ? computeResidualNorm(solution) : 0.0;
    finalizeStats(iterativeSolver.iterations(),
        maxIterations,
        tolerance,
        iterativeSolver.error(),
        residualNorm,
        converged,
        true,
#ifdef SOLVER_USE_MKL
        "PardisoLDLT");
#else
        "SimplicialLDLT");
#endif

    if (!converged) {
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
    reduceSystem(fullK,
        fullF,
        fixedDofs,
        {},
        {},
        reducedK,
        reducedF,
        activeDofs);
}

void LinearSolver::reduceSystem(const Eigen::SparseMatrix<double>& fullK,
    const Eigen::VectorXd& fullF,
    const std::vector<int>& constrainedDofs,
    const std::vector<int>& prescribedDofs,
    const std::vector<double>& prescribedValues,
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

    if (prescribedDofs.size() != prescribedValues.size()) {
        throw std::invalid_argument("Prescribed DOFs and values size mismatch in reduceSystem");
    }

    std::vector<char> isConstrained(totalDof, 0);
    for (int dof : constrainedDofs) {
        if (dof < 0 || dof >= totalDof) {
            throw std::runtime_error("Invalid constrained DOF index: " + std::to_string(dof));
        }
        isConstrained[dof] = 1;
    }

    std::vector<char> isPrescribed(totalDof, 0);
    std::vector<double> prescribedValueByDof(totalDof, 0.0);
    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        int dof = prescribedDofs[i];
        if (dof < 0 || dof >= totalDof) {
            throw std::runtime_error("Invalid prescribed DOF index: " + std::to_string(dof));
        }
        isConstrained[dof] = 1;
        isPrescribed[dof] = 1;
        prescribedValueByDof[dof] = prescribedValues[i];
    }

    activeDofs.clear();
    activeDofs.reserve(totalDof - static_cast<int>(constrainedDofs.size()));
    for (int i = 0; i < totalDof; ++i) {
        if (!isConstrained[i]) {
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

            if (!isConstrained[globalI] && !isConstrained[globalJ]) {
                triplets.emplace_back(globalToActive[globalI], globalToActive[globalJ], it.value());
            }
            else if (!isConstrained[globalI] && isPrescribed[globalJ]) {
                reducedF(globalToActive[globalI]) -= it.value() * prescribedValueByDof[globalJ];
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
