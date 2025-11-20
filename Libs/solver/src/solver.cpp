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

void LinearSolver::reduceSystem(const Eigen::SparseMatrix<double>& fullK,
    const Eigen::VectorXd& fullF,
    const std::vector<int>& fixedDofs,
    Eigen::SparseMatrix<double>& reducedK,
    Eigen::VectorXd& reducedF,
    std::vector<int>& activeDofs) const {
    

    int totalDof = fullK.rows();

    // Проверки безопасности
    if (totalDof == 0) {
        throw std::runtime_error("Empty matrix in reduceSystem");
    }

    if (totalDof != fullF.size()) {
        throw std::runtime_error("Matrix and vector size mismatch in reduceSystem");
    }

    // Определяем активные степени свободы
    activeDofs.clear();
    for (int i = 0; i < totalDof; ++i) {
        if (std::find(fixedDofs.begin(), fixedDofs.end(), i) == fixedDofs.end()) {
            activeDofs.push_back(i);
        }
    }

    int activeDofCount = activeDofs.size();

    // Проверка на пустую систему
    if (activeDofCount == 0) {
        throw std::runtime_error("No active DOFs after applying boundary conditions");
    }

    reducedK.resize(activeDofCount, activeDofCount);
    reducedF.resize(activeDofCount);

    // Строим матрицу перестановок с проверками границ
    std::vector<Eigen::Triplet<double>> triplets;

    for (int i = 0; i < activeDofCount; ++i) {
        int global_i = activeDofs[i];

        // Проверка границ
        if (global_i < 0 || global_i >= totalDof) {
            throw std::runtime_error("Invalid global DOF index: " + std::to_string(global_i));
        }

        reducedF(i) = fullF(global_i);

        for (int j = 0; j < activeDofCount; ++j) {
            int global_j = activeDofs[j];

            if (global_j < 0 || global_j >= totalDof) {
                throw std::runtime_error("Invalid global DOF index: " + std::to_string(global_j));
            }

            double value = fullK.coeff(global_i, global_j);
            triplets.emplace_back(i, j, value);
        }
    }

    reducedK.setFromTriplets(triplets.begin(), triplets.end());
}

void LinearSolver::expandSolution(const Eigen::VectorXd& reducedU,
    const std::vector<int>& fixedDofs,
    const std::vector<int>& activeDofs,
    Eigen::VectorXd& fullU) const {
    int totalDof = fullU.size();
    fullU.setZero();

    // Заполняем активные степени свободы
    for (size_t i = 0; i < activeDofs.size(); ++i) {
        fullU(activeDofs[i]) = reducedU(i);
    }

    // Закрепленные DOF остаются нулевыми (или предписанными значениями)
    // Предписанные перемещения будут установлены отдельно
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

    // Проверки безопасности
    if (totalDof == 0) {
        throw std::runtime_error("Empty matrix in applyPrescribedDisplacements");
    }

    if (totalDof != F.size()) {
        throw std::runtime_error("Matrix and vector size mismatch in applyPrescribedDisplacements");
    }

    reactions.resize(totalDof);
    reactions.setZero();

    // 1. Модифицируем правую часть: F_modified = F - K?? u?
    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        int prescribedDof = prescribedDofs[i];
        double prescribedValue = prescribedValues[i];

        // Проверка границ
        if (prescribedDof < 0 || prescribedDof >= totalDof) {
            throw std::runtime_error("Invalid prescribed DOF index: " + std::to_string(prescribedDof));
        }

        for (int j = 0; j < totalDof; ++j) {
            // Если j не предписанный DOF, модифицируем F[j]
            if (std::find(prescribedDofs.begin(), prescribedDofs.end(), j) == prescribedDofs.end()) {
                // Проверяем, что индексы в пределах
                if (j >= 0 && j < totalDof && prescribedDof >= 0 && prescribedDof < totalDof) {
                    F(j) -= K.coeff(j, prescribedDof) * prescribedValue;
                }
            }
        }
    }

    // 2. Обнуляем строки и столбцы для предписанных DOF, ставим 1 на диагонали
    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        int prescribedDof = prescribedDofs[i];
        double prescribedValue = prescribedValues[i];

        // Проверка границ
        if (prescribedDof < 0 || prescribedDof >= totalDof) {
            throw std::runtime_error("Invalid prescribed DOF index: " + std::to_string(prescribedDof));
        }

        // Обнуляем строку и столбец
        for (int j = 0; j < totalDof; ++j) {
            if (j >= 0 && j < totalDof) {
                K.coeffRef(prescribedDof, j) = 0.0;
                K.coeffRef(j, prescribedDof) = 0.0;
            }
        }

        // Ставим 1 на диагонали
        if (prescribedDof >= 0 && prescribedDof < totalDof) {
            K.coeffRef(prescribedDof, prescribedDof) = 1.0;
            F(prescribedDof) = prescribedValue;
        }
    }
}