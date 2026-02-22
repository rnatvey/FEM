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

//Eigen::VectorXd LinearSolver::solveSystem(
//    const Eigen::SparseMatrix<double>& systemMatrix,
//    const Eigen::VectorXd& rightHandSide)
//{
//   // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
//    Eigen::PardisoLLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
//    solver.compute(systemMatrix);
//    if (solver.info() != Eigen::Success) {
//        throw std::runtime_error("Matrix decomposition failed");
//    }
//    return solver.solve(rightHandSide);
//}


Eigen::VectorXd LinearSolver::solveSystem(
    const Eigen::SparseMatrix<double>& systemMatrix,
    const Eigen::VectorXd& rightHandSide)
{
    // Итеративный решатель ConjugateGradient с предобуславливателем IncompleteCholesky
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
        Eigen::Lower,
        Eigen::IncompleteCholesky<double>> solver;

    // Настройка параметров (опционально)
    solver.setTolerance(1e-10);        // Желаемая точность
    solver.setMaxIterations(2000);    // Максимум итераций

    // Вычисление
    solver.compute(systemMatrix);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }

    // Проверка решения
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
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
            {pos,  -pos, w * w},
            { pos, pos, w * w},
            { -pos,  pos, w * w}
        };
    }

    if (order == 3) {

        
       const auto p1 = -0.774596669241483;  // -√(3/5)
       const auto p2 = 0.0;                 // 0
       const auto p3 = 0.774596669241483;    // √(3/5)
        
        
        const auto w1 = 5.0 / 9.0;            // 5/9
        const auto w2 = 8.0 / 9.0;            // 8/9  
        const auto w3 = 5.0 / 9.0;             // 5/9
        
  

        points = {
            {p1, p1, w1 * w1},
            {p1,  p2, w1 * w2},
            { p1, p3, w1 * w3},
            { p2, p1, w2 * w1},
            { p2,  p2, w2 * w2},
            { p2,  p3, w2 * w3},
            { p3, p1, w3 * w1},
            { p3,  p2, w3 * w2},
            { p3,  p3, w3 * w3}

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

    // Проверки безопасности
    if (totalDof == 0) {
        throw std::runtime_error("Empty matrix in reduceSystem");
    }

    if (totalDof != fullF.size()) {
        throw std::runtime_error("Matrix and vector size mismatch in reduceSystem");
    }

    // Оптимизация 1: Создаем множество закрепленных DOF для быстрого поиска O(1)
    std::vector<char> isFixed(totalDof, 0);
    for (int dof : fixedDofs) {
        if (dof < 0 || dof >= totalDof) {
            throw std::runtime_error("Invalid fixed DOF index: " + std::to_string(dof));
        }
        isFixed[dof] = 1;
    }

    // Определяем активные степени свободы (один проход O(n))
    activeDofs.clear();
    activeDofs.reserve(totalDof - fixedDofs.size());
    for (int i = 0; i < totalDof; ++i) {
        if (!isFixed[i]) {
            activeDofs.push_back(i);
        }
    }

    int activeDofCount = activeDofs.size();

    // Проверка на пустую систему
    if (activeDofCount == 0) {
        throw std::runtime_error("No active DOFs after applying boundary conditions");
    }

    // Оптимизация 2: Создаем обратное отображение global -> active
    std::vector<int> globalToActive(totalDof, -1);
    for (int i = 0; i < activeDofCount; ++i) {
        globalToActive[activeDofs[i]] = i;
    }

    // Оптимизация 3: Заполняем вектор сил
    reducedF.resize(activeDofCount);
    for (int i = 0; i < activeDofCount; ++i) {
        reducedF(i) = fullF[activeDofs[i]];
    }

    // Оптимизация 4: Проходим только по ненулевым элементам разреженной матрицы
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(fullK.nonZeros()); // Резервируем память

    // Проходим по всей матрице, но берем только элементы, где оба индекса активны
    for (int k = 0; k < fullK.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(fullK, k); it; ++it) {
            int global_i = it.row();
            int global_j = it.col();

            // Проверяем, активны ли оба DOF
            if (!isFixed[global_i] && !isFixed[global_j]) {
                int active_i = globalToActive[global_i];
                int active_j = globalToActive[global_j];
                triplets.emplace_back(active_i, active_j, it.value());
            }
        }
    }

    // Создаем редуцированную матрицу
    reducedK.resize(activeDofCount, activeDofCount);
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

    
    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        int prescribedDof = prescribedDofs[i];
        double prescribedValue = prescribedValues[i];

        // Проверка границ
        if (prescribedDof < 0 || prescribedDof >= totalDof) {
            throw std::runtime_error("Invalid prescribed DOF index: " + std::to_string(prescribedDof));
        }

        for (int j = 0; j < totalDof; ++j) {
            
            if (std::find(prescribedDofs.begin(), prescribedDofs.end(), j) == prescribedDofs.end()) {
                
                if (j >= 0 && j < totalDof && prescribedDof >= 0 && prescribedDof < totalDof) {
                    F(j) -= K.coeff(j, prescribedDof) * prescribedValue;
                }
            }
        }
    }

    
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