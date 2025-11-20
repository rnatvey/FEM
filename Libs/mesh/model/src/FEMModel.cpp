#include "FEMModel.h"
#include "planeisometric/Planeisoparametric.h"
#include <chrono>
#include <iostream>
#include <cmath>

FEModel::FEModel() {
    std::cout << "FEModel constructor started..." << std::endl;

    try {
        solver_ = std::make_unique<LinearSolver>();
        std::cout << "LinearSolver created successfully" << std::endl;

        // Инициализируем векторы нулевого размера
        displacements_ = Eigen::VectorXd::Zero(0);
        reactionForces_ = Eigen::VectorXd::Zero(0);
        contactForces_ = Eigen::VectorXd::Zero(0);

        std::cout << "FEModel constructor completed" << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "ERROR in FEModel constructor: " << e.what() << std::endl;
        throw; // Перебрасываем исключение дальше
    }
}

void FEModel::setAssembly(std::shared_ptr<Assembly> assembly) {
    assembly_ = assembly;
}

bool FEModel::solve() {
    if (!assembly_ || !validate()) {
        std::cerr << "FEModel: Invalid assembly or model configuration" << std::endl;
        return false;
    }

    auto startTime = std::chrono::high_resolution_clock::now();

    try {
        // 1. Сборка глобальной матрицы жесткости
        assembly_->assembleGlobalStiffnessMatrix(globalK_);

        // 2. Сборка вектора нагрузок
        //Eigen::VectorXd bodyForces = Eigen::VectorXd::Zero(2); // Объемные силы
        //assembly_->assembleGlobalForceVector(globalF_, bodyForces);

        // 3. Добавляем сосредоточенные силы
        assembly_->assembleConcentratedForces(globalF_);

        // 4. Применение граничных условий
        assembly_->applyBoundaryConditions(globalK_, globalF_);

        // 5. Решение системы
        if (!solveLinearSystem()) {
            std::cerr << "FEModel: Linear system solution failed" << std::endl;
            return false;
        }

        // 6. Расчет реакций
       calculateReactionForces();

        auto endTime = std::chrono::high_resolution_clock::now();
        solutionTime_ = std::chrono::duration<double>(endTime - startTime).count();

        std::cout << "FEModel: Solution completed in " << solutionTime_ << " seconds" << std::endl;
        return true;

    }
    catch (const std::exception& e) {
        std::cerr << "FEModel: Solution error - " << e.what() << std::endl;
        return false;
    }
}

bool FEModel::solveContact() {
    if (!assembly_ || !validate()) {
        std::cerr << "FEModel: Invalid assembly for contact solution" << std::endl;
        return false;
    }

    auto startTime = std::chrono::high_resolution_clock::now();

    try {
        return solveContactIterative();

    }
    catch (const std::exception& e) {
        std::cerr << "FEModel: Contact solution error - " << e.what() << std::endl;
        return false;
    }
}

bool FEModel::solveLinearSystem() {
    try {
        if (!solver_) {
            solver_ = std::make_unique<LinearSolver>();
        }

        // Решаем сокращенную систему
        Eigen::VectorXd reducedDisplacements = solver_->solveSystem(globalK_, globalF_);

        std::cout << "=== Solve Linear System ===" << std::endl;
        std::cout << "Reduced system: " << globalK_.rows() << "x" << globalK_.cols() << std::endl;
        std::cout << "Reduced displacements size: " << reducedDisplacements.size() << std::endl;

        // Восстанавливаем полный вектор перемещений
        int totalDof = assembly_->getTotalDofCount();
        displacements_.resize(totalDof);
        displacements_.setZero();

        const auto& mapping = assembly_->getDofMapping();

        std::cout << "Total DOF: " << totalDof << std::endl;
        std::cout << "Active DOF count: " << mapping.reducedToFull.size() << std::endl;
        std::cout << "Prescribed DOF count: " << mapping.prescribedDofs.size() << std::endl;

        // Заполняем активные DOF
        for (size_t i = 0; i < mapping.reducedToFull.size(); ++i) {
            if (i < reducedDisplacements.size()) {
                int fullDof = mapping.reducedToFull[i];
                if (fullDof >= 0 && fullDof < totalDof) {
                    displacements_(fullDof) = reducedDisplacements(i);
                }
            }
        }

        // Заполняем предписанные перемещения
        for (size_t i = 0; i < mapping.prescribedDofs.size(); ++i) {
            int prescribedDof = mapping.prescribedDofs[i];
            if (prescribedDof >= 0 && prescribedDof < totalDof && i < mapping.prescribedValues.size()) {
                displacements_(prescribedDof) = mapping.prescribedValues[i];
            }
        }

        std::cout << "Full displacements size: " << displacements_.size() << std::endl;
        std::cout << "Displacements norm: " << displacements_.norm() << std::endl;

        return true;

    }
    catch (const std::exception& e) {
        std::cerr << "Linear system solution failed: " << e.what() << std::endl;
        return false;
    }
}

bool FEModel::solveContactIterative() {
    std::cout << "FEModel: Starting iterative contact solution..." << std::endl;

    // Начальное решение без учета контакта
    if (!solve()) {
        return false;
    }

    Eigen::VectorXd prevDisplacements = displacements_;
    bool contactConverged = false;

    for (iterationCount_ = 1; iterationCount_ <= maxIterations_; ++iterationCount_) {
        std::cout << "Contact iteration " << iterationCount_ << std::endl;

        // 1. Определение активных контактных пар
        applyContactConditions();

        // 2. Пересборка системы с учетом контактных жесткостей
        assembly_->assembleGlobalStiffnessMatrix(globalK_);

        // Добавление контактных жесткостей
        // (здесь нужно найти все ContactPlaneElement и добавить их штрафные матрицы)

        Eigen::VectorXd bodyForces = Eigen::VectorXd::Zero(2);
        assembly_->assembleGlobalForceVector(globalF_, bodyForces);

        // Добавление контактных сил
        globalF_ += contactForces_;

        // 3. Применение граничных условий
        assembly_->applyBoundaryConditions(globalK_, globalF_);

        // 4. Решение системы
        if (!solveLinearSystem()) {
            return false;
        }

        // 5. Проверка сходимости
        double displacementNorm = (displacements_ - prevDisplacements).norm();
        double relativeError = displacementNorm / (displacements_.norm() + 1.0e-15);

        std::cout << "Iteration " << iterationCount_ << ": relative error = " << relativeError << std::endl;

        if (relativeError < tolerance_) {
            contactConverged = true;
            break;
        }

        prevDisplacements = displacements_;
    }

    if (contactConverged) {
        calculateReactionForces();
        auto endTime = std::chrono::high_resolution_clock::now();
       // solutionTime_ = std::chrono::duration<double>(endTime - startTime).count();

        std::cout << "FEModel: Contact solution converged in " << iterationCount_
            << " iterations, time = " << solutionTime_ << " seconds" << std::endl;
        return true;
    }
    else {
        std::cerr << "FEModel: Contact solution did not converge after "
            << maxIterations_ << " iterations" << std::endl;
        return false;
    }
}

void FEModel::calculateReactionForces() {
    int totalDof = assembly_->getTotalDofCount();
    reactionForces_.resize(totalDof);
    reactionForces_.setZero();

    Eigen::SparseMatrix<double> fullK;
    Eigen::VectorXd fullF;
    assembly_->assembleGlobalStiffnessMatrix(fullK);

    // Собираем вектор сил БЕЗ граничных условий
    Eigen::VectorXd bodyForces = Eigen::VectorXd::Zero(2);
    assembly_->assembleGlobalForceVector(fullF, bodyForces);
    assembly_->assembleConcentratedForces(fullF);

    std::cout << "Full K size: " << fullK.rows() << "x" << fullK.cols() << std::endl;
    std::cout << "Full F size: " << fullF.size() << std::endl;

    //  Вычисляем реакции: R = K * u - F
    if (fullK.rows() == displacements_.size() && fullK.cols() == displacements_.size()) {
        reactionForces_ = fullK * displacements_ - fullF;
        std::cout << "Reactions calculated successfully" << std::endl;
    }
    else {
        std::cout << "ERROR: Dimension mismatch for reaction calculation" << std::endl;
        std::cout << "K: " << fullK.rows() << "x" << fullK.cols()
            << ", u: " << displacements_.size()
            << ", F: " << fullF.size() << std::endl;
    }

    //  Выводим реакции в закрепленных узлах
    std::cout << "Reaction forces:" << std::endl;
    auto nodes = assembly_->getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        int nodeId = nodes[i]->getId();
        int dofX = i * 2;
        int dofY = i * 2 + 1;

        if (dofX < reactionForces_.size() && dofY < reactionForces_.size()) {
            std::cout << "Node " << nodeId << ": Rx = " << reactionForces_(dofX)
                << ", Ry = " << reactionForces_(dofY) << std::endl;
        }
    }

}

void FEModel::applyContactConditions() {
    // Временная реализация - будет расширена
    contactForces_.resize(globalF_.size());
    contactForces_.setZero();

    // Здесь нужно:
    // 1. Найти все ContactPlaneElement
    // 2. Для каждого вычислить штрафные силы
    // 3. Добавить в contactForces_
}

Eigen::Vector3d FEModel::getElementStress(int elementId, double xi, double eta) const {
    auto element = assembly_->getElement(elementId);
    auto material = assembly_->getMaterial(element->getMaterialId());

    if (element && material) {
        // Используем ПОЛНЫЕ индексы DOF
        auto dofIndices = assembly_->getElementDofIndices(elementId);
        Eigen::VectorXd elementDisplacements = Eigen::VectorXd::Zero(dofIndices.size());

        std::cout << "=== getElementStress with FULL indices ===" << std::endl;
        std::cout << "Element " << elementId << " DOF indices: ";
        for (int idx : dofIndices) {
            std::cout << idx << " ";
        }
        std::cout << std::endl;

        std::cout << "Global displacements size: " << displacements_.size() << std::endl;

        // Заполняем перемещения элемента из полного вектора перемещений
        for (size_t i = 0; i < dofIndices.size(); ++i) {
            if (dofIndices[i] >= 0 && dofIndices[i] < displacements_.size()) {
                elementDisplacements(i) = displacements_(dofIndices[i]);
                std::cout << "DOF " << dofIndices[i] << " = " << displacements_(dofIndices[i]) << std::endl;
            }
            else {
                std::cout << "WARNING: Invalid DOF index " << dofIndices[i] << std::endl;
            }
        }

        std::cout << "Element displacements norm: " << elementDisplacements.norm() << std::endl;

        Eigen::Vector3d stress = element->computeStress(xi, eta, elementDisplacements,
            assembly_->getNodes(), material);
        std::cout << "Computed stress: " << stress.transpose() << std::endl;

        return stress;
    }

    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d FEModel::getElementStrain(int elementId, double xi, double eta) const {
    auto element = assembly_->getElement(elementId);
    auto material = assembly_->getMaterial(element->getMaterialId());

    if (element && material) {
        auto dofIndices = assembly_->getElementDofIndices(elementId);
        Eigen::VectorXd elementDisplacements(dofIndices.size());

        for (size_t i = 0; i < dofIndices.size(); ++i) {
            if (dofIndices[i] >= 0) {
                elementDisplacements(i) = displacements_(dofIndices[i]);
            }
        }

        return element->computeStrain(xi, eta, elementDisplacements,
            assembly_->getNodes(), material);
    }

    return Eigen::Vector3d::Zero();
}

bool FEModel::validate() const {
    if (!assembly_) {
        std::cerr << "FEModel: No assembly set" << std::endl;
        return false;
    }

    if (!assembly_->validate()) {
        std::cerr << "FEModel: Assembly validation failed" << std::endl;
        return false;
    }

    return true;
}

std::vector<Eigen::Vector3d> FEModel::getNodalStresses() const {
    calculateNodalAverages();
    return nodalStresses_;
}

std::vector<Eigen::Vector2d> FEModel::getNodalDisplacements() const {
    calculateNodalAverages();
    return nodalDisplacements_;
}

std::vector<Eigen::Vector3d> FEModel::getNodalStrains() const {
    calculateNodalAverages();
    return nodalStrains_;
}

Eigen::Vector3d FEModel::getNodeStress(int nodeId) const {
    calculateNodalAverages();
    auto nodes = assembly_->getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i]->getId() == nodeId) {
            return nodalStresses_[i];
        }
    }
    return Eigen::Vector3d::Zero();
}

Eigen::Vector2d FEModel::getNodeDisplacement(int nodeId) const {
    calculateNodalAverages();
    auto nodes = assembly_->getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i]->getId() == nodeId) {
            return nodalDisplacements_[i];
        }
    }
    return Eigen::Vector2d::Zero();
}

Eigen::Vector3d FEModel::getNodeStrain(int nodeId) const {
    calculateNodalAverages();
    auto nodes = assembly_->getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i]->getId() == nodeId) {
            return nodalStrains_[i];
        }
    }
    return Eigen::Vector3d::Zero();
}

void FEModel::calculateNodalAverages() const {
    if (nodalDataCalculated_ || !assembly_) {
        return;
    }

    auto nodes = assembly_->getNodes();
    int nodeCount = nodes.size();

    // Инициализация
    nodalStresses_.resize(nodeCount, Eigen::Vector3d::Zero());
    nodalDisplacements_.resize(nodeCount, Eigen::Vector2d::Zero());
    nodalStrains_.resize(nodeCount, Eigen::Vector3d::Zero());

    std::vector<int> stressCount(nodeCount, 0);
    std::vector<int> strainCount(nodeCount, 0);

    // Собираем перемещения из глобального вектора
    for (int i = 0; i < nodeCount; ++i) {
        int nodeId = nodes[i]->getId();
        int dofX = assembly_->getGlobalDofIndex(nodeId, 0);
        int dofY = assembly_->getGlobalDofIndex(nodeId, 1);

        if (dofX >= 0 && dofX < displacements_.size()) {
            nodalDisplacements_[i].x() = displacements_[dofX];
        }
        if (dofY >= 0 && dofY < displacements_.size()) {
            nodalDisplacements_[i].y() = displacements_[dofY];
        }
    }

    // Усреднение напряжений и деформаций по элементам
    auto elements = assembly_->getElements();

    for (const auto& element : elements) {
        auto material = assembly_->getMaterial(element->getMaterialId());
        if (!material) continue;

        // Получаем перемещения элемента
        auto dofIndices = assembly_->getElementDofIndices(element->getId());
        Eigen::VectorXd elementDisplacements = Eigen::VectorXd::Zero(dofIndices.size());

        for (size_t i = 0; i < dofIndices.size(); ++i) {
            if (dofIndices[i] >= 0 && dofIndices[i] < displacements_.size()) {
                elementDisplacements(i) = displacements_[dofIndices[i]];
            }
        }

        // Точки интегрирования для усреднения (углы элемента)
        std::vector<std::pair<double, double>> integrationPoints = {
            {-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}  // Приближенно к узлам
        };

        for (size_t pt = 0; pt < integrationPoints.size(); ++pt) {
            double xi = integrationPoints[pt].first;
            double eta = integrationPoints[pt].second;

            // Вычисляем напряжения и деформации в точке
            Eigen::Vector3d stress = element->computeStress(xi, eta, elementDisplacements,
                nodes, material);
            Eigen::Vector3d strain = element->computeStrain(xi, eta, elementDisplacements,
                nodes, material);

            // Находим соответствующий узел элемента
            auto nodeIds = element->getNodeIds();
            if (pt < nodeIds.size()) {
                int nodeId = nodeIds[pt];

                // Находим индекс узла в глобальном списке
                for (int i = 0; i < nodeCount; ++i) {
                    if (nodes[i]->getId() == nodeId) {
                        nodalStresses_[i] += stress;
                        nodalStrains_[i] += strain;
                        stressCount[i]++;
                        strainCount[i]++;
                        break;
                    }
                }
            }
        }
    }

    // Усреднение по количеству элементов, содержащих узел
    for (int i = 0; i < nodeCount; ++i) {
        if (stressCount[i] > 0) {
            nodalStresses_[i] /= stressCount[i];
        }
        if (strainCount[i] > 0) {
            nodalStrains_[i] /= strainCount[i];
        }
    }

    nodalDataCalculated_ = true;
}