#include "assembly.h"
#include "solver.h"
#include <stdexcept>
#include <iostream>


Assembly::Assembly() {
}

void Assembly::addNode(std::shared_ptr<Node> node) {
    if (!node) {
        throw std::invalid_argument("Cannot add null node");
    }

    // Проверяем уникальность ID
    for (const auto& existingNode : nodes_) {
        if (existingNode->getId() == node->getId()) {
            throw std::invalid_argument("Node with ID " + std::to_string(node->getId()) + " already exists");
        }
    }

    nodes_.push_back(node);
    buildNodeIndexMap(); // Обновляем карту индексов
}

void Assembly::addNodes(const std::vector<std::shared_ptr<Node>>& nodes) {
    for (const auto& node : nodes) {
        addNode(node);
    }
}

void Assembly::addLineofNodes(size_t startId, vector2 start_point,  vector2 end_point, int num_of_nodes ) 
{
    vector2 direction = end_point - start_point;
    const double step = Norm(direction)/num_of_nodes;
    direction.UniVector();
    
    for (int i = 0; i < num_of_nodes;i++)
    {
        auto current_position = start_point + (direction * step * i);
        addNode(std::make_shared<Node>(startId + i, current_position.X, current_position.Y));
    }
    return;
}

std::shared_ptr<Node> Assembly::getNode(int id) const {
    auto it = nodeIdToIndex_.find(id);
    if (it != nodeIdToIndex_.end()) {
        return nodes_[it->second];
    }
    return nullptr;
}

void Assembly::addMaterial(std::shared_ptr<Material> material) {
    if (!material) {
        throw std::invalid_argument("Cannot add null material");
    }

    int materialId = material->getId();
    if (materials_.find(materialId) != materials_.end()) {
        throw std::invalid_argument("Material with ID " + std::to_string(materialId) + " already exists");
    }

    materials_[materialId] = material;
}

std::shared_ptr<Material> Assembly::getMaterial(int id) const {
    auto it = materials_.find(id);
    if (it != materials_.end()) {
        return it->second;
    }
    return nullptr;
}

void Assembly::addElement(std::shared_ptr<BaseElement> element) {
    if (!element) {
        throw std::invalid_argument("Cannot add null element");
    }

    // Проверяем уникальность ID
    for (const auto& existingElement : elements_) {
        if (existingElement->getId() == element->getId()) {
            throw std::invalid_argument("Element with ID " + std::to_string(element->getId()) + " already exists");
        }
    }

    // Проверяем существование узлов
    for (int nodeId : element->getNodeIds()) {
        if (!getNode(nodeId)) {
            throw std::invalid_argument("Node with ID " + std::to_string(nodeId) + " not found for element " + std::to_string(element->getId()));
        }
    }

    // Проверяем существование материала
    if (!getMaterial(element->getMaterialId())) {
        throw std::invalid_argument("Material with ID " + std::to_string(element->getMaterialId()) + " not found for element " + std::to_string(element->getId()));
    }

    elements_.push_back(element);
}

void Assembly::addElements(const std::vector<std::shared_ptr<BaseElement>>& elements) {
    for (const auto& element : elements) {
        addElement(element);
    }
}

std::shared_ptr<BaseElement> Assembly::getElement(int id) const {
    for (const auto& element : elements_) {
        if (element->getId() == id) {
            return element;
        }
    }
    return nullptr;
}

void Assembly::assembleGlobalStiffnessMatrix(Eigen::SparseMatrix<double>& globalK) const {
    int totalDof = getTotalDofCount();
    globalK.resize(totalDof, totalDof);
    globalK.setZero();

    std::vector<Eigen::Triplet<double>> triplets;

    for (const auto& element : elements_) {
        auto material = getMaterial(element->getMaterialId());
        if (!material) {
            continue;
        }

        std::vector<std::shared_ptr<Node>> elementNodes;
        for (int nodeId : element->getNodeIds()) {
            auto node = getNode(nodeId);
            if (node) elementNodes.push_back(node);
        }

        // Вычисляем матрицу жесткости элемента
        Eigen::MatrixXd ke = element->computeStiffnessMatrix(elementNodes, material);

        // Используем ПОЛНЫЕ индексы DOF (до граничных условий)
        std::vector<int> dofIndices = getElementFullDofIndices(element->getId());

        // Добавляем в глобальную матрицу
        for (int i = 0; i < dofIndices.size(); ++i) {
            for (int j = 0; j < dofIndices.size(); ++j) {
                if (dofIndices[i] >= 0 && dofIndices[j] >= 0) {
                    triplets.emplace_back(dofIndices[i], dofIndices[j], ke(i, j));
                }
            }
        }
    }

    globalK.setFromTriplets(triplets.begin(), triplets.end());
}

void Assembly::assembleGlobalForceVector(Eigen::VectorXd& globalF, const Eigen::VectorXd& bodyForces) const {
    int totalDof = getTotalDofCount();
    globalF.resize(totalDof);
    globalF.setZero();

    for (const auto& element : elements_) {
        auto material = getMaterial(element->getMaterialId());
        if (!material) continue;

        std::vector<std::shared_ptr<Node>> elementNodes;
        for (int nodeId : element->getNodeIds()) {
            auto node = getNode(nodeId);
            if (node) elementNodes.push_back(node);
        }

        Eigen::VectorXd fe = element->computeEquivalentNodalForces(bodyForces, elementNodes, material);

        // Используем ПОЛНЫЕ индексы
        auto dofIndices = getElementFullDofIndices(element->getId());

        for (int i = 0; i < dofIndices.size(); ++i) {
            if (dofIndices[i] >= 0) {
                globalF(dofIndices[i]) += fe(i);
            }
        }
    }
}

void Assembly::addFixedNode(int nodeId, bool fixX, bool fixY) {
    BoundaryCondition bc;
    bc.nodeId = nodeId;
    bc.fixX = fixX;
    bc.fixY = fixY;
    bc.hasPrescribedDisplacement = false;
    bc.prescribedDx = 0.0;
    bc.prescribedDy = 0.0;
    boundaryConditions_.push_back(bc);
}

void Assembly::addPrescribedDisplacement(int nodeId, double dx, double dy) {
    BoundaryCondition bc;
    bc.nodeId = nodeId;
    bc.fixX = true;  // Для предписанных перемещений фиксируем обе степени свободы
    bc.fixY = true;
    bc.hasPrescribedDisplacement = true;
    bc.prescribedDx = dx;
    bc.prescribedDy = dy;
    boundaryConditions_.push_back(bc);
}

void Assembly::addPrescribedDisplacementX(int nodeId, double dx) {
    BoundaryCondition bc;
    bc.nodeId = nodeId;
    bc.fixX = true;  // Для предписанных перемещений фиксируем обе степени свободы
    bc.hasPrescribedDisplacement = true;
    bc.prescribedDx = dx;
    boundaryConditions_.push_back(bc);
}

void Assembly::addPrescribedDisplacementY(int nodeId,  double dy) {
    BoundaryCondition bc;
    bc.nodeId = nodeId;
    bc.fixY = true;
    bc.hasPrescribedDisplacement = true;
    bc.prescribedDy = dy;
    boundaryConditions_.push_back(bc);
}

void Assembly::applyBoundaryConditions(Eigen::SparseMatrix<double>& globalK,
    Eigen::VectorXd& globalF) const {
    LinearSolver solver;

    // Сбрасываем mapping
    dofMapping_.fullToReduced.clear();
    dofMapping_.reducedToFull.clear();
    dofMapping_.prescribedDofs.clear();
    dofMapping_.prescribedValues.clear();

    int totalDof = globalK.rows();
    dofMapping_.fullToReduced.resize(totalDof, -1); // Инициализируем -1

    std::vector<int> fixedDofs;
    std::vector<int> prescribedDofs;
    std::vector<double> prescribedValues;

    // Собираем граничные условия
    for (const auto& bc : boundaryConditions_) {
        auto node = getNode(bc.nodeId);
        if (!node) continue;

        if (bc.fixX) {
            int dofX = getGlobalDofIndex(bc.nodeId, 0);
            fixedDofs.push_back(dofX);

            if (bc.hasPrescribedDisplacement) {
                prescribedDofs.push_back(dofX);
                prescribedValues.push_back(bc.prescribedDx);
            }
        }

        if (bc.fixY) {
            int dofY = getGlobalDofIndex(bc.nodeId, 1);
            fixedDofs.push_back(dofY);

            if (bc.hasPrescribedDisplacement) {
                prescribedDofs.push_back(dofY);
                prescribedValues.push_back(bc.prescribedDy);
            }
        }
    }

    // Строим mapping для активных DOF
    int reducedIndex = 0;
    for (int i = 0; i < totalDof; ++i) {
        // Если DOF не закреплен и не имеет предписанного перемещения
        if (std::find(fixedDofs.begin(), fixedDofs.end(), i) == fixedDofs.end() &&
            std::find(prescribedDofs.begin(), prescribedDofs.end(), i) == prescribedDofs.end()) {
            dofMapping_.fullToReduced[i] = reducedIndex;
            dofMapping_.reducedToFull.push_back(i);
            reducedIndex++;
        }
    }

    // Сохраняем предписанные перемещения
    dofMapping_.prescribedDofs = prescribedDofs;
    dofMapping_.prescribedValues = prescribedValues;

    // Исключаем закрепленные DOF из системы
    if (!fixedDofs.empty()) {
        Eigen::SparseMatrix<double> reducedK;
        Eigen::VectorXd reducedF;
        std::vector<int> activeDofs;

        solver.reduceSystem(globalK, globalF, fixedDofs, reducedK, reducedF, activeDofs);
        globalK = reducedK;
        globalF = reducedF;
    }

    // Применяем предписанные перемещения
    if (!prescribedDofs.empty()) {
        Eigen::VectorXd reactions;
        solver.applyPrescribedDisplacements(globalK, globalF, prescribedDofs, prescribedValues, reactions);
    }
}

std::vector<int> Assembly::getElementDofIndices(int elementId) const {
    auto element = getElement(elementId);
    if (!element) {
        throw std::invalid_argument("Element with ID " + std::to_string(elementId) + " not found");
    }

    std::vector<int> fullDofIndices;

    // Получаем полные индексы DOF
    for (int nodeId : element->getNodeIds()) {
        auto node = getNode(nodeId);
        if (node) {
            int nodeIndex = nodeIdToIndex_.at(nodeId);
            fullDofIndices.push_back(nodeIndex * 2);     // DOF X
            fullDofIndices.push_back(nodeIndex * 2 + 1); // DOF Y
        }
    }

    return fullDofIndices; // ВОЗВРАЩАЕМ ПОЛНЫЕ ИНДЕКСЫ!
}

bool Assembly::validate() const {
    if (nodes_.empty()) {
        std::cerr << "Assembly validation failed: No nodes defined" << std::endl;
        return false;
    }

    if (elements_.empty()) {
        std::cerr << "Assembly validation failed: No elements defined" << std::endl;
        return false;
    }

    if (materials_.empty()) {
        std::cerr << "Assembly validation failed: No materials defined" << std::endl;
        return false;
    }

    // Проверяем корректность элементов - передаем ТОЛЬКО узлы элемента
    for (const auto& element : elements_) {
        // Получаем узлы конкретного элемента
        std::vector<std::shared_ptr<Node>> elementNodes;
        for (int nodeId : element->getNodeIds()) {
            auto node = getNode(nodeId);
            if (node) {
                elementNodes.push_back(node);
            }
            else {
                std::cerr << "Assembly validation failed: Element " << element->getId()
                    << " references non-existent node " << nodeId << std::endl;
                return false;
            }
        }

        // Проверяем элемент с его собственными узлами
        if (!element->isValid(elementNodes)) {
            std::cerr << "Assembly validation failed: Element " << element->getId()
                << " is invalid (negative Jacobian or other issue)" << std::endl;
            return false;
        }
    }

    // Проверяем, что все элементы ссылаются на существующие материалы
    for (const auto& element : elements_) {
        if (!getMaterial(element->getMaterialId())) {
            std::cerr << "Assembly validation failed: Element " << element->getId()
                << " references non-existent material " << element->getMaterialId() << std::endl;
            return false;
        }
    }

    return true;
}

int Assembly::getGlobalDofIndex(int nodeId, int direction) const {
    auto it = nodeIdToIndex_.find(nodeId);
    if (it == nodeIdToIndex_.end()) {
        throw std::invalid_argument("Node with ID " + std::to_string(nodeId) + " not found");
    }

    if (direction < 0 || direction > 1) {
        throw std::invalid_argument("Direction must be 0 (X) or 1 (Y)");
    }

    return it->second * 2 + direction;
}

void Assembly::buildNodeIndexMap() {
    nodeIdToIndex_.clear();
    for (size_t i = 0; i < nodes_.size(); ++i) {
        nodeIdToIndex_[nodes_[i]->getId()] = i;
    }
}

std::vector<int> Assembly::getElementFullDofIndices(int elementId) const {
    // Используется при сборке матрицы жесткости - ДО граничных условий
    return getElementDofIndicesInternal(elementId);
}

std::vector<int> Assembly::getElementReducedDofIndices(int elementId) const {
    // Используется после применения граничных условий
    std::vector<int> fullDofIndices = getElementDofIndicesInternal(elementId);
    std::vector<int> reducedDofIndices;

    for (int fullDof : fullDofIndices) {
        if (fullDof >= 0 && fullDof < dofMapping_.fullToReduced.size()) {
            int reducedDof = dofMapping_.fullToReduced[fullDof];
            reducedDofIndices.push_back(reducedDof);
        }
        else {
            reducedDofIndices.push_back(-1); // Некорректный индекс
        }
    }

    return reducedDofIndices;
}

std::vector<int> Assembly::getElementDofIndicesInternal(int elementId) const {
    auto element = getElement(elementId);
    if (!element) {
        throw std::invalid_argument("Element with ID " + std::to_string(elementId) + " not found");
    }

    std::vector<int> fullDofIndices;

    // Получаем полные индексы DOF (всегда от 0 до totalDof-1)
    for (int nodeId : element->getNodeIds()) {
        auto node = getNode(nodeId);
        if (node) {
            auto it = nodeIdToIndex_.find(nodeId);
            if (it != nodeIdToIndex_.end()) {
                int nodeIndex = it->second;
                int dofX = nodeIndex * 2;     // DOF X
                int dofY = nodeIndex * 2 + 1; // DOF Y

                fullDofIndices.push_back(dofX);
                fullDofIndices.push_back(dofY);
            }
        }
    }

    return fullDofIndices;
}

void Assembly::assembleConcentratedForces(Eigen::VectorXd& globalF) const {
    int totalDof = getTotalDofCount();
    globalF.resize(totalDof);
    globalF.setZero();

    for (const auto& force : concentratedForces_) {
        int nodeId = force->getNodeId();


        try {
            // Используем ПОЛНЫЕ индексы
            int dofX = getGlobalDofIndex(nodeId, 0);
            int dofY = getGlobalDofIndex(nodeId, 1);

            if (dofX >= 0 && dofX < globalF.size()) {
                globalF(dofX) += force->getForceX();
            }
            if (dofY >= 0 && dofY < globalF.size()) {
                globalF(dofY) += force->getForceY();
            }
        }
        catch (const std::exception& e) {
            std::cout << "WARNING: Failed to apply force to node " << nodeId << ": " << e.what() << std::endl;
        }
    }
}

void Assembly::addConcentratedForce(std::shared_ptr<ConcentratedForce> force) {
    if (!force) {
        throw std::invalid_argument("Cannot add null force");
    }

    // Проверяем существование узла
    if (!getNode(force->getNodeId())) {
        throw std::invalid_argument("Node with ID " + std::to_string(force->getNodeId()) + " not found for force");
    }

    concentratedForces_.push_back(force);
}

void Assembly::addConcentratedForces(const std::vector<std::shared_ptr<ConcentratedForce>>& forces) {
    for (const auto& force : forces) {
        addConcentratedForce(force);
    }
}

