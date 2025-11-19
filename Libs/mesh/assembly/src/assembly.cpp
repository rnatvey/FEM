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

    // Предварительно резервируем память
    std::vector<Eigen::Triplet<double>> triplets;

    for (const auto& element : elements_) {
        auto material = getMaterial(element->getMaterialId());
        if (!material) {
            throw std::runtime_error("Material not found for element " + std::to_string(element->getId()));
        }

        // Вычисляем матрицу жесткости элемента
        Eigen::MatrixXd ke = element->computeStiffnessMatrix(nodes_, material);

        // Получаем глобальные индексы степеней свободы
        std::vector<int> dofIndices = getElementDofIndices(element->getId());

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
        if (!material) {
            throw std::runtime_error("Material not found for element " + std::to_string(element->getId()));
        }

        // Вычисляем эквивалентные узловые силы
        Eigen::VectorXd fe = element->computeEquivalentNodalForces(bodyForces, nodes_, material);

        // Получаем глобальные индексы степеней свободы
        std::vector<int> dofIndices = getElementDofIndices(element->getId());

        // Добавляем в глобальный вектор
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

void Assembly::applyBoundaryConditions(Eigen::SparseMatrix<double>& globalK,
    Eigen::VectorXd& globalF) const {
    LinearSolver solver;

    for (const auto& bc : boundaryConditions_) {
        auto node = getNode(bc.nodeId);
        if (!node) continue;

        int nodeIndex = nodeIdToIndex_.at(bc.nodeId);

        if (bc.fixX) {
            int dofX = getGlobalDofIndex(bc.nodeId, 0);
            solver.applyBoundaryConditions(globalK, globalF, { dofX });
            if (bc.hasPrescribedDisplacement) {
                globalF(dofX) = bc.prescribedDx * globalK.coeff(dofX, dofX);
            }
        }

        if (bc.fixY) {
            int dofY = getGlobalDofIndex(bc.nodeId, 1);
            solver.applyBoundaryConditions(globalK, globalF, { dofY });
            if (bc.hasPrescribedDisplacement) {
                globalF(dofY) = bc.prescribedDy * globalK.coeff(dofY, dofY);
            }
        }
    }
}

std::vector<int> Assembly::getElementDofIndices(int elementId) const {
    auto element = getElement(elementId);
    if (!element) {
        throw std::invalid_argument("Element with ID " + std::to_string(elementId) + " not found");
    }

    std::vector<int> dofIndices;
    for (int nodeId : element->getNodeIds()) {
        auto node = getNode(nodeId);
        if (node) {
            int nodeIndex = nodeIdToIndex_.at(nodeId);
            dofIndices.push_back(nodeIndex * 2);     // DOF X
            dofIndices.push_back(nodeIndex * 2 + 1); // DOF Y
        }
    }
    return dofIndices;
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

    // Проверяем корректность элементов
    for (const auto& element : elements_) {
        if (!element->isValid(nodes_)) {
            std::cerr << "Assembly validation failed: Element " << element->getId() << " is invalid" << std::endl;
            return false;
        }
    }

    // Проверяем, что все элементы ссылаются на существующие узлы и материалы
    for (const auto& element : elements_) {
        for (int nodeId : element->getNodeIds()) {
            if (!getNode(nodeId)) {
                std::cerr << "Assembly validation failed: Element " << element->getId()
                    << " references non-existent node " << nodeId << std::endl;
                return false;
            }
        }

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