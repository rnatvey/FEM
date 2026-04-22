#include "assembly.h"
#include "solver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include "planeisometric/Planeisoparametric.h"

Assembly::Assembly() {
}

void Assembly::addNode(std::shared_ptr<Node> node) {
    if (!node) {
        throw std::invalid_argument("Cannot add null node");
    }

    // РџСЂРѕРІРµСЂСЏРµРј СѓРЅРёРєР°Р»СЊРЅРѕСЃС‚СЊ ID
    for (const auto& existingNode : nodes_) {
        if (existingNode->getId() == node->getId()) {
            throw std::invalid_argument("Node with ID " + std::to_string(node->getId()) + " already exists");
        }
    }

    nodes_.push_back(node);
    buildNodeIndexMap(); // РћР±РЅРѕРІР»СЏРµРј РєР°СЂС‚Сѓ РёРЅРґРµРєСЃРѕРІ
}

void Assembly::addNodes(const std::vector<std::shared_ptr<Node>>& nodes) {
    for (const auto& node : nodes) {
        addNode(node);
    }
}

void Assembly::addLineofNodes(size_t startId, vector2 start_point,  vector2 end_point, int num_of_nodes ) 
{
    vector2 direction = end_point - start_point;
    const double step = Norm(direction)/(num_of_nodes-1);
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

void Assembly::addFiniteStrainMaterial(std::shared_ptr<FiniteStrainMaterial> material) {
    if (!material) {
        throw std::invalid_argument("Cannot add null finite-strain material");
    }

    const int materialId = material->getId();
    if (finiteStrainMaterials_.find(materialId) != finiteStrainMaterials_.end()) {
        throw std::invalid_argument(
            "Finite-strain material with ID " + std::to_string(materialId) +
            " already exists");
    }

    finiteStrainMaterials_[materialId] = std::move(material);
}

std::shared_ptr<FiniteStrainMaterial> Assembly::getFiniteStrainMaterial(int id) const {
    auto it = finiteStrainMaterials_.find(id);
    if (it != finiteStrainMaterials_.end()) {
        return it->second;
    }
    return nullptr;
}

void Assembly::addElement(std::shared_ptr<BaseElement> element) {
    if (!element) {
        throw std::invalid_argument("Cannot add null element");
    }

    // РџСЂРѕРІРµСЂСЏРµРј СѓРЅРёРєР°Р»СЊРЅРѕСЃС‚СЊ ID
    for (const auto& existingElement : elements_) {
        if (existingElement->getId() == element->getId()) {
            throw std::invalid_argument("Element with ID " + std::to_string(element->getId()) + " already exists");
        }
    }

    // РџСЂРѕРІРµСЂСЏРµРј СЃСѓС‰РµСЃС‚РІРѕРІР°РЅРёРµ СѓР·Р»РѕРІ
    for (int nodeId : element->getNodeIds()) {
        if (!getNode(nodeId)) {
            throw std::invalid_argument("Node with ID " + std::to_string(nodeId) + " not found for element " + std::to_string(element->getId()));
        }
    }

    // РџСЂРѕРІРµСЂСЏРµРј СЃСѓС‰РµСЃС‚РІРѕРІР°РЅРёРµ РјР°С‚РµСЂРёР°Р»Р°
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

void Assembly::addFiniteStrainElement(std::shared_ptr<FiniteStrainQ4Element> element) {
    if (!element) {
        throw std::invalid_argument("Cannot add null finite-strain element");
    }

    for (const auto& existingElement : finiteStrainElements_) {
        if (existingElement->getId() == element->getId()) {
            throw std::invalid_argument(
                "Finite-strain element with ID " + std::to_string(element->getId()) +
                " already exists");
        }
    }

    for (int nodeId : element->getNodeIds()) {
        if (!getNode(nodeId)) {
            throw std::invalid_argument(
                "Node with ID " + std::to_string(nodeId) +
                " not found for finite-strain element " +
                std::to_string(element->getId()));
        }
    }

    if (!getFiniteStrainMaterial(element->getMaterialId())) {
        throw std::invalid_argument(
            "Finite-strain material with ID " +
            std::to_string(element->getMaterialId()) +
            " not found for element " + std::to_string(element->getId()));
    }

    finiteStrainElements_.push_back(std::move(element));
}

void Assembly::addFiniteStrainElements(
    const std::vector<std::shared_ptr<FiniteStrainQ4Element>>& elements) {
    for (const auto& element : elements) {
        addFiniteStrainElement(element);
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

std::shared_ptr<FiniteStrainQ4Element> Assembly::getFiniteStrainElement(int id) const {
    for (const auto& element : finiteStrainElements_) {
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
    triplets.reserve(elements_.size() * 64);

    for (const auto& element : elements_) {
        auto material = getMaterial(element->getMaterialId());
        if (!material) {
            continue;
        }

        std::vector<std::shared_ptr<Node>> elementNodes;
        elementNodes.reserve(element->getNodeIds().size());
        for (int nodeId : element->getNodeIds()) {
            auto node = getNode(nodeId);
            if (node) elementNodes.push_back(node);
        }

        // Р’С‹С‡РёСЃР»СЏРµРј РјР°С‚СЂРёС†Сѓ Р¶РµСЃС‚РєРѕСЃС‚Рё СЌР»РµРјРµРЅС‚Р°
        Eigen::MatrixXd ke = element->computeStiffnessMatrix(elementNodes, material);

        
        std::vector<int> dofIndices = getElementFullDofIndices(element->getId());

        // Р”РѕР±Р°РІР»СЏРµРј РІ РіР»РѕР±Р°Р»СЊРЅСѓСЋ РјР°С‚СЂРёС†Сѓ
        for (int i = 0; i < dofIndices.size(); ++i) {
            for (int j = 0; j < dofIndices.size(); ++j) {
                if (dofIndices[i] >= 0 && dofIndices[j] >= 0) {
                    triplets.emplace_back(dofIndices[i], dofIndices[j], ke(i, j));
                }
            }
        }
    }

    globalK.setFromTriplets(triplets.begin(), triplets.end());
    globalK.makeCompressed();
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

        // РСЃРїРѕР»СЊР·СѓРµРј РџРћР›РќР«Р• РёРЅРґРµРєСЃС‹
        auto dofIndices = getElementFullDofIndices(element->getId());

        for (int i = 0; i < dofIndices.size(); ++i) {
            if (dofIndices[i] >= 0) {
                globalF(dofIndices[i]) += fe(i);
            }
        }
    }
}

void Assembly::assembleFiniteStrainSystem(const Eigen::VectorXd& fullDisplacements,
    Eigen::SparseMatrix<double>& globalTangent,
    Eigen::VectorXd& globalInternalForce,
    double& totalStrainEnergy,
    std::vector<FiniteStrainElementResponse>* elementResponses,
    bool computeTangent) const {
    const int totalDof = getTotalDofCount();
    if (fullDisplacements.size() != totalDof) {
        throw std::invalid_argument(
            "Full displacement vector size mismatch in assembleFiniteStrainSystem");
    }

    globalInternalForce = Eigen::VectorXd::Zero(totalDof);
    totalStrainEnergy = 0.0;

    if (elementResponses) {
        elementResponses->clear();
        elementResponses->reserve(finiteStrainElements_.size());
    }

    std::vector<Eigen::Triplet<double>> triplets;
    if (computeTangent) {
        globalTangent.resize(totalDof, totalDof);
        globalTangent.setZero();
        triplets.reserve(finiteStrainElements_.size() * 64);
    }

    for (const auto& element : finiteStrainElements_) {
        auto material = getFiniteStrainMaterial(element->getMaterialId());
        if (!material) {
            throw std::runtime_error(
                "Finite-strain material " + std::to_string(element->getMaterialId()) +
                " not found during assembly");
        }

        std::vector<std::shared_ptr<Node>> elementNodes;
        elementNodes.reserve(element->getNodeIds().size());
        for (int nodeId : element->getNodeIds()) {
            auto node = getNode(nodeId);
            if (!node) {
                throw std::runtime_error(
                    "Finite-strain element " + std::to_string(element->getId()) +
                    " references missing node " + std::to_string(nodeId));
            }
            elementNodes.push_back(node);
        }

        std::vector<int> dofIndices;
        dofIndices.reserve(element->getNodeIds().size() * 2);
        for (int nodeId : element->getNodeIds()) {
            dofIndices.push_back(getGlobalDofIndex(nodeId, 0));
            dofIndices.push_back(getGlobalDofIndex(nodeId, 1));
        }

        Eigen::VectorXd elementDisplacements = Eigen::VectorXd::Zero(dofIndices.size());
        for (int i = 0; i < static_cast<int>(dofIndices.size()); ++i) {
            elementDisplacements(i) = fullDisplacements(dofIndices[static_cast<size_t>(i)]);
        }

        FiniteStrainElementResponse elementResponse =
            element->evaluateResponse(
                elementNodes,
                elementDisplacements,
                *material,
                computeTangent);

        for (int i = 0; i < static_cast<int>(dofIndices.size()); ++i) {
            globalInternalForce(dofIndices[static_cast<size_t>(i)]) +=
                elementResponse.internalForce(i);
        }

        if (computeTangent) {
            for (int i = 0; i < static_cast<int>(dofIndices.size()); ++i) {
                for (int j = 0; j < static_cast<int>(dofIndices.size()); ++j) {
                    triplets.emplace_back(dofIndices[static_cast<size_t>(i)],
                        dofIndices[static_cast<size_t>(j)],
                        elementResponse.tangent(i, j));
                }
            }
        }

        totalStrainEnergy += elementResponse.totalStrainEnergy;
        if (elementResponses) {
            elementResponses->push_back(std::move(elementResponse));
        }
    }

    if (computeTangent) {
        globalTangent.setFromTriplets(triplets.begin(), triplets.end());
        globalTangent.makeCompressed();
    }
    else {
        globalTangent.resize(totalDof, totalDof);
        globalTangent.setZero();
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
    bc.fixX = true;  // Р”Р»СЏ РїСЂРµРґРїРёСЃР°РЅРЅС‹С… РїРµСЂРµРјРµС‰РµРЅРёР№ С„РёРєСЃРёСЂСѓРµРј РѕР±Рµ СЃС‚РµРїРµРЅРё СЃРІРѕР±РѕРґС‹
    bc.fixY = true;
    bc.hasPrescribedDisplacement = true;
    bc.prescribedDx = dx;
    bc.prescribedDy = dy;
    boundaryConditions_.push_back(bc);
}

void Assembly::addPrescribedDisplacementX(int nodeId, double dx) {
    BoundaryCondition bc;
    bc.nodeId = nodeId;
    bc.fixX = true;
    bc.fixY = false;
    bc.hasPrescribedDisplacement = true;
    bc.prescribedDx = dx;
    bc.prescribedDy = 0.0;
    boundaryConditions_.push_back(bc);
}

void Assembly::addPrescribedDisplacementY(int nodeId,  double dy) {
    BoundaryCondition bc;
    bc.nodeId = nodeId;
    bc.fixX = false;
    bc.fixY = true;
    bc.hasPrescribedDisplacement = true;
    bc.prescribedDx = 0.0;
    bc.prescribedDy = dy;
    boundaryConditions_.push_back(bc);
}

void Assembly::applyBoundaryConditions(Eigen::SparseMatrix<double>& globalK,
    Eigen::VectorXd& globalF) const {
    LinearSolver solver;

    const ConstraintData constraintData = buildConstraintData();
    const int totalDof = globalK.rows();
    populateDofMapping(totalDof, constraintData);

    // Apply all constraints through sparse reduction to avoid mutating the
    // matrix structure for prescribed displacements.
    if (!constraintData.constrainedDofs.empty()) {
        Eigen::SparseMatrix<double> reducedK;
        Eigen::VectorXd reducedF;
        std::vector<int> activeDofs;

        solver.reduceSystem(globalK,
            globalF,
            constraintData.constrainedDofs,
            constraintData.prescribedDofs,
            constraintData.prescribedValues,
            reducedK,
            reducedF,
            activeDofs);
        globalK = std::move(reducedK);
        globalF = std::move(reducedF);
    }
}

Assembly::ConstraintData Assembly::getConstraintData() const {
    return buildConstraintData();
}

void Assembly::initializeDofMapping(int totalDof) const {
    populateDofMapping(totalDof, buildConstraintData());
}

Assembly::ConstraintData Assembly::buildConstraintData() const {
    ConstraintData constraintData;

    auto appendUniqueDof = [](std::vector<int>& dofs, int dof) {
        if (std::find(dofs.begin(), dofs.end(), dof) == dofs.end()) {
            dofs.push_back(dof);
        }
    };
    auto appendOrValidatePrescribedDof =
        [](std::vector<int>& dofs, std::vector<double>& values, int dof, double value) {
            auto dofIt = std::find(dofs.begin(), dofs.end(), dof);
            if (dofIt == dofs.end()) {
                dofs.push_back(dof);
                values.push_back(value);
                return;
            }

            const size_t index = static_cast<size_t>(std::distance(dofs.begin(), dofIt));
            if (index >= values.size()) {
                throw std::runtime_error("Internal prescribed DOF bookkeeping mismatch");
            }

            if (std::abs(values[index] - value) > 1.0e-12) {
                throw std::runtime_error("Conflicting prescribed displacement for DOF " +
                    std::to_string(dof));
            }
        };

    for (const auto& bc : boundaryConditions_) {
        auto node = getNode(bc.nodeId);
        if (!node) {
            continue;
        }

        if (bc.fixX) {
            const int dofX = getGlobalDofIndex(bc.nodeId, 0);
            if (bc.hasPrescribedDisplacement) {
                appendOrValidatePrescribedDof(
                    constraintData.prescribedDofs,
                    constraintData.prescribedValues,
                    dofX,
                    bc.prescribedDx);
            }
            else {
                appendUniqueDof(constraintData.fixedZeroDofs, dofX);
            }
        }

        if (bc.fixY) {
            const int dofY = getGlobalDofIndex(bc.nodeId, 1);
            if (bc.hasPrescribedDisplacement) {
                appendOrValidatePrescribedDof(
                    constraintData.prescribedDofs,
                    constraintData.prescribedValues,
                    dofY,
                    bc.prescribedDy);
            }
            else {
                appendUniqueDof(constraintData.fixedZeroDofs, dofY);
            }
        }
    }

    constraintData.constrainedDofs = constraintData.fixedZeroDofs;
    for (int prescribedDof : constraintData.prescribedDofs) {
        appendUniqueDof(constraintData.constrainedDofs, prescribedDof);
    }

    return constraintData;
}

void Assembly::populateDofMapping(int totalDof,
    const ConstraintData& constraintData) const {
    dofMapping_.fullToReduced.clear();
    dofMapping_.reducedToFull.clear();
    dofMapping_.prescribedDofs.clear();
    dofMapping_.prescribedValues.clear();

    dofMapping_.fullToReduced.resize(totalDof, -1);

    std::vector<char> isConstrained(totalDof, 0);
    for (int dof : constraintData.constrainedDofs) {
        if (dof >= 0 && dof < totalDof) {
            isConstrained[dof] = 1;
        }
    }

    int reducedIndex = 0;
    for (int dof = 0; dof < totalDof; ++dof) {
        if (!isConstrained[dof]) {
            dofMapping_.fullToReduced[dof] = reducedIndex++;
            dofMapping_.reducedToFull.push_back(dof);
        }
    }

    dofMapping_.prescribedDofs = constraintData.prescribedDofs;
    dofMapping_.prescribedValues = constraintData.prescribedValues;
}

std::vector<int> Assembly::getElementDofIndices(int elementId) const {
    auto element = getElement(elementId);
    if (!element) {
        throw std::invalid_argument("Element with ID " + std::to_string(elementId) + " not found");
    }

    std::vector<int> fullDofIndices;

    // РџРѕР»СѓС‡Р°РµРј РїРѕР»РЅС‹Рµ РёРЅРґРµРєСЃС‹ DOF
    for (int nodeId : element->getNodeIds()) {
        auto node = getNode(nodeId);
        if (node) {
            int nodeIndex = nodeIdToIndex_.at(nodeId);
            fullDofIndices.push_back(nodeIndex * 2);     // DOF X
            fullDofIndices.push_back(nodeIndex * 2 + 1); // DOF Y
        }
    }

    return fullDofIndices; // Р’РћР—Р’Р РђР©РђР•Рњ РџРћР›РќР«Р• РРќР”Р•РљРЎР«!
}

bool Assembly::validate() const {
    if (nodes_.empty()) {
        std::cerr << "Assembly validation failed: No nodes defined" << std::endl;
        return false;
    }

    const bool hasLinearElements = !elements_.empty();
    const bool hasFiniteStrainElements = !finiteStrainElements_.empty();
    if (!hasLinearElements && !hasFiniteStrainElements) {
        std::cerr << "Assembly validation failed: No elements defined" << std::endl;
        return false;
    }

    const bool hasLinearMaterials = !materials_.empty();
    const bool hasFiniteStrainMaterials = !finiteStrainMaterials_.empty();
    if (!hasLinearMaterials && !hasFiniteStrainMaterials) {
        std::cerr << "Assembly validation failed: No materials defined" << std::endl;
        return false;
    }

    if (hasLinearElements) {
        if (!hasLinearMaterials) {
            std::cerr << "Assembly validation failed: Linear elements require linear materials"
                      << std::endl;
            return false;
        }

        for (const auto& element : elements_) {
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

            if (!element->isValid(elementNodes)) {
                std::cerr << "Assembly validation failed: Element " << element->getId()
                          << " is invalid (negative Jacobian or other issue)" << std::endl;
                return false;
            }

            if (!getMaterial(element->getMaterialId())) {
                std::cerr << "Assembly validation failed: Element " << element->getId()
                          << " references non-existent material "
                          << element->getMaterialId() << std::endl;
                return false;
            }
        }
    }

    if (hasFiniteStrainElements) {
        if (!hasFiniteStrainMaterials) {
            std::cerr << "Assembly validation failed: Finite-strain elements require "
                      << "finite-strain materials" << std::endl;
            return false;
        }

        for (const auto& element : finiteStrainElements_) {
            std::vector<std::shared_ptr<Node>> elementNodes;
            for (int nodeId : element->getNodeIds()) {
                auto node = getNode(nodeId);
                if (!node) {
                    std::cerr << "Assembly validation failed: Finite-strain element "
                              << element->getId() << " references non-existent node "
                              << nodeId << std::endl;
                    return false;
                }
                elementNodes.push_back(node);
            }

            for (double xi : {-1.0, 1.0}) {
                for (double eta : {-1.0, 1.0}) {
                    if (element->referenceJacobian(xi, eta, elementNodes).determinant() <= 0.0) {
                        std::cerr << "Assembly validation failed: Finite-strain element "
                                  << element->getId()
                                  << " has a non-positive reference Jacobian" << std::endl;
                        return false;
                    }
                }
            }

            if (!getFiniteStrainMaterial(element->getMaterialId())) {
                std::cerr << "Assembly validation failed: Finite-strain element "
                          << element->getId() << " references non-existent material "
                          << element->getMaterialId() << std::endl;
                return false;
            }
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
    // РСЃРїРѕР»СЊР·СѓРµС‚СЃСЏ РїСЂРё СЃР±РѕСЂРєРµ РјР°С‚СЂРёС†С‹ Р¶РµСЃС‚РєРѕСЃС‚Рё - Р”Рћ РіСЂР°РЅРёС‡РЅС‹С… СѓСЃР»РѕРІРёР№
    return getElementDofIndicesInternal(elementId);
}

std::vector<int> Assembly::getElementReducedDofIndices(int elementId) const {
    // РСЃРїРѕР»СЊР·СѓРµС‚СЃСЏ РїРѕСЃР»Рµ РїСЂРёРјРµРЅРµРЅРёСЏ РіСЂР°РЅРёС‡РЅС‹С… СѓСЃР»РѕРІРёР№
    std::vector<int> fullDofIndices = getElementDofIndicesInternal(elementId);
    std::vector<int> reducedDofIndices;

    for (int fullDof : fullDofIndices) {
        if (fullDof >= 0 && fullDof < dofMapping_.fullToReduced.size()) {
            int reducedDof = dofMapping_.fullToReduced[fullDof];
            reducedDofIndices.push_back(reducedDof);
        }
        else {
            reducedDofIndices.push_back(-1); // РќРµРєРѕСЂСЂРµРєС‚РЅС‹Р№ РёРЅРґРµРєСЃ
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

    // РџРѕР»СѓС‡Р°РµРј РїРѕР»РЅС‹Рµ РёРЅРґРµРєСЃС‹ DOF (РІСЃРµРіРґР° РѕС‚ 0 РґРѕ totalDof-1)
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
            // РСЃРїРѕР»СЊР·СѓРµРј РџРћР›РќР«Р• РёРЅРґРµРєСЃС‹
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

    // РџСЂРѕРІРµСЂСЏРµРј СЃСѓС‰РµСЃС‚РІРѕРІР°РЅРёРµ СѓР·Р»Р°
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

void Assembly::addSurfaceLoad(int elementId, int surfaceIndex,
    std::shared_ptr<LoadFunction> loadFunction) {
    if (!loadFunction) {
        throw std::invalid_argument("Cannot add null load function");
    }

    // РџСЂРѕРІРµСЂСЏРµРј СЃСѓС‰РµСЃС‚РІРѕРІР°РЅРёРµ СЌР»РµРјРµРЅС‚Р°
    if (!getElement(elementId)) {
        throw std::invalid_argument("Element with ID " + std::to_string(elementId) + " not found");
    }

    // РџСЂРѕРІРµСЂСЏРµРј РєРѕСЂСЂРµРєС‚РЅРѕСЃС‚СЊ РЅРѕРјРµСЂР° РїРѕРІРµСЂС…РЅРѕСЃС‚Рё
    if (surfaceIndex < 0 || surfaceIndex > 3) {
        throw std::invalid_argument("Surface index must be between 0 and 3");
    }

    // РСЃРїРѕР»СЊР·СѓРµРј emplace_back СЃ РєРѕРЅСЃС‚СЂСѓРєС‚РѕСЂРѕРј
    surfaceLoads_.emplace_back(elementId, surfaceIndex, std::move(loadFunction));
}

void Assembly::assembleSurfaceLoads(Eigen::VectorXd& globalF) const {
    for (const auto& surfaceLoad : surfaceLoads_) {
        try {
            if (!surfaceLoad.loadFunction) continue;

            // РџСЂРёРјРµРЅСЏРµРј РЅР°РіСЂСѓР·РєСѓ Рє РїРѕРІРµСЂС…РЅРѕСЃС‚Рё СЌР»РµРјРµРЅС‚Р°
            Eigen::VectorXd surfaceForces = surfaceLoad.loadFunction->applyToElementSurface(
                surfaceLoad.elementId,
                surfaceLoad.surfaceIndex,
                shared_from_this()
            );

            // РџРѕР»СѓС‡Р°РµРј РёРЅРґРµРєСЃС‹ DOF СЌР»РµРјРµРЅС‚Р°
            auto dofIndices = getElementDofIndices(surfaceLoad.elementId);

            // Р”РѕР±Р°РІР»СЏРµРј РІ РіР»РѕР±Р°Р»СЊРЅС‹Р№ РІРµРєС‚РѕСЂ СЃРёР»
            for (size_t i = 0; i < dofIndices.size(); ++i) {
                if (dofIndices[i] >= 0 && dofIndices[i] < globalF.size()) {
                    globalF(dofIndices[i]) += surfaceForces(i);
                }
            }

        }
        catch (const std::exception& e) {
            std::cerr << "Error applying surface load to element "
                << surfaceLoad.elementId << ": " << e.what() << std::endl;
        }
    }
}
