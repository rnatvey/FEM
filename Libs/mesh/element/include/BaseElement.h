#pragma once
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include "node.h"
#include "material.h"

class BaseElement {
public:
    BaseElement(int id, const std::vector<int>& nodeIds, int materialId);
    virtual ~BaseElement() = default;

    // === Чисто виртуальные функции ===
    virtual Eigen::MatrixXd computeStiffnessMatrix(const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const = 0;

    virtual Eigen::MatrixXd shapeFunctions(double xi, double eta) const = 0;
    virtual Eigen::MatrixXd shapeFunctionsDerivatives(double xi, double eta) const = 0;
    virtual Eigen::Matrix2d jacobian(double xi, double eta, const std::vector<std::shared_ptr<Node>>& nodes) const = 0;
    virtual Eigen::MatrixXd strainDisplacementMatrix(double xi, double eta, const std::vector<std::shared_ptr<Node>>&, const std::shared_ptr<Material>& material) const = 0;

    // === Информационные функции ===
    int getId() const { return id_; }
    int getNodeCount() const { return nodeIds_.size(); }
    int getDofCount() const { return nodeIds_.size() * 2; }
    const std::vector<int>& getNodeIds() const { return nodeIds_; }
    int getMaterialId() const { return materialId_; }

    virtual std::string type() const = 0;
    virtual bool isValid(const std::vector<std::shared_ptr<Node>>& nodes) const = 0;

    // === Расчетные функции ===
    virtual Eigen::Vector3d computeStress(double xi, double eta,
        const Eigen::VectorXd& displacements,
        const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const = 0;

    virtual Eigen::Vector3d computeStrain(double xi, double eta,
        const Eigen::VectorXd& displacements,
        const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const = 0;

    virtual Eigen::VectorXd computeEquivalentNodalForces(const Eigen::VectorXd& bodyForces,
        const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const = 0;

    virtual std::vector<std::shared_ptr<Node>> getElementNodes(
        const std::vector<std::shared_ptr<Node>>& allNodes) const {

        std::vector<std::shared_ptr<Node>> elementNodes;

        // Ищем узлы по их ID
        for (int nodeId : nodeIds_) {
            bool found = false;
            for (const auto& node : allNodes) {
                if (node->getId() == nodeId) {
                    elementNodes.push_back(node);
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Node " + std::to_string(nodeId) +
                    " not found for element " +
                    std::to_string(id_));
            }
        }

        return elementNodes;
    }

protected:
    int id_;
    std::vector<int> nodeIds_;  // Только ID узлов, не координаты!
    int materialId_;
};