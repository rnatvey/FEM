#pragma once
#include "BaseElement.h"
#include <Eigen/Dense>

class PlaneIsoparametricElement : public BaseElement {
public:
    PlaneIsoparametricElement(int id, const std::vector<int>& nodeIds, int materialId);
    ~PlaneIsoparametricElement() override = default;

    // Реализация виртуальных функций
    Eigen::MatrixXd computeStiffnessMatrix(const std::vector<Node>& nodes, 
                                          const Material& material) const override;
    
    Eigen::MatrixXd shapeFunctions(double xi, double eta) const override;
    Eigen::MatrixXd shapeFunctionsDerivatives(double xi, double eta) const override;
    Eigen::Matrix2d jacobian(double xi, double eta, const std::vector<Node>& nodes) const override;
    Eigen::MatrixXd strainDisplacementMatrix(double xi, double eta, 
                                            const std::vector<Node>& nodes,
                                            const Material& material) const override;

    Eigen::Vector3d computeStress(double xi, double eta, 
                                 const Eigen::VectorXd& displacements,
                                 const std::vector<Node>& nodes,
                                 const Material& material) const override;
    
    Eigen::Vector3d computeStrain(double xi, double eta,
                                 const Eigen::VectorXd& displacements,
                                 const std::vector<Node>& nodes,
                                 const Material& material) const override;

    Eigen::VectorXd computeEquivalentNodalForces(const Eigen::VectorXd& bodyForces,
                                                const std::vector<Node>& nodes,
                                                const Material& material) const override;

    bool isValid(const std::vector<Node>& nodes) const override;
    std::string type() const override { return "PlaneIsoparametric"; }

    // Дополнительные методы для контакта
    Eigen::Vector2d getSurfaceNormal(int surfaceIndex, const std::vector<Node>& nodes) const;
    std::vector<Eigen::Vector2d> getSurfacePoints(int surfaceIndex, const std::vector<Node>& nodes) const;

private:
    // Локальные вспомогательные методы
    Eigen::Vector4d shapeFunctionsLocal(double xi, double eta) const;
    Eigen::MatrixXd shapeFunctionsDerivativesLocal(double xi, double eta) const;
    
    // Получение координат узлов элемента
    Eigen::MatrixXd getNodalCoordinates(const std::vector<Node>& nodes) const;
};