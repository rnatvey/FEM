#pragma once

#include "Planeisoparametric.h"


class ContactPlaneElement : public PlaneIsoparametricElement {
public:
    ContactPlaneElement(int id, const std::vector<int>& nodeIds, int materialId,
        double penaltyParameter = 1.0e6, int contactSurface = -1);

    // Штрафная матрица жесткости для контакта
    Eigen::MatrixXd computePenaltyStiffnessMatrix(const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const;

    // Штрафные силы для контакта
    Eigen::VectorXd computePenaltyForces(const Eigen::VectorXd& displacements,
        const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const;

    // Проверка контакта для поверхности
    bool checkContact(int surfaceIndex, const Eigen::VectorXd& displacements,
        const std::vector<std::shared_ptr<Node>>& nodes) const;

    // Установка параметров контакта
    void setPenaltyParameter(double penalty) { penaltyParameter_ = penalty; }
    void setContactSurface(int surface) { contactSurface_ = surface; }
    double getPenaltyParameter() const { return penaltyParameter_; }
    int getContactSurface() const { return contactSurface_; }

    std::string type() const override { return "ContactPlaneElement"; }

private:
    double penaltyParameter_;  // Штрафной параметр
    int contactSurface_;       // Номер контактной поверхности (0-3)

    // Методы для расчета контактных взаимодействий
    double calculateGap(int surfaceIndex, const Eigen::VectorXd& displacements,
        const std::vector<std::shared_ptr<Node>>& nodes) const;
    Eigen::MatrixXd computeSurfaceStiffness(int surfaceIndex, const std::vector<std::shared_ptr<Node>>& nodes,
        const std::shared_ptr<Material>& material) const;
};