
#include "ContactPlane.h"
#include "solver.h"
#include <stdexcept>
#include <cmath>

LinearSolver solver;

ContactPlaneElement::ContactPlaneElement(int id, const std::vector<int>& nodeIds, int materialId,
    double penaltyParameter, int contactSurface)
    : PlaneIsoparametricElement(id, nodeIds, materialId)
    , penaltyParameter_(penaltyParameter)
    , contactSurface_(contactSurface)
{
}

Eigen::MatrixXd ContactPlaneElement::computePenaltyStiffnessMatrix(const std::vector<Node>& nodes,
    const Material& material) const
{
    if (contactSurface_ < 0 || contactSurface_ > 3) {
        throw std::runtime_error("Contact surface not properly set");
    }

    LinearSolver solver;

    auto surfaceIntegrand = [this, &nodes, &material](double param, double) -> Eigen::MatrixXd {
        double local_xi, local_eta;

        switch (contactSurface_) {
        case 0: local_xi = param; local_eta = -1.0; break; // нижняя
        case 1: local_xi = 1.0; local_eta = param; break;  // правая  
        case 2: local_xi = param; local_eta = 1.0; break;  // верхняя
        case 3: local_xi = -1.0; local_eta = param; break; // левая
        default: throw std::runtime_error("Invalid contact surface");
        }

        // Матрица формы
        Eigen::Vector4d N = shapeFunctionsLocal(local_xi, local_eta);

        // Нормаль к поверхности
        Eigen::Vector2d normal = getSurfaceNormal(contactSurface_, nodes);
        Eigen::Matrix2d nnT = normal * normal.transpose();

        // Матрица проекции
        Eigen::MatrixXd projection = Eigen::MatrixXd::Zero(8, 8);
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                projection.block<2, 2>(2 * j, 2 * k) = N(j) * N(k) * nnT;
            }
        }

        // Длина сегмента
        Eigen::Matrix2d J = jacobian(local_xi, local_eta, nodes);
        double surfaceJacobian;
        switch (contactSurface_) {
        case 0: case 2: surfaceJacobian = std::sqrt(J(0, 0) * J(0, 0) + J(0, 1) * J(0, 1)); break;
        case 1: case 3: surfaceJacobian = std::sqrt(J(1, 0) * J(1, 0) + J(1, 1) * J(1, 1)); break;
        default: surfaceJacobian = 1.0;
        }

        return projection * penaltyParameter_ * surfaceJacobian * material.getThickness();
        };

    return solver.computeGaussIntegral(surfaceIntegrand, 2);
}

Eigen::VectorXd ContactPlaneElement::computePenaltyForces(const Eigen::VectorXd& displacements,
    const std::vector<Node>& nodes,
    const Material& material) const
{
    if (contactSurface_ < 0 || contactSurface_ > 3) {
        throw std::runtime_error("Contact surface not properly set");
    }

    Eigen::VectorXd fe_contact = Eigen::VectorXd::Zero(8);

    // Вычисляем зазор
    double gap = calculateGap(contactSurface_, displacements, nodes);

    // Если зазор положительный - контакта нет
    if (gap >= 0) {
        return fe_contact;
    }

    auto forceIntegrand = [this, &nodes, &material, &displacements, gap](double param, double) -> Eigen::MatrixXd {
        double local_xi, local_eta;

        switch (contactSurface_) {
        case 0: local_xi = param; local_eta = -1.0; break;
        case 1: local_xi = 1.0; local_eta = param; break;
        case 2: local_xi = param; local_eta = 1.0; break;
        case 3: local_xi = -1.0; local_eta = param; break;
        }

        // Матрица формы на поверхности
        Eigen::Vector4d N = shapeFunctionsLocal(local_xi, local_eta);

        // Нормаль к поверхности
        Eigen::Vector2d normal = getSurfaceNormal(contactSurface_, nodes);

        // Длина сегмента
        Eigen::Matrix2d J = jacobian(local_xi, local_eta, nodes);
        double surfaceJacobian;
        switch (contactSurface_) {
        case 0: case 2: surfaceJacobian = std::sqrt(J(0, 0) * J(0, 0) + J(0, 1) * J(0, 1)); break;
        case 1: case 3: surfaceJacobian = std::sqrt(J(1, 0) * J(1, 0) + J(1, 1) * J(1, 1)); break;
        default: surfaceJacobian = 1.0;
        }

        // Штрафные силы: fe = -penalty * gap * N^T * n
        Eigen::Vector2d normalForce = -penaltyParameter_ * gap * normal;
        Eigen::VectorXd fe = Eigen::VectorXd::Zero(8);

        for (int j = 0; j < 4; ++j) {
            fe.segment<2>(2 * j) = N(j) * normalForce * surfaceJacobian * material.getThickness();
        }

        return fe;
        };

    Eigen::MatrixXd forceMatrix = solver.computeGaussIntegral(forceIntegrand, 2);
    return forceMatrix.col(0);
}

bool ContactPlaneElement::checkContact(int surfaceIndex, const Eigen::VectorXd& displacements,
    const std::vector<Node>& nodes) const
{
    double gap = calculateGap(surfaceIndex, displacements, nodes);
    return gap < 0.0; // Контакт если зазор отрицательный (проникание)
}

double ContactPlaneElement::calculateGap(int surfaceIndex, const Eigen::VectorXd& displacements,
    const std::vector<Node>& nodes) const
{
    Eigen::Vector2d normal = getSurfaceNormal(surfaceIndex, nodes);

    // Определяем узлы на поверхности
    std::vector<int> surfaceNodes;
    switch (surfaceIndex) {
    case 0: surfaceNodes = { 0, 1 }; break; // нижняя
    case 1: surfaceNodes = { 1, 2 }; break; // правая
    case 2: surfaceNodes = { 2, 3 }; break; // верхняя
    case 3: surfaceNodes = { 3, 0 }; break; // левая
    default: return 0.0;
    }

    // Среднее нормальное перемещение поверхности
    double avgNormalDisp = 0.0;
    for (int localNodeIndex : surfaceNodes) {
        int globalDofIndex = localNodeIndex * 2; // Предполагаем последовательную нумерацию DOF
        Eigen::Vector2d nodeDisp = displacements.segment<2>(globalDofIndex);
        avgNormalDisp += nodeDisp.dot(normal);
    }
    avgNormalDisp /= surfaceNodes.size();

    return avgNormalDisp; // Упрощенный расчет зазора
}

Eigen::MatrixXd ContactPlaneElement::computeSurfaceStiffness(int surfaceIndex, const std::vector<Node>& nodes,
    const Material& material) const
{
    // Временно меняем поверхность для расчета
    int oldSurface = contactSurface_;
    const_cast<ContactPlaneElement*>(this)->contactSurface_ = surfaceIndex;

    auto stiffness = computePenaltyStiffnessMatrix(nodes, material);

    const_cast<ContactPlaneElement*>(this)->contactSurface_ = oldSurface;
    return stiffness;
}