#pragma once
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "assembly.h"
#include "solver.h"

class FEModel {
public:
    FEModel();
    ~FEModel() = default;

    // === Управление моделью ===
    void setAssembly(std::shared_ptr<Assembly> assembly);
    std::shared_ptr<Assembly> getAssembly() const { return assembly_; }

    // === Параметры решения ===
    void setSolverTolerance(double tolerance) { tolerance_ = tolerance; }
    void setMaxIterations(int maxIter) { maxIterations_ = maxIter; }
    void setPenaltyParameter(double penalty) { penaltyParameter_ = penalty; }

    // === Основные методы решения ===
    bool solve();  // Основной метод решения
    bool solveContact();  // Решение контактной задачи (итерационный метод)

    // === Результаты ===
    const Eigen::VectorXd& getDisplacements() const { return displacements_; }
    const Eigen::VectorXd& getReactionForces() const { return reactionForces_; }
    double getSolutionTime() const { return solutionTime_; }
    int getIterationCount() const { return iterationCount_; }

    // === Постпроцессинг ===
    Eigen::Vector3d getElementStress(int elementId, double xi, double eta) const;
    Eigen::Vector3d getElementStrain(int elementId, double xi, double eta) const;

    // Новые методы для узловых величин
    std::vector<Eigen::Vector3d> getNodalStresses() const;
    std::vector<Eigen::Vector2d> getNodalDisplacements() const;
    std::vector<Eigen::Vector3d> getNodalStrains() const;

    // Усредненные напряжения в узлах (для визуализации)
    Eigen::Vector3d getNodeStress(int nodeId) const;
    Eigen::Vector2d getNodeDisplacement(int nodeId) const;
    Eigen::Vector3d getNodeStrain(int nodeId) const;


    // === Валидация ===
    bool validate() const;

private:
    std::shared_ptr<Assembly> assembly_;
    std::unique_ptr<LinearSolver> solver_;

    // Результаты
    Eigen::VectorXd displacements_;
    Eigen::VectorXd reactionForces_;

    // Параметры решения
    double tolerance_ = 1.0e-8;
    int maxIterations_ = 100;
    double penaltyParameter_ = 1.0e6;
    double solutionTime_ = 0.0;
    int iterationCount_ = 0;

    // Вспомогательные методы
    bool solveLinearSystem();  // Решение линейной системы
    bool solveContactIterative();  // Итерационное решение контакта
    void calculateReactionForces();  // Расчет реакций
    void applyContactConditions();  // Применение контактных условий

    // Временные матрицы для контактного решения
    Eigen::SparseMatrix<double> globalK_;
    Eigen::VectorXd globalF_;
    Eigen::VectorXd contactForces_;

    void calculateNodalAverages() const;
    mutable std::vector<Eigen::Vector3d> nodalStresses_;  // Кэш узловых напряжений
    mutable std::vector<Eigen::Vector2d> nodalDisplacements_; // Кэш узловых перемещений
    mutable std::vector<Eigen::Vector3d> nodalStrains_;   // Кэш узловых деформаций
    mutable bool nodalDataCalculated_ = false;

};