#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <string>
#include <vector>

#include "ContactTypes.h"
#include "IRigidPlaneContactSolver.h"
#include "assembly.h"
#include "solver.h"

class FEModel {
public:
    struct PerformanceMetrics {
        int linearIterations = 0;
        int nonlinearIterations = 0;
        int activeSetSize = 0;
        double assemblyTimeSeconds = 0.0;
        double solveTimeSeconds = 0.0;
        double totalTimeSeconds = 0.0;
        double maxPenetration = 0.0;
        double linearResidualEstimate = 0.0;
        double linearResidualNorm = 0.0;
        double equilibriumResidualNorm = 0.0;
        double contactForceNorm = 0.0;
        double contactStateUpdateNorm = 0.0;
        double contactStateRelativeUpdateNorm = 0.0;
        double maxNormalContactMultiplier = 0.0;
        double meanNormalContactMultiplier = 0.0;
        int activeContactGaussPoints = 0;
        bool linearSolveConverged = false;
        bool usedDirectLinearSolver = false;
        std::string linearSolverBackend = "uninitialized";
        std::string contactMethod = "none";
        Eigen::Index matrixNonZeros = 0;
    };

    FEModel();
    ~FEModel() = default;

    void setAssembly(std::shared_ptr<Assembly> assembly);
    std::shared_ptr<Assembly> getAssembly() const { return assembly_; }

    void setSolverTolerance(double tolerance) { tolerance_ = tolerance; }
    void setMaxIterations(int maxIter) { maxIterations_ = maxIter; }
    void setPenaltyParameter(double penalty) { penaltyParameter_ = penalty; }

    void setContactSolver(std::unique_ptr<IRigidPlaneContactSolver> contactSolver);
    void configureRigidPlaneContact(const RigidPlane2D& plane,
        const std::vector<ContactFacet>& facets,
        double penaltyParameter);
    void configureRigidPlanePenaltyContact(const RigidPlane2D& plane,
        const std::vector<ContactFacet>& facets,
        double penaltyParameter);
    void configureRigidPlaneAugmentedLagrangianContact(const RigidPlane2D& plane,
        const std::vector<ContactFacet>& facets,
        const AugmentedLagrangianSettings& settings);
    bool hasContactSolver() const { return contactSolver_ != nullptr; }

    bool solve();
    bool solveContact();

    const Eigen::VectorXd& getDisplacements() const { return displacements_; }
    const Eigen::VectorXd& getReactionForces() const { return reactionForces_; }
    const Eigen::VectorXd& getContactForces() const { return contactForces_; }
    double getSolutionTime() const { return performanceMetrics_.totalTimeSeconds; }
    int getIterationCount() const { return iterationCount_; }
    const PerformanceMetrics& getPerformanceMetrics() const { return performanceMetrics_; }
    RigidPlane2D getContactPlane() const;
    std::vector<ContactFacet> getContactFacets() const;
    std::vector<ContactFacetResult> getContactFacetResults() const;
    ContactState evaluateCurrentContactState() const;

    Eigen::Vector3d getElementStress(int elementId, double xi, double eta) const;
    Eigen::Vector3d getElementStrain(int elementId, double xi, double eta) const;

    std::vector<Eigen::Vector3d> getNodalStresses() const;
    std::vector<Eigen::Vector2d> getNodalDisplacements() const;
    std::vector<Eigen::Vector2d> getNodalReactionForces() const;
    std::vector<Eigen::Vector2d> getNodalContactForces() const;
    std::vector<double> getNodalContactSignedDistances() const;
    std::vector<double> getNodalContactPenetrations() const;
    std::vector<Eigen::Vector3d> getNodalStrains() const;

    Eigen::Vector3d getNodeStress(int nodeId) const;
    Eigen::Vector2d getNodeDisplacement(int nodeId) const;
    Eigen::Vector3d getNodeStrain(int nodeId) const;

    bool validate() const;

private:
    bool solveLinearSystem();
    bool solveContactIterative();
    void calculateReactionForces();
    void applyContactConditions(const Eigen::VectorXd& trialDisplacements,
        ContactIterationInfo& iterationInfo);
    void assembleExternalForces(Eigen::VectorXd& globalF) const;
    std::vector<std::shared_ptr<Node>> getElementNodes(int elementId) const;

    void calculateNodalAverages() const;

    std::shared_ptr<Assembly> assembly_;
    std::unique_ptr<LinearSolver> solver_;
    std::unique_ptr<IRigidPlaneContactSolver> contactSolver_;

    Eigen::VectorXd displacements_;
    Eigen::VectorXd reactionForces_;

    double tolerance_ = 1.0e-8;
    int maxIterations_ = 100;
    double penaltyParameter_ = 1.0e6;
    AugmentedLagrangianSettings augmentedLagrangianSettings_;
    int iterationCount_ = 0;

    Eigen::SparseMatrix<double> globalK_;
    Eigen::VectorXd globalF_;
    Eigen::VectorXd contactForces_;

    PerformanceMetrics performanceMetrics_;
    mutable std::vector<Eigen::Vector3d> nodalStresses_;
    mutable std::vector<Eigen::Vector2d> nodalDisplacements_;
    mutable std::vector<Eigen::Vector3d> nodalStrains_;
    mutable bool nodalDataCalculated_ = false;
};
