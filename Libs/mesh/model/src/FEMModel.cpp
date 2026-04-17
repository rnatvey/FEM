#include "FEMModel.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "ContactTypes.h"
#include "RigidPlaneAugmentedLagrangianContactSolver.h"
#include "RigidPlaneContactSolver.h"
#include "planeisometric/Planeisoparametric.h"

FEModel::FEModel()
    : solver_(std::make_unique<LinearSolver>()),
      displacements_(Eigen::VectorXd::Zero(0)),
      reactionForces_(Eigen::VectorXd::Zero(0)),
      contactForces_(Eigen::VectorXd::Zero(0)) {
}

void FEModel::setAssembly(std::shared_ptr<Assembly> assembly) {
    assembly_ = std::move(assembly);
    nodalDataCalculated_ = false;
}

void FEModel::setContactSolver(std::unique_ptr<IRigidPlaneContactSolver> contactSolver) {
    contactSolver_ = std::move(contactSolver);
}

void FEModel::configureRigidPlaneContact(const RigidPlane2D& plane,
    const std::vector<ContactFacet>& facets,
    double penaltyParameter) {
    configureRigidPlanePenaltyContact(plane, facets, penaltyParameter);
}

void FEModel::configureRigidPlanePenaltyContact(const RigidPlane2D& plane,
    const std::vector<ContactFacet>& facets,
    double penaltyParameter) {
    penaltyParameter_ = penaltyParameter;
    contactSolver_ = std::make_unique<RigidPlaneContactSolver>(assembly_, plane, penaltyParameter_);
    contactSolver_->setContactFacets(facets);
}

void FEModel::configureRigidPlaneAugmentedLagrangianContact(const RigidPlane2D& plane,
    const std::vector<ContactFacet>& facets,
    const AugmentedLagrangianSettings& settings) {
    augmentedLagrangianSettings_ = settings;
    contactSolver_ = std::make_unique<RigidPlaneAugmentedLagrangianContactSolver>(
        assembly_, plane, augmentedLagrangianSettings_);
    contactSolver_->setContactFacets(facets);
}

bool FEModel::solve() {
    if (!assembly_ || !validate()) {
        std::cerr << "FEModel: Invalid assembly or model configuration" << std::endl;
        return false;
    }

    nodalDataCalculated_ = false;
    performanceMetrics_ = {};
    iterationCount_ = 0;

    const auto totalStartTime = std::chrono::high_resolution_clock::now();

    try {
        const auto assemblyStartTime = std::chrono::high_resolution_clock::now();
        assembly_->assembleGlobalStiffnessMatrix(globalK_);
        assembleExternalForces(globalF_);
        assembly_->applyBoundaryConditions(globalK_, globalF_);
        const auto assemblyEndTime = std::chrono::high_resolution_clock::now();

        performanceMetrics_.assemblyTimeSeconds =
            std::chrono::duration<double>(assemblyEndTime - assemblyStartTime).count();
        performanceMetrics_.matrixNonZeros = globalK_.nonZeros();

        if (!solveLinearSystem()) {
            std::cerr << "FEModel: Linear system solution failed" << std::endl;
            return false;
        }

        calculateReactionForces();

        const auto totalEndTime = std::chrono::high_resolution_clock::now();
        performanceMetrics_.totalTimeSeconds =
            std::chrono::duration<double>(totalEndTime - totalStartTime).count();
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "FEModel: Solution error - " << e.what() << std::endl;
        return false;
    }
}

bool FEModel::solveContact() {
    if (!contactSolver_) {
        return solve();
    }

    if (!assembly_ || !validate()) {
        std::cerr << "FEModel: Invalid assembly for contact solution" << std::endl;
        return false;
    }

    nodalDataCalculated_ = false;
    performanceMetrics_ = {};
    performanceMetrics_.contactMethod = std::string(contactSolver_->getMethodName());
    iterationCount_ = 0;

    try {
        contactSolver_->resetState();
        return solveContactIterative();
    }
    catch (const std::exception& e) {
        std::cerr << "FEModel: Contact solution error - " << e.what() << std::endl;
        return false;
    }
}

RigidPlane2D FEModel::getContactPlane() const {
    if (!contactSolver_) {
        throw std::runtime_error("Contact plane requested, but no contact solver is configured");
    }
    return contactSolver_->getPlane();
}

std::vector<ContactFacet> FEModel::getContactFacets() const {
    if (!contactSolver_) {
        return {};
    }
    return contactSolver_->getContactFacets();
}

std::vector<ContactFacetResult> FEModel::getContactFacetResults() const {
    return evaluateCurrentContactState().facetResults;
}

ContactState FEModel::evaluateCurrentContactState() const {
    if (!contactSolver_ || !assembly_ ||
        displacements_.size() != assembly_->getTotalDofCount()) {
        return {};
    }

    Eigen::SparseMatrix<double> contactK;
    Eigen::VectorXd contactF;
    ContactState state;
    contactSolver_->assembleContact(displacements_, contactK, contactF, state);
    return state;
}

bool FEModel::solveLinearSystem() {
    try {
        if (!solver_) {
            solver_ = std::make_unique<LinearSolver>();
        }

        Eigen::VectorXd reducedDisplacements = solver_->solveSystem(globalK_, globalF_);
        const auto& solveStats = solver_->getLastSolveStats();

        performanceMetrics_.linearIterations = solveStats.iterations;
        performanceMetrics_.solveTimeSeconds += solveStats.solveTimeSeconds;
        performanceMetrics_.linearResidualEstimate = solveStats.estimatedError;
        performanceMetrics_.linearResidualNorm = solveStats.relativeResidualNorm;
        performanceMetrics_.linearSolveConverged = solveStats.converged;
        performanceMetrics_.usedDirectLinearSolver = solveStats.usedDirectSolver;
        performanceMetrics_.linearSolverBackend = solveStats.backendName;

        int totalDof = assembly_->getTotalDofCount();
        displacements_.resize(totalDof);
        displacements_.setZero();

        const auto& mapping = assembly_->getDofMapping();
        for (size_t i = 0; i < mapping.reducedToFull.size(); ++i) {
            if (static_cast<Eigen::Index>(i) < reducedDisplacements.size()) {
                int fullDof = mapping.reducedToFull[i];
                if (fullDof >= 0 && fullDof < totalDof) {
                    displacements_(fullDof) = reducedDisplacements(static_cast<Eigen::Index>(i));
                }
            }
        }

        for (size_t i = 0; i < mapping.prescribedDofs.size(); ++i) {
            if (i < mapping.prescribedValues.size()) {
                int prescribedDof = mapping.prescribedDofs[i];
                if (prescribedDof >= 0 && prescribedDof < totalDof) {
                    displacements_(prescribedDof) = mapping.prescribedValues[i];
                }
            }
        }

        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Linear system solution failed: " << e.what() << std::endl;
        return false;
    }
}

bool FEModel::solveContactIterative() {
    if (!contactSolver_) {
        return solve();
    }

    const auto totalStartTime = std::chrono::high_resolution_clock::now();
    const auto assemblyStartTime = std::chrono::high_resolution_clock::now();
    assembly_->assembleGlobalStiffnessMatrix(globalK_);
    assembleExternalForces(globalF_);
    assembly_->applyBoundaryConditions(globalK_, globalF_);
    const auto assemblyEndTime = std::chrono::high_resolution_clock::now();

    performanceMetrics_.assemblyTimeSeconds =
        std::chrono::duration<double>(assemblyEndTime - assemblyStartTime).count();
    performanceMetrics_.matrixNonZeros = globalK_.nonZeros();

    const Eigen::SparseMatrix<double> structuralReducedK = globalK_;
    const Eigen::VectorXd structuralReducedF = globalF_;

    if (!solveLinearSystem()) {
        return false;
    }

    int accumulatedLinearIterations = performanceMetrics_.linearIterations;
    double accumulatedAssemblyTime = performanceMetrics_.assemblyTimeSeconds;
    const bool useNestedAugmentedLagrangianLoop =
        contactSolver_->getMethod() == RigidPlaneContactMethod::AugmentedLagrangian;
    bool converged = false;

    if (!useNestedAugmentedLagrangianLoop) {
        Eigen::VectorXd previousDisplacements = displacements_;

        for (iterationCount_ = 1; iterationCount_ <= maxIterations_; ++iterationCount_) {
            const auto iterationAssemblyStartTime = std::chrono::high_resolution_clock::now();

            ContactIterationInfo contactInfo;
            applyContactConditions(previousDisplacements, contactInfo);

            Eigen::SparseMatrix<double> reducedContactK = contactInfo.stiffness;
            Eigen::VectorXd reducedContactF =
                contactInfo.force + contactInfo.stiffness * previousDisplacements;
            assembly_->applyBoundaryConditions(reducedContactK, reducedContactF);

            globalK_ = structuralReducedK;
            globalK_ += reducedContactK;
            globalK_.makeCompressed();
            globalF_ = structuralReducedF + reducedContactF;

            const auto iterationAssemblyEndTime = std::chrono::high_resolution_clock::now();
            accumulatedAssemblyTime +=
                std::chrono::duration<double>(
                    iterationAssemblyEndTime - iterationAssemblyStartTime).count();
            performanceMetrics_.assemblyTimeSeconds = accumulatedAssemblyTime;
            performanceMetrics_.matrixNonZeros = globalK_.nonZeros();
            performanceMetrics_.nonlinearIterations = iterationCount_;

            if (!solveLinearSystem()) {
                const auto failureTime = std::chrono::high_resolution_clock::now();
                performanceMetrics_.totalTimeSeconds =
                    std::chrono::duration<double>(failureTime - totalStartTime).count();
                return false;
            }
            accumulatedLinearIterations += solver_->getLastSolveStats().iterations;
            performanceMetrics_.linearIterations = accumulatedLinearIterations;

            ContactIterationInfo convergedStateInfo;
            convergedStateInfo.updateInfo = contactSolver_->updateState(displacements_);
            applyContactConditions(displacements_, convergedStateInfo);

            double relativeError =
                (displacements_ - previousDisplacements).norm() /
                (displacements_.norm() + 1.0e-15);
            bool activeSetStable =
                (convergedStateInfo.state.activeFacetIds == contactInfo.state.activeFacetIds);
            const bool contactStateConverged = convergedStateInfo.updateInfo.converged;

            performanceMetrics_.activeSetSize =
                static_cast<int>(convergedStateInfo.state.activeFacetIds.size());
            performanceMetrics_.activeContactGaussPoints =
                convergedStateInfo.updateInfo.activeGaussPointCount;
            performanceMetrics_.maxPenetration = convergedStateInfo.state.maxPenetration;
            performanceMetrics_.contactForceNorm = convergedStateInfo.state.contactForceNorm;
            performanceMetrics_.contactStateUpdateNorm =
                convergedStateInfo.updateInfo.stateUpdateNorm;
            performanceMetrics_.contactStateRelativeUpdateNorm =
                convergedStateInfo.updateInfo.relativeStateUpdateNorm;
            performanceMetrics_.maxNormalContactMultiplier =
                convergedStateInfo.updateInfo.maxNormalMultiplier;
            performanceMetrics_.meanNormalContactMultiplier =
                convergedStateInfo.updateInfo.meanNormalMultiplier;

            previousDisplacements = displacements_;

            if (relativeError < tolerance_ && activeSetStable && contactStateConverged) {
                converged = true;
                break;
            }
        }
    }
    else {
        const int maxEquilibriumIterations =
            std::max(6, std::min(24, maxIterations_));

        for (iterationCount_ = 1; iterationCount_ <= maxIterations_; ++iterationCount_) {
            bool equilibriumConverged = false;
            Eigen::VectorXd equilibriumPreviousDisplacements = displacements_;
            ContactIterationInfo convergedStateInfo;

            for (int equilibriumIteration = 1;
                equilibriumIteration <= maxEquilibriumIterations;
                ++equilibriumIteration) {
                const auto iterationAssemblyStartTime = std::chrono::high_resolution_clock::now();

                ContactIterationInfo contactInfo;
                applyContactConditions(equilibriumPreviousDisplacements, contactInfo);

                Eigen::SparseMatrix<double> reducedContactK = contactInfo.stiffness;
                Eigen::VectorXd reducedContactF =
                    contactInfo.force + contactInfo.stiffness * equilibriumPreviousDisplacements;
                assembly_->applyBoundaryConditions(reducedContactK, reducedContactF);

                globalK_ = structuralReducedK;
                globalK_ += reducedContactK;
                globalK_.makeCompressed();
                globalF_ = structuralReducedF + reducedContactF;

                const auto iterationAssemblyEndTime = std::chrono::high_resolution_clock::now();
                accumulatedAssemblyTime +=
                    std::chrono::duration<double>(
                        iterationAssemblyEndTime - iterationAssemblyStartTime).count();
                performanceMetrics_.assemblyTimeSeconds = accumulatedAssemblyTime;
                performanceMetrics_.matrixNonZeros = globalK_.nonZeros();

                if (!solveLinearSystem()) {
                    const auto failureTime = std::chrono::high_resolution_clock::now();
                    performanceMetrics_.totalTimeSeconds =
                        std::chrono::duration<double>(failureTime - totalStartTime).count();
                    return false;
                }
                accumulatedLinearIterations += solver_->getLastSolveStats().iterations;
                performanceMetrics_.linearIterations = accumulatedLinearIterations;

                applyContactConditions(displacements_, convergedStateInfo);

                const double relativeError =
                    (displacements_ - equilibriumPreviousDisplacements).norm() /
                    (displacements_.norm() + 1.0e-15);
                const bool activeSetStable =
                    (convergedStateInfo.state.activeFacetIds == contactInfo.state.activeFacetIds);

                performanceMetrics_.nonlinearIterations = iterationCount_;
                performanceMetrics_.activeSetSize =
                    static_cast<int>(convergedStateInfo.state.activeFacetIds.size());
                performanceMetrics_.activeContactGaussPoints =
                    convergedStateInfo.state.activeGaussPointCount;
                performanceMetrics_.maxPenetration = convergedStateInfo.state.maxPenetration;
                performanceMetrics_.contactForceNorm = convergedStateInfo.state.contactForceNorm;

                if (relativeError < tolerance_ && activeSetStable) {
                    equilibriumConverged = true;
                    break;
                }

                equilibriumPreviousDisplacements = displacements_;
            }

            if (!equilibriumConverged) {
                std::cerr << "FEModel: Augmented Lagrangian equilibrium solve did not converge "
                          << "during outer iteration " << iterationCount_ << std::endl;
                const auto failureTime = std::chrono::high_resolution_clock::now();
                performanceMetrics_.totalTimeSeconds =
                    std::chrono::duration<double>(failureTime - totalStartTime).count();
                return false;
            }

            convergedStateInfo.updateInfo = contactSolver_->updateState(displacements_);
            applyContactConditions(displacements_, convergedStateInfo);

            performanceMetrics_.nonlinearIterations = iterationCount_;
            performanceMetrics_.activeSetSize =
                static_cast<int>(convergedStateInfo.state.activeFacetIds.size());
            performanceMetrics_.activeContactGaussPoints =
                convergedStateInfo.updateInfo.activeGaussPointCount;
            performanceMetrics_.maxPenetration = convergedStateInfo.state.maxPenetration;
            performanceMetrics_.contactForceNorm = convergedStateInfo.state.contactForceNorm;
            performanceMetrics_.contactStateUpdateNorm =
                convergedStateInfo.updateInfo.stateUpdateNorm;
            performanceMetrics_.contactStateRelativeUpdateNorm =
                convergedStateInfo.updateInfo.relativeStateUpdateNorm;
            performanceMetrics_.maxNormalContactMultiplier =
                convergedStateInfo.updateInfo.maxNormalMultiplier;
            performanceMetrics_.meanNormalContactMultiplier =
                convergedStateInfo.updateInfo.meanNormalMultiplier;

            if (convergedStateInfo.updateInfo.converged) {
                converged = true;
                break;
            }
        }
    }

    if (!converged) {
        std::cerr << "FEModel: Contact solution did not converge after "
                  << maxIterations_ << " iterations" << std::endl;
        const auto failureTime = std::chrono::high_resolution_clock::now();
        performanceMetrics_.totalTimeSeconds =
            std::chrono::duration<double>(failureTime - totalStartTime).count();
        return false;
    }

    calculateReactionForces();
    const auto totalEndTime = std::chrono::high_resolution_clock::now();
    performanceMetrics_.totalTimeSeconds =
        std::chrono::duration<double>(totalEndTime - totalStartTime).count();
    return true;
}

void FEModel::calculateReactionForces() {
    int totalDof = assembly_->getTotalDofCount();
    reactionForces_.resize(totalDof);
    reactionForces_.setZero();

    Eigen::SparseMatrix<double> structuralK;
    Eigen::VectorXd fullF;
    assembly_->assembleGlobalStiffnessMatrix(structuralK);
    assembleExternalForces(fullF);
    contactForces_ = Eigen::VectorXd::Zero(totalDof);

    if (contactSolver_ && performanceMetrics_.nonlinearIterations > 0 &&
        displacements_.size() == totalDof) {
        Eigen::SparseMatrix<double> contactK;
        Eigen::VectorXd rawContactF;
        ContactState contactState;
        contactSolver_->assembleContact(displacements_, contactK, rawContactF, contactState);
        contactForces_ = -rawContactF;
        performanceMetrics_.contactForceNorm = contactState.contactForceNorm;
        reactionForces_ = structuralK * displacements_ - fullF - rawContactF;
    }
    else if (structuralK.rows() == displacements_.size() &&
        structuralK.cols() == displacements_.size()) {
        reactionForces_ = structuralK * displacements_ - fullF;
    }

    const auto& mapping = assembly_->getDofMapping();
    double freeDofResidualSquaredNorm = 0.0;
    for (int fullDof : mapping.reducedToFull) {
        if (fullDof >= 0 && fullDof < reactionForces_.size()) {
            freeDofResidualSquaredNorm += reactionForces_(fullDof) * reactionForces_(fullDof);
        }
    }
    performanceMetrics_.equilibriumResidualNorm = std::sqrt(freeDofResidualSquaredNorm);
}

void FEModel::applyContactConditions(const Eigen::VectorXd& trialDisplacements,
    ContactIterationInfo& iterationInfo) {
    if (!contactSolver_) {
        iterationInfo.stiffness.resize(assembly_->getTotalDofCount(), assembly_->getTotalDofCount());
        iterationInfo.stiffness.setZero();
        iterationInfo.force = Eigen::VectorXd::Zero(assembly_->getTotalDofCount());
        iterationInfo.state = {};
        return;
    }

    contactSolver_->assembleContact(trialDisplacements,
        iterationInfo.stiffness,
        iterationInfo.force,
        iterationInfo.state);
}

void FEModel::assembleExternalForces(Eigen::VectorXd& globalF) const {
    assembly_->assembleConcentratedForces(globalF);
    assembly_->assembleSurfaceLoads(globalF);
}

std::vector<std::shared_ptr<Node>> FEModel::getElementNodes(int elementId) const {
    auto element = assembly_->getElement(elementId);
    if (!element) {
        throw std::invalid_argument("Element not found: " + std::to_string(elementId));
    }

    std::vector<std::shared_ptr<Node>> elementNodes;
    elementNodes.reserve(element->getNodeIds().size());
    for (int nodeId : element->getNodeIds()) {
        auto node = assembly_->getNode(nodeId);
        if (!node) {
            throw std::invalid_argument("Node not found: " + std::to_string(nodeId));
        }
        elementNodes.push_back(node);
    }

    return elementNodes;
}

Eigen::Vector3d FEModel::getElementStress(int elementId, double xi, double eta) const {
    auto element = assembly_->getElement(elementId);
    if (!element) {
        return Eigen::Vector3d::Zero();
    }

    auto material = assembly_->getMaterial(element->getMaterialId());
    if (!material) {
        return Eigen::Vector3d::Zero();
    }

    auto dofIndices = assembly_->getElementDofIndices(elementId);
    Eigen::VectorXd elementDisplacements = Eigen::VectorXd::Zero(dofIndices.size());
    for (size_t i = 0; i < dofIndices.size(); ++i) {
        if (dofIndices[i] >= 0 && dofIndices[i] < displacements_.size()) {
            elementDisplacements(static_cast<Eigen::Index>(i)) = displacements_(dofIndices[i]);
        }
    }

    return element->computeStress(xi, eta, elementDisplacements, getElementNodes(elementId), material);
}

Eigen::Vector3d FEModel::getElementStrain(int elementId, double xi, double eta) const {
    auto element = assembly_->getElement(elementId);
    if (!element) {
        return Eigen::Vector3d::Zero();
    }

    auto material = assembly_->getMaterial(element->getMaterialId());
    if (!material) {
        return Eigen::Vector3d::Zero();
    }

    auto dofIndices = assembly_->getElementDofIndices(elementId);
    Eigen::VectorXd elementDisplacements = Eigen::VectorXd::Zero(dofIndices.size());
    for (size_t i = 0; i < dofIndices.size(); ++i) {
        if (dofIndices[i] >= 0 && dofIndices[i] < displacements_.size()) {
            elementDisplacements(static_cast<Eigen::Index>(i)) = displacements_(dofIndices[i]);
        }
    }

    return element->computeStrain(xi, eta, elementDisplacements, getElementNodes(elementId), material);
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

std::vector<Eigen::Vector2d> FEModel::getNodalReactionForces() const {
    if (!assembly_) {
        return {};
    }

    const auto& nodes = assembly_->getNodes();
    std::vector<Eigen::Vector2d> nodalReactionForces(
        nodes.size(), Eigen::Vector2d::Zero());

    for (size_t i = 0; i < nodes.size(); ++i) {
        const int nodeId = nodes[i]->getId();
        const int dofX = assembly_->getGlobalDofIndex(nodeId, 0);
        const int dofY = assembly_->getGlobalDofIndex(nodeId, 1);
        if (dofX >= 0 && dofX < reactionForces_.size()) {
            nodalReactionForces[i].x() = reactionForces_(dofX);
        }
        if (dofY >= 0 && dofY < reactionForces_.size()) {
            nodalReactionForces[i].y() = reactionForces_(dofY);
        }
    }

    return nodalReactionForces;
}

std::vector<Eigen::Vector2d> FEModel::getNodalContactForces() const {
    if (!assembly_) {
        return {};
    }

    const auto& nodes = assembly_->getNodes();
    std::vector<Eigen::Vector2d> nodalContactForces(
        nodes.size(), Eigen::Vector2d::Zero());

    for (size_t i = 0; i < nodes.size(); ++i) {
        const int nodeId = nodes[i]->getId();
        const int dofX = assembly_->getGlobalDofIndex(nodeId, 0);
        const int dofY = assembly_->getGlobalDofIndex(nodeId, 1);
        if (dofX >= 0 && dofX < contactForces_.size()) {
            nodalContactForces[i].x() = contactForces_(dofX);
        }
        if (dofY >= 0 && dofY < contactForces_.size()) {
            nodalContactForces[i].y() = contactForces_(dofY);
        }
    }

    return nodalContactForces;
}

std::vector<double> FEModel::getNodalContactSignedDistances() const {
    if (!assembly_) {
        return {};
    }

    const auto& nodes = assembly_->getNodes();
    std::vector<double> signedDistances(nodes.size(), 0.0);
    if (!contactSolver_) {
        return signedDistances;
    }

    const auto nodalDisplacements = getNodalDisplacements();
    const RigidPlane2D plane = contactSolver_->getPlane();
    for (size_t i = 0; i < nodes.size() && i < nodalDisplacements.size(); ++i) {
        signedDistances[i] =
            plane.signedDistance(nodes[i]->getCoordinates() + nodalDisplacements[i]);
    }

    return signedDistances;
}

std::vector<double> FEModel::getNodalContactPenetrations() const {
    const auto signedDistances = getNodalContactSignedDistances();
    std::vector<double> penetrations(signedDistances.size(), 0.0);
    for (size_t i = 0; i < signedDistances.size(); ++i) {
        penetrations[i] = std::max(0.0, -signedDistances[i]);
    }
    return penetrations;
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
    int nodeCount = static_cast<int>(nodes.size());

    nodalStresses_.assign(nodeCount, Eigen::Vector3d::Zero());
    nodalDisplacements_.assign(nodeCount, Eigen::Vector2d::Zero());
    nodalStrains_.assign(nodeCount, Eigen::Vector3d::Zero());

    std::vector<int> stressCount(nodeCount, 0);
    std::vector<int> strainCount(nodeCount, 0);

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

    auto elements = assembly_->getElements();
    for (const auto& element : elements) {
        auto material = assembly_->getMaterial(element->getMaterialId());
        if (!material) {
            continue;
        }

        std::vector<std::shared_ptr<Node>> elementNodes;
        for (int nodeId : element->getNodeIds()) {
            auto node = assembly_->getNode(nodeId);
            if (!node) {
                elementNodes.clear();
                break;
            }
            elementNodes.push_back(node);
        }

        if (elementNodes.empty()) {
            continue;
        }

        auto dofIndices = assembly_->getElementDofIndices(element->getId());
        Eigen::VectorXd elementDisplacements = Eigen::VectorXd::Zero(dofIndices.size());
        for (size_t i = 0; i < dofIndices.size(); ++i) {
            if (dofIndices[i] >= 0 && dofIndices[i] < displacements_.size()) {
                elementDisplacements(static_cast<Eigen::Index>(i)) = displacements_[dofIndices[i]];
            }
        }

        const double pos = 0.0;
        const std::vector<std::pair<double, double>> integrationPoints = {
            {-pos, -pos}, {pos, -pos}, {pos, pos}, {-pos, pos}
        };

        const auto& nodeIds = element->getNodeIds();
        for (size_t pt = 0; pt < integrationPoints.size() && pt < nodeIds.size(); ++pt) {
            Eigen::Vector3d stress = element->computeStress(
                integrationPoints[pt].first,
                integrationPoints[pt].second,
                elementDisplacements,
                elementNodes,
                material);
            Eigen::Vector3d strain = element->computeStrain(
                integrationPoints[pt].first,
                integrationPoints[pt].second,
                elementDisplacements,
                elementNodes,
                material);

            int targetNodeId = nodeIds[pt];
            for (int globalIdx = 0; globalIdx < nodeCount; ++globalIdx) {
                if (nodes[globalIdx]->getId() == targetNodeId) {
                    nodalStresses_[globalIdx] += stress;
                    nodalStrains_[globalIdx] += strain;
                    stressCount[globalIdx]++;
                    strainCount[globalIdx]++;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < nodeCount; ++i) {
        if (stressCount[i] > 0) {
            nodalStresses_[i] /= static_cast<double>(stressCount[i]);
        }
        if (strainCount[i] > 0) {
            nodalStrains_[i] /= static_cast<double>(strainCount[i]);
        }
    }

    nodalDataCalculated_ = true;
}
