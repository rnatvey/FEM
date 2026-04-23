#include "FEMModel.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "ContactTypes.h"
#include "NeoHookeanMaterial.h"
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
        const auto structuralAssemblyStartTime = std::chrono::high_resolution_clock::now();
        assembly_->assembleGlobalStiffnessMatrix(globalK_);
        const Eigen::SparseMatrix<double> structuralFullK = globalK_;
        const auto structuralAssemblyEndTime = std::chrono::high_resolution_clock::now();
        performanceMetrics_.structuralAssemblyTimeSeconds =
            std::chrono::duration<double>(
                structuralAssemblyEndTime - structuralAssemblyStartTime).count();
        performanceMetrics_.structuralMatrixNonZeros = globalK_.nonZeros();
        assembleExternalForces(globalF_);
        const Eigen::VectorXd structuralFullF = globalF_;
        Assembly::BoundaryConditionApplicationStats boundaryConditionStats;
        assembly_->applyBoundaryConditions(globalK_, globalF_, &boundaryConditionStats);
        const auto assemblyEndTime = std::chrono::high_resolution_clock::now();

        performanceMetrics_.assemblyTimeSeconds =
            std::chrono::duration<double>(assemblyEndTime - assemblyStartTime).count();
        performanceMetrics_.boundaryConditionTimeSeconds =
            boundaryConditionStats.totalTimeSeconds;
        performanceMetrics_.reductionTimeSeconds =
            boundaryConditionStats.reduceSystemTimeSeconds;
        performanceMetrics_.matrixNonZeros = globalK_.nonZeros();

        if (!solveLinearSystem()) {
            std::cerr << "FEModel: Linear system solution failed" << std::endl;
            return false;
        }

        calculateReactionForcesFromInternalForce(structuralFullK * displacements_, structuralFullF);

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

bool FEModel::solveHyperelastic() {
    if (!assembly_ || !validate()) {
        std::cerr << "FEModel: Invalid assembly or model configuration" << std::endl;
        return false;
    }

    if (!assembly_->hasFiniteStrainModel()) {
        std::cerr << "FEModel: No finite-strain model configured for hyperelastic solve"
                  << std::endl;
        return false;
    }

    const bool useContact = (contactSolver_ != nullptr);
    const bool usePenaltyContact =
        useContact &&
        contactSolver_->getMethod() == RigidPlaneContactMethod::Penalty;
    const bool useAugmentedLagrangianContact =
        useContact &&
        contactSolver_->getMethod() == RigidPlaneContactMethod::AugmentedLagrangian;

    nodalDataCalculated_ = false;
    performanceMetrics_ = {};
    performanceMetrics_.loadSteps = std::max(1, hyperelasticLoadSteps_);
    if (useContact) {
        performanceMetrics_.contactMethod = std::string(contactSolver_->getMethodName());
        contactSolver_->resetState();
    }
    iterationCount_ = 0;

    const int totalDof = assembly_->getTotalDofCount();
    const Assembly::ConstraintData constraintData = assembly_->getConstraintData();
    assembly_->initializeDofMapping(totalDof);

    displacements_ = Eigen::VectorXd::Zero(totalDof);
    reactionForces_ = Eigen::VectorXd::Zero(totalDof);
    contactForces_ = Eigen::VectorXd::Zero(totalDof);

    for (const auto& [materialId, material] : assembly_->getFiniteStrainMaterials()) {
        (void)materialId;
        const auto neoHookeanMaterial =
            std::dynamic_pointer_cast<NeoHookeanMaterial>(material);
        if (!neoHookeanMaterial) {
            continue;
        }

        const double effectivePoissonsRatio =
            neoHookeanMaterial->getEffectivePoissonsRatio();
        const double bulkToShearRatio = neoHookeanMaterial->getBulkToShearRatio();

        performanceMetrics_.maxFiniteStrainEffectivePoissonsRatio =
            std::max(performanceMetrics_.maxFiniteStrainEffectivePoissonsRatio,
                effectivePoissonsRatio);
        performanceMetrics_.maxFiniteStrainBulkToShearRatio =
            std::max(performanceMetrics_.maxFiniteStrainBulkToShearRatio,
                bulkToShearRatio);

        if (neoHookeanMaterial->isNearlyIncompressible()) {
            performanceMetrics_.hasNearIncompressibleFiniteStrainMaterial = true;
        }
    }

    if (performanceMetrics_.hasNearIncompressibleFiniteStrainMaterial) {
        std::cerr
            << "FEModel: Hyperelastic baseline uses displacement-only fully integrated Q4 "
            << "with nearly incompressible Neo-Hookean material (effective nu up to "
            << performanceMetrics_.maxFiniteStrainEffectivePoissonsRatio
            << ", K/mu up to "
            << performanceMetrics_.maxFiniteStrainBulkToShearRatio
            << "). Volumetric locking may occur. "
            << "SRI and mixed u/p are not implemented in this baseline."
            << std::endl;
    }

    Eigen::VectorXd baseExternalForce;
    assembleExternalForces(baseExternalForce);

    const double residualTolerance = std::max(tolerance_, 1.0e-8);
    const double incrementTolerance = std::max(std::sqrt(tolerance_), 1.0e-6);

    Eigen::VectorXd internalForce = Eigen::VectorXd::Zero(totalDof);
    Eigen::VectorXd currentContactForce = Eigen::VectorXd::Zero(totalDof);
    Eigen::VectorXd finalExternalForce = baseExternalForce;
    bool converged = true;
    int cumulativeNonlinearIterations = 0;
    const auto totalStartTime = std::chrono::high_resolution_clock::now();

    const auto buildScaledConstraintData =
        [&constraintData](double loadFactor) -> Assembly::ConstraintData {
            Assembly::ConstraintData scaledConstraintData = constraintData;
            for (double& value : scaledConstraintData.prescribedValues) {
                value *= loadFactor;
            }
            return scaledConstraintData;
        };
    const auto applyBacktrackedFiniteStrainIncrement =
        [this](const Assembly::ConstraintData& activeConstraintData,
            const Eigen::VectorXd& displacementIncrement,
            const Eigen::VectorXd& stepExternalForce,
            double currentResidualNorm) -> bool {
            constexpr int kMaxBacktrackingIterations = 10;
            constexpr double kBacktrackingFactor = 0.5;

            double scale = 1.0;
            for (int lineSearchIteration = 0;
                lineSearchIteration <= kMaxBacktrackingIterations;
                ++lineSearchIteration) {
                Eigen::VectorXd trialDisplacements =
                    displacements_ + scale * displacementIncrement;
                applyPrescribedDisplacements(activeConstraintData, trialDisplacements);

                try {
                    const auto trialAssemblyStartTime =
                        std::chrono::high_resolution_clock::now();
                    Eigen::SparseMatrix<double> trialTangent;
                    Eigen::VectorXd trialInternalForce;
                    double trialStrainEnergy = 0.0;
                    assembly_->assembleFiniteStrainSystem(
                        trialDisplacements,
                        trialTangent,
                        trialInternalForce,
                        trialStrainEnergy,
                        nullptr,
                        false);
                    const auto trialAssemblyEndTime =
                        std::chrono::high_resolution_clock::now();
                    const double trialAssemblyTimeSeconds =
                        std::chrono::duration<double>(
                            trialAssemblyEndTime - trialAssemblyStartTime).count();
                    performanceMetrics_.structuralAssemblyTimeSeconds +=
                        trialAssemblyTimeSeconds;
                    performanceMetrics_.assemblyTimeSeconds += trialAssemblyTimeSeconds;

                    Eigen::VectorXd trialResidual = stepExternalForce - trialInternalForce;
                    if (contactSolver_) {
                        const auto trialContactAssemblyStartTime =
                            std::chrono::high_resolution_clock::now();
                        Eigen::SparseMatrix<double> trialContactStiffness;
                        Eigen::VectorXd trialContactForce;
                        ContactState trialContactState;
                        contactSolver_->assembleContact(
                            trialDisplacements,
                            trialContactStiffness,
                            trialContactForce,
                            trialContactState);
                        const auto trialContactAssemblyEndTime =
                            std::chrono::high_resolution_clock::now();
                        const double trialContactAssemblyTimeSeconds =
                            std::chrono::duration<double>(
                                trialContactAssemblyEndTime -
                                trialContactAssemblyStartTime).count();
                        performanceMetrics_.contactAssemblyTimeSeconds +=
                            trialContactAssemblyTimeSeconds;
                        performanceMetrics_.assemblyTimeSeconds +=
                            trialContactAssemblyTimeSeconds;
                        trialResidual += trialContactForce;
                    }

                    const double trialResidualNorm =
                        computeFreeDofNorm(trialResidual, activeConstraintData);
                    if (!(std::isfinite(trialResidualNorm) &&
                            trialResidualNorm < currentResidualNorm)) {
                        scale *= kBacktrackingFactor;
                        continue;
                    }

                    displacements_ = std::move(trialDisplacements);
                    return true;
                }
                catch (const std::exception& exception) {
                    const std::string_view message(exception.what());
                    const bool invalidDeformationGradient =
                        message.find("positive determinant") != std::string_view::npos ||
                        message.find("positive Jacobian") != std::string_view::npos;
                    if (!invalidDeformationGradient) {
                        throw;
                    }
                }

                scale *= kBacktrackingFactor;
            }

            return false;
        };

    std::vector<double> targetLoadFactors;
    targetLoadFactors.reserve(static_cast<size_t>(performanceMetrics_.loadSteps));
    for (int loadStep = 1; loadStep <= performanceMetrics_.loadSteps; ++loadStep) {
        targetLoadFactors.push_back(
            static_cast<double>(loadStep) / static_cast<double>(performanceMetrics_.loadSteps));
    }

    double convergedLoadFactor = 0.0;
    size_t targetLoadFactorIndex = 0;
    while (targetLoadFactorIndex < targetLoadFactors.size()) {
        const int loadStep = static_cast<int>(targetLoadFactorIndex) + 1;
        const double loadFactor = targetLoadFactors[targetLoadFactorIndex];
        const Assembly::ConstraintData stepConstraintData =
            buildScaledConstraintData(loadFactor);
        const Eigen::VectorXd stepExternalForce = baseExternalForce * loadFactor;
        finalExternalForce = stepExternalForce;
        const double referenceForceNorm =
            std::max(computeFreeDofNorm(stepExternalForce, stepConstraintData), 1.0);

        const Eigen::VectorXd convergedStepDisplacements = displacements_;
        applyPrescribedDisplacements(stepConstraintData, displacements_);

        bool stepConverged = false;
        bool retryableStepFailure = false;
        if (!useAugmentedLagrangianContact) {
            double lastRelativeIncrement = std::numeric_limits<double>::infinity();
            std::vector<int> previousActiveFacetIds;

            for (int stepIteration = 1; stepIteration <= maxIterations_; ++stepIteration) {
                iterationCount_ = ++cumulativeNonlinearIterations;
                const auto assemblyStartTime = std::chrono::high_resolution_clock::now();
                double totalStrainEnergy = 0.0;
                Eigen::SparseMatrix<double> structuralTangent;
                const auto structuralAssemblyStartTime = std::chrono::high_resolution_clock::now();
                try {
                    assembly_->assembleFiniteStrainSystem(
                        displacements_,
                        structuralTangent,
                        internalForce,
                        totalStrainEnergy);
                }
                catch (const std::exception& exception) {
                    const std::string_view message(exception.what());
                    const bool invalidDeformationGradient =
                        message.find("positive determinant") != std::string_view::npos ||
                        message.find("positive Jacobian") != std::string_view::npos;
                    if (!invalidDeformationGradient) {
                        throw;
                    }

                    retryableStepFailure = true;
                    break;
                }
                const auto structuralAssemblyEndTime = std::chrono::high_resolution_clock::now();
                performanceMetrics_.structuralAssemblyTimeSeconds +=
                    std::chrono::duration<double>(
                        structuralAssemblyEndTime - structuralAssemblyStartTime).count();
                performanceMetrics_.structuralMatrixNonZeros = structuralTangent.nonZeros();

                ContactIterationInfo contactInfo;
                bool activeSetStable = true;
                if (usePenaltyContact) {
                    applyContactConditions(displacements_, contactInfo);
                    currentContactForce = contactInfo.force;
                    globalK_ = structuralTangent;
                    globalK_ += contactInfo.stiffness;
                    globalK_.makeCompressed();
                    globalF_ = stepExternalForce + currentContactForce - internalForce;

                    activeSetStable =
                        stepIteration <= 1 ||
                        contactInfo.state.activeFacetIds == previousActiveFacetIds;
                    previousActiveFacetIds = contactInfo.state.activeFacetIds;

                    performanceMetrics_.activeSetSize =
                        static_cast<int>(contactInfo.state.activeFacetIds.size());
                    performanceMetrics_.activeContactGaussPoints =
                        contactInfo.state.activeGaussPointCount;
                    performanceMetrics_.maxPenetration = contactInfo.state.maxPenetration;
                    performanceMetrics_.contactForceNorm = contactInfo.state.contactForceNorm;
                    performanceMetrics_.contactStateUpdateNorm = 0.0;
                    performanceMetrics_.contactStateRelativeUpdateNorm = 0.0;
                    performanceMetrics_.contactStateMaxUpdateMagnitude = 0.0;
                    performanceMetrics_.contactStateRelativeMaxUpdate = 0.0;
                    performanceMetrics_.maxNormalContactMultiplier = 0.0;
                    performanceMetrics_.meanNormalContactMultiplier = 0.0;
                }
                else {
                    currentContactForce.setZero();
                    globalK_ = structuralTangent;
                    globalF_ = stepExternalForce - internalForce;
                    performanceMetrics_.activeSetSize = 0;
                    performanceMetrics_.activeContactGaussPoints = 0;
                    performanceMetrics_.maxPenetration = 0.0;
                    performanceMetrics_.contactForceNorm = 0.0;
                    performanceMetrics_.contactMatrixNonZeros = 0;
                }
                const auto assemblyEndTime = std::chrono::high_resolution_clock::now();

                performanceMetrics_.assemblyTimeSeconds +=
                    std::chrono::duration<double>(assemblyEndTime - assemblyStartTime).count();
                performanceMetrics_.matrixNonZeros = globalK_.nonZeros();
                performanceMetrics_.nonlinearIterations = iterationCount_;
                performanceMetrics_.strainEnergy = totalStrainEnergy;

                const double freeResidualNorm =
                    computeFreeDofNorm(globalF_, stepConstraintData);
                const double relativeResidualNorm = freeResidualNorm / referenceForceNorm;
                performanceMetrics_.equilibriumResidualNorm = freeResidualNorm;

                if (relativeResidualNorm < residualTolerance &&
                    lastRelativeIncrement < incrementTolerance &&
                    (!usePenaltyContact || activeSetStable)) {
                    stepConverged = true;
                    break;
                }

                Eigen::VectorXd displacementIncrement;
                if (!solveReducedIncrementSystem(
                        globalK_, globalF_, stepConstraintData, displacementIncrement)) {
                    retryableStepFailure = true;
                    break;
                }

                const double freeDisplacementNorm =
                    std::max(computeFreeDofNorm(displacements_, stepConstraintData), 1.0);
                const double freeIncrementNorm =
                    computeFreeDofNorm(displacementIncrement, stepConstraintData);
                lastRelativeIncrement = freeIncrementNorm / freeDisplacementNorm;

                if (!applyBacktrackedFiniteStrainIncrement(
                        stepConstraintData,
                        displacementIncrement,
                        stepExternalForce,
                        freeResidualNorm)) {
                    retryableStepFailure = true;
                    break;
                }
            }
        }
        else {
            const int maxEquilibriumIterations =
                std::max(6, std::min(24, maxIterations_));

            for (int outerIteration = 1; outerIteration <= maxIterations_; ++outerIteration) {
                iterationCount_ = ++cumulativeNonlinearIterations;
                bool equilibriumConverged = false;
                double equilibriumLastRelativeIncrement =
                    std::numeric_limits<double>::infinity();
                std::vector<int> previousActiveFacetIds;
                ContactIterationInfo convergedStateInfo;

                for (int equilibriumIteration = 1;
                    equilibriumIteration <= maxEquilibriumIterations;
                    ++equilibriumIteration) {
                    const auto assemblyStartTime = std::chrono::high_resolution_clock::now();
                    double totalStrainEnergy = 0.0;
                    Eigen::SparseMatrix<double> structuralTangent;
                    const auto structuralAssemblyStartTime =
                        std::chrono::high_resolution_clock::now();
                    try {
                        assembly_->assembleFiniteStrainSystem(
                            displacements_,
                            structuralTangent,
                            internalForce,
                            totalStrainEnergy);
                    }
                    catch (const std::exception& exception) {
                        const std::string_view message(exception.what());
                        const bool invalidDeformationGradient =
                            message.find("positive determinant") != std::string_view::npos ||
                            message.find("positive Jacobian") != std::string_view::npos;
                        if (!invalidDeformationGradient) {
                            throw;
                        }

                        retryableStepFailure = true;
                        break;
                    }
                    const auto structuralAssemblyEndTime =
                        std::chrono::high_resolution_clock::now();
                    performanceMetrics_.structuralAssemblyTimeSeconds +=
                        std::chrono::duration<double>(
                            structuralAssemblyEndTime - structuralAssemblyStartTime).count();
                    performanceMetrics_.structuralMatrixNonZeros = structuralTangent.nonZeros();

                    ContactIterationInfo contactInfo;
                    applyContactConditions(displacements_, contactInfo);
                    currentContactForce = contactInfo.force;
                    globalK_ = structuralTangent;
                    globalK_ += contactInfo.stiffness;
                    globalK_.makeCompressed();
                    globalF_ = stepExternalForce + currentContactForce - internalForce;
                    const auto assemblyEndTime = std::chrono::high_resolution_clock::now();

                    performanceMetrics_.assemblyTimeSeconds +=
                        std::chrono::duration<double>(assemblyEndTime - assemblyStartTime)
                            .count();
                    performanceMetrics_.matrixNonZeros = globalK_.nonZeros();
                    performanceMetrics_.nonlinearIterations = iterationCount_;
                    performanceMetrics_.strainEnergy = totalStrainEnergy;
                    performanceMetrics_.activeSetSize =
                        static_cast<int>(contactInfo.state.activeFacetIds.size());
                    performanceMetrics_.activeContactGaussPoints =
                        contactInfo.state.activeGaussPointCount;
                    performanceMetrics_.maxPenetration = contactInfo.state.maxPenetration;
                    performanceMetrics_.contactForceNorm = contactInfo.state.contactForceNorm;

                    const bool activeSetStable =
                        equilibriumIteration <= 1 ||
                        contactInfo.state.activeFacetIds == previousActiveFacetIds;
                    previousActiveFacetIds = contactInfo.state.activeFacetIds;

                    const double freeResidualNorm =
                        computeFreeDofNorm(globalF_, stepConstraintData);
                    const double relativeResidualNorm = freeResidualNorm / referenceForceNorm;
                    performanceMetrics_.equilibriumResidualNorm = freeResidualNorm;

                    if (relativeResidualNorm < residualTolerance &&
                        equilibriumLastRelativeIncrement < incrementTolerance &&
                        activeSetStable) {
                        equilibriumConverged = true;
                        break;
                    }

                    Eigen::VectorXd displacementIncrement;
                    if (!solveReducedIncrementSystem(
                            globalK_, globalF_, stepConstraintData, displacementIncrement)) {
                        retryableStepFailure = true;
                        break;
                    }

                    const double freeDisplacementNorm =
                        std::max(computeFreeDofNorm(displacements_, stepConstraintData), 1.0);
                    const double freeIncrementNorm =
                        computeFreeDofNorm(displacementIncrement, stepConstraintData);
                    equilibriumLastRelativeIncrement =
                        freeIncrementNorm / freeDisplacementNorm;

                    if (!applyBacktrackedFiniteStrainIncrement(
                            stepConstraintData,
                            displacementIncrement,
                            stepExternalForce,
                            freeResidualNorm)) {
                        retryableStepFailure = true;
                        break;
                    }
                }

                if (retryableStepFailure) {
                    break;
                }

                if (!equilibriumConverged) {
                    retryableStepFailure = true;
                    break;
                }

                convergedStateInfo.updateInfo = contactSolver_->updateState(displacements_);
                applyContactConditions(displacements_, convergedStateInfo);

                performanceMetrics_.nonlinearIterations = iterationCount_;
                performanceMetrics_.activeSetSize =
                    static_cast<int>(convergedStateInfo.state.activeFacetIds.size());
                performanceMetrics_.activeContactGaussPoints =
                    convergedStateInfo.updateInfo.activeGaussPointCount;
                performanceMetrics_.maxPenetration = convergedStateInfo.state.maxPenetration;
                performanceMetrics_.contactForceNorm =
                    convergedStateInfo.state.contactForceNorm;
                performanceMetrics_.contactStateUpdateNorm =
                    convergedStateInfo.updateInfo.stateUpdateNorm;
                performanceMetrics_.contactStateRelativeUpdateNorm =
                    convergedStateInfo.updateInfo.relativeStateUpdateNorm;
                performanceMetrics_.contactStateMaxUpdateMagnitude =
                    convergedStateInfo.updateInfo.maxStateUpdateMagnitude;
                performanceMetrics_.contactStateRelativeMaxUpdate =
                    convergedStateInfo.updateInfo.relativeMaxStateUpdate;
                performanceMetrics_.maxNormalContactMultiplier =
                    convergedStateInfo.updateInfo.maxNormalMultiplier;
                performanceMetrics_.meanNormalContactMultiplier =
                    convergedStateInfo.updateInfo.meanNormalMultiplier;

                if (convergedStateInfo.updateInfo.converged) {
                    stepConverged = true;
                    break;
                }
            }
        }

        if (!stepConverged) {
            retryableStepFailure = true;
            displacements_ = convergedStepDisplacements;
            const Assembly::ConstraintData convergedConstraintData =
                buildScaledConstraintData(convergedLoadFactor);
            applyPrescribedDisplacements(convergedConstraintData, displacements_);

            const double stepSize = loadFactor - convergedLoadFactor;
            const bool canSplitStep =
                retryableStepFailure &&
                stepSize > (1.0 / static_cast<double>(maxAdaptiveHyperelasticLoadSteps_)) &&
                targetLoadFactors.size() < static_cast<size_t>(maxAdaptiveHyperelasticLoadSteps_);

            if (!canSplitStep) {
                converged = false;
                break;
            }

            const double splitLoadFactor = convergedLoadFactor + 0.5 * stepSize;
            targetLoadFactors.insert(
                targetLoadFactors.begin() + static_cast<std::ptrdiff_t>(targetLoadFactorIndex),
                splitLoadFactor);
            performanceMetrics_.loadSteps = static_cast<int>(targetLoadFactors.size());

            std::cerr << "FEModel: Hyperelastic load step to factor "
                      << loadFactor << " failed; retrying with adaptive split at "
                      << splitLoadFactor << std::endl;
            continue;
        }

        convergedLoadFactor = loadFactor;
        performanceMetrics_.convergedLoadSteps = static_cast<int>(targetLoadFactorIndex) + 1;
        ++targetLoadFactorIndex;
    }

    if (!converged) {
        std::cerr << "FEModel: Hyperelastic "
                  << (useAugmentedLagrangianContact ? "augmented Lagrangian" : "Newton")
                  << " solve did not converge after "
                  << performanceMetrics_.convergedLoadSteps << " of "
                  << performanceMetrics_.loadSteps << " load steps"
                  << std::endl;
        if (performanceMetrics_.directLinearSolveCount > 0) {
            std::cerr << "FEModel: Direct sparse fallback was used in "
                      << performanceMetrics_.directLinearSolveCount
                      << " of " << performanceMetrics_.linearSolveCount
                      << " Newton linear solves before failure." << std::endl;
        }
        const auto failureTime = std::chrono::high_resolution_clock::now();
        performanceMetrics_.totalTimeSeconds =
            std::chrono::duration<double>(failureTime - totalStartTime).count();
        return false;
    }

    if (useContact) {
        Eigen::SparseMatrix<double> finalStructuralTangent;
        double finalStrainEnergy = 0.0;
        assembly_->assembleFiniteStrainSystem(
            displacements_,
            finalStructuralTangent,
            internalForce,
            finalStrainEnergy);

        ContactIterationInfo finalContactInfo;
        applyContactConditions(displacements_, finalContactInfo);
        currentContactForce = finalContactInfo.force;
        performanceMetrics_.strainEnergy = finalStrainEnergy;
        performanceMetrics_.activeSetSize =
            static_cast<int>(finalContactInfo.state.activeFacetIds.size());
        performanceMetrics_.activeContactGaussPoints =
            finalContactInfo.state.activeGaussPointCount;
        performanceMetrics_.maxPenetration = finalContactInfo.state.maxPenetration;
        performanceMetrics_.contactForceNorm = finalContactInfo.state.contactForceNorm;
        calculateReactionForcesFromInternalForce(
            internalForce, finalExternalForce, currentContactForce);
    }
    else {
        calculateReactionForcesFromInternalForce(internalForce, finalExternalForce);
    }
    if (performanceMetrics_.directLinearSolveCount > 0) {
        std::cerr << "FEModel: Hyperelastic solve used direct sparse fallback in "
                  << performanceMetrics_.directLinearSolveCount
                  << " of " << performanceMetrics_.linearSolveCount
                  << " Newton linear solves"
                  << (performanceMetrics_.hasNearIncompressibleFiniteStrainMaterial
                          ? " (consistent with near-incompressible conditioning)."
                          : ".")
                  << std::endl;
    }
    const auto totalEndTime = std::chrono::high_resolution_clock::now();
    performanceMetrics_.totalTimeSeconds =
        std::chrono::duration<double>(totalEndTime - totalStartTime).count();
    return true;
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

std::vector<FiniteStrainElementResponse> FEModel::evaluateCurrentFiniteStrainElementResponses() const {
    if (!assembly_ || !assembly_->hasFiniteStrainModel() ||
        displacements_.size() != assembly_->getTotalDofCount()) {
        return {};
    }

    Eigen::SparseMatrix<double> tangent;
    Eigen::VectorXd internalForce;
    double totalStrainEnergy = 0.0;
    std::vector<FiniteStrainElementResponse> responses;
    assembly_->assembleFiniteStrainSystem(
        displacements_,
        tangent,
        internalForce,
        totalStrainEnergy,
        &responses,
        false);
    return responses;
}

void FEModel::accumulateLinearSolveStats(const LinearSolver::SolveStats& solveStats) {
    performanceMetrics_.linearSolveCount += 1;
    performanceMetrics_.linearIterations += solveStats.iterations;
    performanceMetrics_.solveTimeSeconds += solveStats.solveTimeSeconds;
    performanceMetrics_.linearResidualEstimate = solveStats.estimatedError;
    performanceMetrics_.linearResidualNorm = solveStats.relativeResidualNorm;
    performanceMetrics_.linearSolveConverged = solveStats.converged;

    if (solveStats.usedDirectSolver) {
        performanceMetrics_.usedDirectLinearSolver = true;
        performanceMetrics_.directLinearSolveCount += 1;
    }

    if (performanceMetrics_.linearSolverBackend == "uninitialized") {
        performanceMetrics_.linearSolverBackend = solveStats.backendName;
    }
    else if (performanceMetrics_.linearSolverBackend != solveStats.backendName) {
        performanceMetrics_.linearSolverBackend = "mixed";
    }
}

bool FEModel::solveLinearSystem() {
    try {
        if (!solver_) {
            solver_ = std::make_unique<LinearSolver>();
        }

        Eigen::VectorXd reducedDisplacements = solver_->solveSystem(globalK_, globalF_);
        const auto& solveStats = solver_->getLastSolveStats();
        accumulateLinearSolveStats(solveStats);

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

bool FEModel::solveReducedIncrementSystem(const Eigen::SparseMatrix<double>& tangent,
    const Eigen::VectorXd& residual,
    const Assembly::ConstraintData& constraintData,
    Eigen::VectorXd& displacementIncrement) {
    try {
        if (!solver_) {
            solver_ = std::make_unique<LinearSolver>();
        }

        Eigen::SparseMatrix<double> reducedTangent;
        Eigen::VectorXd reducedResidual;
        std::vector<int> activeDofs;
        const auto reduceStartTime = std::chrono::high_resolution_clock::now();
        solver_->reduceSystem(tangent,
            residual,
            constraintData.constrainedDofs,
            reducedTangent,
            reducedResidual,
            activeDofs);
        const auto reduceEndTime = std::chrono::high_resolution_clock::now();
        const double reductionTimeSeconds =
            std::chrono::duration<double>(reduceEndTime - reduceStartTime).count();
        performanceMetrics_.reductionTimeSeconds += reductionTimeSeconds;
        performanceMetrics_.boundaryConditionTimeSeconds += reductionTimeSeconds;

        const Eigen::VectorXd reducedIncrement =
            solver_->solveSystem(reducedTangent, reducedResidual);
        const auto& solveStats = solver_->getLastSolveStats();
        accumulateLinearSolveStats(solveStats);

        displacementIncrement = Eigen::VectorXd::Zero(tangent.rows());
        solver_->expandSolution(
            reducedIncrement,
            constraintData.constrainedDofs,
            activeDofs,
            displacementIncrement);
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Hyperelastic Newton increment solve failed: " << e.what() << std::endl;
        return false;
    }
}

bool FEModel::solveContactIterative() {
    if (!contactSolver_) {
        return solve();
    }

    const auto totalStartTime = std::chrono::high_resolution_clock::now();
    const auto assemblyStartTime = std::chrono::high_resolution_clock::now();
    const auto structuralAssemblyStartTime = std::chrono::high_resolution_clock::now();
    assembly_->assembleGlobalStiffnessMatrix(globalK_);
    const Eigen::SparseMatrix<double> structuralFullK = globalK_;
    const auto structuralAssemblyEndTime = std::chrono::high_resolution_clock::now();
    performanceMetrics_.structuralAssemblyTimeSeconds =
        std::chrono::duration<double>(
            structuralAssemblyEndTime - structuralAssemblyStartTime).count();
    performanceMetrics_.structuralMatrixNonZeros = globalK_.nonZeros();
    assembleExternalForces(globalF_);
    const Eigen::VectorXd structuralFullF = globalF_;
    Assembly::BoundaryConditionApplicationStats boundaryConditionStats;
    assembly_->applyBoundaryConditions(globalK_, globalF_, &boundaryConditionStats);
    const auto assemblyEndTime = std::chrono::high_resolution_clock::now();

    performanceMetrics_.assemblyTimeSeconds =
        std::chrono::duration<double>(assemblyEndTime - assemblyStartTime).count();
    performanceMetrics_.boundaryConditionTimeSeconds =
        boundaryConditionStats.totalTimeSeconds;
    performanceMetrics_.reductionTimeSeconds =
        boundaryConditionStats.reduceSystemTimeSeconds;
    performanceMetrics_.matrixNonZeros = globalK_.nonZeros();

    const Eigen::SparseMatrix<double> structuralReducedK = globalK_;
    const Eigen::VectorXd structuralReducedF = globalF_;

    if (!solveLinearSystem()) {
        return false;
    }

    double accumulatedAssemblyTime = performanceMetrics_.assemblyTimeSeconds;
    const bool useNestedAugmentedLagrangianLoop =
        contactSolver_->getMethod() == RigidPlaneContactMethod::AugmentedLagrangian;
    bool converged = false;
    ContactIterationInfo finalContactInfo;
    bool haveFinalContactInfo = false;

    if (!useNestedAugmentedLagrangianLoop) {
        Eigen::VectorXd previousDisplacements = displacements_;

        for (iterationCount_ = 1; iterationCount_ <= maxIterations_; ++iterationCount_) {
            const auto iterationAssemblyStartTime = std::chrono::high_resolution_clock::now();

            ContactIterationInfo contactInfo;
            applyContactConditions(previousDisplacements, contactInfo);

            Eigen::SparseMatrix<double> reducedContactK = contactInfo.stiffness;
            Eigen::VectorXd reducedContactF =
                contactInfo.force + contactInfo.stiffness * previousDisplacements;
            Assembly::BoundaryConditionApplicationStats contactBoundaryConditionStats;
            assembly_->applyBoundaryConditions(
                reducedContactK, reducedContactF, &contactBoundaryConditionStats);
            performanceMetrics_.boundaryConditionTimeSeconds +=
                contactBoundaryConditionStats.totalTimeSeconds;
            performanceMetrics_.reductionTimeSeconds +=
                contactBoundaryConditionStats.reduceSystemTimeSeconds;
            performanceMetrics_.contactMatrixNonZeros = reducedContactK.nonZeros();

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
            performanceMetrics_.contactStateMaxUpdateMagnitude =
                convergedStateInfo.updateInfo.maxStateUpdateMagnitude;
            performanceMetrics_.contactStateRelativeMaxUpdate =
                convergedStateInfo.updateInfo.relativeMaxStateUpdate;
            performanceMetrics_.maxNormalContactMultiplier =
                convergedStateInfo.updateInfo.maxNormalMultiplier;
            performanceMetrics_.meanNormalContactMultiplier =
                convergedStateInfo.updateInfo.meanNormalMultiplier;
            finalContactInfo = convergedStateInfo;
            haveFinalContactInfo = true;

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
                Assembly::BoundaryConditionApplicationStats contactBoundaryConditionStats;
                assembly_->applyBoundaryConditions(
                    reducedContactK, reducedContactF, &contactBoundaryConditionStats);
                performanceMetrics_.boundaryConditionTimeSeconds +=
                    contactBoundaryConditionStats.totalTimeSeconds;
                performanceMetrics_.reductionTimeSeconds +=
                    contactBoundaryConditionStats.reduceSystemTimeSeconds;
                performanceMetrics_.contactMatrixNonZeros = reducedContactK.nonZeros();

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
            performanceMetrics_.contactStateMaxUpdateMagnitude =
                convergedStateInfo.updateInfo.maxStateUpdateMagnitude;
            performanceMetrics_.contactStateRelativeMaxUpdate =
                convergedStateInfo.updateInfo.relativeMaxStateUpdate;
            performanceMetrics_.maxNormalContactMultiplier =
                convergedStateInfo.updateInfo.maxNormalMultiplier;
            performanceMetrics_.meanNormalContactMultiplier =
                convergedStateInfo.updateInfo.meanNormalMultiplier;
            finalContactInfo = convergedStateInfo;
            haveFinalContactInfo = true;

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

    if (haveFinalContactInfo) {
        calculateReactionForcesFromInternalForce(
            structuralFullK * displacements_, structuralFullF, finalContactInfo.force);
    }
    else {
        calculateReactionForcesFromInternalForce(structuralFullK * displacements_, structuralFullF);
    }
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

void FEModel::calculateReactionForcesFromInternalForce(const Eigen::VectorXd& internalForce,
    const Eigen::VectorXd& externalForce) {
    reactionForces_ = internalForce - externalForce;
    contactForces_ = Eigen::VectorXd::Zero(reactionForces_.size());

    const auto& mapping = assembly_->getDofMapping();
    double freeDofResidualSquaredNorm = 0.0;
    for (int fullDof : mapping.reducedToFull) {
        if (fullDof >= 0 && fullDof < reactionForces_.size()) {
            freeDofResidualSquaredNorm += reactionForces_(fullDof) * reactionForces_(fullDof);
        }
    }
    performanceMetrics_.equilibriumResidualNorm = std::sqrt(freeDofResidualSquaredNorm);
}

void FEModel::calculateReactionForcesFromInternalForce(const Eigen::VectorXd& internalForce,
    const Eigen::VectorXd& externalForce,
    const Eigen::VectorXd& contactForce) {
    reactionForces_ = internalForce - externalForce - contactForce;
    contactForces_ = -contactForce;

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
        performanceMetrics_.contactMatrixNonZeros = 0;
        return;
    }

    const auto contactAssemblyStartTime = std::chrono::high_resolution_clock::now();
    contactSolver_->assembleContact(trialDisplacements,
        iterationInfo.stiffness,
        iterationInfo.force,
        iterationInfo.state);
    const auto contactAssemblyEndTime = std::chrono::high_resolution_clock::now();
    performanceMetrics_.contactAssemblyTimeSeconds +=
        std::chrono::duration<double>(
            contactAssemblyEndTime - contactAssemblyStartTime).count();
    performanceMetrics_.contactMatrixNonZeros = iterationInfo.stiffness.nonZeros();
}

void FEModel::assembleExternalForces(Eigen::VectorXd& globalF) const {
    assembly_->assembleConcentratedForces(globalF);
    assembly_->assembleSurfaceLoads(globalF);
}

void FEModel::applyPrescribedDisplacements(
    const Assembly::ConstraintData& constraintData,
    Eigen::VectorXd& fullDisplacements) const {
    for (int dof : constraintData.fixedZeroDofs) {
        if (dof >= 0 && dof < fullDisplacements.size()) {
            fullDisplacements(dof) = 0.0;
        }
    }

    for (size_t i = 0; i < constraintData.prescribedDofs.size(); ++i) {
        if (i >= constraintData.prescribedValues.size()) {
            break;
        }

        const int dof = constraintData.prescribedDofs[i];
        if (dof >= 0 && dof < fullDisplacements.size()) {
            fullDisplacements(dof) = constraintData.prescribedValues[i];
        }
    }
}

double FEModel::computeFreeDofNorm(const Eigen::VectorXd& fullVector,
    const Assembly::ConstraintData& constraintData) const {
    if (fullVector.size() == 0) {
        return 0.0;
    }

    std::vector<char> isConstrained(static_cast<size_t>(fullVector.size()), 0);
    for (int dof : constraintData.constrainedDofs) {
        if (dof >= 0 && dof < fullVector.size()) {
            isConstrained[static_cast<size_t>(dof)] = 1;
        }
    }

    double squaredNorm = 0.0;
    for (Eigen::Index dof = 0; dof < fullVector.size(); ++dof) {
        if (!isConstrained[static_cast<size_t>(dof)]) {
            squaredNorm += fullVector(dof) * fullVector(dof);
        }
    }

    return std::sqrt(squaredNorm);
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
    std::unordered_map<int, int> nodeIdToGlobalIndex;
    nodeIdToGlobalIndex.reserve(static_cast<size_t>(nodeCount));
    for (int i = 0; i < nodeCount; ++i) {
        nodeIdToGlobalIndex[nodes[static_cast<size_t>(i)]->getId()] = i;
    }

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
        elementNodes.reserve(element->getNodeIds().size());
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

        // For the linear Q4 path these values are true nodal evaluations in the
        // element parent space, not duplicated center-point samples.
        constexpr std::array<std::pair<double, double>, 4> nodalEvaluationPoints = {{
            {1.0, 1.0},
            {-1.0, 1.0},
            {-1.0, -1.0},
            {1.0, -1.0}
        }};

        const auto& nodeIds = element->getNodeIds();
        for (size_t pt = 0; pt < nodalEvaluationPoints.size() && pt < nodeIds.size(); ++pt) {
            Eigen::Vector3d stress = element->computeStress(
                nodalEvaluationPoints[pt].first,
                nodalEvaluationPoints[pt].second,
                elementDisplacements,
                elementNodes,
                material);
            Eigen::Vector3d strain = element->computeStrain(
                nodalEvaluationPoints[pt].first,
                nodalEvaluationPoints[pt].second,
                elementDisplacements,
                elementNodes,
                material);

            const int targetNodeId = nodeIds[pt];
            const auto targetNodeIt = nodeIdToGlobalIndex.find(targetNodeId);
            if (targetNodeIt == nodeIdToGlobalIndex.end()) {
                continue;
            }
            const int globalIdx = targetNodeIt->second;
            nodalStresses_[globalIdx] += stress;
            nodalStrains_[globalIdx] += strain;
            stressCount[globalIdx]++;
            strainCount[globalIdx]++;
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
