#include <iostream>

#include "HyperelasticScenarioSupport.h"

namespace {

constexpr int kNodesX = 4;
constexpr int kNodesY = 4;
constexpr double kYoungsModulus = 100.0;
constexpr double kPoissonsRatio = 0.30;
constexpr double kThickness = 1.0;
constexpr double kPrescribedTopEdgeDy = -0.02;
constexpr double kAugmentationParameter = 1.0e5;

} // namespace

int main() {
    try {
        const Eigen::Vector2d lowerCorner(-0.5, 0.05);
        const Eigen::Vector2d upperCorner(0.5, 1.05);

        const auto hyperelasticCase = HyperelasticScenarioSupport::buildFiniteStrainBlockCase(
            kNodesX,
            kNodesY,
            lowerCorner,
            upperCorner,
            kYoungsModulus,
            kPoissonsRatio,
            kThickness,
            kPrescribedTopEdgeDy);

        FEModel model;
        model.setAssembly(hyperelasticCase.assembly);
        model.setSolverTolerance(1.0e-8);
        model.setMaxIterations(10);
        AugmentedLagrangianSettings settings;
        settings.augmentationParameter = kAugmentationParameter;
        model.configureRigidPlaneAugmentedLagrangianContact(
            RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0},
            {ContactFacet{1, 0}},
            settings);

        if (model.solveHyperelastic()) {
            std::cerr << "Hyperelastic contact guard failed: solveHyperelastic() unexpectedly "
                         "accepted augmented-Lagrangian contact configuration"
                      << std::endl;
            return 1;
        }

        std::cout << "Hyperelastic contact guard" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "expected_behavior=hyperelastic_augmented_lagrangian_contact_not_supported_yet"
                  << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "hyperelastic_contact_guard_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
