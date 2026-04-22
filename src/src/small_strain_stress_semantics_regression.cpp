#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>

#include "ConcentratedForce.h"
#include "FEMModel.h"
#include "assembly.h"
#include "material.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kStressMatchTolerance = 1.0e-10;
constexpr double kNonUniformStressTolerance = 1.0e-4;

std::shared_ptr<Assembly> buildSingleElementStressGradientCase() {
    auto assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<Material>(1, 2.1e5, 0.30, 1.0);
    assembly->addMaterial(material);

    assembly->addNode(std::make_shared<Node>(1, 1.10, 1.00));  // top-right
    assembly->addNode(std::make_shared<Node>(2, -0.15, 0.95)); // top-left
    assembly->addNode(std::make_shared<Node>(3, 0.00, 0.00));  // bottom-left
    assembly->addNode(std::make_shared<Node>(4, 1.00, -0.05)); // bottom-right

    assembly->addElement(std::make_shared<PlaneIsoparametricElement>(
        1,
        std::vector<int>{1, 2, 3, 4},
        material->getId()));

    assembly->addFixedNode(3, true, true);
    assembly->addFixedNode(4, false, true);
    assembly->addFixedNode(2, true, false);

    assembly->addConcentratedForce(
        std::make_shared<ConcentratedForce>(1, 1250.0, -750.0));
    assembly->addConcentratedForce(
        std::make_shared<ConcentratedForce>(2, 0.0, 300.0));

    return assembly;
}

double stressMismatchNorm(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs) {
    return (lhs - rhs).norm();
}

} // namespace

int main() {
    try {
        auto assembly = buildSingleElementStressGradientCase();

        FEModel model;
        model.setAssembly(assembly);
        model.setSolverTolerance(1.0e-10);

        if (!model.solve()) {
            std::cerr << "Small-strain stress semantics regression solve failed" << std::endl;
            return 1;
        }

        const auto nodalStresses = model.getNodalStresses();
        if (nodalStresses.size() != 4) {
            std::cerr << "Unexpected nodal stress count: " << nodalStresses.size() << std::endl;
            return 1;
        }

        const auto element = assembly->getElement(1);
        const auto material = assembly->getMaterial(element->getMaterialId());
        const auto elementNodes = element->getElementNodes(assembly->getNodes());
        const auto& fullDisplacements = model.getDisplacements();

        Eigen::VectorXd elementDisplacements(8);
        elementDisplacements.setZero();
        for (size_t nodeIndex = 0; nodeIndex < element->getNodeIds().size(); ++nodeIndex) {
            const int nodeId = element->getNodeIds()[nodeIndex];
            const int dofX = assembly->getGlobalDofIndex(nodeId, 0);
            const int dofY = assembly->getGlobalDofIndex(nodeId, 1);
            elementDisplacements(static_cast<Eigen::Index>(2 * nodeIndex)) = fullDisplacements[dofX];
            elementDisplacements(static_cast<Eigen::Index>(2 * nodeIndex + 1)) = fullDisplacements[dofY];
        }

        constexpr std::array<std::pair<double, double>, 4> localNodeCoordinates = {{
            {1.0, 1.0},
            {-1.0, 1.0},
            {-1.0, -1.0},
            {1.0, -1.0}
        }};

        std::array<Eigen::Vector3d, 4> directNodalStresses;
        for (size_t i = 0; i < localNodeCoordinates.size(); ++i) {
            directNodalStresses[i] = element->computeStress(
                localNodeCoordinates[i].first,
                localNodeCoordinates[i].second,
                elementDisplacements,
                elementNodes,
                material);

            const double mismatch = stressMismatchNorm(nodalStresses[i], directNodalStresses[i]);
            const double referenceNorm = std::max(1.0, directNodalStresses[i].norm());
            if (mismatch > kStressMatchTolerance * referenceNorm) {
                std::cerr << "Recovered nodal stress does not match true nodal evaluation at local node "
                          << i << ": mismatch=" << mismatch << std::endl;
                return 1;
            }
        }

        double maxPairwiseDifference = 0.0;
        for (size_t i = 0; i < directNodalStresses.size(); ++i) {
            for (size_t j = i + 1; j < directNodalStresses.size(); ++j) {
                maxPairwiseDifference =
                    std::max(maxPairwiseDifference,
                        stressMismatchNorm(directNodalStresses[i], directNodalStresses[j]));
            }
        }

        if (maxPairwiseDifference <= kNonUniformStressTolerance) {
            std::cerr << "Regression case is not stress-varying enough; max pairwise difference="
                      << maxPairwiseDifference << std::endl;
            return 1;
        }

        const Eigen::Vector3d centerStress = model.getElementStress(1, 0.0, 0.0);
        double maxCenterDifference = 0.0;
        for (const auto& nodalStress : nodalStresses) {
            maxCenterDifference = std::max(
                maxCenterDifference,
                stressMismatchNorm(nodalStress, centerStress));
        }

        if (maxCenterDifference <= kNonUniformStressTolerance) {
            std::cerr << "Nodal stresses still collapse to the element center value; max center difference="
                      << maxCenterDifference << std::endl;
            return 1;
        }

        std::cout << "Small-strain stress semantics regression" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "max_pairwise_nodal_stress_difference=" << maxPairwiseDifference << std::endl;
        std::cout << "max_center_to_nodal_stress_difference=" << maxCenterDifference << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "small_strain_stress_semantics_regression_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
