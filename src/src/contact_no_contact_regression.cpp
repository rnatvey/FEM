#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

#include "FEMModel.h"
#include "assembly.h"
#include "material.h"
#include "meshgenerator.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kCoordinateTolerance = 1.0e-9;
constexpr double kSolutionTolerance = 1.0e-10;

std::shared_ptr<Assembly> buildBlockAssembly() {
    auto assembly = std::make_shared<Assembly>();
    auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
    assembly->addMaterial(material);

    constexpr int nodesX = 11;
    constexpr int nodesY = 11;
    const Eigen::Vector2d lowerCorner(-0.5, 0.05);
    const Eigen::Vector2d upperCorner(0.5, 1.05);

    std::vector<std::vector<int>> nodeGrid(nodesY, std::vector<int>(nodesX, 0));
    int nextNodeId = 1;
    for (int j = 0; j < nodesY; ++j) {
        const double eta = static_cast<double>(j) / static_cast<double>(nodesY - 1);
        const double y = lowerCorner.y() + eta * (upperCorner.y() - lowerCorner.y());
        for (int i = 0; i < nodesX; ++i) {
            const double xi = static_cast<double>(i) / static_cast<double>(nodesX - 1);
            const double x = lowerCorner.x() + xi * (upperCorner.x() - lowerCorner.x());
            assembly->addNode(std::make_shared<Node>(nextNodeId, x, y));
            nodeGrid[j][i] = nextNodeId++;
        }
    }

    int nextElementId = 1;
    for (int j = 0; j < nodesY - 1; ++j) {
        for (int i = 0; i < nodesX - 1; ++i) {
            assembly->addElement(std::make_shared<PlaneIsoparametricElement>(
                nextElementId++,
                std::vector<int>{
                    nodeGrid[j + 1][i + 1],
                    nodeGrid[j + 1][i],
                    nodeGrid[j][i],
                    nodeGrid[j][i + 1]},
                material->getId()));
        }
    }

    const auto& nodes = assembly->getNodes();
    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    int anchorNodeId = -1;
    double smallestAnchorX = std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        const Eigen::Vector2d coordinates = node->getCoordinates();
        if (std::abs(coordinates.y() - maxY) > kCoordinateTolerance) {
            continue;
        }

        assembly->addPrescribedDisplacementY(node->getId(), -0.02);
        if (coordinates.x() < smallestAnchorX) {
            smallestAnchorX = coordinates.x();
            anchorNodeId = node->getId();
        }
    }

    if (anchorNodeId < 0) {
        throw std::runtime_error("Failed to select an anchor node");
    }

    assembly->addFixedNode(anchorNodeId, true, false);
    return assembly;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

} // namespace

int main() {
    try {
        auto baselineAssembly = buildBlockAssembly();
        auto baselineModel = std::make_shared<FEModel>();
        baselineModel->setAssembly(baselineAssembly);
        baselineModel->setSolverTolerance(1.0e-8);
        require(baselineModel->solve(), "Baseline solve() failed");

        auto contactAssembly = buildBlockAssembly();
        auto contactMeshGenerator = std::make_shared<MeshGenerator>(contactAssembly);
        auto contactModel = std::make_shared<FEModel>();
        contactModel->setAssembly(contactAssembly);
        contactModel->setSolverTolerance(1.0e-8);
        contactModel->setMaxIterations(10);

        const auto contactFacets = contactMeshGenerator->collectExteriorFacets(1, true, 1.0e-8);
        require(!contactFacets.empty(), "No candidate contact facets found");

        const RigidPlane2D farPlane{Eigen::Vector2d(0.0, 1.0), -10.0};
        contactModel->configureRigidPlaneContact(farPlane, contactFacets, 1.0e7);
        require(contactModel->solveContact(), "solveContact() failed in no-contact scenario");

        const Eigen::VectorXd displacementDifference =
            baselineModel->getDisplacements() - contactModel->getDisplacements();
        const double differenceNorm = displacementDifference.norm();
        const auto& metrics = contactModel->getPerformanceMetrics();

        require(differenceNorm < kSolutionTolerance,
            "solveContact() does not match solve() when the plane is far away");
        require(metrics.activeSetSize == 0, "Active set must remain empty without contact");
        require(metrics.maxPenetration == 0.0, "Penetration must remain zero without contact");

        std::cout << "contact_no_contact_regression_passed" << std::endl;
        std::cout << "difference_norm=" << differenceNorm << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "active_set_size=" << metrics.activeSetSize << std::endl;
        std::cout << "linear_solver_backend=" << metrics.linearSolverBackend << std::endl;
        std::cout << "equilibrium_residual_norm=" << metrics.equilibriumResidualNorm << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "contact_no_contact_regression_failed: " << exception.what() << std::endl;
        return 1;
    }
}
