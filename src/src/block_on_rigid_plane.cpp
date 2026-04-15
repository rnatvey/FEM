#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include "FEMModel.h"
#include "RigidPlaneContactSolver.h"
#include "assembly.h"
#include "material.h"
#include "meshgenerator.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kCoordinateTolerance = 1.0e-9;

} // namespace

int main() {
    try {
        auto assembly = std::make_shared<Assembly>();
        auto model = std::make_shared<FEModel>();
        auto meshGenerator = std::make_shared<MeshGenerator>(assembly);

        auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
        assembly->addMaterial(material);

        const Eigen::Vector2d lowerCorner(-0.5, 0.05);
        const Eigen::Vector2d upperCorner(0.5, 1.05);
        constexpr int nodesX = 11;
        constexpr int nodesY = 11;
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

            assembly->addPrescribedDisplacementY(node->getId(), -0.10);
            if (coordinates.x() < smallestAnchorX) {
                smallestAnchorX = coordinates.x();
                anchorNodeId = node->getId();
            }
        }

        if (anchorNodeId < 0) {
            std::cerr << "Failed to select an anchor node on the top edge" << std::endl;
            return 1;
        }

        assembly->addFixedNode(anchorNodeId, true, false);

        auto contactFacets = meshGenerator->collectExteriorFacets(1, true, 1.0e-8);
        if (contactFacets.empty()) {
            std::cerr << "No contact facets were found on the lower boundary" << std::endl;
            return 1;
        }

        model->setAssembly(assembly);
        model->setMaxIterations(30);
        model->setSolverTolerance(1.0e-8);

        const bool linearSuccess = model->solve();
        std::cout << "precontact_linear_success=" << std::boolalpha << linearSuccess << std::endl;

        RigidPlaneContactSolver precheckSolver(
            assembly,
            RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0},
            1.0e7);
        precheckSolver.setContactFacets(contactFacets);

        Eigen::SparseMatrix<double> precontactK;
        Eigen::VectorXd precontactF;
        ContactState precontactState;
        precheckSolver.assembleContact(
            model->getDisplacements(),
            precontactK,
            precontactF,
            precontactState);

        std::cout << "precontact_active_set_size=" << precontactState.activeFacetIds.size() << std::endl;
        std::cout << "precontact_max_penetration=" << precontactState.maxPenetration << std::endl;

        model->configureRigidPlaneContact(
            RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0},
            contactFacets,
            1.0e7);

        const bool success = model->solveContact();
        const auto& metrics = model->getPerformanceMetrics();

        std::cout << "Block on rigid plane verification" << std::endl;
        std::cout << "success=" << std::boolalpha << success << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "active_set_size=" << metrics.activeSetSize << std::endl;
        std::cout << "max_penetration=" << metrics.maxPenetration << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        std::cout << "matrix_nonzeros=" << metrics.matrixNonZeros << std::endl;

        const auto nodalDisplacements = model->getNodalDisplacements();
        double minimumSignedDistance = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < nodes.size() && i < nodalDisplacements.size(); ++i) {
            const Eigen::Vector2d deformedPosition =
                nodes[i]->getCoordinates() + nodalDisplacements[i];
            minimumSignedDistance = std::min(minimumSignedDistance, deformedPosition.y());
        }
        std::cout << "minimum_signed_distance_to_plane=" << minimumSignedDistance << std::endl;

        return success ? 0 : 1;
    }
    catch (const std::exception& exception) {
        std::cerr << "Verification case failed: " << exception.what() << std::endl;
        return 1;
    }
}
