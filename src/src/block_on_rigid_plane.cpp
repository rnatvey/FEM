#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>

#include "FEMModel.h"
#include "ResultFileExporter.h"
#include "assembly.h"
#include "material.h"
#include "meshgenerator.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kCoordinateTolerance = 1.0e-9;

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

        assembly->addPrescribedDisplacementY(node->getId(), -0.10);
        if (coordinates.x() < smallestAnchorX) {
            smallestAnchorX = coordinates.x();
            anchorNodeId = node->getId();
        }
    }

    if (anchorNodeId < 0) {
        throw std::runtime_error("Failed to select an anchor node on the top edge");
    }

    assembly->addFixedNode(anchorNodeId, true, false);
    return assembly;
}

double computeMinimumSignedDistanceToPlane(const FEModel& model,
    const std::shared_ptr<Assembly>& assembly,
    const RigidPlane2D& plane) {
    const auto& nodes = assembly->getNodes();
    const auto nodalDisplacements = model.getNodalDisplacements();

    double minimumSignedDistance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < nodes.size() && i < nodalDisplacements.size(); ++i) {
        const Eigen::Vector2d deformedPosition =
            nodes[i]->getCoordinates() + nodalDisplacements[i];
        minimumSignedDistance =
            std::min(minimumSignedDistance, plane.signedDistance(deformedPosition));
    }

    return minimumSignedDistance;
}

} // namespace

int main() {
    try {
        auto assembly = buildBlockAssembly();
        auto model = std::make_shared<FEModel>();
        auto meshGenerator = std::make_shared<MeshGenerator>(assembly);

        const RigidPlane2D rigidPlane{Eigen::Vector2d(0.0, 1.0), 0.0};
        const auto contactFacets = meshGenerator->collectExteriorFacets(1, true, 1.0e-8);
        if (contactFacets.empty()) {
            throw std::runtime_error("No contact facets were found on the lower boundary");
        }

        model->setAssembly(assembly);
        model->setMaxIterations(30);
        model->setSolverTolerance(1.0e-8);
        model->configureRigidPlaneContact(rigidPlane, contactFacets, 1.0e7);

        const bool success = model->solveContact();
        const auto& metrics = model->getPerformanceMetrics();
        const double minimumSignedDistance =
            computeMinimumSignedDistanceToPlane(*model, assembly, rigidPlane);

        std::cout << "Block on rigid plane example" << std::endl;
        std::cout << "success=" << std::boolalpha << success << std::endl;
        std::cout << "contact_facets=" << contactFacets.size() << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "active_set_size=" << metrics.activeSetSize << std::endl;
        std::cout << "max_penetration=" << metrics.maxPenetration << std::endl;
        std::cout << "minimum_signed_distance_to_plane=" << minimumSignedDistance << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        std::cout << "linear_solver_backend=" << metrics.linearSolverBackend << std::endl;
        std::cout << "linear_relative_residual=" << metrics.linearResidualNorm << std::endl;
        std::cout << "equilibrium_residual_norm=" << metrics.equilibriumResidualNorm << std::endl;
        std::cout << "contact_force_norm=" << metrics.contactForceNorm << std::endl;
        std::cout << "assembly_time_seconds=" << metrics.assemblyTimeSeconds << std::endl;
        std::cout << "solve_time_seconds=" << metrics.solveTimeSeconds << std::endl;
        std::cout << "matrix_nonzeros=" << metrics.matrixNonZeros << std::endl;

        if (success) {
            ResultFileExportOptions exportOptions;
            exportOptions.outputDirectory = std::filesystem::path("results") / "block_on_rigid_plane";
            exportOptions.baseName = "solution";
            exportOptions.extraStringMetrics = {
                {"case_name", "block_on_rigid_plane"}
            };
            exportOptions.extraNumericMetrics = {
                {"minimum_signed_distance", minimumSignedDistance},
                {"candidate_contact_facets", static_cast<double>(contactFacets.size())}
            };
            const auto exportArtifacts = ResultFileExporter::exportSolution(*model, exportOptions);
            std::cout << "results_vtu=" << exportArtifacts.vtuPath.string() << std::endl;
            std::cout << "results_metrics_json=" << exportArtifacts.metricsJsonPath.string() << std::endl;
        }

        return success ? 0 : 1;
    }
    catch (const std::exception& exception) {
        std::cerr << "Block on rigid plane example failed: " << exception.what() << std::endl;
        return 1;
    }
}
