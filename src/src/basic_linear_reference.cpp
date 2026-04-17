#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>

#include "FEMModel.h"
#include "ResultFileExporter.h"
#include "assembly.h"
#include "material.h"
#include "node.h"
#include "planeisometric/Planeisoparametric.h"

namespace {

constexpr double kCoordinateTolerance = 1.0e-9;
constexpr double kPrescribedTopEdgeDy = -0.02;
constexpr double kMaxPrescribedDisplacementErrorTolerance = 1.0e-12;
constexpr double kEquilibriumResidualTolerance = 1.0e-5;

struct LinearReferenceCase {
    std::shared_ptr<Assembly> assembly;
    int anchorNodeId = -1;
    int topEdgeNodeCount = 0;
};

LinearReferenceCase buildLinearReferenceCase() {
    LinearReferenceCase linearCase;
    linearCase.assembly = std::make_shared<Assembly>();

    auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
    linearCase.assembly->addMaterial(material);

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
            linearCase.assembly->addNode(std::make_shared<Node>(nextNodeId, x, y));
            nodeGrid[j][i] = nextNodeId++;
        }
    }

    int nextElementId = 1;
    for (int j = 0; j < nodesY - 1; ++j) {
        for (int i = 0; i < nodesX - 1; ++i) {
            linearCase.assembly->addElement(std::make_shared<PlaneIsoparametricElement>(
                nextElementId++,
                std::vector<int>{
                    nodeGrid[j + 1][i + 1],
                    nodeGrid[j + 1][i],
                    nodeGrid[j][i],
                    nodeGrid[j][i + 1]},
                material->getId()));
        }
    }

    const auto& nodes = linearCase.assembly->getNodes();
    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    double smallestAnchorX = std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        const Eigen::Vector2d coordinates = node->getCoordinates();
        if (std::abs(coordinates.y() - maxY) > kCoordinateTolerance) {
            continue;
        }

        linearCase.assembly->addPrescribedDisplacementY(node->getId(), kPrescribedTopEdgeDy);
        ++linearCase.topEdgeNodeCount;
        if (coordinates.x() < smallestAnchorX) {
            smallestAnchorX = coordinates.x();
            linearCase.anchorNodeId = node->getId();
        }
    }

    if (linearCase.anchorNodeId < 0) {
        throw std::runtime_error("Failed to select an anchor node on the top edge");
    }

    linearCase.assembly->addFixedNode(linearCase.anchorNodeId, true, false);
    return linearCase;
}

double computeMaximumTopEdgeDisplacementError(const FEModel& model,
    const std::shared_ptr<Assembly>& assembly) {
    const auto& nodes = assembly->getNodes();
    const auto nodalDisplacements = model.getNodalDisplacements();

    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& node : nodes) {
        maxY = std::max(maxY, node->getCoordinates().y());
    }

    double maxError = 0.0;
    for (size_t i = 0; i < nodes.size() && i < nodalDisplacements.size(); ++i) {
        if (std::abs(nodes[i]->getCoordinates().y() - maxY) > kCoordinateTolerance) {
            continue;
        }
        maxError = std::max(maxError, std::abs(nodalDisplacements[i].y() - kPrescribedTopEdgeDy));
    }

    return maxError;
}

Eigen::Vector2d computeReactionForceSum(const FEModel& model) {
    Eigen::Vector2d reactionSum = Eigen::Vector2d::Zero();
    for (const auto& nodalReaction : model.getNodalReactionForces()) {
        reactionSum += nodalReaction;
    }
    return reactionSum;
}

} // namespace

int main() {
    try {
        const auto linearCase = buildLinearReferenceCase();

        FEModel model;
        model.setAssembly(linearCase.assembly);
        model.setSolverTolerance(1.0e-8);

        const bool success = model.solve();
        if (!success) {
            std::cerr << "Basic linear reference solve failed" << std::endl;
            return 1;
        }

        const auto& metrics = model.getPerformanceMetrics();
        const double maxTopEdgeDisplacementError =
            computeMaximumTopEdgeDisplacementError(model, linearCase.assembly);
        const Eigen::Vector2d reactionForceSum = computeReactionForceSum(model);

        if (maxTopEdgeDisplacementError > kMaxPrescribedDisplacementErrorTolerance) {
            std::cerr << "Top-edge prescribed displacement mismatch: "
                      << maxTopEdgeDisplacementError << std::endl;
            return 1;
        }
        if (metrics.equilibriumResidualNorm > kEquilibriumResidualTolerance) {
            std::cerr << "Equilibrium residual is above the acceptance threshold: "
                      << metrics.equilibriumResidualNorm << std::endl;
            return 1;
        }

        ResultFileExportOptions exportOptions;
        exportOptions.outputDirectory = std::filesystem::path("results") / "basic_linear_reference";
        exportOptions.baseName = "solution";
        exportOptions.extraStringMetrics = {
            {"case_name", "basic_linear_reference"},
            {"verification_category", "linear_baseline"}
        };
        exportOptions.extraNumericMetrics = {
            {"prescribed_top_edge_dy", kPrescribedTopEdgeDy},
            {"top_edge_node_count", static_cast<double>(linearCase.topEdgeNodeCount)},
            {"anchor_node_id", static_cast<double>(linearCase.anchorNodeId)},
            {"maximum_prescribed_displacement_error", maxTopEdgeDisplacementError},
            {"reaction_force_sum_x", reactionForceSum.x()},
            {"reaction_force_sum_y", reactionForceSum.y()}
        };

        const auto exportArtifacts = ResultFileExporter::exportSolution(model, exportOptions);

        std::cout << "Basic linear reference scenario" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "nodes=" << linearCase.assembly->getNodeCount() << std::endl;
        std::cout << "elements=" << linearCase.assembly->getElementCount() << std::endl;
        std::cout << "top_edge_node_count=" << linearCase.topEdgeNodeCount << std::endl;
        std::cout << "maximum_prescribed_displacement_error="
                  << maxTopEdgeDisplacementError << std::endl;
        std::cout << "reaction_force_sum_x=" << reactionForceSum.x() << std::endl;
        std::cout << "reaction_force_sum_y=" << reactionForceSum.y() << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        std::cout << "equilibrium_residual_norm=" << metrics.equilibriumResidualNorm << std::endl;
        std::cout << "assembly_time_seconds=" << metrics.assemblyTimeSeconds << std::endl;
        std::cout << "solve_time_seconds=" << metrics.solveTimeSeconds << std::endl;
        std::cout << "matrix_nonzeros=" << metrics.matrixNonZeros << std::endl;
        std::cout << "results_vtu=" << exportArtifacts.vtuPath.string() << std::endl;
        std::cout << "results_metrics_json=" << exportArtifacts.metricsJsonPath.string() << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "basic_linear_reference_failed: " << exception.what() << std::endl;
        return 1;
    }
}
