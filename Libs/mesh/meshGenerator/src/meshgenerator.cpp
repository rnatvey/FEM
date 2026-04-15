#include "meshgenerator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "ConcentratedForce.h"
#include "constants.h"
#include "loadFunction.h"
#include "planeisometric/Planeisoparametric.h"

using namespace Constants;

MeshGenerator::MeshGenerator(std::shared_ptr<Assembly> assembly)
    : assembly_(std::move(assembly)) {
    if (!assembly_) {
        throw std::invalid_argument("Assembly cannot be null");
    }
}

void MeshGenerator::addBlock(const Block& block) {
    if (block.edges.size() != 4) {
        throw std::invalid_argument("Block must have exactly 4 edges");
    }

    if (block.nodesX < 2 || block.nodesY < 2) {
        throw std::invalid_argument("Block must have at least 2 nodes in each direction");
    }

    generateBlockMesh(block);
}

void MeshGenerator::createRectangle(const Eigen::Vector2d& corner1, const Eigen::Vector2d& corner2,
    int nodesX, int nodesY, int materialId) {
    Eigen::Vector2d p1 = corner1;
    Eigen::Vector2d p2(corner2.x(), corner1.y());
    Eigen::Vector2d p3 = corner2;
    Eigen::Vector2d p4(corner1.x(), corner2.y());

    Block block;
    block.edges = {
        Geometry::createLine(p1, p2),
        Geometry::createLine(p2, p3),
        Geometry::createLine(p3, p4),
        Geometry::createLine(p4, p1)
    };
    block.nodesX = nodesX;
    block.nodesY = nodesY;
    block.materialId = materialId;

    addBlock(block);
}

void MeshGenerator::createAnnulus(const Eigen::Vector2d& center, double innerRadius, double outerRadius,
    double startAngle, double endAngle,
    int radialLayers, int circumferentialNodes, int materialId) {
    Block block;
    block.edges.push_back(Geometry::createArc(center, innerRadius, startAngle, endAngle));
    block.edges.push_back(Geometry::createArc(center, outerRadius, startAngle, endAngle));
    block.edges.push_back(Geometry::createLine(
        center + innerRadius * Eigen::Vector2d(std::cos(startAngle), std::sin(startAngle)),
        center + outerRadius * Eigen::Vector2d(std::cos(startAngle), std::sin(startAngle))));
    block.edges.push_back(Geometry::createLine(
        center + innerRadius * Eigen::Vector2d(std::cos(endAngle), std::sin(endAngle)),
        center + outerRadius * Eigen::Vector2d(std::cos(endAngle), std::sin(endAngle))));
    block.nodesX = circumferentialNodes;
    block.nodesY = radialLayers;
    block.materialId = materialId;

    addBlock(block);
}

void MeshGenerator::createTriangle(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
    const Eigen::Vector2d& p3, int divisions, int materialId) {
    (void)p1;
    (void)p2;
    (void)p3;
    (void)divisions;
    (void)materialId;
    std::cout << "Triangle meshing is not implemented for Q4-only workflow" << std::endl;
}

Eigen::Vector2d MeshGenerator::transfiniteInterpolation(const Block& block, double xi, double eta) const {
    Eigen::Vector2d p1 = block.edges[0](0.0);
    Eigen::Vector2d p2 = block.edges[0](1.0);
    Eigen::Vector2d p3 = block.edges[2](1.0);
    Eigen::Vector2d p4 = block.edges[2](0.0);

    Eigen::Vector2d c1 = block.edges[0](xi);
    Eigen::Vector2d c3 = block.edges[2](xi);
    Eigen::Vector2d c2 = block.edges[1](eta);
    Eigen::Vector2d c4 = block.edges[3](eta);

    return (1.0 - eta) * c1 + eta * c3 +
        (1.0 - xi) * c4 + xi * c2 -
        ((1.0 - xi) * (1.0 - eta) * p1 +
            xi * (1.0 - eta) * p2 +
            xi * eta * p3 +
            (1.0 - xi) * eta * p4);
}

void MeshGenerator::generateBlockMesh(const Block& block) {
    std::vector<std::vector<int>> nodeGrid(block.nodesY, std::vector<int>(block.nodesX, 0));

    std::vector<Eigen::Vector2d> bottomEdgePoints = Geometry::discretizeCurve(block.edges[0], block.nodesX);
    for (int i = 0; i < block.nodesX; ++i) {
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId, bottomEdgePoints[i].x(), bottomEdgePoints[i].y()));
        nodeGrid[0][i] = nodeId;
    }

    std::vector<Eigen::Vector2d> topEdgePoints = Geometry::discretizeCurve(block.edges[2], block.nodesX);
    for (int i = 0; i < block.nodesX; ++i) {
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId, topEdgePoints[i].x(), topEdgePoints[i].y()));
        nodeGrid[block.nodesY - 1][i] = nodeId;
    }

    std::vector<Eigen::Vector2d> leftEdgePoints = Geometry::discretizeCurve(block.edges[3], block.nodesY);
    for (int j = 1; j < block.nodesY - 1; ++j) {
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId, leftEdgePoints[j].x(), leftEdgePoints[j].y()));
        nodeGrid[j][0] = nodeId;
    }

    std::vector<Eigen::Vector2d> rightEdgePoints = Geometry::discretizeCurve(block.edges[1], block.nodesY);
    for (int j = 1; j < block.nodesY - 1; ++j) {
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId, rightEdgePoints[j].x(), rightEdgePoints[j].y()));
        nodeGrid[j][block.nodesX - 1] = nodeId;
    }

    for (int j = 1; j < block.nodesY - 1; ++j) {
        double eta = static_cast<double>(j) / static_cast<double>(block.nodesY - 1);
        Eigen::Vector2d leftPoint = block.edges[3](eta);
        Eigen::Vector2d rightPoint = block.edges[1](eta);

        for (int i = 1; i < block.nodesX - 1; ++i) {
            double xi = static_cast<double>(i) / static_cast<double>(block.nodesX - 1);
            Eigen::Vector2d coords = leftPoint * (1.0 - xi) + rightPoint * xi;
            int nodeId = nextNodeId_++;
            assembly_->addNode(std::make_shared<Node>(nodeId, coords.x(), coords.y()));
            nodeGrid[j][i] = nodeId;
        }
    }

    for (int j = 0; j < block.nodesY - 1; ++j) {
        for (int i = 0; i < block.nodesX - 1; ++i) {
            std::vector<int> nodeIds = {
                nodeGrid[j][i],
                nodeGrid[j][i + 1],
                nodeGrid[j + 1][i + 1],
                nodeGrid[j + 1][i]
            };

            assembly_->addElement(std::make_shared<PlaneIsoparametricElement>(
                nextElementId_++, nodeIds, block.materialId));
        }
    }
}

void MeshGenerator::createAnnulusSimple(const Eigen::Vector2d& center,
    double innerRadius, double outerRadius,
    double startAngle, double endAngle,
    int radialLayers, int circumferentialNodes,
    int materialId) {
    createAnnulusGraded(center, innerRadius, outerRadius,
        startAngle, endAngle,
        radialLayers, circumferentialNodes,
        materialId, AnnulusGrading{});
}

void MeshGenerator::createAnnulusGraded(const Eigen::Vector2d& center,
    double innerRadius, double outerRadius,
    double startAngle, double endAngle,
    int radialLayers, int circumferentialNodes,
    int materialId,
    const AnnulusGrading& grading) {
    if (radialLayers < 2 || circumferentialNodes < 2) {
        throw std::invalid_argument("Invalid annulus mesh size");
    }
    if (innerRadius >= outerRadius) {
        throw std::invalid_argument("Inner radius must be smaller than outer radius");
    }

    double startRad = (std::abs(startAngle) > TWO_PI) ? startAngle * DEG_TO_RAD : startAngle;
    double endRad = (std::abs(endAngle) > TWO_PI) ? endAngle * DEG_TO_RAD : endAngle;
    if (endRad <= startRad) {
        endRad += TWO_PI;
    }

    std::vector<double> angularParams;
    if (grading.useAngularBias && grading.contactHalfAngle > 0.0) {
        double centerAngle = normalizeAngleToSweep(grading.contactCenterAngle, startRad, endRad);
        double halfAngle = std::abs(grading.contactHalfAngle);
        angularParams = buildDensityMappedParameters(
            circumferentialNodes,
            [=](double s) {
                double angle = startRad + s * (endRad - startRad);
                double distance = std::abs(angle - centerAngle);
                double normalizedDistance = distance / std::max(halfAngle, 1.0e-12);
                return 1.0 + grading.angularBiasStrength * std::exp(-4.0 * normalizedDistance * normalizedDistance);
            });
    }
    else {
        angularParams = buildDensityMappedParameters(circumferentialNodes, [](double) { return 1.0; });
    }

    std::vector<double> radialParams;
    if (grading.useRadialBias) {
        radialParams = buildDensityMappedParameters(
            radialLayers,
            [=](double s) {
                return 1.0 + grading.radialBiasToOuterStrength * s * s;
            });
    }
    else {
        radialParams = buildDensityMappedParameters(radialLayers, [](double) { return 1.0; });
    }

    std::vector<double> angles(circumferentialNodes);
    for (int i = 0; i < circumferentialNodes; ++i) {
        angles[i] = startRad + angularParams[i] * (endRad - startRad);
    }

    std::vector<double> radii(radialLayers);
    for (int i = 0; i < radialLayers; ++i) {
        radii[i] = innerRadius + radialParams[i] * (outerRadius - innerRadius);
    }

    std::vector<std::vector<int>> nodeGrid(radialLayers, std::vector<int>(circumferentialNodes, 0));
    for (int layer = 0; layer < radialLayers; ++layer) {
        for (int node = 0; node < circumferentialNodes; ++node) {
            double x = center.x() + radii[layer] * std::cos(angles[node]);
            double y = center.y() + radii[layer] * std::sin(angles[node]);

            int nodeId = nextNodeId_++;
            assembly_->addNode(std::make_shared<Node>(nodeId, x, y));
            nodeGrid[layer][node] = nodeId;
        }
    }

    for (int layer = 0; layer < radialLayers - 1; ++layer) {
        for (int segment = 0; segment < circumferentialNodes - 1; ++segment) {
            std::vector<int> nodeIds = {
                nodeGrid[layer][segment + 1],
                nodeGrid[layer][segment],
                nodeGrid[layer + 1][segment],
                nodeGrid[layer + 1][segment + 1]
            };

            assembly_->addElement(std::make_shared<PlaneIsoparametricElement>(
                nextElementId_++, nodeIds, materialId));
        }
    }
}

std::vector<int> MeshGenerator::findContactNodes(double contactCenterX,
    double contactHalfWidth,
    double maxYtolerance) const {
    std::vector<int> contactNodeIds;
    const auto& nodes = assembly_->getNodes();
    if (nodes.empty()) {
        return contactNodeIds;
    }

    double minY = nodes.front()->getCoordinates().y();
    for (const auto& node : nodes) {
        minY = std::min(minY, node->getCoordinates().y());
    }

    for (const auto& node : nodes) {
        Eigen::Vector2d coords = node->getCoordinates();
        if (std::abs(coords.x() - contactCenterX) <= contactHalfWidth &&
            std::abs(coords.y() - minY) <= maxYtolerance) {
            contactNodeIds.push_back(node->getId());
        }
    }

    return contactNodeIds;
}

void MeshGenerator::applyParabolicContactToNodes(double maxPressure,
    double contactHalfWidth,
    double contactCenterX,
    double totalForce) {
    auto parabolicLoad = LoadFunction::parabolicPressure(maxPressure, contactHalfWidth, contactCenterX);
    auto contactNodeIds = findContactNodes(contactCenterX, contactHalfWidth, 0.01);

    if (contactNodeIds.empty()) {
        return;
    }

    std::vector<double> nodalPressures;
    nodalPressures.reserve(contactNodeIds.size());
    double integratedPressure = 0.0;

    for (int nodeId : contactNodeIds) {
        auto node = assembly_->getNode(nodeId);
        Eigen::Vector2d pressure = parabolicLoad.distribution_(node->getCoordinates().x(),
            node->getCoordinates().y(), Eigen::Vector2d(0.0, -1.0));
        double scalarPressure = std::abs(pressure.y());
        nodalPressures.push_back(scalarPressure);
        integratedPressure += scalarPressure;
    }

    if (integratedPressure <= 0.0) {
        return;
    }

    double scaleFactor = totalForce / integratedPressure;
    for (size_t i = 0; i < contactNodeIds.size(); ++i) {
        assembly_->addConcentratedForce(std::make_shared<ConcentratedForce>(
            contactNodeIds[i], 0.0, -nodalPressures[i] * scaleFactor));
    }
}

std::vector<double> MeshGenerator::buildDensityMappedParameters(int pointCount,
    const std::function<double(double)>& densityFunction) {
    if (pointCount < 2) {
        throw std::invalid_argument("At least two points are required");
    }

    constexpr int sampleCount = 2048;
    const double step = 1.0 / static_cast<double>(sampleCount - 1);

    std::vector<double> cumulative(sampleCount, 0.0);
    double previousDensity = std::max(1.0e-12, densityFunction(0.0));
    for (int i = 1; i < sampleCount; ++i) {
        double s = i * step;
        double density = std::max(1.0e-12, densityFunction(s));
        cumulative[i] = cumulative[i - 1] + 0.5 * (previousDensity + density) * step;
        previousDensity = density;
    }

    const double totalWeight = cumulative.back();
    for (double& value : cumulative) {
        value /= totalWeight;
    }

    std::vector<double> parameters(pointCount, 0.0);
    for (int pointIndex = 1; pointIndex < pointCount - 1; ++pointIndex) {
        double target = static_cast<double>(pointIndex) / static_cast<double>(pointCount - 1);
        auto upper = std::lower_bound(cumulative.begin(), cumulative.end(), target);
        int upperIndex = static_cast<int>(std::distance(cumulative.begin(), upper));
        int lowerIndex = std::max(0, upperIndex - 1);

        double lowerCdf = cumulative[lowerIndex];
        double upperCdf = cumulative[upperIndex];
        double alpha = (upperCdf > lowerCdf) ? (target - lowerCdf) / (upperCdf - lowerCdf) : 0.0;
        parameters[pointIndex] = (static_cast<double>(lowerIndex) + alpha) * step;
    }

    parameters.front() = 0.0;
    parameters.back() = 1.0;
    return parameters;
}

double MeshGenerator::normalizeAngleToSweep(double angle, double startAngle, double endAngle) {
    double normalized = angle;
    while (normalized < startAngle) {
        normalized += TWO_PI;
    }
    while (normalized > endAngle) {
        normalized -= TWO_PI;
    }
    if (normalized < startAngle) {
        normalized = startAngle;
    }
    if (normalized > endAngle) {
        normalized = endAngle;
    }
    return normalized;
}
