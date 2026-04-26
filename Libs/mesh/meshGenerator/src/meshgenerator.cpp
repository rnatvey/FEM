#include "meshgenerator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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

    const std::vector<double> uniformAngularParams =
        buildDensityMappedParameters(circumferentialNodes, [](double) { return 1.0; });
    std::vector<double> biasedAngularParams;
    if (grading.useAngularBias && grading.contactHalfAngle > 0.0) {
        double centerAngle = normalizeAngleToSweep(grading.contactCenterAngle, startRad, endRad);
        double halfAngle = std::abs(grading.contactHalfAngle);
        biasedAngularParams = buildDensityMappedParameters(
            circumferentialNodes,
            [=](double s) {
                double angle = startRad + s * (endRad - startRad);
                double distance = std::abs(angle - centerAngle);
                double normalizedDistance = distance / std::max(halfAngle, 1.0e-12);
                return 1.0 + grading.angularBiasStrength * std::exp(-4.0 * normalizedDistance * normalizedDistance);
            });
    }
    else {
        biasedAngularParams = uniformAngularParams;
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

    std::vector<double> radii(radialLayers);
    for (int i = 0; i < radialLayers; ++i) {
        radii[i] = innerRadius + radialParams[i] * (outerRadius - innerRadius);
    }

    std::vector<std::vector<double>> anglesByLayer(
        radialLayers,
        std::vector<double>(circumferentialNodes, startRad));
    for (int layer = 0; layer < radialLayers; ++layer) {
        const double physicalRadialCoordinate =
            (radii[layer] - innerRadius) / std::max(outerRadius - innerRadius, 1.0e-15);
        const double angularBiasBlend = grading.localizeAngularBiasToOuterSurface
            ? std::clamp(
                std::pow(std::clamp(physicalRadialCoordinate, 0.0, 1.0),
                    std::max(grading.angularBiasOuterLocalizationPower, 1.0)),
                0.0, 1.0)
            : 1.0;

        for (int node = 0; node < circumferentialNodes; ++node) {
            const double blendedParameter =
                uniformAngularParams[node] +
                angularBiasBlend * (biasedAngularParams[node] - uniformAngularParams[node]);
            anglesByLayer[layer][node] = startRad + blendedParameter * (endRad - startRad);
        }
    }

    std::vector<std::vector<int>> nodeGrid(radialLayers, std::vector<int>(circumferentialNodes, 0));
    for (int layer = 0; layer < radialLayers; ++layer) {
        for (int node = 0; node < circumferentialNodes; ++node) {
            const double angle = anglesByLayer[layer][node];
            double x = center.x() + radii[layer] * std::cos(angle);
            double y = center.y() + radii[layer] * std::sin(angle);

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

MeshGenerator::RingMeshDiagnostics MeshGenerator::generateTireRingGraded(
    const Eigen::Vector2d& center,
    double innerRadius,
    double outerRadius,
    const RingMeshControl& control) {
    AnnulusGrading grading;
    grading.useAngularBias = control.useAngularBias;
    grading.useRadialBias = control.useRadialBias;
    grading.localizeAngularBiasToOuterSurface = control.localizeAngularBiasToOuterSurface;
    grading.contactCenterAngle = control.contactCenterAngle;
    grading.contactHalfAngle = control.contactHalfAngle;
    grading.angularBiasStrength = control.angularBiasStrength;
    grading.radialBiasToOuterStrength = control.radialBiasToOuterStrength;
    grading.angularBiasOuterLocalizationPower = control.angularBiasOuterLocalizationPower;

    createAnnulusGraded(center,
        innerRadius,
        outerRadius,
        control.startAngle,
        control.endAngle,
        control.radialLayers,
        control.circumferentialNodes,
        control.materialId,
        grading);

    double startRad = (std::abs(control.startAngle) > TWO_PI)
        ? control.startAngle * DEG_TO_RAD
        : control.startAngle;
    double endRad = (std::abs(control.endAngle) > TWO_PI)
        ? control.endAngle * DEG_TO_RAD
        : control.endAngle;
    if (endRad <= startRad) {
        endRad += TWO_PI;
    }

    const std::vector<double> uniformAngularParams =
        buildDensityMappedParameters(control.circumferentialNodes, [](double) { return 1.0; });
    std::vector<double> biasedAngularParams;
    if (grading.useAngularBias && grading.contactHalfAngle > 0.0) {
        const double centerAngle =
            normalizeAngleToSweep(grading.contactCenterAngle, startRad, endRad);
        const double halfAngle = std::abs(grading.contactHalfAngle);
        biasedAngularParams = buildDensityMappedParameters(
            control.circumferentialNodes,
            [=](double s) {
                const double angle = startRad + s * (endRad - startRad);
                const double distance = std::abs(angle - centerAngle);
                const double normalizedDistance = distance / std::max(halfAngle, 1.0e-12);
                return 1.0 + grading.angularBiasStrength *
                    std::exp(-4.0 * normalizedDistance * normalizedDistance);
            });
    }
    else {
        biasedAngularParams = uniformAngularParams;
    }

    std::vector<double> radialParams;
    if (grading.useRadialBias) {
        radialParams = buildDensityMappedParameters(
            control.radialLayers,
            [=](double s) {
                return 1.0 + grading.radialBiasToOuterStrength * s * s;
            });
    }
    else {
        radialParams = buildDensityMappedParameters(control.radialLayers, [](double) { return 1.0; });
    }

    std::vector<double> radii(control.radialLayers, 0.0);
    for (int i = 0; i < control.radialLayers; ++i) {
        radii[i] = innerRadius + radialParams[i] * (outerRadius - innerRadius);
    }

    std::vector<std::vector<double>> anglesByLayer(
        control.radialLayers,
        std::vector<double>(control.circumferentialNodes, startRad));
    for (int layer = 0; layer < control.radialLayers; ++layer) {
        const double physicalRadialCoordinate =
            (radii[layer] - innerRadius) / std::max(outerRadius - innerRadius, 1.0e-15);
        const double angularBiasBlend = grading.localizeAngularBiasToOuterSurface
            ? std::clamp(
                std::pow(std::clamp(physicalRadialCoordinate, 0.0, 1.0),
                    std::max(grading.angularBiasOuterLocalizationPower, 1.0)),
                0.0, 1.0)
            : 1.0;

        for (int node = 0; node < control.circumferentialNodes; ++node) {
            const double blendedParameter =
                uniformAngularParams[node] +
                angularBiasBlend * (biasedAngularParams[node] - uniformAngularParams[node]);
            anglesByLayer[layer][node] = startRad + blendedParameter * (endRad - startRad);
        }
    }

    return buildRingMeshDiagnostics(radii, anglesByLayer, outerRadius);
}

MeshGenerator::TireContactMeshResult MeshGenerator::generateTireContactRingMesh(
    const TireContactMeshControl& control) {
    if (control.radialLayers < 2 || control.circumferentialNodes < 2) {
        throw std::invalid_argument("Invalid tire-contact mesh size");
    }
    if (control.innerRadius >= control.outerRadius) {
        throw std::invalid_argument("Inner radius must be smaller than outer radius");
    }
    if (control.outerRadiusTolerance < 0.0) {
        throw std::invalid_argument("Outer radius tolerance must be non-negative");
    }
    if (control.candidateFacetWindowScale <= 0.0) {
        throw std::invalid_argument("Candidate facet window scale must be positive");
    }

    RingMeshControl ringControl;
    ringControl.startAngle = control.startAngle;
    ringControl.endAngle = control.endAngle;
    ringControl.radialLayers = control.radialLayers;
    ringControl.circumferentialNodes = control.circumferentialNodes;
    ringControl.materialId = control.materialId;
    ringControl.useAngularBias = control.refineCircumferentiallyNearContact;
    ringControl.useRadialBias = control.refineRadiallyToOuterSurface;
    ringControl.localizeAngularBiasToOuterSurface =
        control.localizeCircumferentialRefinementToOuterSurface;
    ringControl.contactCenterAngle = control.expectedContactCenterAngle;
    ringControl.contactHalfAngle = control.expectedContactHalfAngle;
    ringControl.angularBiasStrength = control.circumferentialRefinementStrength;
    ringControl.radialBiasToOuterStrength = control.radialRefinementStrength;
    ringControl.angularBiasOuterLocalizationPower = control.circumferentialLocalizationPower;

    TireContactMeshResult result;
    result.diagnostics = generateTireRingGraded(
        control.center,
        control.innerRadius,
        control.outerRadius,
        ringControl);

    double startRad = (std::abs(control.startAngle) > TWO_PI)
        ? control.startAngle * DEG_TO_RAD
        : control.startAngle;
    double endRad = (std::abs(control.endAngle) > TWO_PI)
        ? control.endAngle * DEG_TO_RAD
        : control.endAngle;
    if (endRad <= startRad) {
        endRad += TWO_PI;
    }

    result.normalizedStartAngle = startRad;
    result.normalizedEndAngle = endRad;
    result.normalizedContactCenterAngle =
        normalizeAngleToSweep(control.expectedContactCenterAngle, startRad, endRad);

    if (control.expectedContactHalfAngle <= 0.0) {
        result.candidateFacetWindowStartAngle = startRad;
        result.candidateFacetWindowEndAngle = endRad;
    }
    else {
        const double maxHalfWindow = 0.5 * (endRad - startRad);
        const double candidateHalfWindow = std::min(
            maxHalfWindow,
            std::abs(control.expectedContactHalfAngle) * control.candidateFacetWindowScale);

        result.candidateFacetWindowStartAngle = normalizeAngleToSweep(
            result.normalizedContactCenterAngle - candidateHalfWindow, startRad, endRad);
        result.candidateFacetWindowEndAngle = normalizeAngleToSweep(
            result.normalizedContactCenterAngle + candidateHalfWindow, startRad, endRad);
    }

    result.candidateContactFacets = collectBoundaryFacetsByPolarWindow(
        control.center,
        control.outerRadius,
        control.outerRadiusTolerance,
        result.candidateFacetWindowStartAngle,
        result.candidateFacetWindowEndAngle);

    if (result.candidateContactFacets.empty()) {
        result.candidateFacetWindowStartAngle = startRad;
        result.candidateFacetWindowEndAngle = endRad;
        result.candidateContactFacets = collectBoundaryFacetsByPolarWindow(
            control.center,
            control.outerRadius,
            control.outerRadiusTolerance,
            startRad,
            endRad);
    }

    return result;
}

MeshGenerator::TireContactAnalysisSetup MeshGenerator::generateTireContactAnalysisSetup(
    const TireContactAnalysisControl& control) {
    if (control.innerRadiusTolerance < 0.0) {
        throw std::invalid_argument("Inner radius tolerance must be non-negative");
    }
    if (control.addInnerBoundaryAnchor &&
        !control.anchorFixX &&
        !control.anchorFixY) {
        throw std::invalid_argument("Anchor must constrain at least one DOF");
    }
    if (control.prescribeInnerBoundaryX && control.addInnerBoundaryAnchor && control.anchorFixX) {
        throw std::invalid_argument(
            "Anchor X constraint conflicts with prescribed inner-boundary X displacement");
    }
    if (control.prescribeInnerBoundaryY && control.addInnerBoundaryAnchor && control.anchorFixY) {
        throw std::invalid_argument(
            "Anchor Y constraint conflicts with prescribed inner-boundary Y displacement");
    }

    TireContactAnalysisSetup setup;
    setup.mesh = generateTireContactRingMesh(control.mesh);
    setup.rigidPlane = control.rigidPlane;

    double anchorSelectorCoordinate = control.anchorSelectMinimumX
        ? std::numeric_limits<double>::infinity()
        : -std::numeric_limits<double>::infinity();
    double anchorAngleDistance = std::numeric_limits<double>::infinity();

    for (const auto& node : assembly_->getNodes()) {
        const Eigen::Vector2d coordinates = node->getCoordinates();
        const double radius = (coordinates - control.mesh.center).norm();
        if (std::abs(radius - control.mesh.innerRadius) > control.innerRadiusTolerance) {
            continue;
        }

        setup.innerBoundaryNodeIds.push_back(node->getId());

        if (control.prescribeInnerBoundaryX && control.prescribeInnerBoundaryY) {
            assembly_->addPrescribedDisplacement(node->getId(),
                control.prescribedInnerBoundaryDx,
                control.prescribedInnerBoundaryDy);
        }
        else if (control.prescribeInnerBoundaryX) {
            assembly_->addPrescribedDisplacementX(node->getId(), control.prescribedInnerBoundaryDx);
        }
        else if (control.prescribeInnerBoundaryY) {
            assembly_->addPrescribedDisplacementY(node->getId(), control.prescribedInnerBoundaryDy);
        }

        if (!control.addInnerBoundaryAnchor) {
            continue;
        }

        const double selectorCoordinate = coordinates.x();
        bool updateAnchor = false;
        if (control.anchorSelectContactCenterAngle) {
            const Eigen::Vector2d relativeCoordinates = coordinates - control.mesh.center;
            const double nodeAngle =
                std::atan2(relativeCoordinates.y(), relativeCoordinates.x());
            const double angleDifference =
                std::atan2(std::sin(nodeAngle - control.mesh.expectedContactCenterAngle),
                    std::cos(nodeAngle - control.mesh.expectedContactCenterAngle));
            const double angleDistance = std::abs(angleDifference);
            updateAnchor = angleDistance < anchorAngleDistance;
            if (updateAnchor) {
                anchorAngleDistance = angleDistance;
            }
        }
        else {
            updateAnchor = control.anchorSelectMinimumX
                ? (selectorCoordinate < anchorSelectorCoordinate)
                : (selectorCoordinate > anchorSelectorCoordinate);
        }
        if (updateAnchor) {
            anchorSelectorCoordinate = selectorCoordinate;
            setup.anchorNodeId = node->getId();
        }
    }

    if (setup.innerBoundaryNodeIds.empty()) {
        throw std::runtime_error("No inner-boundary nodes found for tire-contact setup");
    }

    if (control.addInnerBoundaryAnchor) {
        if (setup.anchorNodeId < 0) {
            throw std::runtime_error("Failed to select an inner-boundary anchor node");
        }
        assembly_->addFixedNode(setup.anchorNodeId, control.anchorFixX, control.anchorFixY);
    }

    return setup;
}

std::vector<ContactFacet> MeshGenerator::collectBoundaryFacetsByCoordinate(int axis,
    double coordinateValue,
    double tolerance) const {
    if (axis < 0 || axis > 1) {
        throw std::invalid_argument("Axis must be 0 (x) or 1 (y)");
    }
    if (tolerance < 0.0) {
        throw std::invalid_argument("Tolerance must be non-negative");
    }

    struct OrderedFacet {
        ContactFacet facet;
        double tangentCoordinate = 0.0;
    };

    std::vector<OrderedFacet> orderedFacets;
    for (const auto& element : assembly_->getElements()) {
        auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
        if (!planeElement) {
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

        if (elementNodes.size() != 4) {
            continue;
        }

        static constexpr int surfaceNodePairs[4][2] = {
            {2, 3}, // eta = -1
            {0, 3}, // xi = 1
            {0, 1}, // eta = 1
            {1, 2}  // xi = -1
        };

        for (int surfaceIndex = 0; surfaceIndex < 4; ++surfaceIndex) {
            const Eigen::Vector2d firstPoint =
                elementNodes[surfaceNodePairs[surfaceIndex][0]]->getCoordinates();
            const Eigen::Vector2d secondPoint =
                elementNodes[surfaceNodePairs[surfaceIndex][1]]->getCoordinates();

            const bool matchesBoundary =
                std::abs(firstPoint(axis) - coordinateValue) <= tolerance &&
                std::abs(secondPoint(axis) - coordinateValue) <= tolerance;

            if (!matchesBoundary) {
                continue;
            }

            OrderedFacet orderedFacet;
            orderedFacet.facet = ContactFacet{element->getId(), surfaceIndex};
            orderedFacet.tangentCoordinate =
                0.5 * (firstPoint(1 - axis) + secondPoint(1 - axis));
            orderedFacets.push_back(orderedFacet);
        }
    }

    std::sort(orderedFacets.begin(), orderedFacets.end(),
        [](const OrderedFacet& lhs, const OrderedFacet& rhs) {
            if (lhs.tangentCoordinate == rhs.tangentCoordinate) {
                if (lhs.facet.elementId == rhs.facet.elementId) {
                    return lhs.facet.surfaceIndex < rhs.facet.surfaceIndex;
                }
                return lhs.facet.elementId < rhs.facet.elementId;
            }
            return lhs.tangentCoordinate < rhs.tangentCoordinate;
        });

    std::vector<ContactFacet> facets;
    facets.reserve(orderedFacets.size());
    for (const auto& orderedFacet : orderedFacets) {
        facets.push_back(orderedFacet.facet);
    }
    return facets;
}

std::vector<ContactFacet> MeshGenerator::collectExteriorFacets(int axis,
    bool selectMinimum,
    double tolerance) const {
    if (assembly_->getNodes().empty()) {
        return {};
    }

    double boundaryCoordinate = assembly_->getNodes().front()->getCoordinates()(axis);
    for (const auto& node : assembly_->getNodes()) {
        boundaryCoordinate = selectMinimum
            ? std::min(boundaryCoordinate, node->getCoordinates()(axis))
            : std::max(boundaryCoordinate, node->getCoordinates()(axis));
    }

    return collectBoundaryFacetsByCoordinate(axis, boundaryCoordinate, tolerance);
}

std::vector<ContactFacet> MeshGenerator::collectBoundaryFacetsByPolarWindow(
    const Eigen::Vector2d& center,
    double targetRadius,
    double radiusTolerance,
    double startAngle,
    double endAngle) const {
    if (radiusTolerance < 0.0) {
        throw std::invalid_argument("Radius tolerance must be non-negative");
    }

    struct OrderedFacet {
        ContactFacet facet;
        double angle = 0.0;
    };

    std::vector<OrderedFacet> orderedFacets;
    for (const auto& element : assembly_->getElements()) {
        auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
        if (!planeElement) {
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

        if (elementNodes.size() != 4) {
            continue;
        }

        static constexpr int surfaceNodePairs[4][2] = {
            {2, 3},
            {0, 3},
            {0, 1},
            {1, 2}
        };

        for (int surfaceIndex = 0; surfaceIndex < 4; ++surfaceIndex) {
            const Eigen::Vector2d p0 =
                elementNodes[surfaceNodePairs[surfaceIndex][0]]->getCoordinates();
            const Eigen::Vector2d p1 =
                elementNodes[surfaceNodePairs[surfaceIndex][1]]->getCoordinates();

            const double r0 = (p0 - center).norm();
            const double r1 = (p1 - center).norm();
            if (std::abs(r0 - targetRadius) > radiusTolerance ||
                std::abs(r1 - targetRadius) > radiusTolerance) {
                continue;
            }

            const Eigen::Vector2d midpoint = 0.5 * (p0 + p1);
            const double angle = std::atan2(midpoint.y() - center.y(), midpoint.x() - center.x());
            if (!isAngleWithinSweep(angle, startAngle, endAngle)) {
                continue;
            }

            orderedFacets.push_back({ContactFacet{element->getId(), surfaceIndex}, angle});
        }
    }

    std::sort(orderedFacets.begin(), orderedFacets.end(),
        [](const OrderedFacet& lhs, const OrderedFacet& rhs) {
            if (lhs.angle == rhs.angle) {
                if (lhs.facet.elementId == rhs.facet.elementId) {
                    return lhs.facet.surfaceIndex < rhs.facet.surfaceIndex;
                }
                return lhs.facet.elementId < rhs.facet.elementId;
            }
            return lhs.angle < rhs.angle;
        });

    std::vector<ContactFacet> facets;
    facets.reserve(orderedFacets.size());
    for (const auto& orderedFacet : orderedFacets) {
        facets.push_back(orderedFacet.facet);
    }
    return facets;
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

MeshGenerator::RingMeshDiagnostics MeshGenerator::buildRingMeshDiagnostics(
    const std::vector<double>& radii,
    const std::vector<std::vector<double>>& anglesByLayer,
    double outerRadius) {
    RingMeshDiagnostics diagnostics;
    if (radii.size() < 2 || anglesByLayer.size() != radii.size() || anglesByLayer.front().size() < 2) {
        return diagnostics;
    }

    diagnostics.minRadialStep = std::numeric_limits<double>::infinity();
    diagnostics.minAngularStep = std::numeric_limits<double>::infinity();
    diagnostics.minOuterArcStep = std::numeric_limits<double>::infinity();
    diagnostics.minAspectRatio = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i + 1 < radii.size(); ++i) {
        const double radialStep = std::abs(radii[i + 1] - radii[i]);
        diagnostics.minRadialStep = std::min(diagnostics.minRadialStep, radialStep);
        diagnostics.maxRadialStep = std::max(diagnostics.maxRadialStep, radialStep);
    }

    const auto& outerAngles = anglesByLayer.back();
    for (size_t i = 0; i + 1 < outerAngles.size(); ++i) {
        const double angularStep = std::abs(outerAngles[i + 1] - outerAngles[i]);
        const double outerArcStep = outerRadius * angularStep;
        diagnostics.minOuterArcStep = std::min(diagnostics.minOuterArcStep, outerArcStep);
        diagnostics.maxOuterArcStep = std::max(diagnostics.maxOuterArcStep, outerArcStep);
    }

    for (size_t radialIndex = 0; radialIndex + 1 < radii.size(); ++radialIndex) {
        const double radialStep = std::abs(radii[radialIndex + 1] - radii[radialIndex]);
        const double meanRadius = 0.5 * (radii[radialIndex + 1] + radii[radialIndex]);
        const auto& lowerAngles = anglesByLayer[radialIndex];
        const auto& upperAngles = anglesByLayer[radialIndex + 1];
        for (size_t angleIndex = 0; angleIndex + 1 < lowerAngles.size(); ++angleIndex) {
            const double lowerAngularStep =
                std::abs(lowerAngles[angleIndex + 1] - lowerAngles[angleIndex]);
            const double upperAngularStep =
                std::abs(upperAngles[angleIndex + 1] - upperAngles[angleIndex]);
            const double representativeAngularStep = 0.5 * (lowerAngularStep + upperAngularStep);
            diagnostics.minAngularStep =
                std::min(diagnostics.minAngularStep, representativeAngularStep);
            diagnostics.maxAngularStep =
                std::max(diagnostics.maxAngularStep, representativeAngularStep);

            const double tangentialStep = meanRadius * representativeAngularStep;
            const double longerStep = std::max(radialStep, tangentialStep);
            const double shorterStep = std::max(1.0e-15, std::min(radialStep, tangentialStep));
            const double aspectRatio = longerStep / shorterStep;
            diagnostics.minAspectRatio = std::min(diagnostics.minAspectRatio, aspectRatio);
            diagnostics.maxAspectRatio = std::max(diagnostics.maxAspectRatio, aspectRatio);
        }
    }

    return diagnostics;
}

bool MeshGenerator::isAngleWithinSweep(double angle, double startAngle, double endAngle) {
    double normalizedStart = startAngle;
    double normalizedEnd = endAngle;
    while (normalizedEnd <= normalizedStart) {
        normalizedEnd += TWO_PI;
    }

    double normalizedAngle = angle;
    while (normalizedAngle < normalizedStart) {
        normalizedAngle += TWO_PI;
    }
    while (normalizedAngle > normalizedEnd) {
        normalizedAngle -= TWO_PI;
    }

    return normalizedAngle >= normalizedStart && normalizedAngle <= normalizedEnd;
}
