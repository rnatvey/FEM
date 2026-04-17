#pragma once

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <vector>

#include "ContactTypes.h"
#include "assembly.h"
#include "geometry.h"

class MeshGenerator {
public:
    struct Block {
        std::vector<Geometry::ParametricCurve> edges;
        int nodesX = 0;
        int nodesY = 0;
        int materialId = 0;
    };

    struct AnnulusGrading {
        bool useAngularBias = false;
        bool useRadialBias = false;
        bool localizeAngularBiasToOuterSurface = false;
        double contactCenterAngle = 0.0;
        double contactHalfAngle = 0.0;
        double angularBiasStrength = 4.0;
        double radialBiasToOuterStrength = 2.0;
        double angularBiasOuterLocalizationPower = 2.0;
    };

    struct RingMeshControl {
        double startAngle = 0.0;
        double endAngle = 2.0 * EIGEN_PI;
        int radialLayers = 0;
        int circumferentialNodes = 0;
        int materialId = 0;
        bool useAngularBias = true;
        bool useRadialBias = true;
        bool localizeAngularBiasToOuterSurface = false;
        double contactCenterAngle = -0.5 * EIGEN_PI;
        double contactHalfAngle = 0.25 * EIGEN_PI;
        double angularBiasStrength = 4.0;
        double radialBiasToOuterStrength = 2.0;
        double angularBiasOuterLocalizationPower = 2.0;
    };

    struct RingMeshDiagnostics {
        double minRadialStep = 0.0;
        double maxRadialStep = 0.0;
        double minAngularStep = 0.0;
        double maxAngularStep = 0.0;
        double minOuterArcStep = 0.0;
        double maxOuterArcStep = 0.0;
        double minAspectRatio = 0.0;
        double maxAspectRatio = 0.0;
    };

    struct TireContactMeshControl {
        Eigen::Vector2d center = Eigen::Vector2d::Zero();
        double innerRadius = 0.0;
        double outerRadius = 0.0;
        double startAngle = 0.0;
        double endAngle = 2.0 * EIGEN_PI;
        int radialLayers = 0;
        int circumferentialNodes = 0;
        int materialId = 0;

        bool refineCircumferentiallyNearContact = true;
        bool refineRadiallyToOuterSurface = true;
        bool localizeCircumferentialRefinementToOuterSurface = true;
        double expectedContactCenterAngle = -0.5 * EIGEN_PI;
        double expectedContactHalfAngle = 0.25 * EIGEN_PI;
        double circumferentialRefinementStrength = 4.0;
        double radialRefinementStrength = 2.0;
        double circumferentialLocalizationPower = 2.0;

        // Expand the candidate contact window beyond the expected patch
        // so the active set can grow without remeshing.
        double candidateFacetWindowScale = 3.0;
        double outerRadiusTolerance = 1.0e-6;
    };

    struct TireContactMeshResult {
        RingMeshDiagnostics diagnostics;
        std::vector<ContactFacet> candidateContactFacets;
        double normalizedStartAngle = 0.0;
        double normalizedEndAngle = 0.0;
        double normalizedContactCenterAngle = 0.0;
        double candidateFacetWindowStartAngle = 0.0;
        double candidateFacetWindowEndAngle = 0.0;
    };

    struct TireContactAnalysisControl {
        TireContactMeshControl mesh;
        RigidPlane2D rigidPlane = RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0};
        double innerRadiusTolerance = 1.0e-6;
        bool prescribeInnerBoundaryX = false;
        bool prescribeInnerBoundaryY = true;
        double prescribedInnerBoundaryDx = 0.0;
        double prescribedInnerBoundaryDy = 0.0;
        bool addInnerBoundaryAnchor = true;
        bool anchorFixX = true;
        bool anchorFixY = false;
        bool anchorSelectMinimumX = true;
    };

    struct TireContactAnalysisSetup {
        TireContactMeshResult mesh;
        RigidPlane2D rigidPlane;
        std::vector<int> innerBoundaryNodeIds;
        int anchorNodeId = -1;
    };

    explicit MeshGenerator(std::shared_ptr<Assembly> assembly);

    void addBlock(const Block& block);

    void createRectangle(const Eigen::Vector2d& corner1, const Eigen::Vector2d& corner2,
        int nodesX, int nodesY, int materialId);
    void createAnnulus(const Eigen::Vector2d& center, double innerRadius, double outerRadius,
        double startAngle, double endAngle,
        int radialLayers, int circumferentialNodes, int materialId);
    void createTriangle(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3,
        int divisions, int materialId);

    void setNodeIdStart(int startId) { nextNodeId_ = startId; }
    void setElementIdStart(int startId) { nextElementId_ = startId; }

    void createAnnulusSimple(const Eigen::Vector2d& center,
        double innerRadius, double outerRadius,
        double startAngle, double endAngle,
        int radialLayers, int circumferentialNodes,
        int materialId);

    void createAnnulusGraded(const Eigen::Vector2d& center,
        double innerRadius, double outerRadius,
        double startAngle, double endAngle,
        int radialLayers, int circumferentialNodes,
        int materialId,
        const AnnulusGrading& grading);
    RingMeshDiagnostics generateTireRingGraded(const Eigen::Vector2d& center,
        double innerRadius,
        double outerRadius,
        const RingMeshControl& control);
    TireContactMeshResult generateTireContactRingMesh(const TireContactMeshControl& control);
    TireContactAnalysisSetup generateTireContactAnalysisSetup(
        const TireContactAnalysisControl& control);

    std::vector<ContactFacet> collectBoundaryFacetsByCoordinate(int axis,
        double coordinateValue,
        double tolerance) const;
    std::vector<ContactFacet> collectExteriorFacets(int axis,
        bool selectMinimum,
        double tolerance) const;
    std::vector<ContactFacet> collectBoundaryFacetsByPolarWindow(const Eigen::Vector2d& center,
        double targetRadius,
        double radiusTolerance,
        double startAngle,
        double endAngle) const;

    std::vector<int> findContactNodes(double contactCenterX,
        double contactHalfWidth,
        double maxYtolerance) const;
    void applyParabolicContactToNodes(double maxPressure,
        double contactHalfWidth,
        double contactCenterX,
        double totalForce);

private:
    void generateBlockMesh(const Block& block);
    Eigen::Vector2d transfiniteInterpolation(const Block& block, double xi, double eta) const;

    static std::vector<double> buildDensityMappedParameters(int pointCount,
        const std::function<double(double)>& densityFunction);
    static double normalizeAngleToSweep(double angle, double startAngle, double endAngle);
    static RingMeshDiagnostics buildRingMeshDiagnostics(const std::vector<double>& radii,
        const std::vector<std::vector<double>>& anglesByLayer,
        double outerRadius);
    static bool isAngleWithinSweep(double angle, double startAngle, double endAngle);

    std::shared_ptr<Assembly> assembly_;
    int nextNodeId_ = 1;
    int nextElementId_ = 1;
};
