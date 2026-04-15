#pragma once

#include <Eigen/Dense>
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
        double contactCenterAngle = 0.0;
        double contactHalfAngle = 0.0;
        double angularBiasStrength = 4.0;
        double radialBiasToOuterStrength = 2.0;
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

    std::vector<ContactFacet> collectBoundaryFacetsByCoordinate(int axis,
        double coordinateValue,
        double tolerance) const;
    std::vector<ContactFacet> collectExteriorFacets(int axis,
        bool selectMinimum,
        double tolerance) const;

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

    std::shared_ptr<Assembly> assembly_;
    int nextNodeId_ = 1;
    int nextElementId_ = 1;
};
