#include "geometry.h"
#include <cmath>
#include <stdexcept>

// Вспомогательная функция для нормализации угла
double normalizeAngle(double angle) {
    while (angle > 2 * M_PI) angle -= 2 * M_PI;
    while (angle < 0) angle += 2 * M_PI;
    return angle;
}

Geometry::ParametricCurve Geometry::createLine(const Eigen::Vector2d& start, const Eigen::Vector2d& end) {
    return [start, end](double t) -> Eigen::Vector2d {
        if (t < 0 || t > 1) {
            throw std::invalid_argument("Parameter t must be in [0, 1] for line");
        }
        return start + t * (end - start);
        };
}

Geometry::ParametricCurve Geometry::createArc(const Eigen::Vector2d& center, double radius,
    double startAngle, double endAngle) {
    // Нормализуем углы
    startAngle = normalizeAngle(startAngle);
    endAngle = normalizeAngle(endAngle);

    // Если endAngle меньше startAngle, добавляем 2π
    if (endAngle < startAngle) {
        endAngle += 2 * M_PI;
    }

    double totalAngle = endAngle - startAngle;

    return [center, radius, startAngle, totalAngle](double t) -> Eigen::Vector2d {
        if (t < 0 || t > 1) {
            throw std::invalid_argument("Parameter t must be in [0, 1] for arc");
        }
        double angle = startAngle + t * totalAngle;
        return center + radius * Eigen::Vector2d(std::cos(angle), std::sin(angle));
        };
}

Geometry::ParametricCurve Geometry::createEllipticArc(const Eigen::Vector2d& center,
    double rx, double ry,
    double startAngle, double endAngle) {
    // Нормализуем углы
    startAngle = normalizeAngle(startAngle);
    endAngle = normalizeAngle(endAngle);

    if (endAngle < startAngle) {
        endAngle += 2 * M_PI;
    }

    double totalAngle = endAngle - startAngle;

    return [center, rx, ry, startAngle, totalAngle](double t) -> Eigen::Vector2d {
        if (t < 0 || t > 1) {
            throw std::invalid_argument("Parameter t must be in [0, 1] for elliptic arc");
        }
        double angle = startAngle + t * totalAngle;
        return center + Eigen::Vector2d(rx * std::cos(angle), ry * std::sin(angle));
        };
}

Geometry::ParametricCurve Geometry::createSpline(const std::vector<Eigen::Vector2d>& controlPoints) {
    if (controlPoints.size() < 2) {
        throw std::invalid_argument("Spline requires at least 2 control points");
    }

    return [controlPoints](double t) -> Eigen::Vector2d {
        if (t < 0 || t > 1) {
            throw std::invalid_argument("Parameter t must be in [0, 1] for spline");
        }

        // Простая линейная интерполяция между точками (позже можно улучшить до кубического сплайна)
        double segmentLength = 1.0 / (controlPoints.size() - 1);
        int segment = static_cast<int>(t / segmentLength);

        if (segment >= controlPoints.size() - 1) {
            segment = controlPoints.size() - 2;
        }

        double localT = (t - segment * segmentLength) / segmentLength;
        return controlPoints[segment] * (1 - localT) + controlPoints[segment + 1] * localT;
        };
}

std::vector<Eigen::Vector2d> Geometry::discretizeCurve(const ParametricCurve& curve,
    int numPoints) {
    if (numPoints < 2) {
        throw std::invalid_argument("At least 2 points required for discretization");
    }

    std::vector<Eigen::Vector2d> points;
    points.reserve(numPoints);

    for (int i = 0; i < numPoints; ++i) {
        double t = static_cast<double>(i) / (numPoints - 1);
        points.push_back(curve(t));
    }

    return points;
}