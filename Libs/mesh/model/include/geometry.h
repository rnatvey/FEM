#pragma once
#include <Eigen/Dense>
#include <functional>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Geometry {
public:
    // Параметрическая кривая: точка = f(t), t ∈ [0, 1]
    using ParametricCurve = std::function<Eigen::Vector2d(double)>;

    // Линейный сегмент
    static ParametricCurve createLine(const Eigen::Vector2d& start, const Eigen::Vector2d& end);

    // Дуга окружности
    static ParametricCurve createArc(const Eigen::Vector2d& center, double radius,
        double startAngle, double endAngle);

    // Эллиптическая дуга
    static ParametricCurve createEllipticArc(const Eigen::Vector2d& center,
        double rx, double ry,
        double startAngle, double endAngle);

    // Сплайн (позже можно добавить NURBS по Бате, стр. 410)
    static ParametricCurve createSpline(const std::vector<Eigen::Vector2d>& controlPoints);

    // Генерация точек вдоль кривой
    static std::vector<Eigen::Vector2d> discretizeCurve(const ParametricCurve& curve,
        int numPoints);
};