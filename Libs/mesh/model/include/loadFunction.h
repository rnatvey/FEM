#pragma once
#include <Eigen/Dense>
#include <functional>
#include <memory>
#include "geometry.h"

class LoadFunction {
public:
    // Функция распределения нагрузки: нагрузка = f(x, y, нормаль)
    using Distribution = std::function<Eigen::Vector2d(double x, double y, const Eigen::Vector2d& normal)>;

    LoadFunction(Distribution func) : distribution_(func) {}

    // Применение к поверхности элемента
    Eigen::VectorXd applyToElementSurface(int elementId, int surfaceIndex,
        const std::shared_ptr<Assembly>& assembly) const;

    // Стандартные распределения
    static LoadFunction constantLoad(double fx, double fy);
    static LoadFunction linearPressure(double p0, double gradient); // Нарастающее давление
    static LoadFunction hertzianPressure(double maxPressure, double contactWidth, double x0);
    static LoadFunction sinusoidalLoad(double amplitude, double wavelength);

private:
    Distribution distribution_;
};