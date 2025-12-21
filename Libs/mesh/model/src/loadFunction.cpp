#include "loadFunction.h"
#include "assembly.h"
#include "planeisometric/Planeisoparametric.h"
#include <cmath>
#include <iostream>

Eigen::VectorXd LoadFunction::applyToElementSurface(int elementId, int surfaceIndex,
    const std::shared_ptr<Assembly>& assembly) const {
    auto element = assembly->getElement(elementId);
    auto material = assembly->getMaterial(element->getMaterialId());
    auto nodes = assembly->getNodes();

    if (!element || !material) {
        throw std::invalid_argument("Invalid element or material");
    }

    // Преобразуем в PlaneIsoparametricElement для доступа к методам поверхности
    auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
    if (!planeElement) {
        throw std::runtime_error("Element is not a PlaneIsoparametricElement");
    }

    // Точки Гаусса для интегрирования по поверхности
    const double gp = 0.577350269189626;
    const std::vector<double> gaussPoints = { -gp, gp };
    const std::vector<double> weights = { 1.0, 1.0 };

    Eigen::VectorXd surfaceForces = Eigen::VectorXd::Zero(element->getDofCount());

    for (int i = 0; i < 2; ++i) {
        double param = gaussPoints[i];
        double weight = weights[i];

        // Определяем локальные координаты на поверхности
        double xi, eta;
        switch (surfaceIndex) {
        case 0: xi = param; eta = -1.0; break; // нижняя
        case 1: xi = 1.0; eta = param; break;  // правая
        case 2: xi = param; eta = 1.0; break;  // верхняя
        case 3: xi = -1.0; eta = param; break; // левая
        default:
            throw std::invalid_argument("Invalid surface index");
        }

        // Вычисляем глобальные координаты точки
        Eigen::MatrixXd N = planeElement->shapeFunctions(xi, eta);
        Eigen::Vector2d point(0, 0);
        auto nodeCoords = planeElement->getNodalCoordinates(nodes);

        for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx) {
            double shapeFunc = N(0, nodeIdx * 2); // N_x или N_y одинаковы
            point += shapeFunc * nodeCoords.row(nodeIdx).transpose();
        }

        // Нормаль к поверхности
        Eigen::Vector2d normal = planeElement->getSurfaceNormal(surfaceIndex, nodes);

        // Вычисляем нагрузку в точке
        Eigen::Vector2d load = distribution_(point.x(), point.y(), normal);

        // Длина сегмента (якобиан преобразования на поверхности)
        Eigen::Matrix2d J = planeElement->jacobian(xi, eta, nodes);
        double surfaceJacobian;
        switch (surfaceIndex) {
        case 0: case 2:
            surfaceJacobian = std::sqrt(J(0, 0) * J(0, 0) + J(0, 1) * J(0, 1));
            break;
        case 1: case 3:
            surfaceJacobian = std::sqrt(J(1, 0) * J(1, 0) + J(1, 1) * J(1, 1));
            break;
        default:
            surfaceJacobian = 1.0;
        }

        // Добавляем вклад в узловые силы
        for (int nodeIdx = 0; nodeIdx < 4; ++nodeIdx) {
            double shapeFunc = N(0, nodeIdx * 2);
            int dofX = nodeIdx * 2;
            int dofY = nodeIdx * 2 + 1;

            surfaceForces(dofX) += shapeFunc * load.x() * weight * surfaceJacobian * material->getThickness();
            surfaceForces(dofY) += shapeFunc * load.y() * weight * surfaceJacobian * material->getThickness();
        }
    }

    return surfaceForces;
}

LoadFunction LoadFunction::constantLoad(double fx, double fy) {
    return LoadFunction([fx, fy](double x, double y, const Eigen::Vector2d& normal) {
        return Eigen::Vector2d(fx, fy);
        });
}

LoadFunction LoadFunction::linearPressure(double p0, double gradient) {
    return LoadFunction([p0, gradient](double x, double y, const Eigen::Vector2d& normal) {
        double pressure = p0 + gradient * x; // Линейная зависимость от x
        return pressure * normal; // Давление действует по нормали
        });
}

LoadFunction LoadFunction::hertzianPressure(double maxPressure, double contactWidth, double x0) {
    return LoadFunction([maxPressure, contactWidth, x0](double x, double y, const Eigen::Vector2d& normal) {
        double dx = x - x0;  // Расстояние от центра контакта

        // Если точка вне зоны контакта
        if (std::abs(dx) > contactWidth) {
            return Eigen::Vector2d(0, 0);
        }

        // Эллиптическое распределение Герца
        double relativeX = dx / contactHalfWidth;  // x/a
        double pressure = maxPressure * std::sqrt(1.0 - relativeX * relativeX);

        // Давление действует по нормали к поверхности
        return pressure * normal;
    });

LoadFunction LoadFunction::sinusoidalLoad(double amplitude, double wavelength) {
    return LoadFunction([amplitude, wavelength](double x, double y, const Eigen::Vector2d& normal) {
        double pressure = amplitude * std::sin(2 * M_PI * x / wavelength);
        return pressure * normal;
        });
}