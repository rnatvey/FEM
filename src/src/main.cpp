//#include <iostream>
//#include "slau.h"
//#include "vector.h"
//#include "matrix.h"
//#include <Eigen/Dense>
//
//int main()
//{
//	std::cout << "111"<<std::endl;
//	std::cout << testing();
//}
#include <iostream>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>

// Пример 1: Интегрирование простой функции
double simple_function(double x) {
    return x * x;  // ?x? dx = x?/3
}

// Пример 2: Гауссова функция
double gaussian(double x) {
    return std::exp(-x * x);
}

// Пример 3: Тригонометрическая функция
double trig_function(double x) {
    return std::sin(x);
}

int main() {
    using namespace boost::math::quadrature;

    // Метод tanh-sinh (очень точный для гладких функций)
    tanh_sinh<double> integrator;

    std::cout << "=== Boost.Math Numerical Integration Examples ===\n" << std::endl;

    // Пример 1: ?x? dx от 0 до 1 = 1/3 ? 0.333...
    try {
        double result1 = integrator.integrate(simple_function, 0.0, 1.0);
        double exact1 = 1.0 / 3.0;
        std::cout << "?x? dx from 0 to 1:" << std::endl;
        std::cout << "  Result: " << result1 << std::endl;
        std::cout << "  Exact:  " << exact1 << std::endl;
        std::cout << "  Error:  " << std::abs(result1 - exact1) << "\n" << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    // Пример 2: ?exp(-x?) dx от -? до ? = ?? ? 1.77245
    try {
        double result2 = integrator.integrate(gaussian, -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());
        double exact2 = boost::math::constants::root_pi<double>();
        std::cout << "?exp(-x?) dx from -? to ?:" << std::endl;
        std::cout << "  Result: " << result2 << std::endl;
        std::cout << "  Exact:  " << exact2 << std::endl;
        std::cout << "  Error:  " << std::abs(result2 - exact2) << "\n" << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    // Пример 3: ?sin(x) dx от 0 до ? = 2
    try {
        double result3 = integrator.integrate(trig_function, 0.0, boost::math::constants::pi<double>());
        double exact3 = 2.0;
        std::cout << "?sin(x) dx from 0 to ?:" << std::endl;
        std::cout << "  Result: " << result3 << std::endl;
        std::cout << "  Exact:  " << exact3 << std::endl;
        std::cout << "  Error:  " << std::abs(result3 - exact3) << "\n" << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}