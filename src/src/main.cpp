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
#include "assembly.h"
#include "node.h"
#include "material.h"
#include "planeisometric/Planeisoparametric.h"

int main() {
    std::cout << "FEM Contact Solver Test" << std::endl;

    try {
        // Создаем сборку
        auto assembly = std::make_shared<Assembly>();

        // Добавляем материалы
        auto material = std::make_shared<Material>(1, 2.0e11, 0.3, 0.1);
        assembly->addMaterial(material);

        // Добавляем узлы
        assembly->addNode(std::make_shared<Node>(1, 0.0, 0.0));
        assembly->addNode(std::make_shared<Node>(2, 1.0, 0.0));
        assembly->addNode(std::make_shared<Node>(3, 1.0, 1.0));
        assembly->addNode(std::make_shared<Node>(4, 0.0, 1.0));

        // Добавляем элемент
        std::vector<int> nodeIds = { 1, 2, 3, 4 };
        auto element = std::make_shared<PlaneIsoparametricElement>(1, nodeIds, 1);
        assembly->addElement(element);

        // Валидация
        if (assembly->validate()) {
            std::cout << "Assembly validation passed!" << std::endl;

            // Сборка матрицы жесткости
            Eigen::SparseMatrix<double> globalK;
            assembly->assembleGlobalStiffnessMatrix(globalK);
            std::cout << "Global stiffness matrix size: " << globalK.rows() << "x" << globalK.cols() << std::endl;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}