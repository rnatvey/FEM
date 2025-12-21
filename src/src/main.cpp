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
#include <iomanip> 
#include <cmath> 
#include "assembly.h"
#include "node.h"
#include "material.h"
#include "planeisometric/Planeisoparametric.h"
#include "FEMModel.h"

void testStressCalculation() {
    std::cout << "=== Stress Calculation Test ===" << std::endl;

    // Простой квадратный элемент 1x1
    Eigen::MatrixXd nodeCoords(4, 2);
    nodeCoords << 0, 0,
        1, 0,
        1, 1,
        0, 1;

    // Создаем узелы
    std::vector<std::shared_ptr<Node>> nodes;
    for (int i = 0; i < 4; ++i) {
        nodes.push_back(std::make_shared<Node>(i + 1, nodeCoords(i, 0), nodeCoords(i, 1)));
    }

    // Материал
    auto material = std::make_shared<Material>(2, 2.0e11, 0.3, 0.1);

    // Элемент
    std::vector<int> nodeIds = { 1, 2, 3, 4 };
    auto element = std::make_shared<PlaneIsoparametricElement>(1, nodeIds, 2);

    // Нулевые перемещения
    Eigen::VectorXd zeroDisplacements = Eigen::VectorXd::Zero(8);

    // Должны получить нулевые напряжения
    Eigen::Vector3d stress = element->computeStress(0, 0, zeroDisplacements, nodes, material);
    std::cout << "Stress with zero displacements: " << stress.transpose() << std::endl;

    // Проверим B-матрицу в центре
    Eigen::MatrixXd B = element->strainDisplacementMatrix(0, 0, nodes, material);
    Eigen::Vector3d strain = B * zeroDisplacements;
    std::cout << "Strain with zero displacements: " << strain.transpose() << std::endl;

    // Проверим матрицу D отдельно
    Eigen::Matrix3d D = material->getElasticityMatrix();
    std::cout << "D * zero_strain: " << (D * strain).transpose() << std::endl;
}
void debugStiffnessAssembly() {
    std::cout << "=== Debug Stiffness Assembly ===" << std::endl;

    auto assembly = std::make_shared<Assembly>();

    // Простой материал
    auto material = std::make_shared<Material>(1, 2.0e11, 0.3, 0.1);
    assembly->addMaterial(material);

    // Простой квадратный элемент
    assembly->addNode(std::make_shared<Node>(1, 0.0, 0.0));
    assembly->addNode(std::make_shared<Node>(2, 1.0, 0.0));
    assembly->addNode(std::make_shared<Node>(3, 1.0, 1.0));
    assembly->addNode(std::make_shared<Node>(4, 0.0, 1.0));

    std::vector<int> nodeIds = { 1, 2, 3, 4 };
    auto element = std::make_shared<PlaneIsoparametricElement>(1, nodeIds, 1);
    assembly->addElement(element);

    // Проверим сборку матрицы жесткости
    Eigen::SparseMatrix<double> globalK;
    assembly->assembleGlobalStiffnessMatrix(globalK);

    std::cout << "=== Final Check ===" << std::endl;
    std::cout << "Global matrix size: " << globalK.rows() << "x" << globalK.cols() << std::endl;
    std::cout << "Global matrix non-zeros: " << globalK.nonZeros() << std::endl;
    std::cout << "Global matrix norm: " << globalK.norm() << std::endl;

    if (globalK.norm() < 1e-10) {
        std::cout << "ERROR: Global stiffness matrix is ZERO!" << std::endl;
    }
    else {
        std::cout << "SUCCESS: Global stiffness matrix is non-zero" << std::endl;
    }
}
int main() {
    std::cout << "FEM Contact Solver Test" << std::endl;

    try {
        //debugStiffnessAssembly();
        // Создаем сборку
        auto assembly = std::make_shared<Assembly>();
        auto model = std::make_shared<FEModel>();
        // Добавляем материалы
        auto material = std::make_shared<Material>(1, 6.0, 0.49, 1);
        assembly->addMaterial(material);

        // Добавляем узлы
        vector2 st1(0.0, 0.0);
        vector2 en1(20.0, 0.0);
        vector2 st2(0.0, 4.0);
        vector2 en2(20.0, 4.0);
        vector2 st3(0.0, 8.0);
        vector2 en3(20.0, 8.0);
        vector2 st4(0.0, 12.0);
        vector2 en4(20.0, 12.0);
        assembly->addLineofNodes(1, st1, en1, 5);
        assembly->addLineofNodes(6, st2, en2, 5);
        assembly->addLineofNodes(11, st3, en3, 5);
        assembly->addLineofNodes(16, st4, en4, 5);
        std::cout << "==========================================" << std::endl;
        // Добавляем элементы
        for (int i = 1; i < 5; i++)
        {
            std::vector<int> nodeIds = { i, i + 1, i + 6, i + 5 };
            auto element = std::make_shared<PlaneIsoparametricElement>(i, nodeIds, 1);
            assembly->addElement(element);
        }
        for (int i = 6; i < 10; i++)
        {
            std::vector<int> nodeIds = { i, i + 1, i + 6, i + 5 };
            auto element = std::make_shared<PlaneIsoparametricElement>(i - 1, nodeIds, 1);
            assembly->addElement(element);
        }
        for (int i = 11; i < 15; i++)
        {
            std::vector<int> nodeIds = { i, i + 1, i + 6, i + 5 };
            auto element = std::make_shared<PlaneIsoparametricElement>(i - 2, nodeIds, 1);
            assembly->addElement(element);
        }
        assembly->addFixedNode(1, 1, 1);
        assembly->addFixedNode(6, 1, 0);
        assembly->addFixedNode(11, 1, 0);
        assembly->addFixedNode(16, 1, 0);

        auto ndforc1 = std::make_shared<ConcentratedForce>(5, 0.0, -50e-3);
        assembly->addConcentratedForce(ndforc1);

        //Валидация
            if (assembly->validate()) {
                std::cout << "Assembly validation passed!" << std::endl;


                // Сборка матрицы жесткости
                Eigen::SparseMatrix<double> globalK;
                Eigen::VectorXd globalF;
               
                //assembly->assembleGlobalForceVector(globalF, Force);


                assembly->assembleGlobalStiffnessMatrix(globalK);
                assembly->assembleConcentratedForces(globalF);
                std::cout << "Global stiffness matrix size: " << globalK.rows() << "x" << globalK.cols() << std::endl;
               // std::cout << globalK << std::endl;
                std::cout << "===========================globalK============================================" << std::endl;
                assembly->applyBoundaryConditions(globalK, globalF);
                //std::cout << globalK << std::endl;
                std::cout << "=========================globalF==============================================" << std::endl;
                //std::cout << globalF << std::endl;
                std::cout << "============================Force===========================================" << std::endl;
               // std::cout << Force << std::endl;
                std::cout << "============================DofCount===========================================" << std::endl;
                std::cout << assembly->getTotalDofCount() << std::endl;
                std::cout << "===============================k*F========================================" << std::endl;
                //Eigen::Vector2d sila (5.0,5.0);
                //std::cout << globalK*sila << std::endl;

            }

        model->setAssembly(assembly);
       // Eigen::MatrixXd B = PlaneIsoparametricElement::strainDisplacementMatrix(xi, eta, nodes, material);
        
        if (model->solve()) {
            std::cout << "Solution successful!" << std::endl;
         
            
            auto nodalDisp = model->getNodalDisplacements();
            std::cout << "Nodal Displacements:" << std::endl;
            for (size_t i = 0; i < nodalDisp.size(); ++i) {
                std::cout << "Node " << (i + 1) << ": ("
                    << nodalDisp[i].x() << ", " << nodalDisp[i].y() << ")" << std::endl;
            }

            // Выводим узловые напряжения
            auto nodalStress = model->getNodalStresses();
            std::cout << "Nodal Stresses (sigma_xx, sigma_yy, tau_xy):" << std::endl;
            for (size_t i = 0; i < nodalStress.size(); ++i) {
                std::cout << "Node " << (i + 1) << ": ("
                    << nodalStress[i].x() << ", " << nodalStress[i].y() << ", "
                    << nodalStress[i].z() << ")" << std::endl;
            }

            std::cout << "Solution time: " << model->getSolutionTime() << " seconds" << std::endl;
        }
        else {
            std::cerr << "Solution failed!" << std::endl;
            return 1;
        }
        
        std::cout << "=======================================================================" << std::endl;
        //std::cout << model->getDisplacements() << std::endl;
        
        //std::cout << model->getDisplacements() << std::endl;
        std::cout << "=======================================================================" << std::endl;
       // std::cout << model->getElementStress(5, 0.0, 0.0) << std::endl;
        std::cout << "=======================================================================" << std::endl;
      //  std::cout << model->getElementStress(6, 0.0, 0.0) << std::endl;
       // testStressCalculation();
        
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}