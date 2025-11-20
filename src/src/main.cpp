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
        auto material = std::make_shared<Material>(1, 6.0e6, 0.49, 1);
        assembly->addMaterial(material);

        // Добавляем узлы
        assembly->addNode(std::make_shared<Node>(1, 0.0, 0.0));
        assembly->addFixedNode(1, 1, 1);
        assembly->addNode(std::make_shared<Node>(2, 20.0, 0.0));
        assembly->addFixedNode(2, 0, 1);
        assembly->addNode(std::make_shared<Node>(3, 20.0, 12.0));
        assembly->addFixedNode(3, 0, 1);
        assembly->addNode(std::make_shared<Node>(4, 0.0, 12.0));
        assembly->addFixedNode(4, 1, 1);

        auto ndforc1 = std::make_shared<ConcentratedForce>(2, 11843663.639952626079, 0.0);
        auto ndforc2 = std::make_shared<ConcentratedForce>(3, 11843663.639952627942, 0.0);
        assembly->addConcentratedForce(ndforc1);
        assembly->addConcentratedForce(ndforc2);

        std::vector<int> nodeIds = { 1, 2, 3, 4 };
        auto element = std::make_shared<PlaneIsoparametricElement>(1, nodeIds, 1);
        std::cout << element->getNodeCount() << std::endl;
        assembly->addElement(element);
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
                std::cout << globalK << std::endl;
                std::cout << "===========================globalK============================================" << std::endl;
                assembly->applyBoundaryConditions(globalK, globalF);
                std::cout << globalK << std::endl;
                std::cout << "=========================globalF==============================================" << std::endl;
                std::cout << globalF << std::endl;
                std::cout << "============================Force===========================================" << std::endl;
               // std::cout << Force << std::endl;
                std::cout << "============================DofCount===========================================" << std::endl;
                std::cout << assembly->getTotalDofCount() << std::endl;
                std::cout << "===============================k*F========================================" << std::endl;
                //Eigen::Vector3d Force{ 0.0, 5.0,5.0 };
                //cout.precision(15);
                //auto sila = globalK * Force;
               // std::cout << std::setprecision(20) << sila<< std::endl;
                //1.18437e+07
            }
        
        model->setAssembly(assembly);

        
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
        
        std::cout << model->getDisplacements() << std::endl;
        std::cout << "=======================================================================" << std::endl;
        std::cout << model->getElementStress(1, 0.0, 0.0) << std::endl;
        std::cout << "=======================================================================" << std::endl;
       // testStressCalculation();
        
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}