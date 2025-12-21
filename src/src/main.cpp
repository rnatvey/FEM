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
#include "meshgenerator.h"
#include "constants.h"
#include "loadFunction.h"
using namespace Constants;

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
    

    try {
    
        // Создаем сборку
        auto assembly = std::make_shared<Assembly>();
        auto model = std::make_shared<FEModel>();
        auto meshGen = std::make_shared<MeshGenerator>(assembly);
        // Добавляем материалы
        auto material = std::make_shared<Material>(1, 6.0, 0.49, 1);
        assembly->addMaterial(material);
        Eigen::Vector2d center(0.0, 0.0);      // Центр колеса
        double innerRadius = 30;              // Внутренний радиус (обод), м
        double outerRadius = 50;              // Внешний радиус (протектор), м
        double contactAngle = 15.0 * DEG_TO_RAD; // Угол контакта (±15° от вертикали)

        // Создаем кольцевое сечение (90° дуга)
        double startAngle = (PI + contactAngle * 3.0);                    // -45°
        double endAngle = startAngle + PI/2.0;                 // -135° (90° дуга)

        int radialLayers = 5;          // Слоев по толщине
        int circumferentialNodes = 10;  // Узлов по окружности

        meshGen->createAnnulusSimple(center, innerRadius, outerRadius,
            startAngle, endAngle,
            radialLayers, circumferentialNodes,
            1);
        std::cout << "   mesh done: " << assembly->getElementCount() << " elements, "
            << assembly->getNodeCount() << " nodes" << std::endl;

       //////////////////////

        auto nodes = assembly->getNodes();
        int innerNodesCount = 0;

        for (const auto& node : nodes) {
            Eigen::Vector2d coords = node->getCoordinates();
            double radius = coords.norm();

            // Узлы на внутреннем радиусе (с допуском)
            if (std::abs(radius - innerRadius) < 0.5) {
                assembly->addFixedNode(node->getId(), true, true); // Закрепляем X и Y
                innerNodesCount++;
                //std::cout << "   node was fixed: " << node->getId() <<"cords:" << node->getCoordinates() << std::endl;

            }
        }
        std::cout << "   nodes fixed: " << innerNodesCount << std::endl;

        double maxContactPressure = 1;     // 1.5 МПа максимальное давление
        double contactHalfWidth = 50*std::sin(contactAngle);        // 4 см полуширина контакта
        double contactCenterX = 0.0;           // Центр контакта

        auto parabolicLoad = LoadFunction::parabolicPressure(
            maxContactPressure, contactHalfWidth, contactCenterX);

        Eigen::Vector2d normal(0, 1);
        double minY = 1e9;

        for (const auto& node : assembly->getNodes()) {
            minY = std::min(minY, node->getCoordinates().y());
        }

        auto contactNodes = meshGen->findContactNodes(
            contactCenterX, contactHalfWidth, 0.01);

        if (contactNodes.empty()) {
            std::cout << "ERROR: No contact nodes found!" << std::endl;
            return -1;
        }

        std::cout << "Applying load to " << contactNodes.size() << " nodes" << std::endl;

        // Вычисляем суммарное давление для нормировки
        double totalRawPressure = 0;
        std::vector<double> rawPressures;

        for (int nodeId : contactNodes) {
            auto node = assembly->getNode(nodeId);
            Eigen::Vector2d coords = node->getCoordinates();

            Eigen::Vector2d pressure = parabolicLoad.distribution_(
                coords.x(), coords.y(), normal);

            double p = std::abs(pressure.y());
            rawPressures.push_back(p);
            totalRawPressure += p;
        }

        // Нормируем
        

        // Прикладываем силы
        for (size_t i = 0; i < contactNodes.size(); ++i) {
            int nodeId = contactNodes[i];
            double finalForce = rawPressures[i]; // Минус = вниз

            auto force = std::make_shared<ConcentratedForce>(nodeId, 0.0, finalForce);
            assembly->addConcentratedForce(force);

            // Выводим информацию для первых 5 узлов
            if (i < 5) {
                auto node = assembly->getNode(nodeId);
                std::cout << "Node " << nodeId << " (x=" << node->getCoordinates().x()
                    << "): force = " << finalForce << " N" << std::endl;
            }
        }

        std::cout << "... and " << (contactNodes.size() - 5) << " more nodes" << std::endl;
        //model->setAssembly(assembly);
        //model->setSolverTolerance(1.0e-6);



        //std::cout <<"+++++++++" << std::endl;
        //auto chek = assembly->getElement(1)->getNodeIds();
        //std::cout << assembly->getNode(chek[0])->getCoordinates()<<"," << assembly->getNode(chek[1])->getCoordinates()
        //    << "," << assembly->getNode(chek[2])->getCoordinates()
        //    << "," << assembly->getNode(chek[3])->getCoordinates() << std::endl;
        //std::cout << "+++++++++" << std::endl;

        //Валидация
            if (assembly->validate()) {
                std::cout << "Assembly validation passed!" << std::endl;


                // Сборка матрицы жесткости
                Eigen::SparseMatrix<double> globalK;
                Eigen::VectorXd globalF;
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
    //    
    //    std::cout << "=======================================================================" << std::endl;
    //    //std::cout << model->getDisplacements() << std::endl;
    //    
    //    //std::cout << model->getDisplacements() << std::endl;
    //    std::cout << "=======================================================================" << std::endl;
    //   // std::cout << model->getElementStress(5, 0.0, 0.0) << std::endl;
    //    std::cout << "=======================================================================" << std::endl;
    //  //  std::cout << model->getElementStress(6, 0.0, 0.0) << std::endl;
    //   // testStressCalculation();
    //    
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}