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
#include <fstream>

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
        auto material = std::make_shared<Material>(1, 12.728, 0.48, 1);
        assembly->addMaterial(material);
        //
        Eigen::Vector2d center(0.0, 0.0);      // Центр колеса
        double innerRadius = 272.75;              // Внутренний радиус (обод), м
        double outerRadius = 312.75;              // Внешний радиус (протектор), м
        double contactAngle = 0.08249400479616308; // Угол контакта (±15° от вертикали)

        // Создаем кольцевое сечение (90° дуга)
        double startAngle = (PI + 60 * DEG_TO_RAD );                    // -45°
        double endAngle = startAngle + 60 * DEG_TO_RAD;                 // -135° (90° дуга)

        int radialLayers = 40;          // Слоев по толщине
        int circumferentialNodes = 101;  // Узлов по окружности

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
            if (std::abs(radius - innerRadius) < 0.005) {
                assembly->addFixedNode(node->getId(), true, true); // Закрепляем X и Y
                innerNodesCount++;
                //std::cout << "   node was fixed: " << node->getId() <<"cords:" << node->getCoordinates() << std::endl;

            }
        }
        std::cout << "   nodes fixed: " << innerNodesCount << std::endl;

        double maxContactPressure = 858716.234090089e-6;     // 1 МПа максимальное давление
        double contactHalfWidth = outerRadius * std::sin(contactAngle);        //  полуширина контакта
        double contactCenterX = 0.0;           // Центр контакта
        auto hertzLoad = std::make_shared<LoadFunction>(LoadFunction::parabolicPressure(maxContactPressure, contactHalfWidth, contactCenterX));
        auto elements = assembly->getElements();
        int loadedElements = 0;

        for (const auto& element : elements) {
            auto planeElement = std::dynamic_pointer_cast<PlaneIsoparametricElement>(element);
            if (!planeElement) continue;

            // Получаем координаты узлов элемента

            std::vector<std::shared_ptr<Node>> elementNodes;
            for (int nodeId : element->getNodeIds()) {
                auto node = assembly->getNode(nodeId);
                if (node) {
                    elementNodes.push_back(node);
                }
                else {
                    std::cerr << "Assembly validation failed: Element " << element->getId()
                        << " references non-existent node " << nodeId << std::endl;
                    return false;
                }
            }
            auto nodeCoords = planeElement->getNodalCoordinates(elementNodes);
            // Проверяем радиусы узлов элемента
            bool isOuterElement = false;
            for (int i = 0; i < 4; ++i) {
                double radius = nodeCoords.row(i).norm();
                if (std::abs(radius - outerRadius) < 0.005) { // Допуск 0.5 м
                    isOuterElement = true;
                    break;
                }
            }
            elementNodes.clear();
            if (isOuterElement) {
                // Для кольцевой сетки нижняя грань обычно surfaceIndex = 0
                // (между узлами 0-1 в локальных координатах)
                assembly->addSurfaceLoad(element->getId(), 0, hertzLoad);
                loadedElements++;
            }
        }

        std::cout << "Pressure applied to " << loadedElements << " outer elements" << std::endl;
        //Валидация
            if (assembly->validate()) {
                std::cout << "Assembly validation passed!" << std::endl;


                 //Сборка матрицы жесткости
               // Eigen::SparseMatrix<double> globalK;
               // Eigen::VectorXd globalF;
               // assembly->assembleGlobalStiffnessMatrix(globalK);
               // assembly->assembleConcentratedForces(globalF);
               // assembly->assembleSurfaceLoads(globalF);
               // std::cout << "Global stiffness matrix size: " << globalK.rows() << "x" << globalK.cols() << std::endl;
               //// std::cout << globalK << std::endl;
               // std::cout << "===========================globalK============================================" << std::endl;
               // assembly->applyBoundaryConditions(globalK, globalF);
               // //std::cout << globalK << std::endl;
               // std::cout << "=========================globalF==============================================" << std::endl;
               // //std::cout << globalF << std::endl;
               // std::cout << "============================Force===========================================" << std::endl;
               //// std::cout << Force << std::endl;
               // std::cout << "============================DofCount===========================================" << std::endl;
               // std::cout << assembly->getTotalDofCount() << std::endl;
               // std::cout << "===============================k*F========================================" << std::endl;
               // //Eigen::Vector2d sila (5.0,5.0);
               // //std::cout << globalK*sila << std::endl;

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

             //Выводим узловые напряжения
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

        // далее попытка  портировать 
      /*  std::ofstream file("tireTest.csv");
        if (!file.is_open()) {
            std::cerr << "cont open this file!" << std::endl;
            return 1;
        }
        file << " nodeID, x, y, ux, uy, sx, sy, sxy\n";
        std::vector<std::tuple<size_t, const double, const double, double, double, double,double,double>> data;
        auto node_out = assembly->getNodes();
        auto nodalDisp_out = model->getNodalDisplacements();
        auto nodalStress_out = model->getNodalStresses();
        for (size_t i = 0; i < node_out.size(); ++i) {
            data.push_back(std::make_tuple(
                i + 1,
                node_out[i]->getCoordinates()(0),
                node_out[i]->getCoordinates()(1),
                nodalDisp_out[i].x(),
                nodalDisp_out[i].y(),
                nodalStress_out[i].x(),
                nodalStress_out[i].y(),
                nodalStress_out[i].z()
            ));
        }
        for (const auto& row : data) {
            file << std::get<0>(row) << ","
                << std::get<1>(row) << ","
                << std::get<2>(row) << ","
                << std::get<3>(row) << ","
                << std::get<4>(row) << ","
                << std::get<5>(row) << ","
                << std::get<6>(row) << ","
                << std::get<7>(row) << "\n";
        }
        file.close();
        std::cout << "data saved as tireTest.csv" << std::endl;*/
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}