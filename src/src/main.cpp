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
int main() {
    

    try {
    
        // Создаем сборку
        auto assembly = std::make_shared<Assembly>();
        auto model = std::make_shared<FEModel>();
        auto meshGen = std::make_shared<MeshGenerator>(assembly);
        // Добавляем материалы
        auto material = std::make_shared<Material>(1, 11.84, 0.48, 1);
        assembly->addMaterial(material);
        //
        Eigen::Vector2d center(0.0, 0.0);      
        double innerRadius = 250;              
        double outerRadius = 300;              
        double contactAngle = 5 * DEG_TO_RAD; 

        // Создаем кольцевое сечение 
        double startAngle = (PI + 60 * DEG_TO_RAD );                    
        double endAngle = startAngle + 60 * DEG_TO_RAD;                 

        int radialLayers = 250;          // Слоев по толщине
        int circumferentialNodes = 251;  // Узлов по окружности

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
                assembly->addFixedNode(node->getId(), true, true); 
                innerNodesCount++;
               

            }
        }
        std::cout << "   nodes fixed: " << innerNodesCount << std::endl;

        double maxContactPressure = 1.0; // 858716.234090089e-6;     //  максимальное давление
        double contactHalfWidth = outerRadius * std::sin(contactAngle);        //  полуширина контакта
        double contactCenterX = 0.0;           
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
            }

       model->setAssembly(assembly);
       
        
        if (model->solve()) {
            std::cout << "Solution successful!" << std::endl;
         
            
            //auto nodalDisp = model->getNodalDisplacements();
            //std::cout << "Nodal Displacements:" << std::endl;
            //for (size_t i = 0; i < nodalDisp.size(); ++i) {
            //    std::cout << "Node " << (i + 1) << ": ("
            //        << nodalDisp[i].x() << ", " << nodalDisp[i].y() << ")" << std::endl;
            //}

            // //Выводим узловые напряжения
            //auto nodalStress = model->getNodalStresses();
            //std::cout << "Nodal Stresses (sigma_xx, sigma_yy, tau_xy):" << std::endl;
            //for (size_t i = 0; i < nodalStress.size(); ++i) {
            //    std::cout << "Node " << (i + 1) << ": ("
            //        << nodalStress[i].x() << ", " << nodalStress[i].y() << ", "
            //        << nodalStress[i].z() << ")" << std::endl;
            //}

            std::cout << "Solution time: " << model->getSolutionTime() << " seconds" << std::endl;
        }
        else {
            std::cerr << "Solution failed!" << std::endl;
            return 1;
        }


        // далее попытка  портировать 
        std::ofstream file("tireREF5.csv");
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
        std::cout << "data saved as tireREF5.csv" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}