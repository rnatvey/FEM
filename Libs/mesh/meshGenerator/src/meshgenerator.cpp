#include "meshgenerator.h"
#include "planeisometric/Planeisoparametric.h"
#include <iostream>
#include <cmath>
#include "constants.h"
#include "loadFunction.h"
using namespace Constants;

MeshGenerator::MeshGenerator(std::shared_ptr<Assembly> assembly)
    : assembly_(assembly), nextNodeId_(1), nextElementId_(1) {
    if (!assembly_) {
        throw std::invalid_argument("Assembly cannot be null");
    }
}

void MeshGenerator::addBlock(const Block& block) {
    if (block.edges.size() != 4) {
        throw std::invalid_argument("Block must have exactly 4 edges");
    }

    if (block.nodesX < 2 || block.nodesY < 2) {
        throw std::invalid_argument("Block must have at least 2 nodes in each direction");
    }

    generateBlockMesh(block);
}

void MeshGenerator::createRectangle(const Eigen::Vector2d& corner1, const Eigen::Vector2d& corner2,
    int nodesX, int nodesY, int materialId) {
    std::cout << "Creating rectangle mesh: "
        << nodesX << "x" << nodesY << " elements" << std::endl;

    // Создаем 4 края прямоугольника
    Eigen::Vector2d p1 = corner1;
    Eigen::Vector2d p2(corner2.x(), corner1.y());
    Eigen::Vector2d p3 = corner2;
    Eigen::Vector2d p4(corner1.x(), corner2.y());

    Block block;
    block.edges = {
        Geometry::createLine(p1, p2),  // нижний край
        Geometry::createLine(p2, p3),  // правый край
        Geometry::createLine(p3, p4),  // верхний край
        Geometry::createLine(p4, p1)   // левый край
    };

    block.nodesX = nodesX;
    block.nodesY = nodesY;
    block.materialId = materialId;

    addBlock(block);
}

void MeshGenerator::createAnnulus(const Eigen::Vector2d& center, double innerRadius, double outerRadius,
    double startAngle, double endAngle,
    int radialLayers, int circumferentialNodes, int materialId) {
    std::cout << "Creating annulus mesh: "
        << radialLayers << "x" << circumferentialNodes << " elements" << std::endl;

    // Создаем 4 края кольцевого сектора
    Block block;

    // Внутренняя дуга
    block.edges.push_back(Geometry::createArc(center, innerRadius, startAngle, endAngle));

    // Внешняя дуга
    block.edges.push_back(Geometry::createArc(center, outerRadius, startAngle, endAngle));

    // Радиальные стороны
    block.edges.push_back(Geometry::createLine(
        center + innerRadius * Eigen::Vector2d(std::cos(startAngle), std::sin(startAngle)),
        center + outerRadius * Eigen::Vector2d(std::cos(startAngle), std::sin(startAngle))
    ));

    block.edges.push_back(Geometry::createLine(
        center + innerRadius * Eigen::Vector2d(std::cos(endAngle), std::sin(endAngle)),
        center + outerRadius * Eigen::Vector2d(std::cos(endAngle), std::sin(endAngle))
    ));

    block.nodesX = circumferentialNodes;
    block.nodesY = radialLayers;
    block.materialId = materialId;

    addBlock(block);
}

void MeshGenerator::createTriangle(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
    const Eigen::Vector2d& p3, int divisions, int materialId) {
    std::cout << "Creating triangle mesh with " << divisions << " divisions" << std::endl;

    // Для треугольника используем два блока или специальную логику
    // Упрощенная реализация: разбиваем на 4-угольники через добавление центральной точки

    Eigen::Vector2d center = (p1 + p2 + p3) / 3.0;

    // Создаем 3 блока (по одному на каждую сторону треугольника)
    // Это упрощенный подход, можно улучшить

    // Блок 1: p1-p2-center
    Block block1;
    block1.edges = {
        Geometry::createLine(p1, p2),
        Geometry::createLine(p2, center),
        Geometry::createLine(center, p1),
        Geometry::createLine(p1, p1)  // Фиктивный четвертый край
    };
    block1.nodesX = divisions;
    block1.nodesY = divisions;
    block1.materialId = materialId;

    // Аналогично для остальных блоков...
    // (это упрощенная реализация, нужна более сложная логика для треугольников)

    std::cout << "WARNING: Triangle meshing is simplified. Consider using quadrilateral elements only." << std::endl;
}

Eigen::Vector2d MeshGenerator::transfiniteInterpolation(const Block& block, double xi, double eta) const {
    // Трансфинитная интерполяция Кунса (Coons patch)
    // r(ξ,η) = (1-η)C1(ξ) + ηC3(ξ) + (1-ξ)C4(η) + ξC2(η) 
    //          - [(1-ξ)(1-η)P1 + ξ(1-η)P2 + ξηP3 + (1-ξ)ηP4]

    // Угловые точки
    Eigen::Vector2d P1 = block.edges[0](0.0);  // (0,0)
    Eigen::Vector2d P2 = block.edges[0](1.0);  // (1,0)
    Eigen::Vector2d P3 = block.edges[2](1.0);  // (1,1)
    Eigen::Vector2d P4 = block.edges[2](0.0);  // (0,1)

    // Кривые в параметрическом виде
    Eigen::Vector2d C1 = block.edges[0](xi);  // нижний край (η=0)
    Eigen::Vector2d C3 = block.edges[2](xi);  // верхний край (η=1)
    Eigen::Vector2d C2 = block.edges[1](eta); // правый край (ξ=1)
    Eigen::Vector2d C4 = block.edges[3](eta); // левый край (ξ=0)

    // Интерполяция Кунса
    Eigen::Vector2d result =
        (1 - eta) * C1 + eta * C3 +
        (1 - xi) * C4 + xi * C2 -
        ((1 - xi) * (1 - eta) * P1 +
            xi * (1 - eta) * P2 +
            xi * eta * P3 +
            (1 - xi) * eta * P4);

    return result;
}

void MeshGenerator::generateBlockMesh(const Block& block) {
    std::cout << "Generating block mesh with "
        << block.nodesX << "x" << block.nodesY << " nodes" << std::endl;

    // Проверяем, что у нас 4 грани
    if (block.edges.size() != 4) {
        throw std::invalid_argument("Block must have exactly 4 edges");
    }

    // 1. Генерация узлов ТОЛЬКО НА ГРАНИЦАХ с последующей прямолинейной интерполяцией
    std::vector<std::vector<int>> nodeGrid(block.nodesY, std::vector<int>(block.nodesX, 0));

    // 1.1. Генерируем узлы на нижней грани (edge 0)
    std::vector<Eigen::Vector2d> bottomEdgePoints =
        Geometry::discretizeCurve(block.edges[0], block.nodesX);
    for (int i = 0; i < block.nodesX; ++i) {
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId,
            bottomEdgePoints[i].x(), bottomEdgePoints[i].y()));
        nodeGrid[0][i] = nodeId;
    }

    // 1.2. Генерируем узлы на верхней грани (edge 2)
    std::vector<Eigen::Vector2d> topEdgePoints =
        Geometry::discretizeCurve(block.edges[2], block.nodesX);
    for (int i = 0; i < block.nodesX; ++i) {
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId,
            topEdgePoints[i].x(), topEdgePoints[i].y()));
        nodeGrid[block.nodesY - 1][i] = nodeId;
    }

    // 1.3. Генерируем узлы на левой грани (edge 3) - кроме уже созданных углов
    std::vector<Eigen::Vector2d> leftEdgePoints =
        Geometry::discretizeCurve(block.edges[3], block.nodesY);
    for (int j = 1; j < block.nodesY - 1; ++j) { // Пропускаем углы
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId,
            leftEdgePoints[j].x(), leftEdgePoints[j].y()));
        nodeGrid[j][0] = nodeId;
    }

    // 1.4. Генерируем узлы на правой грани (edge 1) - кроме уже созданных углов
    std::vector<Eigen::Vector2d> rightEdgePoints =
        Geometry::discretizeCurve(block.edges[1], block.nodesY);
    for (int j = 1; j < block.nodesY - 1; ++j) { // Пропускаем углы
        int nodeId = nextNodeId_++;
        assembly_->addNode(std::make_shared<Node>(nodeId,
            rightEdgePoints[j].x(), rightEdgePoints[j].y()));
        nodeGrid[j][block.nodesX - 1] = nodeId;
    }

    // 1.5. Генерируем внутренние узлы через билинейную интерполяцию
    for (int j = 1; j < block.nodesY - 1; ++j) {
        double eta = static_cast<double>(j) / (block.nodesY - 1);

        // Точки на левой и правой гранях для этого eta
        Eigen::Vector2d leftPoint = block.edges[3](eta);
        Eigen::Vector2d rightPoint = block.edges[1](eta);

        for (int i = 1; i < block.nodesX - 1; ++i) {
            double xi = static_cast<double>(i) / (block.nodesX - 1);

            // Линейная интерполяция между левой и правой гранями
            Eigen::Vector2d coords = leftPoint * (1 - xi) + rightPoint * xi;

            int nodeId = nextNodeId_++;
            assembly_->addNode(std::make_shared<Node>(nodeId, coords.x(), coords.y()));
            nodeGrid[j][i] = nodeId;
        }
    }

    // 2. Генерация элементов
    int elementsGenerated = 0;
    for (int j = 0; j < block.nodesY - 1; ++j) {
        for (int i = 0; i < block.nodesX - 1; ++i) {
            // Проверяем, что все 4 узла элемента существуют
            if (nodeGrid[j][i] > 0 && nodeGrid[j][i + 1] > 0 &&
                nodeGrid[j + 1][i + 1] > 0 && nodeGrid[j + 1][i] > 0) {

                // Четырехугольный элемент
                std::vector<int> nodeIds = {
                    nodeGrid[j][i],        // нижний-левый
                    nodeGrid[j][i + 1],      // нижний-правый
                    nodeGrid[j + 1][i + 1],    // верхний-правый
                    nodeGrid[j + 1][i]       // верхний-левый
                };

                auto element = std::make_shared<PlaneIsoparametricElement>(
                    nextElementId_++, nodeIds, block.materialId);

                assembly_->addElement(element);
                elementsGenerated++;
            }
        }
    }

    std::cout << "Generated " << elementsGenerated << " elements" << std::endl;
    std::cout << "Node IDs from " << nextNodeId_ - block.nodesX * block.nodesY
        << " to " << nextNodeId_ - 1 << std::endl;
}

void MeshGenerator::createAnnulusSimple(const Eigen::Vector2d& center,
    double innerRadius, double outerRadius,
    double startAngle, double endAngle,
    int radialLayers, int circumferentialNodes,
    int materialId) {
    std::cout << "Creating annulus mesh (simple method): "
        << radialLayers << " radial layers, "
        << circumferentialNodes << " nodes per layer" << std::endl;

    

    // Проверка параметров
    if (radialLayers < 1 || circumferentialNodes < 3) {
        throw std::invalid_argument("Invalid mesh parameters");
    }

    if (innerRadius >= outerRadius) {
        throw std::invalid_argument("Inner radius must be less than outer radius");
    }

    // Преобразуем углы в радианы, если нужно
    if (startAngle > TWO_PI || endAngle > TWO_PI) {
        startAngle *= DEG_TO_RAD;
        endAngle *= DEG_TO_RAD;
    }

    // Вычисляем угловой шаг
    double totalAngle = endAngle - startAngle;
    if (totalAngle <= 0) {
        totalAngle += TWO_PI;
    }
    double angularStep = totalAngle / (circumferentialNodes - 1);

    // Вычисляем радиальный шаг
    double radialStep = (outerRadius - innerRadius) / radialLayers;

    // 1. СОЗДАНИЕ УЗЛОВ по слоям (от внутреннего радиуса к внешнему)
    std::vector<std::vector<int>> nodeGrid(radialLayers + 1,
        std::vector<int>(circumferentialNodes, 0));

    for (int layer = 0; layer <= radialLayers-1; ++layer) {
        double currentRadius = innerRadius + layer * radialStep;

        for (int node = 0; node < circumferentialNodes; ++node) {
            double angle = startAngle + node * angularStep;

            // Вычисляем координаты
            double x = center.x() + currentRadius * std::cos(angle);
            double y = center.y() + currentRadius * std::sin(angle);

            // Создаем узел
            int nodeId = nextNodeId_++;
            assembly_->addNode(std::make_shared<Node>(nodeId, x, y));
            nodeGrid[layer][node] = nodeId;

            std::cout << "Node " << nodeId << ": layer=" << layer
                << ", angle=" << angle * RAD_TO_DEG << "°, "
                << "r=" << currentRadius << " -> ("
                << x << ", " << y << ")" << std::endl;
        }
    }

    // 2. СОЗДАНИЕ ЭЛЕМЕНТОВ (четырехугольников между слоями)
    int elementsCreated = 0;

    for (int layer = 0; layer < radialLayers-1; ++layer) {
        for (int segment = 0; segment < circumferentialNodes - 1; ++segment) {
            // Четырехугольный элемент между двумя радиальными слоями
            std::vector<int> nodeIds = {
               nodeGrid[layer + 1][segment],       // Внешний-левый
               nodeGrid[layer + 1][segment + 1],   // Внешний-правый
               nodeGrid[layer][segment + 1],       // Внутренний-правый
               nodeGrid[layer][segment]          // Внутренний-левый
            };

            // Создаем элемент
            auto element = std::make_shared<PlaneIsoparametricElement>(
                nextElementId_++, nodeIds, materialId);

            assembly_->addElement(element);
            elementsCreated++;

            std::cout << "Element " << nextElementId_ - 1 << ": nodes "
                << nodeIds[0] << "-" << nodeIds[1] << "-"
                << nodeIds[2] << "-" << nodeIds[3] << std::endl;
        }
    }

    // 3. ЗАМЫКАНИЕ ПОСЛЕДНЕГО СЕГМЕНТА (если полный круг)
    if (std::abs(totalAngle - TWO_PI) < 1e-6) {
        for (int layer = 0; layer < radialLayers; ++layer) {
            std::vector<int> nodeIds = {
                nodeGrid[layer][circumferentialNodes - 1],  // Последний узел слоя
                nodeGrid[layer][0],                         // Первый узел слоя
                nodeGrid[layer + 1][0],                     // Первый узел след. слоя
                nodeGrid[layer + 1][circumferentialNodes - 1] // Последний узел след. слоя
            };

            auto element = std::make_shared<PlaneIsoparametricElement>(
                nextElementId_++, nodeIds, materialId);

            assembly_->addElement(element);
            elementsCreated++;
        }
    }

    std::cout << "Annulus mesh created: " << elementsCreated << " elements, "
        << (radialLayers ) * circumferentialNodes << " nodes" << std::endl;
}

std::vector<int> MeshGenerator::findContactNodes(double contactCenterX,
    double contactHalfWidth,
    double maxYtolerance) const {
    std::vector<int> contactNodeIds;

    auto nodes = assembly_->getNodes();

    for (const auto& node : nodes) {
        Eigen::Vector2d coords = node->getCoordinates();
        double x = coords.x();
        double y = coords.y();

        // Проверяем, находится ли узел в зоне контакта:
        // 1. По X: в пределах ±contactHalfWidth от центра
        // 2. По Y: близко к минимальной Y-координате (нижние узлы)

        if (std::abs(x - contactCenterX) <= contactHalfWidth) {
            // Находим минимальную Y координату во всей сетке
            static double minY = 1e9;
            if (y < minY) minY = y;

            // Проверяем, близко ли узел к нижней границе
            if (std::abs(y - minY) <= maxYtolerance) {
                contactNodeIds.push_back(node->getId());
            }
        }
    }

    std::cout << "Found " << contactNodeIds.size()
        << " nodes in contact zone (x ∈ ["
        << contactCenterX - contactHalfWidth << ", "
        << contactCenterX + contactHalfWidth << "])" << std::endl;

    return contactNodeIds;
}

void MeshGenerator::applyParabolicContactToNodes(double maxPressure,
    double contactHalfWidth,
    double contactCenterX,
    double totalForce) {
    std::cout << "\n=== Applying PARABOLIC contact to nodes ===" << std::endl;

    // Используем готовую функцию параболического распределения
    auto parabolicLoad = LoadFunction::parabolicPressure(maxPressure,
        contactHalfWidth,
        contactCenterX);

    // 1. Находим нижние узлы в зоне контакта
    auto contactNodeIds = findContactNodes(contactCenterX,
        contactHalfWidth,
        0.01); // 1 см допуск

    if (contactNodeIds.empty()) {
        std::cout << "WARNING: No contact nodes found!" << std::endl;
        return;
    }

    std::cout << "Found " << contactNodeIds.size()
        << " nodes in contact zone" << std::endl;

    // 2. Вычисляем давление в каждом узле по параболическому закону
    std::vector<double> nodalPressures;
    double integralPressure = 0.0;

    for (int nodeId : contactNodeIds) {
        auto node = assembly_->getNode(nodeId);
        Eigen::Vector2d coords = node->getCoordinates();

        // Нормаль вниз (0, -1) для плоской поверхности
        Eigen::Vector2d normal(0, -1);

        // Вычисляем давление в этой точке по параболическому закону
        Eigen::Vector2d pressureVec = parabolicLoad.distribution_(
            coords.x(), coords.y(), normal);

        double pressure = std::abs(pressureVec.y()); // Вертикальная компонента

        nodalPressures.push_back(pressure);
        integralPressure += pressure;

        std::cout << "Node " << nodeId << " (x=" << coords.x()
            << "): pressure = " << pressure / 1e6 << " MPa" << std::endl;
    }

    // 3. Нормируем так, чтобы суммарная сила была равна totalForce
    if (integralPressure > 0) {
        double scaleFactor = totalForce / integralPressure;
        std::cout << "Scale factor for normalization: " << scaleFactor << std::endl;

        // 4. Прикладываем силы к узлам
        for (size_t i = 0; i < contactNodeIds.size(); ++i) {
            double scaledPressure = nodalPressures[i] * scaleFactor;
            int nodeId = contactNodeIds[i];

            // Создаем сосредоточенную силу
            auto force = std::make_shared<ConcentratedForce>(
                nodeId, 0.0, -scaledPressure); // Вертикально вниз

            assembly_->addConcentratedForce(force);

            std::cout << "Applied force to node " << nodeId
                << ": " << -scaledPressure << " N" << std::endl;
        }

        // 5. Проверяем баланс сил
        double appliedForceSum = 0;
        for (size_t i = 0; i < contactNodeIds.size(); ++i) {
            appliedForceSum += nodalPressures[i] * scaleFactor;
        }

        std::cout << "Total applied force: " << appliedForceSum << " N" << std::endl;
        std::cout << "Target force: " << totalForce << " N" << std::endl;
        std::cout << "Difference: " << std::abs(appliedForceSum - totalForce)
            << " N ("
            << std::abs(appliedForceSum - totalForce) / totalForce * 100
            << "%)" << std::endl;
    }
}