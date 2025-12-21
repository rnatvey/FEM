#include "meshgenerator.h"
#include "planeisometric/Planeisoparametric.h"
#include <iostream>
#include <cmath>

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

    // 1. Генерация узлов
    std::vector<std::vector<int>> nodeGrid(block.nodesY, std::vector<int>(block.nodesX, 0));

    for (int j = 0; j < block.nodesY; ++j) {
        double eta = static_cast<double>(j) / (block.nodesY - 1);

        for (int i = 0; i < block.nodesX; ++i) {
            double xi = static_cast<double>(i) / (block.nodesX - 1);

            // Вычисляем координаты через трансфинитную интерполяцию
            Eigen::Vector2d coords = transfiniteInterpolation(block, xi, eta);

            // Добавляем узел
            int nodeId = nextNodeId_++;
            assembly_->addNode(std::make_shared<Node>(nodeId, coords.x(), coords.y()));
            nodeGrid[j][i] = nodeId;
        }
    }

    // 2. Генерация элементов
    for (int j = 0; j < block.nodesY - 1; ++j) {
        for (int i = 0; i < block.nodesX - 1; ++i) {
            // Четырехугольный элемент
            std::vector<int> nodeIds = {
                nodeGrid[j][i],      // нижний-левый
                nodeGrid[j][i + 1],    // нижний-правый
                nodeGrid[j + 1][i + 1],  // верхний-правый
                nodeGrid[j + 1][i]     // верхний-левый
            };

            auto element = std::make_shared<PlaneIsoparametricElement>(
                nextElementId_++, nodeIds, block.materialId);

            assembly_->addElement(element);
        }
    }

    std::cout << "Generated " << (block.nodesX - 1) * (block.nodesY - 1) << " elements" << std::endl;
}