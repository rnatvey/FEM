#pragma once
#include "assembly.h"
#include "geometry.h"

class MeshGenerator {
public:
    MeshGenerator(std::shared_ptr<Assembly> assembly);

    // === Блочное построение сетки (метод блоков по Бате, стр. 392) ===

    // Создание блока (подобласти) с трансфинитной интерполяцией
    struct Block {
        std::vector<Geometry::ParametricCurve> edges; // 4 кривые для квадрилатерального блока
        int nodesX, nodesY; // Количество узлов по параметрам
        int materialId;
    };

    void addBlock(const Block& block);

    // === Автоматическая генерация типовых областей ===

    // Прямоугольная область с регулярной сеткой
    void createRectangle(const Eigen::Vector2d& corner1, const Eigen::Vector2d& corner2,
        int nodesX, int nodesY, int materialId);

    // Кольцевая область (для шины)
    void createAnnulus(const Eigen::Vector2d& center, double innerRadius, double outerRadius,
        double startAngle, double endAngle,
        int radialLayers, int circumferentialNodes, int materialId);

    // Треугольная область (для адаптивных сеток)
    void createTriangle(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3,
        int divisions, int materialId);

    // === Управление нумерацией ===
    void setNodeIdStart(int startId) { nextNodeId_ = startId; }
    void setElementIdStart(int startId) { nextElementId_ = startId; }

private:
    std::shared_ptr<Assembly> assembly_;
    int nextNodeId_ = 1;
    int nextElementId_ = 1;

    // Генерация сетки в блоке методом трансфинитной интерполяции
    void generateBlockMesh(const Block& block);

    // Создание изопараметрического отображения (по Бате, стр. 395)
    Eigen::Vector2d transfiniteInterpolation(const Block& block, double xi, double eta) const;
};