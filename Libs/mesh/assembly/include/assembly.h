#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "node.h"
#include "FiniteStrainElementResponse.h"
#include "FiniteStrainMaterial.h"
#include "material.h"
#include "BaseElement.h"
#include "planeisometric/FiniteStrainQ4Element.h"
#include "vector.h"
#include "ConcentratedForce.h"
#include "loadFunction.h"
//class ConcentratedForce;
class vector2;
class LoadFunction;
class Assembly : public std::enable_shared_from_this<Assembly> {
public:
    struct ConstraintData {
        std::vector<int> fixedZeroDofs;
        std::vector<int> constrainedDofs;
        std::vector<int> prescribedDofs;
        std::vector<double> prescribedValues;
    };

    struct BoundaryConditionApplicationStats {
        double totalTimeSeconds = 0.0;
        double reduceSystemTimeSeconds = 0.0;
        Eigen::Index inputMatrixNonZeros = 0;
        Eigen::Index outputMatrixNonZeros = 0;
        int constrainedDofCount = 0;
        int activeDofCount = 0;
    };

    Assembly();
    ~Assembly() = default;

    struct SurfaceLoad
    {
        int elementId;
        int surfaceIndex; // 0-нижняя, 1-правая, 2-верхняя, 3-левая
        std::shared_ptr<LoadFunction> loadFunction;
        SurfaceLoad(int elemId, int surfIdx, std::shared_ptr<LoadFunction> loadFunc)
            : elementId(elemId), surfaceIndex(surfIdx), loadFunction(std::move(loadFunc)) {
        }
    };

    // === Управление узлами ===
    void addNode(std::shared_ptr<Node> node);
    void addNodes(const std::vector<std::shared_ptr<Node>>& nodes);
    void addLineofNodes(size_t startId,  vector2 start_point, vector2 end_point, int num_of_nodes);
    std::shared_ptr<Node> getNode(int id) const;
    const std::vector<std::shared_ptr<Node>>& getNodes() const { return nodes_; }
    int getNodeCount() const { return nodes_.size(); }


    // === Управление материалами ===
    void addMaterial(std::shared_ptr<Material> material);
    std::shared_ptr<Material> getMaterial(int id) const;
    const std::unordered_map<int, std::shared_ptr<Material>>& getMaterials() const { return materials_; }

    void addFiniteStrainMaterial(std::shared_ptr<FiniteStrainMaterial> material);
    std::shared_ptr<FiniteStrainMaterial> getFiniteStrainMaterial(int id) const;
    const std::unordered_map<int, std::shared_ptr<FiniteStrainMaterial>>&
        getFiniteStrainMaterials() const { return finiteStrainMaterials_; }

    // === Управление элементами ===
    void addElement(std::shared_ptr<BaseElement> element);
    void addElements(const std::vector<std::shared_ptr<BaseElement>>& elements);
    std::shared_ptr<BaseElement> getElement(int id) const;
    const std::vector<std::shared_ptr<BaseElement>>& getElements() const { return elements_; }
    int getElementCount() const { return elements_.size(); }

    void addFiniteStrainElement(std::shared_ptr<FiniteStrainQ4Element> element);
    void addFiniteStrainElements(
        const std::vector<std::shared_ptr<FiniteStrainQ4Element>>& elements);
    std::shared_ptr<FiniteStrainQ4Element> getFiniteStrainElement(int id) const;
    const std::vector<std::shared_ptr<FiniteStrainQ4Element>>&
        getFiniteStrainElements() const { return finiteStrainElements_; }
    bool hasFiniteStrainModel() const { return !finiteStrainElements_.empty(); }

    // === Сборка глобальной системы ===
    void assembleGlobalStiffnessMatrix(Eigen::SparseMatrix<double>& globalK) const;
    void assembleGlobalForceVector(Eigen::VectorXd& globalF, const Eigen::VectorXd& bodyForces) const;
    void assembleFiniteStrainSystem(const Eigen::VectorXd& fullDisplacements,
        Eigen::SparseMatrix<double>& globalTangent,
        Eigen::VectorXd& globalInternalForce,
        double& totalStrainEnergy,
        std::vector<FiniteStrainElementResponse>* elementResponses = nullptr,
        bool computeTangent = true) const;

    // === Граничные условия ===
    void addFixedNode(int nodeId, bool fixX = true, bool fixY = true);
    void addPrescribedDisplacement(int nodeId, double dx, double dy);
    void addPrescribedDisplacementX(int nodeId, double dx);
    void addPrescribedDisplacementY(int nodeId, double dy);
    void applyBoundaryConditions(Eigen::SparseMatrix<double>& globalK,
        Eigen::VectorXd& globalF,
        BoundaryConditionApplicationStats* stats = nullptr) const;

    // === Система уравнений ===
    int getTotalDofCount() const { return nodes_.size() * 2; }
    std::vector<int> getElementDofIndices(int elementId) const;

    // Для сборки матрицы жесткости (до граничных условий) - использует полные индексы
    std::vector<int> getElementFullDofIndices(int elementId) const;

    // Для работы с решениями (после граничных условий) - использует сокращенные индексы
    std::vector<int> getElementReducedDofIndices(int elementId) const;


    // === Валидация ===
    bool validate() const;

    int getGlobalDofIndex(int nodeId, int direction) const; // 0-x, 1-y

    struct DofMapping {
        std::vector<int> fullToReduced;  // fullDof -> reducedDof (-1 если закреплен)
        std::vector<int> reducedToFull;  // reducedDof -> fullDof
        std::vector<int> prescribedDofs; // Индексы DOF с предписанными перемещениями
        std::vector<double> prescribedValues; // Значения предписанных перемещений
    };
    // Силы
    void addConcentratedForce(std::shared_ptr<ConcentratedForce> force);
    void addConcentratedForces(const std::vector<std::shared_ptr<ConcentratedForce>>& forces);
    const std::vector<std::shared_ptr<ConcentratedForce>>& getConcentratedForces() const { return concentratedForces_; }
    void addSurfaceLoad(int elementId, int surfaceIndex, std::shared_ptr<LoadFunction> loadFunction);
    const std::vector<SurfaceLoad>& getSurfaceLoads() const { return surfaceLoads_; }
    // Сборка вектора сосредоточенных сил
    void assembleConcentratedForces(Eigen::VectorXd& globalF) const;
    void assembleSurfaceLoads(Eigen::VectorXd& globalF) const;

    //=========================================

    const DofMapping& getDofMapping() const { return dofMapping_; }
    ConstraintData getConstraintData() const;
    void initializeDofMapping(int totalDof) const;

private:
    std::vector<std::shared_ptr<Node>> nodes_;
    std::unordered_map<int, std::shared_ptr<Material>> materials_;
    std::unordered_map<int, std::shared_ptr<FiniteStrainMaterial>> finiteStrainMaterials_;
    std::vector<std::shared_ptr<BaseElement>> elements_;
    std::vector<std::shared_ptr<FiniteStrainQ4Element>> finiteStrainElements_;
    std::vector<std::shared_ptr<ConcentratedForce>> concentratedForces_;
    std::vector<int> getElementDofIndicesInternal(int elementId) const;
    // Граничные условия
    struct BoundaryCondition {
        int nodeId;
        bool fixX, fixY;
        double prescribedDx, prescribedDy;
        bool hasPrescribedDisplacement;
    };
    std::vector<BoundaryCondition> boundaryConditions_;
    std::vector<SurfaceLoad> surfaceLoads_;
    // Вспомогательные методы
   
    void buildNodeIndexMap();

    std::unordered_map<int, int> nodeIdToIndex_; // Быстрый поиск узлов по ID

    mutable DofMapping dofMapping_;
    ConstraintData buildConstraintData() const;
    void populateDofMapping(int totalDof, const ConstraintData& constraintData) const;

};

