#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>

#include "AppRuntimeSupport.h"
#include "HyperelasticScenarioSupport.h"
#include "ResultFileExporter.h"

namespace {

constexpr int kNodesX = 7;
constexpr int kNodesY = 7;
constexpr double kYoungsModulus = 2.0e5;
constexpr double kPoissonsRatio = 0.30;
constexpr double kThickness = 1.0;
constexpr double kPrescribedTopEdgeDy = -0.005;
constexpr double kMaximumPrescribedDisplacementErrorTolerance = 1.0e-12;
constexpr double kRelativeDisplacementDifferenceTolerance = 5.0e-3;
constexpr double kMinimumJacobianDeterminant = 0.95;
constexpr double kEquilibriumResidualTolerance = 1.0e-5;

} // namespace

int main() {
    try {
        const Eigen::Vector2d lowerCorner(-0.5, 0.05);
        const Eigen::Vector2d upperCorner(0.5, 1.05);

        const auto linearCase = HyperelasticScenarioSupport::buildLinearBlockCase(
            kNodesX,
            kNodesY,
            lowerCorner,
            upperCorner,
            kYoungsModulus,
            kPoissonsRatio,
            kThickness,
            kPrescribedTopEdgeDy);
        const auto hyperelasticCase = HyperelasticScenarioSupport::buildFiniteStrainBlockCase(
            kNodesX,
            kNodesY,
            lowerCorner,
            upperCorner,
            kYoungsModulus,
            kPoissonsRatio,
            kThickness,
            kPrescribedTopEdgeDy);

        FEModel linearModel;
        linearModel.setAssembly(linearCase.assembly);
        linearModel.setSolverTolerance(1.0e-8);
        if (!linearModel.solve()) {
            std::cerr << "Linear baseline solve failed in hyperelastic sanity scenario"
                      << std::endl;
            return 1;
        }

        FEModel hyperelasticModel;
        hyperelasticModel.setAssembly(hyperelasticCase.assembly);
        hyperelasticModel.setSolverTolerance(1.0e-8);
        hyperelasticModel.setMaxIterations(20);
        if (!hyperelasticModel.solveHyperelastic()) {
            std::cerr << "Hyperelastic no-contact sanity solve failed" << std::endl;
            return 1;
        }

        const double maxTopEdgeDisplacementError =
            HyperelasticScenarioSupport::computeMaximumTopEdgeDisplacementError(
                hyperelasticModel, hyperelasticCase.assembly, kPrescribedTopEdgeDy);
        const Eigen::VectorXd displacementDifference =
            hyperelasticModel.getDisplacements() - linearModel.getDisplacements();
        const double relativeDisplacementDifference =
            displacementDifference.norm() /
            std::max(linearModel.getDisplacements().norm(), 1.0);
        const auto diagnostics =
            HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(hyperelasticModel);
        const auto reactionForceSum =
            HyperelasticScenarioSupport::computeReactionForceSum(hyperelasticModel);
        const auto& metrics = hyperelasticModel.getPerformanceMetrics();

        if (maxTopEdgeDisplacementError > kMaximumPrescribedDisplacementErrorTolerance) {
            std::cerr << "Hyperelastic top-edge prescribed displacement mismatch: "
                      << maxTopEdgeDisplacementError << std::endl;
            return 1;
        }
        if (relativeDisplacementDifference > kRelativeDisplacementDifferenceTolerance) {
            std::cerr << "Hyperelastic solution deviates too much from small-strain baseline: "
                      << relativeDisplacementDifference << std::endl;
            return 1;
        }
        if (diagnostics.gaussPointCount <= 0 ||
            diagnostics.minJ <= kMinimumJacobianDeterminant) {
            std::cerr << "Hyperelastic sanity scenario has unacceptable det(F): "
                      << diagnostics.minJ << std::endl;
            return 1;
        }
        if (metrics.equilibriumResidualNorm > kEquilibriumResidualTolerance) {
            std::cerr << "Hyperelastic equilibrium residual is above tolerance: "
                      << metrics.equilibriumResidualNorm << std::endl;
            return 1;
        }

        ResultFileExportOptions exportOptions;
        exportOptions.outputDirectory =
            AppRuntimeSupport::caseOutputDirectory("hyperelastic_no_contact_sanity");
        exportOptions.baseName = "solution";
        exportOptions.extraStringMetrics = {
            {"case_name", "hyperelastic_no_contact_sanity"},
            {"verification_category", "hyperelastic_baseline"},
            {"reference_solution", "small_strain_linear"}
        };
        exportOptions.extraNumericMetrics = {
            {"prescribed_top_edge_dy", kPrescribedTopEdgeDy},
            {"top_edge_node_count", static_cast<double>(hyperelasticCase.topEdgeNodeCount)},
            {"anchor_node_id", static_cast<double>(hyperelasticCase.anchorNodeId)},
            {"relative_displacement_difference_to_linear", relativeDisplacementDifference},
            {"maximum_prescribed_displacement_error", maxTopEdgeDisplacementError},
            {"minimum_jacobian_determinant", diagnostics.minJ},
            {"maximum_jacobian_determinant", diagnostics.maxJ},
            {"maximum_strain_energy_density", diagnostics.maxStrainEnergyDensity},
            {"reaction_force_sum_x", reactionForceSum.x()},
            {"reaction_force_sum_y", reactionForceSum.y()}
        };

        const auto exportArtifacts =
            ResultFileExporter::exportSolution(hyperelasticModel, exportOptions);
        AppRuntimeSupport::runPostprocessorIfEnabled(exportOptions.outputDirectory);

        std::cout << "Hyperelastic no-contact sanity scenario" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "nodes=" << hyperelasticCase.assembly->getNodeCount() << std::endl;
        std::cout << "elements=" << hyperelasticCase.assembly->getFiniteStrainElements().size()
                  << std::endl;
        std::cout << "relative_displacement_difference_to_linear="
                  << relativeDisplacementDifference << std::endl;
        std::cout << "maximum_prescribed_displacement_error="
                  << maxTopEdgeDisplacementError << std::endl;
        std::cout << "minimum_jacobian_determinant=" << diagnostics.minJ << std::endl;
        std::cout << "maximum_jacobian_determinant=" << diagnostics.maxJ << std::endl;
        std::cout << "strain_energy=" << metrics.strainEnergy << std::endl;
        std::cout << "nonlinear_iterations=" << metrics.nonlinearIterations << std::endl;
        std::cout << "linear_iterations=" << metrics.linearIterations << std::endl;
        std::cout << "equilibrium_residual_norm=" << metrics.equilibriumResidualNorm
                  << std::endl;
        std::cout << "results_vtu=" << exportArtifacts.vtuPath.string() << std::endl;
        std::cout << "results_metrics_json=" << exportArtifacts.metricsJsonPath.string()
                  << std::endl;
        return 0;
    }
    catch (const std::exception& exception) {
        std::cerr << "hyperelastic_no_contact_sanity_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
