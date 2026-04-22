#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>

#include "AppRuntimeSupport.h"
#include "HyperelasticScenarioSupport.h"
#include "ResultFileExporter.h"

namespace {

constexpr int kNodesX = 6;
constexpr int kNodesY = 6;
constexpr double kYoungsModulus = 25.0;
constexpr double kPoissonsRatio = 0.48;
constexpr double kThickness = 1.0;
constexpr double kPrescribedTopEdgeDy = -0.10;
constexpr double kMaximumPrescribedDisplacementErrorTolerance = 1.0e-12;
constexpr double kMinimumJacobianDeterminant = 0.20;
constexpr double kMaximumJacobianDeterminant = 2.50;
constexpr double kEquilibriumResidualTolerance = 1.0e-5;

} // namespace

int main() {
    try {
        const Eigen::Vector2d lowerCorner(-0.5, 0.05);
        const Eigen::Vector2d upperCorner(0.5, 1.05);

        const auto hyperelasticCase = HyperelasticScenarioSupport::buildFiniteStrainBlockCase(
            kNodesX,
            kNodesY,
            lowerCorner,
            upperCorner,
            kYoungsModulus,
            kPoissonsRatio,
            kThickness,
            kPrescribedTopEdgeDy);

        FEModel model;
        model.setAssembly(hyperelasticCase.assembly);
        model.setSolverTolerance(1.0e-8);
        model.setMaxIterations(35);
        if (!model.solveHyperelastic()) {
            std::cerr << "Hyperelastic large-deformation no-contact solve failed"
                      << std::endl;
            return 1;
        }

        const double maxTopEdgeDisplacementError =
            HyperelasticScenarioSupport::computeMaximumTopEdgeDisplacementError(
                model, hyperelasticCase.assembly, kPrescribedTopEdgeDy);
        const auto diagnostics =
            HyperelasticScenarioSupport::evaluateFiniteStrainDiagnostics(model);
        const auto reactionForceSum =
            HyperelasticScenarioSupport::computeReactionForceSum(model);
        const auto& metrics = model.getPerformanceMetrics();

        if (maxTopEdgeDisplacementError > kMaximumPrescribedDisplacementErrorTolerance) {
            std::cerr << "Large-deformation top-edge prescribed displacement mismatch: "
                      << maxTopEdgeDisplacementError << std::endl;
            return 1;
        }
        if (diagnostics.gaussPointCount <= 0 ||
            diagnostics.minJ <= kMinimumJacobianDeterminant ||
            diagnostics.maxJ >= kMaximumJacobianDeterminant) {
            std::cerr << "Large-deformation scenario has unacceptable det(F) range: ["
                      << diagnostics.minJ << ", " << diagnostics.maxJ << "]" << std::endl;
            return 1;
        }
        if (diagnostics.maxStrainEnergyDensity <= 0.0 ||
            diagnostics.maxSecondPiolaStressNorm <= 0.0) {
            std::cerr << "Large-deformation scenario produced non-physical finite-strain response"
                      << std::endl;
            return 1;
        }
        if (!metrics.hasNearIncompressibleFiniteStrainMaterial) {
            std::cerr << "Large-deformation scenario did not flag nearly incompressible material"
                      << std::endl;
            return 1;
        }
        if (metrics.equilibriumResidualNorm > kEquilibriumResidualTolerance) {
            std::cerr << "Large-deformation equilibrium residual is above tolerance: "
                      << metrics.equilibriumResidualNorm << std::endl;
            return 1;
        }

        ResultFileExportOptions exportOptions;
        exportOptions.outputDirectory =
            AppRuntimeSupport::caseOutputDirectory(
                "hyperelastic_large_deformation_no_contact");
        exportOptions.baseName = "solution";
        exportOptions.extraStringMetrics = {
            {"case_name", "hyperelastic_large_deformation_no_contact"},
            {"verification_category", "hyperelastic_large_deformation"}
        };
        exportOptions.extraNumericMetrics = {
            {"prescribed_top_edge_dy", kPrescribedTopEdgeDy},
            {"top_edge_node_count", static_cast<double>(hyperelasticCase.topEdgeNodeCount)},
            {"anchor_node_id", static_cast<double>(hyperelasticCase.anchorNodeId)},
            {"minimum_jacobian_determinant", diagnostics.minJ},
            {"maximum_jacobian_determinant", diagnostics.maxJ},
            {"maximum_strain_energy_density", diagnostics.maxStrainEnergyDensity},
            {"maximum_second_piola_stress_norm", diagnostics.maxSecondPiolaStressNorm},
            {"reaction_force_sum_x", reactionForceSum.x()},
            {"reaction_force_sum_y", reactionForceSum.y()}
        };

        const auto exportArtifacts = ResultFileExporter::exportSolution(model, exportOptions);
        AppRuntimeSupport::runPostprocessorIfEnabled(exportOptions.outputDirectory);

        std::cout << "Hyperelastic large-deformation no-contact scenario" << std::endl;
        std::cout << "success=true" << std::endl;
        std::cout << "nodes=" << hyperelasticCase.assembly->getNodeCount() << std::endl;
        std::cout << "elements=" << hyperelasticCase.assembly->getFiniteStrainElements().size()
                  << std::endl;
        std::cout << "maximum_prescribed_displacement_error="
                  << maxTopEdgeDisplacementError << std::endl;
        std::cout << "minimum_jacobian_determinant=" << diagnostics.minJ << std::endl;
        std::cout << "maximum_jacobian_determinant=" << diagnostics.maxJ << std::endl;
        std::cout << "maximum_strain_energy_density="
                  << diagnostics.maxStrainEnergyDensity << std::endl;
        std::cout << "maximum_second_piola_stress_norm="
                  << diagnostics.maxSecondPiolaStressNorm << std::endl;
        std::cout << "strain_energy=" << metrics.strainEnergy << std::endl;
        std::cout << "near_incompressible_material="
                  << (metrics.hasNearIncompressibleFiniteStrainMaterial ? 1 : 0)
                  << std::endl;
        std::cout << "direct_linear_solve_count=" << metrics.directLinearSolveCount
                  << std::endl;
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
        std::cerr << "hyperelastic_large_deformation_no_contact_failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
