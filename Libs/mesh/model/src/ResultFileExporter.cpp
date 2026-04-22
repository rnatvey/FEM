#include "ResultFileExporter.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

#include "FEMModel.h"
#include "NeoHookeanMaterial.h"
#include "assembly.h"

namespace {

std::string formatNumber(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

std::string formatInteger(long long value) {
    return std::to_string(value);
}

std::string jsonEscape(std::string_view text) {
    std::string escaped;
    escaped.reserve(text.size() + 8);

    for (char character : text) {
        switch (character) {
        case '\\':
            escaped += "\\\\";
            break;
        case '\"':
            escaped += "\\\"";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            escaped += character;
            break;
        }
    }

    return escaped;
}

std::string xmlEscape(std::string_view text) {
    std::string escaped;
    escaped.reserve(text.size() + 8);

    for (char character : text) {
        switch (character) {
        case '&':
            escaped += "&amp;";
            break;
        case '<':
            escaped += "&lt;";
            break;
        case '>':
            escaped += "&gt;";
            break;
        case '\"':
            escaped += "&quot;";
            break;
        case '\'':
            escaped += "&apos;";
            break;
        default:
            escaped += character;
            break;
        }
    }

    return escaped;
}

double vectorMagnitude(const Eigen::Vector2d& vector) {
    return vector.norm();
}

double resolveRepresentativePoissonsRatio(const Assembly& assembly) {
    const auto& elements = assembly.getElements();
    if (!elements.empty()) {
        double poissonsRatioSum = 0.0;
        int contributingElementCount = 0;
        for (const auto& element : elements) {
            const auto material = assembly.getMaterial(element->getMaterialId());
            if (!material) {
                continue;
            }

            poissonsRatioSum += material->getPoissonsRatio();
            ++contributingElementCount;
        }

        if (contributingElementCount > 0) {
            return poissonsRatioSum / static_cast<double>(contributingElementCount);
        }
    }

    const auto& materials = assembly.getMaterials();
    if (!materials.empty()) {
        double poissonsRatioSum = 0.0;
        int materialCount = 0;
        for (const auto& [materialId, material] : materials) {
            (void)materialId;
            if (!material) {
                continue;
            }

            poissonsRatioSum += material->getPoissonsRatio();
            ++materialCount;
        }

        if (materialCount > 0) {
            return poissonsRatioSum / static_cast<double>(materialCount);
        }
    }

    return 0.0;
}

double computePlaneStrainSigmaZZ(const Eigen::Vector3d& stress2d, double poissonsRatio) {
    return poissonsRatio * (stress2d.x() + stress2d.y());
}

double computeVonMisesStress(const Eigen::Vector3d& stress2d, double sigmaZZ) {
    const double sxx = stress2d.x();
    const double syy = stress2d.y();
    const double txy = stress2d.z();
    return std::sqrt(std::max(
        0.0,
        sxx * sxx - sxx * syy + syy * syy
            + sigmaZZ * sigmaZZ - syy * sigmaZZ - sxx * sigmaZZ
            + 3.0 * txy * txy));
}

Eigen::Matrix2d voigtStressToTensor(const Eigen::Vector3d& stressVoigt) {
    Eigen::Matrix2d stress = Eigen::Matrix2d::Zero();
    stress(0, 0) = stressVoigt.x();
    stress(1, 1) = stressVoigt.y();
    stress(0, 1) = stressVoigt.z();
    stress(1, 0) = stressVoigt.z();
    return stress;
}

struct FiniteStrainCellSummary {
    int elementId = -1;
    int materialId = -1;
    double meanJ = 0.0;
    double minJ = 1.0;
    double maxJ = 1.0;
    double meanStrainEnergyDensity = 0.0;
    double meanSigmaZZ = 0.0;
    Eigen::Vector3d meanGreenLagrangeStrain = Eigen::Vector3d::Zero();
    Eigen::Vector3d meanSecondPiolaStress = Eigen::Vector3d::Zero();
    Eigen::Vector3d meanCauchyStress = Eigen::Vector3d::Zero();
};

double finiteStrainSigmaZZ(const FiniteStrainGaussPointResult& gaussPoint,
    const FiniteStrainMaterial& material) {
    const auto* neoHookeanMaterial = dynamic_cast<const NeoHookeanMaterial*>(&material);
    if (neoHookeanMaterial == nullptr) {
        return 0.0;
    }

    const double J = gaussPoint.jacobianDeterminant;
    if (!(J > 0.0)) {
        return 0.0;
    }

    return neoHookeanMaterial->getLameFirstParameter() * std::log(J) / J;
}

FiniteStrainCellSummary summarizeFiniteStrainElementResponse(
    const FiniteStrainElementResponse& response,
    const FiniteStrainMaterial& material,
    int elementId,
    int materialId) {
    FiniteStrainCellSummary summary;
    summary.elementId = elementId;
    summary.materialId = materialId;

    if (response.gaussPointResults.empty()) {
        return summary;
    }

    summary.minJ = response.gaussPointResults.front().jacobianDeterminant;
    summary.maxJ = summary.minJ;

    for (const auto& gaussPoint : response.gaussPointResults) {
        const double J = gaussPoint.jacobianDeterminant;
        const Eigen::Matrix2d F = gaussPoint.deformationGradient;
        const Eigen::Matrix2d secondPiolaStress =
            voigtStressToTensor(gaussPoint.secondPiolaStressVoigt);
        Eigen::Matrix2d cauchyStress = Eigen::Matrix2d::Zero();
        if (J > 0.0) {
            cauchyStress = (1.0 / J) * F * secondPiolaStress * F.transpose();
        }

        summary.meanJ += J;
        summary.minJ = std::min(summary.minJ, J);
        summary.maxJ = std::max(summary.maxJ, J);
        summary.meanStrainEnergyDensity += gaussPoint.strainEnergyDensity;
        summary.meanGreenLagrangeStrain += gaussPoint.greenLagrangeStrainVoigt;
        summary.meanSecondPiolaStress += gaussPoint.secondPiolaStressVoigt;
        summary.meanCauchyStress += Eigen::Vector3d(
            cauchyStress(0, 0),
            cauchyStress(1, 1),
            cauchyStress(0, 1));
        summary.meanSigmaZZ += finiteStrainSigmaZZ(gaussPoint, material);
    }

    const double gaussPointCount =
        static_cast<double>(response.gaussPointResults.size());
    summary.meanJ /= gaussPointCount;
    summary.meanStrainEnergyDensity /= gaussPointCount;
    summary.meanGreenLagrangeStrain /= gaussPointCount;
    summary.meanSecondPiolaStress /= gaussPointCount;
    summary.meanCauchyStress /= gaussPointCount;
    summary.meanSigmaZZ /= gaussPointCount;
    return summary;
}

int vtkCellTypeForNodeCount(int nodeCount) {
    switch (nodeCount) {
    case 2:
        return 3;  // VTK_LINE
    case 3:
        return 5;  // VTK_TRIANGLE
    case 4:
        return 9;  // VTK_QUAD
    default:
        return 7;  // VTK_POLYGON
    }
}

template <typename ValueFormatter>
void writeScalarDataArray(std::ostream& stream,
    std::string_view name,
    const int count,
    ValueFormatter formatter) {
    stream << "        <DataArray type=\"Float64\" Name=\""
           << xmlEscape(name)
           << "\" format=\"ascii\">\n";
    stream << "          ";
    for (int i = 0; i < count; ++i) {
        if (i > 0) {
            stream << ' ';
        }
        stream << formatter(i);
    }
    stream << "\n";
    stream << "        </DataArray>\n";
}

template <typename ValueFormatter>
void writeIntDataArray(std::ostream& stream,
    std::string_view name,
    const int count,
    ValueFormatter formatter) {
    stream << "        <DataArray type=\"Int32\" Name=\""
           << xmlEscape(name)
           << "\" format=\"ascii\">\n";
    stream << "          ";
    for (int i = 0; i < count; ++i) {
        if (i > 0) {
            stream << ' ';
        }
        stream << formatter(i);
    }
    stream << "\n";
    stream << "        </DataArray>\n";
}

template <typename ValueFormatter>
void writeVector3DataArray(std::ostream& stream,
    std::string_view name,
    const int count,
    ValueFormatter formatter) {
    stream << "        <DataArray type=\"Float64\" Name=\""
           << xmlEscape(name)
           << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    stream << "          ";
    for (int i = 0; i < count; ++i) {
        const Eigen::Vector3d value = formatter(i);
        if (i > 0) {
            stream << ' ';
        }
        stream << formatNumber(value.x()) << ' '
               << formatNumber(value.y()) << ' '
               << formatNumber(value.z());
    }
    stream << "\n";
    stream << "        </DataArray>\n";
}

} // namespace

ResultFileExportArtifacts ResultFileExporter::exportSolution(
    const FEModel& model,
    const ResultFileExportOptions& options) {
    const auto assembly = model.getAssembly();
    if (!assembly) {
        throw std::runtime_error("Cannot export results without an assembly");
    }

    if (options.outputDirectory.empty()) {
        throw std::invalid_argument("Result export output directory cannot be empty");
    }

    std::filesystem::create_directories(options.outputDirectory);

    const std::string baseName =
        options.baseName.empty() ? std::string("solution") : options.baseName;
    const std::filesystem::path vtuPath = options.outputDirectory / (baseName + ".vtu");
    const std::filesystem::path metricsJsonPath = options.outputDirectory / "metrics.json";
    const std::filesystem::path contactFacetCsvPath =
        options.outputDirectory / "contact_facets.csv";

    const auto& nodes = assembly->getNodes();
    const auto& linearElements = assembly->getElements();
    const auto& finiteStrainElements = assembly->getFiniteStrainElements();
    const bool exportFiniteStrain =
        model.hasFiniteStrainModel() && !finiteStrainElements.empty();
    const auto finiteStrainResponses =
        exportFiniteStrain ? model.evaluateCurrentFiniteStrainElementResponses()
                           : std::vector<FiniteStrainElementResponse>{};
    const auto nodalDisplacements = model.getNodalDisplacements();
    const auto nodalReactionForces = model.getNodalReactionForces();
    const auto nodalContactForces = model.getNodalContactForces();
    const auto nodalSignedDistances = model.getNodalContactSignedDistances();
    const auto nodalPenetrations = model.getNodalContactPenetrations();
    const auto performanceMetrics = model.getPerformanceMetrics();
    const double representativePoissonsRatio = resolveRepresentativePoissonsRatio(*assembly);

    std::vector<Eigen::Vector3d> nodalStresses;
    std::vector<Eigen::Vector3d> nodalStrains;
    std::vector<double> nodalSigmaZZ(nodes.size(), 0.0);
    std::vector<double> nodalJacobianDeterminant(nodes.size(), 1.0);
    std::vector<double> nodalStrainEnergyDensity(nodes.size(), 0.0);
    std::vector<FiniteStrainCellSummary> finiteStrainCellSummaries;

    if (exportFiniteStrain) {
        nodalStresses.assign(nodes.size(), Eigen::Vector3d::Zero());
        nodalStrains.assign(nodes.size(), Eigen::Vector3d::Zero());
        nodalJacobianDeterminant.assign(nodes.size(), 0.0);
        nodalStrainEnergyDensity.assign(nodes.size(), 0.0);
        std::vector<int> nodalContributionCount(nodes.size(), 0);

        finiteStrainCellSummaries.reserve(finiteStrainElements.size());
        for (size_t elementIndex = 0;
            elementIndex < finiteStrainElements.size() &&
            elementIndex < finiteStrainResponses.size();
            ++elementIndex) {
            const auto& element = finiteStrainElements[elementIndex];
            const auto material = assembly->getFiniteStrainMaterial(element->getMaterialId());
            if (!material) {
                continue;
            }

            finiteStrainCellSummaries.push_back(
                summarizeFiniteStrainElementResponse(
                    finiteStrainResponses[elementIndex],
                    *material,
                    element->getId(),
                    element->getMaterialId()));

            const auto& summary = finiteStrainCellSummaries.back();
            for (int nodeId : element->getNodeIds()) {
                auto nodeIt = std::find_if(nodes.begin(),
                    nodes.end(),
                    [nodeId](const std::shared_ptr<Node>& node) {
                        return node && node->getId() == nodeId;
                    });
                if (nodeIt == nodes.end()) {
                    continue;
                }
                const size_t nodeIndex =
                    static_cast<size_t>(std::distance(nodes.begin(), nodeIt));
                nodalStresses[nodeIndex] += summary.meanCauchyStress;
                nodalStrains[nodeIndex] += summary.meanGreenLagrangeStrain;
                nodalSigmaZZ[nodeIndex] += summary.meanSigmaZZ;
                nodalJacobianDeterminant[nodeIndex] += summary.meanJ;
                nodalStrainEnergyDensity[nodeIndex] += summary.meanStrainEnergyDensity;
                nodalContributionCount[nodeIndex] += 1;
            }
        }

        if (finiteStrainCellSummaries.size() < finiteStrainElements.size()) {
            const size_t previousSize = finiteStrainCellSummaries.size();
            finiteStrainCellSummaries.resize(finiteStrainElements.size());
            for (size_t elementIndex = previousSize;
                elementIndex < finiteStrainElements.size();
                ++elementIndex) {
                finiteStrainCellSummaries[elementIndex].elementId =
                    finiteStrainElements[elementIndex]->getId();
                finiteStrainCellSummaries[elementIndex].materialId =
                    finiteStrainElements[elementIndex]->getMaterialId();
            }
        }

        for (size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) {
            if (nodalContributionCount[nodeIndex] <= 0) {
                nodalJacobianDeterminant[nodeIndex] = 1.0;
                continue;
            }
            const double count = static_cast<double>(nodalContributionCount[nodeIndex]);
            nodalStresses[nodeIndex] /= count;
            nodalStrains[nodeIndex] /= count;
            nodalSigmaZZ[nodeIndex] /= count;
            nodalJacobianDeterminant[nodeIndex] /= count;
            nodalStrainEnergyDensity[nodeIndex] /= count;
        }
    }
    else {
        nodalStresses = model.getNodalStresses();
        nodalStrains = model.getNodalStrains();
        for (size_t nodeIndex = 0; nodeIndex < nodalStresses.size() &&
            nodeIndex < nodalSigmaZZ.size(); ++nodeIndex) {
            nodalSigmaZZ[nodeIndex] =
                computePlaneStrainSigmaZZ(nodalStresses[nodeIndex], representativePoissonsRatio);
        }
    }

    std::unordered_map<int, int> nodeIdToPointIndex;
    nodeIdToPointIndex.reserve(nodes.size());
    for (size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) {
        nodeIdToPointIndex[nodes[nodeIndex]->getId()] = static_cast<int>(nodeIndex);
    }

    std::vector<int> connectivity;
    std::vector<int> offsets;
    std::vector<int> cellTypes;
    const size_t cellCount = exportFiniteStrain
        ? finiteStrainElements.size()
        : linearElements.size();
    connectivity.reserve(cellCount * 4);
    offsets.reserve(cellCount);
    cellTypes.reserve(cellCount);

    std::unordered_map<int, int> elementIdToCellIndex;
    elementIdToCellIndex.reserve(cellCount);

    int offset = 0;
    if (exportFiniteStrain) {
        for (size_t elementIndex = 0; elementIndex < finiteStrainElements.size(); ++elementIndex) {
            const auto& element = finiteStrainElements[elementIndex];
            elementIdToCellIndex[element->getId()] = static_cast<int>(elementIndex);

            for (int nodeId : element->getNodeIds()) {
                auto nodeIt = nodeIdToPointIndex.find(nodeId);
                if (nodeIt == nodeIdToPointIndex.end()) {
                    throw std::runtime_error(
                        "Finite-strain element references a node that is not present in the assembly");
                }
                connectivity.push_back(nodeIt->second);
                ++offset;
            }

            offsets.push_back(offset);
            cellTypes.push_back(vtkCellTypeForNodeCount(
                static_cast<int>(element->getNodeIds().size())));
        }
    }
    else {
        for (size_t elementIndex = 0; elementIndex < linearElements.size(); ++elementIndex) {
            const auto& element = linearElements[elementIndex];
            elementIdToCellIndex[element->getId()] = static_cast<int>(elementIndex);

            for (int nodeId : element->getNodeIds()) {
                auto nodeIt = nodeIdToPointIndex.find(nodeId);
                if (nodeIt == nodeIdToPointIndex.end()) {
                    throw std::runtime_error(
                        "Element references a node that is not present in the assembly");
                }
                connectivity.push_back(nodeIt->second);
                ++offset;
            }

            offsets.push_back(offset);
            cellTypes.push_back(vtkCellTypeForNodeCount(element->getNodeCount()));
        }
    }

    std::vector<int> elementIds(cellCount, 0);
    std::vector<int> materialIds(cellCount, 0);
    std::vector<int> candidateContactFacet(cellCount, 0);
    std::vector<int> candidateContactSurfaceIndex(cellCount, -1);
    std::vector<int> activeContactFacet(cellCount, 0);
    std::vector<int> activeContactSurfaceIndex(cellCount, -1);

    if (exportFiniteStrain) {
        for (size_t elementIndex = 0; elementIndex < finiteStrainElements.size(); ++elementIndex) {
            elementIds[elementIndex] = finiteStrainElements[elementIndex]->getId();
            materialIds[elementIndex] = finiteStrainElements[elementIndex]->getMaterialId();
        }
    }
    else {
        for (size_t elementIndex = 0; elementIndex < linearElements.size(); ++elementIndex) {
            elementIds[elementIndex] = linearElements[elementIndex]->getId();
            materialIds[elementIndex] = linearElements[elementIndex]->getMaterialId();
        }
    }

    ContactState contactState;
    std::vector<ContactFacetResult> contactFacetResults;
    if (model.hasContactSolver()) {
        const auto contactFacets = model.getContactFacets();
        for (const auto& facet : contactFacets) {
            auto elementIt = elementIdToCellIndex.find(facet.elementId);
            if (elementIt == elementIdToCellIndex.end()) {
                continue;
            }
            const int cellIndex = elementIt->second;
            candidateContactFacet[cellIndex] = 1;
            candidateContactSurfaceIndex[cellIndex] = facet.surfaceIndex;
        }

        contactState = model.evaluateCurrentContactState();
        contactFacetResults = contactState.facetResults;
        for (int activeFacetId : contactState.activeFacetIds) {
            if (activeFacetId < 0 || activeFacetId >= static_cast<int>(contactFacets.size())) {
                continue;
            }
            const auto& facet = contactFacets[activeFacetId];
            auto elementIt = elementIdToCellIndex.find(facet.elementId);
            if (elementIt == elementIdToCellIndex.end()) {
                continue;
            }
            const int cellIndex = elementIt->second;
            activeContactFacet[cellIndex] = 1;
            activeContactSurfaceIndex[cellIndex] = facet.surfaceIndex;
        }
    }

    std::ofstream vtuStream(vtuPath);
    if (!vtuStream.is_open()) {
        throw std::runtime_error("Failed to open VTU output file: " + vtuPath.string());
    }

    vtuStream << "<?xml version=\"1.0\"?>\n";
    vtuStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtuStream << "  <UnstructuredGrid>\n";
    vtuStream << "    <Piece NumberOfPoints=\"" << nodes.size()
              << "\" NumberOfCells=\"" << cellCount << "\">\n";
    vtuStream << "      <PointData Vectors=\"displacement\">\n";

    writeVector3DataArray(vtuStream, "displacement", static_cast<int>(nodes.size()),
        [&](int index) {
            const Eigen::Vector2d displacement =
                index < static_cast<int>(nodalDisplacements.size())
                ? nodalDisplacements[index]
                : Eigen::Vector2d::Zero();
            return Eigen::Vector3d(displacement.x(), displacement.y(), 0.0);
        });
    writeVector3DataArray(vtuStream, "deformed_position", static_cast<int>(nodes.size()),
        [&](int index) {
            const Eigen::Vector2d coordinates = nodes[static_cast<size_t>(index)]->getCoordinates();
            const Eigen::Vector2d displacement =
                index < static_cast<int>(nodalDisplacements.size())
                ? nodalDisplacements[index]
                : Eigen::Vector2d::Zero();
            const Eigen::Vector2d deformedPosition = coordinates + displacement;
            return Eigen::Vector3d(deformedPosition.x(), deformedPosition.y(), 0.0);
        });
    writeScalarDataArray(vtuStream, "displacement_magnitude", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalDisplacements.size())
                ? vectorMagnitude(nodalDisplacements[index])
                : 0.0);
        });

    writeVector3DataArray(vtuStream, "stress_2d", static_cast<int>(nodes.size()),
        [&](int index) {
            const Eigen::Vector3d stress =
                index < static_cast<int>(nodalStresses.size())
                ? nodalStresses[index]
                : Eigen::Vector3d::Zero();
            return stress;
        });
    writeScalarDataArray(vtuStream, "sigma_xx", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStresses.size())
                ? nodalStresses[index].x()
                : 0.0);
        });
    writeScalarDataArray(vtuStream, "sigma_yy", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStresses.size())
                ? nodalStresses[index].y()
                : 0.0);
        });
    writeScalarDataArray(vtuStream, "tau_xy", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStresses.size())
                ? nodalStresses[index].z()
                : 0.0);
        });
    writeScalarDataArray(vtuStream, "sigma_zz", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalSigmaZZ.size())
                ? nodalSigmaZZ[index]
                : 0.0);
        });
    writeScalarDataArray(vtuStream, "von_mises_stress", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStresses.size())
                ? computeVonMisesStress(
                    nodalStresses[index],
                    index < static_cast<int>(nodalSigmaZZ.size()) ? nodalSigmaZZ[index] : 0.0)
                : 0.0);
        });
    if (exportFiniteStrain) {
        writeVector3DataArray(vtuStream, "cauchy_stress_2d", static_cast<int>(nodes.size()),
            [&](int index) {
                return index < static_cast<int>(nodalStresses.size())
                    ? nodalStresses[index]
                    : Eigen::Vector3d::Zero();
            });
        writeVector3DataArray(vtuStream, "green_lagrange_strain_2d", static_cast<int>(nodes.size()),
            [&](int index) {
                return index < static_cast<int>(nodalStrains.size())
                    ? nodalStrains[index]
                    : Eigen::Vector3d::Zero();
            });
        writeScalarDataArray(vtuStream, "jacobian_determinant", static_cast<int>(nodes.size()),
            [&](int index) {
                return formatNumber(index < static_cast<int>(nodalJacobianDeterminant.size())
                    ? nodalJacobianDeterminant[index]
                    : 1.0);
            });
        writeScalarDataArray(vtuStream, "strain_energy_density", static_cast<int>(nodes.size()),
            [&](int index) {
                return formatNumber(index < static_cast<int>(nodalStrainEnergyDensity.size())
                    ? nodalStrainEnergyDensity[index]
                    : 0.0);
            });
    }

    writeVector3DataArray(vtuStream, "strain_2d", static_cast<int>(nodes.size()),
        [&](int index) {
            const Eigen::Vector3d strain =
                index < static_cast<int>(nodalStrains.size())
                ? nodalStrains[index]
                : Eigen::Vector3d::Zero();
            return strain;
        });
    writeScalarDataArray(vtuStream, "strain_xx", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStrains.size())
                ? nodalStrains[index].x()
                : 0.0);
        });
    writeScalarDataArray(vtuStream, "strain_yy", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStrains.size())
                ? nodalStrains[index].y()
                : 0.0);
        });
    writeScalarDataArray(vtuStream, "gamma_xy", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalStrains.size())
                ? nodalStrains[index].z()
                : 0.0);
        });

    writeVector3DataArray(vtuStream, "reaction_force", static_cast<int>(nodes.size()),
        [&](int index) {
            const Eigen::Vector2d reaction =
                index < static_cast<int>(nodalReactionForces.size())
                ? nodalReactionForces[index]
                : Eigen::Vector2d::Zero();
            return Eigen::Vector3d(reaction.x(), reaction.y(), 0.0);
        });
    writeScalarDataArray(vtuStream, "reaction_force_magnitude", static_cast<int>(nodes.size()),
        [&](int index) {
            return formatNumber(index < static_cast<int>(nodalReactionForces.size())
                ? vectorMagnitude(nodalReactionForces[index])
                : 0.0);
        });

    if (model.hasContactSolver()) {
        writeVector3DataArray(vtuStream, "contact_force", static_cast<int>(nodes.size()),
            [&](int index) {
                const Eigen::Vector2d contactForce =
                    index < static_cast<int>(nodalContactForces.size())
                    ? nodalContactForces[index]
                    : Eigen::Vector2d::Zero();
                return Eigen::Vector3d(contactForce.x(), contactForce.y(), 0.0);
            });
        writeScalarDataArray(vtuStream, "contact_force_magnitude", static_cast<int>(nodes.size()),
            [&](int index) {
                return formatNumber(index < static_cast<int>(nodalContactForces.size())
                    ? vectorMagnitude(nodalContactForces[index])
                    : 0.0);
            });
        writeScalarDataArray(vtuStream, "rigid_plane_signed_distance", static_cast<int>(nodes.size()),
            [&](int index) {
                return formatNumber(index < static_cast<int>(nodalSignedDistances.size())
                    ? nodalSignedDistances[index]
                    : 0.0);
            });
        writeScalarDataArray(vtuStream, "rigid_plane_penetration", static_cast<int>(nodes.size()),
            [&](int index) {
                return formatNumber(index < static_cast<int>(nodalPenetrations.size())
                    ? nodalPenetrations[index]
                    : 0.0);
            });
    }

    vtuStream << "      </PointData>\n";
    vtuStream << "      <CellData>\n";
    writeIntDataArray(vtuStream, "element_id", static_cast<int>(cellCount),
        [&](int index) { return formatInteger(elementIds[index]); });
    writeIntDataArray(vtuStream, "material_id", static_cast<int>(cellCount),
        [&](int index) { return formatInteger(materialIds[index]); });

    if (exportFiniteStrain) {
        writeScalarDataArray(vtuStream, "cell_jacobian_determinant", static_cast<int>(cellCount),
            [&](int index) { return formatNumber(finiteStrainCellSummaries[static_cast<size_t>(index)].meanJ); });
        writeScalarDataArray(vtuStream, "cell_jacobian_determinant_min", static_cast<int>(cellCount),
            [&](int index) { return formatNumber(finiteStrainCellSummaries[static_cast<size_t>(index)].minJ); });
        writeScalarDataArray(vtuStream, "cell_jacobian_determinant_max", static_cast<int>(cellCount),
            [&](int index) { return formatNumber(finiteStrainCellSummaries[static_cast<size_t>(index)].maxJ); });
        writeScalarDataArray(vtuStream, "cell_strain_energy_density", static_cast<int>(cellCount),
            [&](int index) { return formatNumber(finiteStrainCellSummaries[static_cast<size_t>(index)].meanStrainEnergyDensity); });
        writeVector3DataArray(vtuStream, "cell_cauchy_stress_2d", static_cast<int>(cellCount),
            [&](int index) { return finiteStrainCellSummaries[static_cast<size_t>(index)].meanCauchyStress; });
        writeScalarDataArray(vtuStream, "cell_sigma_zz", static_cast<int>(cellCount),
            [&](int index) { return formatNumber(finiteStrainCellSummaries[static_cast<size_t>(index)].meanSigmaZZ); });
        writeScalarDataArray(vtuStream, "cell_von_mises_stress", static_cast<int>(cellCount),
            [&](int index) {
                const auto& summary = finiteStrainCellSummaries[static_cast<size_t>(index)];
                return formatNumber(computeVonMisesStress(summary.meanCauchyStress, summary.meanSigmaZZ));
            });
        writeVector3DataArray(vtuStream, "cell_second_piola_stress_2d", static_cast<int>(cellCount),
            [&](int index) { return finiteStrainCellSummaries[static_cast<size_t>(index)].meanSecondPiolaStress; });
        writeVector3DataArray(vtuStream, "cell_green_lagrange_strain_2d", static_cast<int>(cellCount),
            [&](int index) { return finiteStrainCellSummaries[static_cast<size_t>(index)].meanGreenLagrangeStrain; });
    }

    if (model.hasContactSolver()) {
        writeIntDataArray(vtuStream, "candidate_contact_facet", static_cast<int>(cellCount),
            [&](int index) { return formatInteger(candidateContactFacet[index]); });
        writeIntDataArray(vtuStream, "candidate_contact_surface_index", static_cast<int>(cellCount),
            [&](int index) { return formatInteger(candidateContactSurfaceIndex[index]); });
        writeIntDataArray(vtuStream, "active_contact_facet", static_cast<int>(cellCount),
            [&](int index) { return formatInteger(activeContactFacet[index]); });
        writeIntDataArray(vtuStream, "active_contact_surface_index", static_cast<int>(cellCount),
            [&](int index) { return formatInteger(activeContactSurfaceIndex[index]); });
    }

    vtuStream << "      </CellData>\n";
    vtuStream << "      <Points>\n";
    vtuStream << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    vtuStream << "          ";
    for (size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) {
        if (nodeIndex > 0) {
            vtuStream << ' ';
        }
        const Eigen::Vector2d coordinates = nodes[nodeIndex]->getCoordinates();
        vtuStream << formatNumber(coordinates.x()) << ' '
                  << formatNumber(coordinates.y()) << ' '
                  << formatNumber(0.0);
    }
    vtuStream << "\n";
    vtuStream << "        </DataArray>\n";
    vtuStream << "      </Points>\n";
    vtuStream << "      <Cells>\n";
    vtuStream << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    vtuStream << "          ";
    for (size_t i = 0; i < connectivity.size(); ++i) {
        if (i > 0) {
            vtuStream << ' ';
        }
        vtuStream << connectivity[i];
    }
    vtuStream << "\n";
    vtuStream << "        </DataArray>\n";
    vtuStream << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    vtuStream << "          ";
    for (size_t i = 0; i < offsets.size(); ++i) {
        if (i > 0) {
            vtuStream << ' ';
        }
        vtuStream << offsets[i];
    }
    vtuStream << "\n";
    vtuStream << "        </DataArray>\n";
    vtuStream << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    vtuStream << "          ";
    for (size_t i = 0; i < cellTypes.size(); ++i) {
        if (i > 0) {
            vtuStream << ' ';
        }
        vtuStream << cellTypes[i];
    }
    vtuStream << "\n";
    vtuStream << "        </DataArray>\n";
    vtuStream << "      </Cells>\n";
    vtuStream << "    </Piece>\n";
    vtuStream << "  </UnstructuredGrid>\n";
    vtuStream << "</VTKFile>\n";
    vtuStream.close();

    double maxDisplacementMagnitude = 0.0;
    for (const auto& displacement : nodalDisplacements) {
        maxDisplacementMagnitude = std::max(maxDisplacementMagnitude, displacement.norm());
    }

    double maxReactionForceMagnitude = 0.0;
    for (const auto& reactionForce : nodalReactionForces) {
        maxReactionForceMagnitude = std::max(maxReactionForceMagnitude, reactionForce.norm());
    }

    double maxContactForceMagnitude = 0.0;
    for (const auto& contactForce : nodalContactForces) {
        maxContactForceMagnitude = std::max(maxContactForceMagnitude, contactForce.norm());
    }

    double maxNodalPenetration = 0.0;
    double minSignedDistance = 0.0;
    if (!nodalPenetrations.empty()) {
        maxNodalPenetration = *std::max_element(nodalPenetrations.begin(), nodalPenetrations.end());
    }
    if (!nodalSignedDistances.empty()) {
        minSignedDistance = *std::min_element(nodalSignedDistances.begin(), nodalSignedDistances.end());
    }

    double totalNormalForce = 0.0;
    double activeContactArea = 0.0;
    double contactPatchLength = 0.0;
    double maxFacetAveragePressure = 0.0;
    double centerOfPressureX = 0.0;
    double centerOfPressureY = 0.0;
    int activeFacetCount = 0;
    for (const auto& facetResult : contactFacetResults) {
        totalNormalForce += facetResult.integratedNormalForce;
        activeContactArea += facetResult.activeArea;
        contactPatchLength += facetResult.activeLength;
        maxFacetAveragePressure = std::max(maxFacetAveragePressure, facetResult.averagePressure);
        if (facetResult.active) {
            ++activeFacetCount;
        }
        centerOfPressureX +=
            facetResult.integratedNormalForce * facetResult.referenceMidpoint.x();
        centerOfPressureY +=
            facetResult.integratedNormalForce * facetResult.referenceMidpoint.y();
    }
    if (totalNormalForce > 0.0) {
        centerOfPressureX /= totalNormalForce;
        centerOfPressureY /= totalNormalForce;
    }
    const double meanActivePressure =
        activeContactArea > 0.0 ? totalNormalForce / activeContactArea : 0.0;

    const auto& dofMapping = assembly->getDofMapping();
    const long long nodeCount = static_cast<long long>(nodes.size());
    const long long elementCount = static_cast<long long>(cellCount);
    const long long totalDof = static_cast<long long>(assembly->getTotalDofCount());
    const long long freeDof = static_cast<long long>(dofMapping.reducedToFull.size());
    const long long prescribedDof = static_cast<long long>(dofMapping.prescribedDofs.size());
    const long long constrainedDof = totalDof - freeDof;
    const long long candidateFacetCount = static_cast<long long>(contactFacetResults.size());
    const long long activeFacetCountLongLong = static_cast<long long>(activeFacetCount);

    double minCellJacobianDeterminant = 1.0;
    double maxCellJacobianDeterminant = 1.0;
    double meanCellJacobianDeterminant = 1.0;
    double minCellStrainEnergyDensity = 0.0;
    double maxCellStrainEnergyDensity = 0.0;
    double meanCellStrainEnergyDensity = 0.0;
    if (!finiteStrainCellSummaries.empty()) {
        minCellJacobianDeterminant =
            finiteStrainCellSummaries.front().minJ;
        maxCellJacobianDeterminant =
            finiteStrainCellSummaries.front().maxJ;
        meanCellJacobianDeterminant = 0.0;
        minCellStrainEnergyDensity =
            finiteStrainCellSummaries.front().meanStrainEnergyDensity;
        maxCellStrainEnergyDensity =
            finiteStrainCellSummaries.front().meanStrainEnergyDensity;
        meanCellStrainEnergyDensity = 0.0;

        for (const auto& summary : finiteStrainCellSummaries) {
            minCellJacobianDeterminant =
                std::min(minCellJacobianDeterminant, summary.minJ);
            maxCellJacobianDeterminant =
                std::max(maxCellJacobianDeterminant, summary.maxJ);
            meanCellJacobianDeterminant += summary.meanJ;

            minCellStrainEnergyDensity =
                std::min(minCellStrainEnergyDensity, summary.meanStrainEnergyDensity);
            maxCellStrainEnergyDensity =
                std::max(maxCellStrainEnergyDensity, summary.meanStrainEnergyDensity);
            meanCellStrainEnergyDensity += summary.meanStrainEnergyDensity;
        }

        meanCellJacobianDeterminant /=
            static_cast<double>(finiteStrainCellSummaries.size());
        meanCellStrainEnergyDensity /=
            static_cast<double>(finiteStrainCellSummaries.size());
    }

    if (model.hasContactSolver()) {
        std::ofstream contactFacetStream(contactFacetCsvPath);
        if (!contactFacetStream.is_open()) {
            throw std::runtime_error(
                "Failed to open contact facet CSV output file: " + contactFacetCsvPath.string());
        }

        contactFacetStream << "facet_id,element_id,surface_index,active,active_gauss_points,"
                           << "reference_midpoint_x,reference_midpoint_y,"
                           << "deformed_midpoint_x,deformed_midpoint_y,"
                           << "normal_x,normal_y,thickness,facet_length,active_length,"
                           << "integrated_area,active_area,average_gap,average_penetration,"
                           << "maximum_penetration,integrated_normal_force,average_pressure\n";

        for (const auto& facetResult : contactFacetResults) {
            contactFacetStream << facetResult.facetId << ','
                               << facetResult.elementId << ','
                               << facetResult.surfaceIndex << ','
                               << (facetResult.active ? 1 : 0) << ','
                               << facetResult.activeGaussPointCount << ','
                               << formatNumber(facetResult.referenceMidpoint.x()) << ','
                               << formatNumber(facetResult.referenceMidpoint.y()) << ','
                               << formatNumber(facetResult.deformedMidpoint.x()) << ','
                               << formatNumber(facetResult.deformedMidpoint.y()) << ','
                               << formatNumber(facetResult.normal.x()) << ','
                               << formatNumber(facetResult.normal.y()) << ','
                               << formatNumber(facetResult.thickness) << ','
                               << formatNumber(facetResult.facetLength) << ','
                               << formatNumber(facetResult.activeLength) << ','
                               << formatNumber(facetResult.integratedArea) << ','
                               << formatNumber(facetResult.activeArea) << ','
                               << formatNumber(facetResult.averageGap) << ','
                               << formatNumber(facetResult.averagePenetration) << ','
                               << formatNumber(facetResult.maximumPenetration) << ','
                               << formatNumber(facetResult.integratedNormalForce) << ','
                               << formatNumber(facetResult.averagePressure) << '\n';
        }
    }

    std::ofstream metricsStream(metricsJsonPath);
    if (!metricsStream.is_open()) {
        throw std::runtime_error("Failed to open metrics JSON output file: " + metricsJsonPath.string());
    }

    metricsStream << "{\n";
    metricsStream << "  \"base_name\": \"" << jsonEscape(baseName) << "\",\n";
    metricsStream << "  \"files\": {\n";
    metricsStream << "    \"vtu\": \"" << jsonEscape(vtuPath.filename().string()) << "\",\n";
    metricsStream << "    \"metrics\": \"" << jsonEscape(metricsJsonPath.filename().string()) << "\"";
    if (model.hasContactSolver()) {
        metricsStream << ",\n";
        metricsStream << "    \"contact_facets\": \""
                      << jsonEscape(contactFacetCsvPath.filename().string()) << "\"\n";
    }
    else {
        metricsStream << "\n";
    }
    metricsStream << "  },\n";
    metricsStream << "  \"units\": {\n";
    metricsStream << "    \"length\": \"mm\",\n";
    metricsStream << "    \"stress\": \"MPa\",\n";
    metricsStream << "    \"pressure\": \"MPa\",\n";
    metricsStream << "    \"elastic_modulus\": \"MPa\",\n";
    metricsStream << "    \"strain\": \"1\",\n";
    metricsStream << "    \"force_internal\": \"N\",\n";
    metricsStream << "    \"force_display\": \"kN\"\n";
    metricsStream << "  },\n";
    metricsStream << "  \"counts\": {\n";
    metricsStream << "    \"nodes\": " << nodeCount << ",\n";
    metricsStream << "    \"elements\": " << elementCount << ",\n";
    metricsStream << "    \"total_dofs\": " << totalDof << ",\n";
    metricsStream << "    \"free_dofs\": " << freeDof << ",\n";
    metricsStream << "    \"constrained_dofs\": " << constrainedDof << ",\n";
    metricsStream << "    \"prescribed_dofs\": " << prescribedDof << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"timings\": {\n";
    metricsStream << "    \"assembly_time_seconds\": " << formatNumber(performanceMetrics.assemblyTimeSeconds) << ",\n";
    metricsStream << "    \"solve_time_seconds\": " << formatNumber(performanceMetrics.solveTimeSeconds) << ",\n";
    metricsStream << "    \"total_time_seconds\": " << formatNumber(performanceMetrics.totalTimeSeconds) << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"formulation\": {\n";
    metricsStream << "    \"finite_strain\": " << (exportFiniteStrain ? "true" : "false") << ",\n";
    metricsStream << "    \"point_stress_measure\": \""
                  << (exportFiniteStrain ? "cauchy" : "small_strain_constitutive")
                  << "\",\n";
    metricsStream << "    \"point_strain_measure\": \""
                  << (exportFiniteStrain ? "green_lagrange" : "small_strain")
                  << "\"\n";
    metricsStream << "  },\n";
    metricsStream << "  \"iterations\": {\n";
    metricsStream << "    \"linear_solves\": " << performanceMetrics.linearSolveCount << ",\n";
    metricsStream << "    \"direct_linear_solves\": "
                  << performanceMetrics.directLinearSolveCount << ",\n";
    metricsStream << "    \"nonlinear_iterations\": " << performanceMetrics.nonlinearIterations << ",\n";
    metricsStream << "    \"linear_iterations\": " << performanceMetrics.linearIterations << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"solver\": {\n";
    metricsStream << "    \"backend\": \"" << jsonEscape(performanceMetrics.linearSolverBackend) << "\",\n";
    metricsStream << "    \"converged\": " << (performanceMetrics.linearSolveConverged ? "true" : "false") << ",\n";
    metricsStream << "    \"used_direct_solver\": " << (performanceMetrics.usedDirectLinearSolver ? "true" : "false") << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"residuals\": {\n";
    metricsStream << "    \"linear_residual_estimate\": " << formatNumber(performanceMetrics.linearResidualEstimate) << ",\n";
    metricsStream << "    \"linear_relative_residual_norm\": " << formatNumber(performanceMetrics.linearResidualNorm) << ",\n";
    metricsStream << "    \"equilibrium_residual_norm\": " << formatNumber(performanceMetrics.equilibriumResidualNorm) << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"matrix\": {\n";
    metricsStream << "    \"nnz\": " << performanceMetrics.matrixNonZeros << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"finite_strain\": {\n";
    metricsStream << "    \"configured\": " << (exportFiniteStrain ? "true" : "false") << ",\n";
    metricsStream << "    \"strain_energy\": " << formatNumber(performanceMetrics.strainEnergy) << ",\n";
    metricsStream << "    \"has_near_incompressible_material\": "
                  << (performanceMetrics.hasNearIncompressibleFiniteStrainMaterial ? "true" : "false") << ",\n";
    metricsStream << "    \"max_effective_poissons_ratio\": "
                  << formatNumber(performanceMetrics.maxFiniteStrainEffectivePoissonsRatio) << ",\n";
    metricsStream << "    \"max_bulk_to_shear_ratio\": "
                  << formatNumber(performanceMetrics.maxFiniteStrainBulkToShearRatio) << ",\n";
    metricsStream << "    \"min_cell_jacobian_determinant\": "
                  << formatNumber(minCellJacobianDeterminant) << ",\n";
    metricsStream << "    \"mean_cell_jacobian_determinant\": "
                  << formatNumber(meanCellJacobianDeterminant) << ",\n";
    metricsStream << "    \"max_cell_jacobian_determinant\": "
                  << formatNumber(maxCellJacobianDeterminant) << ",\n";
    metricsStream << "    \"min_cell_strain_energy_density\": "
                  << formatNumber(minCellStrainEnergyDensity) << ",\n";
    metricsStream << "    \"mean_cell_strain_energy_density\": "
                  << formatNumber(meanCellStrainEnergyDensity) << ",\n";
    metricsStream << "    \"max_cell_strain_energy_density\": "
                  << formatNumber(maxCellStrainEnergyDensity) << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"contact\": {\n";
    metricsStream << "    \"configured\": " << (model.hasContactSolver() ? "true" : "false") << ",\n";
    metricsStream << "    \"method\": \"" << jsonEscape(performanceMetrics.contactMethod) << "\",\n";
    metricsStream << "    \"candidate_contact_facets\": " << candidateFacetCount << ",\n";
    metricsStream << "    \"active_contact_facets\": " << activeFacetCountLongLong << ",\n";
    metricsStream << "    \"active_contact_gauss_points\": "
                  << performanceMetrics.activeContactGaussPoints << ",\n";
    metricsStream << "    \"contact_force_norm\": " << formatNumber(performanceMetrics.contactForceNorm) << ",\n";
    metricsStream << "    \"contact_state_update_norm\": "
                  << formatNumber(performanceMetrics.contactStateUpdateNorm) << ",\n";
    metricsStream << "    \"contact_state_relative_update_norm\": "
                  << formatNumber(performanceMetrics.contactStateRelativeUpdateNorm) << ",\n";
    metricsStream << "    \"contact_state_max_update_magnitude\": "
                  << formatNumber(performanceMetrics.contactStateMaxUpdateMagnitude) << ",\n";
    metricsStream << "    \"contact_state_relative_max_update\": "
                  << formatNumber(performanceMetrics.contactStateRelativeMaxUpdate) << ",\n";
    metricsStream << "    \"max_normal_multiplier\": "
                  << formatNumber(performanceMetrics.maxNormalContactMultiplier) << ",\n";
    metricsStream << "    \"mean_normal_multiplier\": "
                  << formatNumber(performanceMetrics.meanNormalContactMultiplier) << ",\n";
    metricsStream << "    \"max_penetration\": " << formatNumber(performanceMetrics.maxPenetration) << ",\n";
    metricsStream << "    \"max_nodal_penetration\": " << formatNumber(maxNodalPenetration) << ",\n";
    metricsStream << "    \"minimum_signed_distance\": " << formatNumber(minSignedDistance) << ",\n";
    metricsStream << "    \"total_normal_force\": " << formatNumber(totalNormalForce) << ",\n";
    metricsStream << "    \"active_contact_area\": " << formatNumber(activeContactArea) << ",\n";
    metricsStream << "    \"contact_patch_length\": " << formatNumber(contactPatchLength) << ",\n";
    metricsStream << "    \"max_facet_average_pressure\": " << formatNumber(maxFacetAveragePressure) << ",\n";
    metricsStream << "    \"mean_active_pressure\": " << formatNumber(meanActivePressure) << ",\n";
    metricsStream << "    \"center_of_pressure_x\": " << formatNumber(centerOfPressureX) << ",\n";
    metricsStream << "    \"center_of_pressure_y\": " << formatNumber(centerOfPressureY) << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"extrema\": {\n";
    metricsStream << "    \"max_displacement_magnitude\": " << formatNumber(maxDisplacementMagnitude) << ",\n";
    metricsStream << "    \"max_reaction_force_magnitude\": " << formatNumber(maxReactionForceMagnitude) << ",\n";
    metricsStream << "    \"max_contact_force_magnitude\": " << formatNumber(maxContactForceMagnitude) << "\n";
    metricsStream << "  },\n";
    metricsStream << "  \"extra\": {\n";

    bool wroteExtraMetric = false;
    metricsStream << "    \"plane_strain_poissons_ratio_for_von_mises\": "
                  << formatNumber(representativePoissonsRatio);
    wroteExtraMetric = true;
    for (const auto& [name, value] : options.extraStringMetrics) {
        if (wroteExtraMetric) {
            metricsStream << ",\n";
        }
        metricsStream << "    \"" << jsonEscape(name) << "\": \"" << jsonEscape(value) << "\"";
        wroteExtraMetric = true;
    }
    for (const auto& [name, value] : options.extraNumericMetrics) {
        if (wroteExtraMetric) {
            metricsStream << ",\n";
        }
        metricsStream << "    \"" << jsonEscape(name) << "\": " << formatNumber(value);
        wroteExtraMetric = true;
    }
    if (wroteExtraMetric) {
        metricsStream << '\n';
    }
    metricsStream << "  }\n";
    metricsStream << "}\n";
    metricsStream.close();

    return {vtuPath, metricsJsonPath, model.hasContactSolver() ? contactFacetCsvPath : std::filesystem::path()};
}
