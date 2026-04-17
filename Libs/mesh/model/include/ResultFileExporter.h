#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

class FEModel;

struct ResultFileExportOptions {
    std::filesystem::path outputDirectory;
    std::string baseName = "solution";
    std::vector<std::pair<std::string, double>> extraNumericMetrics;
    std::vector<std::pair<std::string, std::string>> extraStringMetrics;
};

struct ResultFileExportArtifacts {
    std::filesystem::path vtuPath;
    std::filesystem::path metricsJsonPath;
    std::filesystem::path contactFacetCsvPath;
};

class ResultFileExporter {
public:
    static ResultFileExportArtifacts exportSolution(
        const FEModel& model,
        const ResultFileExportOptions& options);
};
