#pragma once

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>

#ifndef FEM_REPO_ROOT
#error FEM_REPO_ROOT must be defined for application targets
#endif

namespace AppRuntimeSupport {

inline std::filesystem::path repoRoot() {
    return std::filesystem::path(FEM_REPO_ROOT);
}

inline std::filesystem::path resultsRoot() {
    return repoRoot() / "results";
}

inline std::filesystem::path caseOutputDirectory(std::string_view caseName) {
    return resultsRoot() / std::filesystem::path(caseName);
}

inline std::filesystem::path postprocessScriptPath() {
    return repoRoot() / "scripts" / "postprocess_results.py";
}

inline std::string quoteCommandArgument(const std::filesystem::path& path) {
    return "\"" + path.string() + "\"";
}

inline void runPostprocessorIfEnabled(const std::filesystem::path& targetPath) {
    if (const char* skipAutoPostprocess = std::getenv("FEM_SKIP_AUTO_POSTPROCESS");
        skipAutoPostprocess != nullptr && std::string_view(skipAutoPostprocess) == "1") {
        return;
    }

    const std::filesystem::path scriptPath = postprocessScriptPath();
    if (!std::filesystem::exists(scriptPath)) {
        std::cerr << "Skipping automatic postprocessing: script not found at "
                  << scriptPath << std::endl;
        return;
    }

    const std::string command =
        "python " + quoteCommandArgument(scriptPath) + " " + quoteCommandArgument(targetPath);
    const int exitCode = std::system(command.c_str());
    if (exitCode != 0) {
        std::cerr << "Automatic postprocessing returned exit code " << exitCode
                  << " for " << targetPath << std::endl;
    }
}

} // namespace AppRuntimeSupport
