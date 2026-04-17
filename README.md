# FEM Tire Contact Solver

Этот репозиторий содержит 2D Q4 FEM-решатель с нелинейным penalty-contact для постановки
`шина -> жесткая плоскость` без трения.

Основной инженерный workflow такой:

1. C++ собирает задачу и решает ее.
2. C++ пишет результаты в файловый pipeline:
   - `solution.vtu`
   - `metrics.json`
   - для contact-case дополнительно `contact_facets.csv`
   - при серийных расчетах дополнительно `study_summary.csv`
3. Python используется только для чтения результатов, графиков и автоматических проверок.

## Быстрый старт

После конфигурации CMake доступны основные исполняемые файлы:

- `FEMBasicLinearReference`
  Базовая линейная verification-задача без контакта.
- `FEM`
  Канонический одиночный расчет контакта шины с жесткой плоскостью.
- `FEMContactFarPlaneScenario`
  Contact-enabled сценарий, в котором плоскость вынесена далеко и контакта быть не должно.
- `FEMBlockOnRigidPlane`
  Маленький sanity-case для контактного контура.
- `FEMRingContactStudy`
  Серия расчетов по штрафному параметру и типу сетки.
- `FEMMainScaleContactStudy`
  Более тяжелая серия contact-расчетов на параметрах из старого `main.cpp`,
  где одновременно отслеживаются влияние `penalty` и влияние contact-focused сетки.
- `FEMContactNoContactRegression`
  Проверка, что `solveContact()` не портит решение при отсутствии контакта.
- `FEMContactBlockRegression`
  Регрессия для блока на жесткой плоскости.

### Сборка

Пример для Windows + Visual Studio:

```powershell
cmake -S . -B build
cmake --build build --config Release --target FEM FEMBlockOnRigidPlane FEMRingContactStudy FEMMainScaleContactStudy
```

Если `cmake` не лежит в `PATH`, можно использовать `cmake.exe` из установки Visual Studio.

### Запуск одиночного tire-contact case

```powershell
.\build\bin\Release\FEM.exe
```

Этот target использует файл [src/src/tire_contact_single_case.cpp](/c:/Users/Admin/Fem_diplom/FEM/src/src/tire_contact_single_case.cpp:1)
как канонический пример.

Примеры и study-сценарии теперь всегда пишут результаты в `results/...` в корне репозитория,
независимо от текущей рабочей директории.
После успешного расчета Python-постпроцессор запускается автоматически, поэтому рядом с
`solution.vtu` и `metrics.json` сразу появляются PNG-графики.
Для `FEM.exe` результаты появятся в:

- [results/tire_contact_single_case/solution.vtu](/c:/Users/Admin/Fem_diplom/FEM/results/tire_contact_single_case/solution.vtu)
- [results/tire_contact_single_case/metrics.json](/c:/Users/Admin/Fem_diplom/FEM/results/tire_contact_single_case/metrics.json)

### Запуск block sanity-case

```powershell
.\build\bin\Release\FEMBlockOnRigidPlane.exe
```

Результаты:

- [results/block_on_rigid_plane/solution.vtu](/c:/Users/Admin/Fem_diplom/FEM/results/block_on_rigid_plane/solution.vtu)
- [results/block_on_rigid_plane/metrics.json](/c:/Users/Admin/Fem_diplom/FEM/results/block_on_rigid_plane/metrics.json)

### Запуск серии tire-contact cases

```powershell
.\build\bin\Release\FEMRingContactStudy.exe
```

### Запуск main-scale contact study

```powershell
.\build\bin\Release\FEMMainScaleContactStudy.exe
```

Этот target использует геометрию и материал из старого [src/src/main.cpp](/c:/Users/Admin/Fem_diplom/FEM/src/src/main.cpp:1),
но решает уже контактную задачу с жесткой плоскостью.
По умолчанию это тяжелый прогон по нескольким уровням сетки и нескольким значениям `penalty`.

Результаты:

- [results/main_scale_contact_study/study_summary.csv](/c:/Users/Admin/Fem_diplom/FEM/results/main_scale_contact_study/study_summary.csv)
- [results/main_scale_contact_study/summary_metrics.png](/c:/Users/Admin/Fem_diplom/FEM/results/main_scale_contact_study/summary_metrics.png)
- [results/main_scale_contact_study/summary_contact_metrics.png](/c:/Users/Admin/Fem_diplom/FEM/results/main_scale_contact_study/summary_contact_metrics.png)

Для быстрой sanity-проверки можно ограничить прогон одной комбинацией `mesh x penalty`:

```powershell
$env:FEM_MAIN_SCALE_CONTACT_QUICK='1'
.\build\bin\Release\FEMMainScaleContactStudy.exe
```

Результаты:

- [results/ring_contact_study/study_summary.csv](/c:/Users/Admin/Fem_diplom/FEM/results/ring_contact_study/study_summary.csv)
- по одной папке на кейс, каждая с `solution.vtu` и `metrics.json`

## Что лежит в файловом pipeline

### `solution.vtu`

Основной файловый формат обмена между C++-решателем и Python-постпроцессором.

Экспортируются:

- `Points`
- `Cells / connectivity / offsets / types`
- nodal `displacement`
- nodal `stress_2d`, `strain_2d`
- nodal `reaction_force`
- для contact-case:
  - `contact_force`
  - `rigid_plane_signed_distance`
  - `rigid_plane_penetration`
  - cell flags `candidate_contact_facet`
  - cell flags `active_contact_facet`

Экспортер реализован в:

- [ResultFileExporter.h](/c:/Users/Admin/Fem_diplom/FEM/Libs/mesh/model/include/ResultFileExporter.h:10)
- [ResultFileExporter.cpp](/c:/Users/Admin/Fem_diplom/FEM/Libs/mesh/model/src/ResultFileExporter.cpp:175)

### `metrics.json`

Содержит:

- число узлов и элементов
- общее число DOF, свободные и закрепленные DOF
- времена `assembly/solve/total`
- `linear_iterations`
- `nonlinear_iterations`
- `matrix nnz`
- residuals
- `max_penetration`
- `contact_force_norm`
- `contact_patch_length`
- `max_facet_average_pressure`
- `total_normal_force`
- полезные экстремумы по displacement / reaction / contact force

### `contact_facets.csv`

Facet-level экспорт для контактного постпроцессинга. Для каждой candidate-фасетки записываются:

- `facet_id`, `element_id`, `surface_index`
- флаг `active`
- midpoint в исходной и деформированной конфигурации
- `facet_length`, `active_length`
- `average_gap`, `average_penetration`, `maximum_penetration`
- `integrated_normal_force`
- `average_pressure`

## Как построить графики и проверки

На текущем этапе ParaView не является частью рабочего сценария. Официальный путь такой:

1. C++ пишет `solution.vtu` и `metrics.json`.
2. Python читает эти файлы.
3. Python сам сохраняет все основные графики в `.png`.

### Python: PyVista + matplotlib

Скрипт:

- [scripts/postprocess_results.py](/c:/Users/Admin/Fem_diplom/FEM/scripts/postprocess_results.py:1)

Установка зависимостей:

```powershell
python -m pip install pyvista vtk matplotlib
```

Пример ручного запуска для одиночного case:

```powershell
python scripts/postprocess_results.py results/tire_contact_single_case
```

Пример ручного запуска для серии расчетов:

```powershell
python scripts/postprocess_results.py results/ring_contact_study
```

Скрипт:

- читает `metrics.json` и `.vtu`
- делает базовые consistency-checks
- сам строит PNG по кейсам
- строит summary plot по серии расчетов
- сохраняет обзорный `case_overview.png` с ключевыми метриками
- использует строгий технический стиль с русскими подписями, единицами измерения и шрифтом `Times New Roman`

Типовые файлы, которые появляются рядом с `solution.vtu`:

- `computational_mesh.png`
- `displacement_magnitude.png`
- `sigma_yy.png`
- `von_mises_stress.png`
- `reaction_force_magnitude.png`
- `contact_force_magnitude.png`, если контактное поле есть
- `penetration.png`, если контактное поле есть
- `active_contact_facet.png`
- `candidate_contact_facet.png`
- `case_overview.png`

Для tire-ring case дополнительно сохраняются графики, ближе к формату РПЗ:

- `ring_contour_stress_profiles.png`
  Профили `sigma_rr`, `sigma_tt`, `tau_rtheta` вдоль внутреннего, срединного и внешнего контуров.
- `ring_radial_section_profiles.png`
  Профили `sigma_rr`, `sigma_tt`, `u_r` по радиальному сечению через центр ожидаемого контакта.
- `contact_patch_profiles.png`
  Профили penetration, интегральной нормальной силы и среднего контактного давления по contact facets.

Для этих же графиков рядом сохраняются CSV:

- `ring_contour_stress_profiles.csv`
- `ring_radial_section_profiles.csv`
- `contact_patch_profiles.csv`

Для серии расчетов в корневой папке дополнительно сохраняется:

- `summary_metrics.png`
- `summary_contact_metrics.png`
- `summary_contact_metrics.csv`

Python не содержит логику решателя и используется только для анализа/визуализации.

## Verification

Полная карта Stage 6 с целями сценариев, критериями успеха и обязательным postprocessing
собрана в [docs/verification_scenarios.md](/c:/Users/Admin/Fem_diplom/FEM/docs/verification_scenarios.md:1).

## Как собрать задачу с нуля в пустом `main`

Ниже минимальная последовательность, которой достаточно, чтобы из пустого `main` собрать и решить
задачу контакта шины с жесткой плоскостью.

### 1. Создать `Assembly` и материал

```cpp
auto assembly = std::make_shared<Assembly>();
auto material = std::make_shared<Material>(1, 2.0e5, 0.30, 1.0);
assembly->addMaterial(material);
```

### 2. Создать structured graded Q4 mesh под tire contact

Используй `MeshGenerator::TireContactAnalysisControl` и
`MeshGenerator::generateTireContactAnalysisSetup(...)`.

Минимальный каркас:

```cpp
MeshGenerator meshGenerator(assembly);

MeshGenerator::TireContactAnalysisControl control;
control.mesh.center = Eigen::Vector2d(0.0, 0.58);
control.mesh.innerRadius = 0.25;
control.mesh.outerRadius = 0.50;
control.mesh.startAngle = -120.0 * DEG_TO_RAD;
control.mesh.endAngle = -60.0 * DEG_TO_RAD;
control.mesh.radialLayers = 9;
control.mesh.circumferentialNodes = 61;
control.mesh.materialId = material->getId();

control.mesh.refineCircumferentiallyNearContact = true;
control.mesh.refineRadiallyToOuterSurface = true;
control.mesh.expectedContactCenterAngle = -90.0 * DEG_TO_RAD;
control.mesh.expectedContactHalfAngle = 12.0 * DEG_TO_RAD;
control.mesh.circumferentialRefinementStrength = 6.0;
control.mesh.radialRefinementStrength = 2.5;
control.mesh.candidateFacetWindowScale = 3.0;

control.rigidPlane = RigidPlane2D{Eigen::Vector2d(0.0, 1.0), 0.0};
control.prescribedInnerBoundaryDy = -0.12;

const auto setup = meshGenerator.generateTireContactAnalysisSetup(control);
```

Что делает этот helper:

- строит structured graded ring mesh без hanging nodes
- собирает `candidateContactFacets`
- возвращает `RigidPlane2D`
- накладывает типовые BC на внутреннюю границу
- выбирает anchor-узел для устранения лишнего rigid-body mode

Смотри:

- [meshgenerator.h](/c:/Users/Admin/Fem_diplom/FEM/Libs/mesh/meshGenerator/include/meshgenerator.h:55)
- [meshgenerator.cpp](/c:/Users/Admin/Fem_diplom/FEM/Libs/mesh/meshGenerator/src/meshgenerator.cpp:339)
- [meshgenerator.cpp](/c:/Users/Admin/Fem_diplom/FEM/Libs/mesh/meshGenerator/src/meshgenerator.cpp:426)

### 3. Создать `FEModel` и настроить контакт

```cpp
FEModel model;
model.setAssembly(assembly);
model.setSolverTolerance(1.0e-8);
model.setMaxIterations(40);
model.configureRigidPlaneContact(
    setup.rigidPlane,
    setup.mesh.candidateContactFacets,
    1.0e7);
```

### 4. Решить

```cpp
const bool success = model.solveContact();
if (!success) {
    return 1;
}
```

### 5. Экспортировать результаты

```cpp
ResultFileExportOptions exportOptions;
exportOptions.outputDirectory = std::filesystem::path("results") / "my_case";
exportOptions.baseName = "solution";

const auto artifacts = ResultFileExporter::exportSolution(model, exportOptions);
```

После этого C++ уже положит в папку:

- `results/my_case/solution.vtu`
- `results/my_case/metrics.json`

Дальше построй графики:

```powershell
python scripts/postprocess_results.py results/my_case
```

После постпроцессинга смотри:

- `results/my_case/displacement_magnitude.png`
- `results/my_case/sigma_yy.png`
- `results/my_case/case_overview.png`

## Где смотреть рабочие примеры

Если нужен краткий и чистый пример, смотри:

- [src/src/tire_contact_single_case.cpp](/c:/Users/Admin/Fem_diplom/FEM/src/src/tire_contact_single_case.cpp:1)

Если нужен серийный прогон нескольких кейсов:

- [src/src/ring_contact_study.cpp](/c:/Users/Admin/Fem_diplom/FEM/src/src/ring_contact_study.cpp:1)

Если нужен маленький sanity-case для контактного контура:

- [src/src/block_on_rigid_plane.cpp](/c:/Users/Admin/Fem_diplom/FEM/src/src/block_on_rigid_plane.cpp:1)

## Практическое замечание

В репозитории еще лежит старый файл `src/src/main.cpp` с историческим сценарием.
Каноническим стартовым примером теперь считается target `FEM`, который собирается из
[src/src/tire_contact_single_case.cpp](/c:/Users/Admin/Fem_diplom/FEM/src/src/tire_contact_single_case.cpp:1).
