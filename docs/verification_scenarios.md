# Verification Scenarios

Этот документ фиксирует Stage 6 для текущего FEM-решателя: какой набор сценариев считается
базовым для верификации и инженерного контроля, что именно проверяет каждый сценарий,
по каким метрикам он считается успешным и какие данные должны сохраняться в postprocessing.

Набор deliberately небольшой и привязан к реальному фокусу проекта:

1. базовая линейная задача без контакта;
2. контактный сценарий `контакта нет`;
3. простой rigid-plane contact regression на блоке;
4. `tire -> rigid plane` как основной инженерный сценарий.

## Сводка

| Scenario | Target / Artifact | Purpose |
| --- | --- | --- |
| Basic linear reference | `FEMBasicLinearReference` | Проверка базовой линейной сборки, BC и postprocessing без контакта |
| Contact, plane far away | `FEMContactFarPlaneScenario`, `FEMContactNoContactRegression` | Проверка, что `solveContact()` не портит обычную задачу, если контакта фактически нет |
| Block on rigid plane | `FEMBlockOnRigidPlane`, `FEMContactBlockRegression` | Быстрая регрессия корректности unilateral penalty contact на простой геометрии |
| Tire to rigid plane | `FEM`, `FEMRingContactStudy` | Основной инженерный сценарий контакта шины с жесткой плоскостью |

## 1. Basic Linear Reference

### Target

- `FEMBasicLinearReference`
- Results:
  - `results/basic_linear_reference/solution.vtu`
  - `results/basic_linear_reference/metrics.json`

### Что проверяется

- корректность линейной сборки `K` и `F`;
- корректность наложения fixed DOF и prescribed displacement;
- корректность восстановления полей `u`, `sigma`, `epsilon`, `reaction`;
- базовая работоспособность sparse solver pipeline без контакта.

### Конфигурация

- 2D Q4 block;
- линейно-упругий материал;
- верхняя грань задается перемещением `uy = -0.02`;
- один anchor-узел фиксируется по `x` для устранения лишнего rigid-body mode;
- контакт не включается.

### Метрики успеха

- `success=true`;
- `maximum_prescribed_displacement_error < 1e-12`;
- `equilibrium_residual_norm <= 1e-5`;
- решение без `NaN/Inf`;
- осмысленные `reaction_force` на закрепленной границе;
- `linear_iterations > 0` и корректно заполнены solver metrics.

### Что сохранять в postprocessing

Поля:

- `displacement`, `displacement_magnitude`;
- `stress_2d`, `sigma_xx`, `sigma_yy`, `tau_xy`, `von_mises_stress`;
- `strain_2d`, `strain_xx`, `strain_yy`, `gamma_xy`;
- `reaction_force`, `reaction_force_magnitude`.

Метрики:

- `counts.*`;
- `timings.*`;
- `iterations.linear_iterations`;
- `solver.*`;
- `residuals.*`;
- `matrix.nnz`;
- `extrema.max_displacement_magnitude`;
- `extrema.max_reaction_force_magnitude`;
- `extra.prescribed_top_edge_dy`;
- `extra.maximum_prescribed_displacement_error`;
- `extra.reaction_force_sum_x`, `extra.reaction_force_sum_y`.

## 2. Contact Scenario: Plane Far Away

### Target

- Engineering scenario: `FEMContactFarPlaneScenario`
- Strict regression: `FEMContactNoContactRegression`
- Results:
  - `results/contact_far_plane_scenario/solution.vtu`
  - `results/contact_far_plane_scenario/metrics.json`
  - `results/contact_far_plane_scenario/contact_facets.csv`

### Что проверяется

- контактный контур не должен менять решение, если плоскость вынесена далеко;
- `solveContact()` должен совпадать с обычным `solve()`;
- active set должен оставаться пустым;
- контактные силы, penetration и patch metrics должны оставаться нулевыми.

### Конфигурация

- та же базовая block-задача, что и в линейном reference;
- включен `RigidPlaneContactSolver`;
- плоскость вынесена на `offset = -10.0`;
- penalty parameter `1e7`.

### Метрики успеха

- `success=true`;
- `difference_norm_to_baseline < 1e-10`;
- `active_contact_facets = 0`;
- `max_penetration = 0`;
- `total_normal_force = 0`;
- `contact_patch_length = 0`;
- `minimum_signed_distance > 0`;
- `equilibrium_residual_norm` остается в линейно-адекватном диапазоне.

### Что сохранять в postprocessing

Поля:

- все поля линейного reference;
- `contact_force`, `contact_force_magnitude`;
- `rigid_plane_signed_distance`;
- `rigid_plane_penetration`;
- cell flags `candidate_contact_facet`, `active_contact_facet`.

Метрики:

- все линейные метрики;
- `contact.configured`;
- `contact.candidate_contact_facets`;
- `contact.active_contact_facets`;
- `contact.max_penetration`;
- `contact.total_normal_force`;
- `contact.contact_patch_length`;
- `contact.max_facet_average_pressure`;
- `extra.difference_norm_to_baseline`;
- `extra.far_plane_offset`;
- facet-level `contact_facets.csv` должен быть полностью неактивным.

## 3. Block On Rigid Plane

### Target

- Engineering scenario: `FEMBlockOnRigidPlane`
- Strict regression: `FEMContactBlockRegression`
- Results:
  - `results/block_on_rigid_plane/solution.vtu`
  - `results/block_on_rigid_plane/metrics.json`
  - `results/block_on_rigid_plane/contact_facets.csv`

### Что проверяется

- корректность unilateral frictionless rigid-plane contact на простой геометрии;
- адекватность active set;
- корректность penalty RHS/stiffness;
- ограниченность penetration после сходимости;
- корректность contact force и reaction force.

### Конфигурация

- прямоугольный block;
- верхняя грань нагружается prescribed displacement;
- нижняя внешняя грань участвует в контакте с плоскостью `y = 0`;
- penalty parameter `1e7`.

### Метрики успеха

Для строгой регрессии используются текущие пороги из кода:

- `precontact active facets = 10`;
- `precontact max penetration > 0.04`;
- `converged active_contact_facets = 10`;
- `nonlinear_iterations <= 10`;
- `0 <= max_penetration < 2e-3`;
- `minimum_signed_distance_to_plane > -2e-3`;
- `contact_force_norm > 0`.

Для инженерной оценки дополнительно полезно контролировать:

- `total_normal_force > 0`;
- `contact_patch_length > 0`;
- отсутствие паразитических активных фасеток вне нижней границы;
- разумный уровень `equilibrium_residual_norm`.

### Что сохранять в postprocessing

Поля:

- все поля linear/reference;
- все контактные поля;
- `contact_facets.csv` как основной источник для pressure/patch postprocessing.

Метрики:

- `contact.active_contact_facets`;
- `contact.max_penetration`;
- `contact.total_normal_force`;
- `contact.contact_patch_length`;
- `contact.max_facet_average_pressure`;
- `contact.mean_active_pressure`;
- `contact.center_of_pressure_x`, `contact.center_of_pressure_y`;
- `residuals.linear_relative_residual_norm`;
- `residuals.equilibrium_residual_norm`.

Графики:

- `contact_patch_profiles.png`;
- `case_overview.png`;
- карты `penetration`, `contact_force_magnitude`, `active_contact_facet`.

## 4. Tire To Rigid Plane

### Target

- Single-case: `FEM`
- Parametric study: `FEMRingContactStudy`
- Results:
  - `results/tire_contact_single_case/*`
  - `results/ring_contact_study/*`

### Что проверяется

- работоспособность основного tire-contact pipeline;
- корректность structured graded ring mesh и candidate contact facets;
- корректность нелинейного active-set contact solve;
- адекватность contact patch, penetration и center of pressure;
- устойчивость к изменению penalty parameter;
- сравнение `uniform` vs `graded` сетки на инженерных метриках.

### Конфигурация

- structured graded Q4 ring;
- внутренний контур нагружается prescribed displacement;
- внешний контур контактирует с rigid plane;
- основная penalty study: `1e6`, `1e7`, `1e8`;
- основной single-case: graded mesh, penalty `1e7`.

### Метрики успеха

Single-case:

- `success=true`;
- `active_contact_facets > 0`;
- `contact_patch_length > 0`;
- `total_normal_force > 0`;
- `max_penetration > 0`, но значительно меньше imposed displacement;
- `center_of_pressure_x` близок к плоскости симметрии для симметричного кейса;
- `equilibrium_residual_norm` остается малым;
- нет неактивных/аномальных facet pressures вне ожидаемого patch.

Penalty study:

- с ростом penalty `max_penetration` должен убывать монотонно;
- `total_normal_force` не должен “гулять” хаотично между penalty-level;
- `graded` mesh должна давать не худшую локализацию patch и не деградировать force balance;
- `summary_contact_metrics.csv` должен позволять сравнение по:
  - `contact_patch_length`;
  - `max_average_pressure`;
  - `mean_active_pressure`;
  - `total_normal_force`;
  - `center_of_pressure_arc`.

### Что сохранять в postprocessing

Поля:

- все стандартные displacement/stress/strain/reaction fields;
- все contact fields;
- facet-level `contact_facets.csv`.

Метрики:

- все timing/solver/matrix metrics;
- `contact.max_penetration`;
- `contact.total_normal_force`;
- `contact.contact_patch_length`;
- `contact.max_facet_average_pressure`;
- `contact.mean_active_pressure`;
- `contact.center_of_pressure_x`, `contact.center_of_pressure_y`;
- `extra.mesh_label`;
- `extra.penalty_parameter`;
- mesh diagnostics:
  - `mesh_min_radial_step`, `mesh_max_radial_step`;
  - `mesh_min_angular_step`, `mesh_max_angular_step`;
  - `mesh_min_aspect_ratio`, `mesh_max_aspect_ratio`.

Графики:

- `ring_contour_stress_profiles.png`;
- `ring_radial_section_profiles.png`;
- `contact_patch_profiles.png`;
- `summary_metrics.png`;
- `summary_contact_metrics.png`.

## 5. Hyperelastic Baseline Verification

### Summary

| Scenario | Target / Artifact | Purpose |
| --- | --- | --- |
| Hyperelastic no-contact sanity | `FEMHyperelasticNoContactSanity` | Проверка finite-strain Neo-Hookean ветки в режиме малых деформаций через сравнение с linear baseline |
| Hyperelastic large deformation, no contact | `FEMHyperelasticLargeDeformationNoContact` | Проверка сходимости, `det(F)` и baseline-поведения при конечных деформациях без контакта |
| Hyperelastic contact guard | `FEMHyperelasticContactGuard` | Фиксация текущего ограничения: `solveHyperelastic()` пока не подключен к rigid-plane contact |

### 5.1 Hyperelastic No-Contact Sanity

Target:

- `FEMHyperelasticNoContactSanity`
- Results:
  - `results/hyperelastic_no_contact_sanity/solution.vtu`
  - `results/hyperelastic_no_contact_sanity/metrics.json`

Что проверяется:

- отдельная ветка `FiniteStrainMaterial -> NeoHookeanMaterial -> FiniteStrainQ4Element -> solveHyperelastic()`;
- корректность Newton solve с prescribed displacements в incremental form;
- согласованность finite-strain baseline с линейным решением в режиме малых деформаций;
- корректность finite-strain export pipeline: `Cauchy`, `Green-Lagrange`, `det(F)`, `strain energy`.

Конфигурация:

- 2D Q4 block;
- compressible Neo-Hookean, параметры из `E = 2e5 МПа`, `nu = 0.30`, `thickness = 1 мм`;
- верхняя грань: prescribed `uy = -0.005 мм`;
- нижняя грань: фиксирована по `y`;
- один нижний узел фиксируется по `x`.

Метрики успеха:

- `success=true`;
- `maximum_prescribed_displacement_error < 1e-12`;
- `relative_displacement_difference_to_linear < 5e-3`;
- `minimum_jacobian_determinant > 0.95`;
- `equilibrium_residual_norm <= 1e-5`;
- `nonlinear_iterations > 0`.

Что сохранять в postprocessing:

Поля:

- стандартные `displacement`, `reaction_force`;
- finite-strain point fields: `deformed_position`, `stress_2d` as Cauchy, `strain_2d` as Green-Lagrange, `sigma_zz`, `jacobian_determinant`, `strain_energy_density`;
- finite-strain cell fields: `cell_jacobian_determinant`, `cell_cauchy_stress_2d`, `cell_second_piola_stress_2d`, `cell_green_lagrange_strain_2d`.

Метрики:

- `formulation.*`;
- `finite_strain.*`;
- `iterations.*`;
- `residuals.*`;
- `extra.relative_displacement_difference_to_linear`;
- `extra.minimum_jacobian_determinant`, `extra.maximum_jacobian_determinant`.

Ограничения:

- это baseline displacement-only Q4;
- сценарий intentionally small-strain-like и нужен именно как sanity-check новой архитектуры, а не как демонстрация преимуществ hyperelasticity.

### 5.2 Hyperelastic Large Deformation Without Contact

Target:

- `FEMHyperelasticLargeDeformationNoContact`
- Results:
  - `results/hyperelastic_large_deformation_no_contact/solution.vtu`
  - `results/hyperelastic_large_deformation_no_contact/metrics.json`

Что проверяется:

- сходимость `solveHyperelastic()` на конечных деформациях;
- экспорт `det(F)` и strain-energy полей;
- поведение baseline Neo-Hookean при nearly incompressible материале;
- корректная диагностика current limitation по locking-risk.

Конфигурация:

- 2D Q4 block;
- compressible / nearly incompressible Neo-Hookean с `E = 25 МПа`, `nu = 0.48`, `thickness = 1 мм`;
- верхняя грань: prescribed `uy = -0.10 мм`;
- нижняя грань: фиксирована по `y`;
- один нижний узел фиксируется по `x`.

Метрики успеха:

- `success=true`;
- `maximum_prescribed_displacement_error < 1e-12`;
- `minimum_jacobian_determinant > 0.20`;
- `maximum_jacobian_determinant < 2.50`;
- `maximum_strain_energy_density > 0`;
- `maximum_second_piola_stress_norm > 0`;
- `has_near_incompressible_material = true`;
- `equilibrium_residual_norm <= 1e-5`.

Что сохранять в postprocessing:

- все finite-strain поля из sanity-case;
- `jacobian_determinant.png`;
- `strain_energy_density.png`;
- `case_overview.png` с finite-strain metrics.

Ограничения:

- это еще не mixed `u/p` и не `SRI`;
- для `nu = 0.48` сценарий показывает only baseline architecture and diagnostics;
- отсутствие direct-solver fallback в конкретном кейсе не означает, что near-incompressible режим всегда будет комфортен для `CG/IC`.

### 5.3 Hyperelastic Block On Rigid Plane

Target:

- `FEMHyperelasticBlockOnRigidPlane`
- Results:
  - `results/hyperelastic_block_on_rigid_plane/solution.vtu`
  - `results/hyperelastic_block_on_rigid_plane/metrics.json`
  - `results/hyperelastic_block_on_rigid_plane/contact_facets.csv`

Что проверяется:

- первый runnable `hyperelastic + penalty rigid-plane contact` path;
- корректная сборка остатка `F_ext + F_contact - F_int(u)` и касательной `K_t + K_contact`;
- работоспособность existing rigid-plane penalty solver на finite-strain Q4 facets;
- совместимость finite-strain export pipeline с contact fields.

Конфигурация:

- 2D Q4 block;
- compressible Neo-Hookean с `E = 1e3 МПа`, `nu = 0.30`, `thickness = 1 мм`;
- верхняя грань: prescribed `uy = -0.08 мм`;
- нижняя грань взаимодействует с rigid plane `y = 0`;
- один верхний узел фиксируется по `x`;
- penalty parameter `1e6`.

Метрики успеха:

- `success=true`;
- `active_contact_facets > 0`;
- `0 < max_penetration < 2e-2`;
- `contact_force_norm > 0`;
- `minimum_jacobian_determinant > 0.90`;
- finite-strain solve converges without breaking contact export.

Что сохранять в postprocessing:

- все finite-strain поля;
- все contact fields;
- `contact_facets.csv`;
- `penetration.png`, `contact_force_magnitude.png`, `active_contact_facet.png`;
- `jacobian_determinant.png`, `strain_energy_density.png`.

Ограничения:

- это только penalty-contact branch;
- `hyperelastic + augmented Lagrangian` пока не подключен;
- contact kinematics пока small-motion-in-contact-surface sense, without a fully consistent large-sliding linearization.

### 5.4 Hyperelastic + Rigid-Plane Contact Guard

Target:

- `FEMHyperelasticContactGuard`

Что проверяется:

- что current code path не пытается silently решать неподключенный
  `hyperelastic + augmented Lagrangian rigid-plane contact`;
- ограничение фиксируется явно и воспроизводимо.

Критерий успеха:

- test passes only if `solveHyperelastic()` returns `false` when
  augmented-Lagrangian contact solver is configured.

Ограничение:

- это не инженерный contact-case, а guard regression на текущее архитектурное состояние.

### 5.5 Tire-Oriented Hyperelastic Scenario

Текущее состояние:

- пока не оформлен как runnable regression/scenario.

Причина:

- hyperelastic penalty-contact path уже есть для простого блока, но tire-oriented
  hyperelastic pipeline еще не подключен и не верифицирован;
- до этого шага tire-oriented hyperelastic case дал бы ложное ощущение готовности.

Что должно появиться позже:

- simplified tire ring with hyperelastic material;
- затем `hyperelastic + rigid-plane contact` for tire-like geometry;
- после этого — уже сравнение linear elastic vs Neo-Hookean tire response.

## Practical Use

Рекомендуемый минимальный Stage 6 workflow:

1. Запустить `FEMBasicLinearReference`.
2. Запустить `FEMContactFarPlaneScenario`.
3. Запустить строгие regressions:
   - `FEMContactNoContactRegression`
   - `FEMContactBlockRegression`
   Быстрый automated run:

   ```powershell
   "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\ctest.exe" --test-dir build/src -C Release --output-on-failure
   ```
4. Запустить инженерные contact scenarios:
   - `FEMBlockOnRigidPlane`
   - `FEM`
   - `FEMRingContactStudy`
5. Прогнать Python postprocessor по `results/...` и убедиться, что:
   - summary-графики построились;
   - facet-level contact metrics согласованы с `metrics.json`;
   - нет аномалий по penetration/pressure/center-of-pressure.
