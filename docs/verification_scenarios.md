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
