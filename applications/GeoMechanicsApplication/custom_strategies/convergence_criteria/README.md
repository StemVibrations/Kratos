# Local error criteria

`GeoLocalErrorCriteria` implements local error criteria. On top of a global (equilibrium) criterion, it
checks the accuracy of the individual stress points:

- inaccurate plastic points of soil elements 
- non-linear inaccurate elastic points of soil elements 
- inaccurate plastic points of interfaces 

## Definition

For every stress point $j$ the equilibrium stress $\sigma_{eq,j}$ (the linearized response) is compared with
the constitutive stress $\sigma_{c,j}$ (the actual material response):

$$\sigma_{eq,j} = \sigma_{c,j-1} + D^e \delta\varepsilon_j$$

$$\text{LocalError}_j = \frac{\left\| \sigma_{eq,j} - \sigma_{c,j} \right\|}{\max\left(\tau_{max,j},\ c,\ \sigma_{ref}\right)}$$

where $\sigma_{ref}$ is the `minimum_reference_stress` (1 kPa by default) for plastic points and $p_{ref}/200$
for elastic points (`reference_pressure` $p_{ref}$ is 100 kPa by default). For interface points, only the
shear tractions are compared ($\left\| \tau_{eq,j} - \tau_{c,j} \right\|$, such that a changing normal traction (e.g. of an interface that opens) does not make a point
inaccurate. A stress point is inaccurate when its local error exceeds the tolerated local error. Convergence requires, for each kind of stress point:

$$N_{inaccurate} < f \cdot N + 3$$

where $f$ is the tolerated fraction of inaccurate points (0.1 by default) and $N$ is the number of soil
plastic points, soil elastic points (including those of inactive elements) or interface plastic points.

The criterion does not store the equilibrium stress itself. Every constitutive law reports the deviation of
its constitutive stress from the elastic predictor of the step,
$P = \sigma_0 + D^e (\varepsilon - \varepsilon_0) - \sigma_c$. Since $\sigma_0$ and $D^e$ are constant within a
step, $\sigma_{eq,j} - \sigma_{c,j} = P_j - P_{j-1}$, with $P_0 = 0$ at the start of each step.

## Usage

In the `solver_settings` of the project parameters:

```json
"convergence_criterion": "residual_criterion",
"use_local_error_criteria": true,
"local_error_criteria_settings": {
    "tolerated_error": 0.01,
    "minimum_reference_stress": 1.0,
    "reference_pressure": 100.0
}
```

The local error criteria are then combined with the chosen global criterion (both need to be satisfied).
All settings (and their defaults) are:

| Setting                                                 | Default | Description                                           |
|---------------------------------------------------------|---------|-------------------------------------------------------|
| `tolerated_error`                                       | 0.01    | default of the three tolerated local errors below     |
| `tolerated_soil_plastic_local_error`                    | 0.01    |                                                       |
| `tolerated_soil_elastic_local_error`                    | 0.01    |                                                       |
| `tolerated_interface_plastic_local_error`               | 0.01    |                                                       |
| `tolerated_inaccurate_soil_plastic_point_fraction`      | 0.1     |                                                       |
| `tolerated_inaccurate_soil_elastic_point_fraction`      | 0.1     |                                                       |
| `tolerated_inaccurate_interface_plastic_point_fraction` | 0.1     |                                                       |
| `number_of_always_tolerated_inaccurate_points`          | 3       |                                                       |
| `minimum_reference_stress`                              | 1.0e3   | 1 kPa (in the stress unit of the model)               |
| `reference_pressure`                                    | 1.0e5   | $p_{ref}$ = 100 kPa (in the stress unit of the model) |
| `echo_level`                                            | 1       | prints the numbers of (inaccurate) stress points      |

The defaults of `minimum_reference_stress` and `reference_pressure` assume that stresses are expressed in Pa.
Use 1.0 and 100.0 for models in kPa.

The evaluation needs the material responses at the latest solution, therefore the criterion asks the
strategy to rebuild the right-hand side before `PostCriteria` is called.

## Supported constitutive laws

Constitutive laws take part by implementing `LocalErrorDataProvider`
(see `custom_constitutive/local_error_data_provider.h`).

| Constitutive law                                   | Stress point | Plastic when                         | $D^e$                          |
|----------------------------------------------------|--------------|--------------------------------------|--------------------------------|
| `GeoMohrCoulombWithTensionCutOff2D/3D`             | soil         | the elastic predictor is inadmissible | elastic matrix                |
| `GeoInterfaceCoulombWithTensionCutOff`             | interface    | the elastic predictor is inadmissible | normal/shear stiffness        |
| `SmallStrainUDSM2DPlaneStrainLaw`, `SmallStrainUDSM3DLaw` | soil  | IPL of the UDSM is non-zero          | UDSM, IDTASK = 6              |
| `SmallStrainUDSM2D/3DInterfaceLaw`                 | interface    | IPL of the UDSM is non-zero          | UDSM, IDTASK = 6              |
| `SmallStrainUMAT3DLaw`, `SmallStrainUMAT2DPlaneStrainLaw` | soil  | the stress deviates from the predictor | UMAT tangent for a zero strain increment |
| `SmallStrainUMAT*InterfaceLaw`                     | interface    | the stress deviates from the predictor | UMAT tangent for a zero strain increment |
| linear elastic soil laws                           | soil         | never                                | -                              |
| linear elastic interface laws                      | interface    | never                                | -                              |
| `LinearElastic2DBeamLaw`, `TrussBackboneConstitutiveLaw` | not considered | -                        | -                              |

Notes:

- Only stress points with a stress-dependent elastic stiffness can be inaccurate *elastic* points. This
  information is only available for UDSMs (the ISTRSDEP attribute).
- UMATs provide neither a plasticity indicator nor an elastic stiffness matrix. Their elastic predictor
  is based on the stiffness that the UMAT returns for a zero strain increment from the state at the start
  of the step, and non-linear elastic UMAT stress points are regarded as plastic points.
- The cohesion is taken from `GEO_COHESION` or from the `UMAT_PARAMETERS` (at `INDEX_OF_UMAT_C_PARAMETER`);
  when neither is available, zero is used.
- The tolerated number of inaccurate interface points is based on the number of *interface* plastic points.

# Global force error criterion

`GeoGlobalForceErrorCriteria` implements a check on force residuals:

$$\text{GlobalErrorForce} = \frac{\left\| r \right\|_2}{CSP \left\| f_{ext}^{inact} \right\|_2 + \left\| f_{int}^{act} \right\|_2} < \text{tolerance}$$

where

- $r$ is the out-of-balance force vector (the residual of the strategy),
- $f_{ext}^{inact}$ are the external forces that were already in equilibrium at the start of the stage,
  i.e. the internal forces at the start of the stage, $f_{int,0}$,
- $f_{int}^{act} = f_{int} - f_{int,0}$ are the internal forces that have developed in the stage so far,
- $CSP$ is the current stiffness parameter (see below).

Only the free displacement degrees of freedom are taken into account. The residuals of other degrees of
freedom (e.g. water pressures) need to be checked by another criterion, and rotations are not checked (a
check on moment residuals is not implemented).

## Reference forces

The criterion needs the reference forces from the strategy, so it can only be used with the Newton-Raphson
strategy of this application (`GeoMechanicsNewtonRaphsonStrategy`, i.e. `ResidualBasedNewtonRaphsonStrategyTwo`).
That strategy applies a load fraction $\lambda$ of the stage unbalance $r_{stage} = f_{ext} - f_{int,0}$ in
each step, so its residual is $r = f_{ext} - f_{int} - (1 - \lambda)\, r_{stage}$. Hence:

$$f_{int}^{act} = f_{int} - f_{int,0} = \lambda\, r_{stage} - r$$

This assumes that the external forces don't change within a stage, other than through the stage unbalance.

The internal forces at the start of the stage are evaluated once per stage, element by element, as
$f_{int,0} = f_{ext,e} - (f_{ext,e} - f_{int,e})$, where $f_{ext,e}$ is the `EXTERNAL_FORCES_VECTOR` of the
element (its body forces) and $f_{ext,e} - f_{int,e}$ is its right-hand side. Elements that don't provide
`EXTERNAL_FORCES_VECTOR` are regarded as unloaded, i.e. their self weight is not part of $f_{ext}^{inact}$.
The following elements provide it:

- the `UPwSmallStrainElement` family (including the FIC, updated Lagrangian and axisymmetric variants),
- `SmallStrainUPwDiffOrderElement`,
- `UPwInterfaceElement`,
- the `LinearTrussElement` and `LinearTimoshenkoBeamElement2D2N` of the StructuralMechanicsApplication.

A criterion that is passed directly to the strategy is linked automatically. When it is part of a combined
criterion (e.g. together with the local error criteria), it is linked by means of
`SetGlobalForceErrorCriteria`, which the Python solver takes care of.

## Current stiffness parameter

$$CSP = \frac{\int_{V_{act}} \Delta\varepsilon \cdot \Delta\sigma \, dV}{\int_{V_{act}} \Delta\varepsilon \cdot D^e \Delta\varepsilon \, dV}$$

i.e. the ratio of the total and the elastic strain energy increments of the soil stress points, where the
increments are taken with respect to the start of the step. The CSP equals one when all stress points behave
elastically and it approaches zero at failure, so the inactive loads no longer contribute to the reference
once failure is approached. It is limited to the range $[0, 1]$.

Note that this ratio is sometimes defined the other way around (elastic over total). That would give a
CSP that grows without bounds at failure, which would make the criterion ever more lenient when failure is
approached, rather than stricter. Therefore, the ratio is taken as total over elastic.

The strain energy increments are provided by the constitutive laws through `LocalErrorDataProvider`
(`TotalStrainEnergyIncrement` and `ElasticStrainEnergyIncrement`): the Mohr-Coulomb model, UDSMs, UMATs and the
incremental linear elastic law. Interface and structural stress points are not taken into account. The volume
of a stress point is the domain size of its element divided by its number of stress points.

## Usage

In the `solver_settings` of the project parameters:

```json
"strategy_type": "newton_raphson",
"convergence_criterion": "global_force_error_criterion",
"global_force_error_criterion_settings": {
    "force_relative_tolerance": 0.01
}
```

Use `"global_force_error_and_water_pressure_criterion"` to check the water pressures as well (by means of the
`water_pressure_relative_tolerance` and `water_pressure_absolute_tolerance`). The criterion can be combined with
the local error criteria (`"use_local_error_criteria": true`). All settings (and their defaults) are:

| Setting                    | Default | Description                                             |
|----------------------------|---------|---------------------------------------------------------|
| `force_relative_tolerance` | 0.01    | tolerated global force error                            |
| `force_absolute_tolerance` | 1.0e-9  | - (a residual norm below this value is always accepted) |
| `echo_level`               | 1       | prints the global force error, CSP and the norms        |
