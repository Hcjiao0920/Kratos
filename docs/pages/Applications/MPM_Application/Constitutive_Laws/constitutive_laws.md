---
title: Constitutive Laws
keywords: mpm constitutive laws
tags: [mpm constitutive laws]
sidebar: mpm_application
summary: 
---

In this section we discuss the constitutive laws that are implemented in the `MPMApplication` and for each of them we list the variables that must be defined in the `ParticleMaterials.json` file.

## Linear Elasticity

The **linear elastic** constitutive law is identified by the following labels:

- `LinearElasticIsotropicPlaneStrain2DLaw`: two-dimensional problem, plane strain formulation;
- `LinearElasticIsotropicPlaneStress2DLaw`: two-dimensional problem, plane stress formulation;
- `LinearElasticIsotropic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`

## Hyperelastic NeoHookean

The **hyperelastic NeoHookean** constitutive law is identified by the following labels:

- `HyperElasticNeoHookeanPlaneStrain2DLaw`: two-dimensional problem, plane strain and irreducible formulation;
- `HyperElasticNeoHookeanPlaneStrainUP2DLaw`: two-dimensional problem, plane strain and mixed formulation;
- `HyperElasticNeoHookean3DLaw`: three-dimensional problem, irreducible formulation;
- `HyperElasticNeoHookeanUP3DLaw`: three-dimensional problem, mixed formulation.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`

## Mohr Coulomb

The **plastic Mohr Coulomb** constitutive law is identified by the following labels:

- `HenckyMCPlasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `HenckyMCPlastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`
- `COHESION`
- `INTERNAL_FRICTION_ANGLE`
- `INTERNAL_DILATANCY_ANGLE`

## Mohr Coulomb Strain Softening

The **Mohr Coulomb** with **Strain Softening** constitutive law is identified by the following labels:

- `HenckyMCStrainSofteningPlasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `HenckyMCStrainSofteningPlastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`
- `COHESION`: cohesion (peak)
- `COHESION_RESIDUAL`: cohesion (residual)
- `INTERNAL_FRICTION_ANGLE`: internal friction angle (peak)
- `INTERNAL_FRICTION_ANGLE_RESIDUAL`: internal friction angle (residual)
- `INTERNAL_DILATANCY_ANGLE`: internal dilatancy angle (peak)
- `INTERNAL_DILATANCY_ANGLE_RESIDUAL`: internal dilatancy angle (residual)
- `SHAPE_FUNCTION_BETA`: exponential softening beta coefficient

## Modified Cam Clay

The **modified cam clay** constitutive law is identified by the following labels:

- `HenckyBorjaCamClayPlasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `HenckyBorjaCamClayPlastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `PRE_CONSOLIDATION_STRESS`: preconsolidation pressure
- `OVER_CONSOLIDATION_RATIO`: over Consolidation Ratio (OCR)
- `SWELLING_SLOPE`: slope of swelling line
- `NORMAL_COMPRESSION_SLOPE`: slope of Normal Consolidation Line (NCL)
- `CRITICAL_STATE_LINE`: slope of Critical State Line (CSL)
- `INITIAL_SHEAR_MODULUS`: initial Shear Modulus
- `ALPHA_SHEAR`: volumetric-deviatoric coupling constant

## Newtonian Fluid

The **displacement-based Newtonian fluid** constitutive law is identified by the following labels:

- `DispNewtonianFluidPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `DispNewtonianFluid3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `BULK_MODULUS`
- `DYNAMIC_VISCOSITY`

## Sand Hypoplasticity

The **sand hypoplasticity** constitutive law implements the rate-type model of von Wolffersdorff (1996) with the intergranular-strain extension of Niemunis & Herle (1997), integrated with an adaptive error-controlled Runge-Kutta-Fehlberg 2(3) substepping scheme (Fellin & Ostermann). It is identified by the following labels:

- `SandHypoplasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `SandHypoplastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `SAND_HYPOPLASTIC_PARAMETERS`: vector of exactly 16 entries, in the order
  1. `phi`: critical-state friction angle (degrees)
  2. `p_t`: tensile shift of the mean pressure (kPa)
  3. `h_s`: granular hardness (kPa)
  4. `n`: exponent of the Bauer compression law
  5. `e_d0`: reference void ratio at maximum density
  6. `e_c0`: reference critical-state void ratio
  7. `e_i0`: reference void ratio at maximum looseness
  8. `alpha`: pyknotropy exponent
  9. `beta`: barotropy exponent
  10. `m_R`: intergranular-strain stiffness multiplier for 180-degree strain-path reversals
  11. `m_T`: intergranular-strain stiffness multiplier for 90-degree strain-path reversals
  12. `R_max`: intergranular-strain reference length
  13. `beta_r`: intergranular-strain evolution exponent
  14. `chi`: intergranular-strain interpolation exponent
  15. `bulk_w`: bulk modulus of the pore fluid (0 for dry conditions)
  16. `e_0`: initial void ratio (values > 10 are interpreted as raw void ratio + 10; values <= 10 are scaled with the Bauer law from the initial mean stress)
- `INITIAL_STRESS_VECTOR` (optional): initial Cauchy stress at the material points, Voigt notation `[s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]`; hypoplasticity is a rate-type model and requires a physically meaningful (compressive) initial stress state.

Setting the intergranular-strain multiplier `m_R < 0.5` disables the intergranular-strain extension.
