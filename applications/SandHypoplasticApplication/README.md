# SandHypoplasticApplication

Project-local Kratos application that hosts the parallel sand-hypoplastic constitutive-law paths outside `MPMApplication`.

`MPMApplication` remains the particle-method host. This application owns the sand-hypoplastic law registrations and shared variables used by Path A/B/C tests.

Currently registered:

- `SandHypoplasticFortranDllLaw` (Path A dynamic-link Fortran bridge)
- `SAND_HYPO_FORTRAN_DLL_PATH`
- `SAND_HYPOPLASTIC_PROPS_16`

`INITIAL_STRESS_VECTOR` is not redefined here; use the Kratos core variable.