# README

A singularity container was used for the `meffil` environment. `plink` and `king` are also available in the singularity container.

The singularitry definition file `meffil.def` includes instructions to build the container.
A pre-built container is hosted at `quay.io` and can be pulled using the following command:

```bash
singularity pull oras://quay.io/datatecnica/meffil
```