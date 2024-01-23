# MicrotubuleSpringModel

This code generates microtubules which are modelled as beads representing tubulin monomers connected with up to four neighbours via bonds.

## Installation

Clone or pull the repository. Once downloaded, start a Julia repl in this directory, then run
```
using Pkg
Pkg.instantiate()
```
Alternatively, you can type `]` to go to the package manager in the repl and simply type `instantiate` there.

**This will not work on the server since GLMakie will not work. You should `Pkg.rm("GLMakie") on the server.Then remove any imports/uses of GLMakie in the MicrotubuleSpringModel package and scripts you wish to run.**

## Configurations

The configurations for running as simulation are specified in the `MicrotubuleConfig` or the `PatchConfig` structs.

The `Configurations.jl` package is used to load these parameters from a `.toml` file stored in the `config/` directory.

The configuration is broken down into lattice parameters, spring constants, bond angles, and iteration parameters. Most of these parameters will not change between simulations/experiments besides the iteration parameters. Using a simple `euler` scheme seems more efficient than the `rk4` scheme since using larger timesteps with `rk4` seems to be unstable for some reason. These deterministic schemes are useful for mechanical tests.

## Simulating

Once a config has been loaded, the MT can be constructed using the `create_lattice` function, which returns the mutable positions/orientations separately from the static information contained in `BeadPars`.

Examples of various simulations for specific tests can be found in `experiments/`.

## Scripts

In the `scripts/` directory, various files can be found to simulate and visualise MTs. Hopefully the context is clear even though they are usually rough workings.

`lengthscales.jl` has the derivation of the values for the damping constants as well as the Brownian motion values and equilibrium spring lengths laterally.

## TODO 

- [ ] ensure lattice can be constructed with different protofilament numbers etc 
- [ ] ensure the lattice is then still stable with the same natural bond angles.
- [ ] Test various kinesin mechanisms 

## Potential extensions

- [ ] In theory the simulations could be faster if the algorithm iterates over the springs and distributes the forces to the beads rather than directly iterating over the beads.
- [ ] 