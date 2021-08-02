# Monte Carlo simulations in the NPT ensemble

This is actually a fork of [mcfort](https://github.com/edwinb-ai/mcfort)
but this version works exclusively in the NPT ensemble.

It was made for the sole purpose of computing the **equation of state**
of fluids using the NPT ensemble method and the Monte Carlo simulation
technique.

It currently has several potentials implemented, but the source code must
be modified by hand to enable them.

It has basic multithreading capabilities to accelerate the energy computation.

## Build

To build and use this code, first create a build directory

```shell
mkdir build
```

then, inside this directory run the following

```shell
cmake -GNinja ..
```

This command expects that you have installed the [ninja](https://ninja-build.org/)
build system.

Then, just run

```shell
ninja
```

In case you only have `make` just run

```shell
cmake ..
```

and then run

```shell
make
```

Either way, an executable `mcnpt` will be created inside the build directory.
This is the main executable.

## Usage

The `mcnpt` executable expects a file called `input.in` with *six values* in it:

- _Packing fraction_, a value between strictly larger than 0 and strictly less than 1.0. Any value in between is acceptable. This value will determine the density of the system.
- _Reduced temperature_.
- _Reduced pressure_.
- _Volume displacement_, logarithm based. This is just an initial value, as the code already handles the volume displacement automatically to set it between 15-20% acceptance.
- _Number of particles_, an integer. The larger this value, the longer it will take to run a simulation.
- _Number of Monte Carlo cycles_, this is the value to use in order to equilibrate and accumulate results. The default is to use half this number for equilibration and half for averaging results.

## Observables

The only observable that this computes is the average density given a value for the reduced
pressure. It also provides the standard error for this average.

With both the pressure and the density, one can plot the equation of state for a given
potential.

If you are looking for a way to compute the radial distribution function or the structure
factor, please do look into [mcfort](https://github.com/edwinb-ai/mcfort).
