# Monte Carlo simulations of pseudo hard spheres

This code can run Monte Carlo simulations using a pseudo hard sphere
potential.
It is written from the ground up to be very fast and efficient.

## Build

To build and use this code, first create a build directory

```shell
mkdir build
```

then, inside this directory

```shell
cmake -GNinja ..
```

This command expects that you have installed the [ninja](https://ninja-build.org/)
build system.

In case you only have `make` just run

```shell
cmake ..
```

Either way, an executable `mcfort` will be created inside the build directory.
This is the main executable.

## Usage

The `mcfort` executable expects a file called `input.in` with *four values* in it:

- Packing fraction, a value between (0, 0.49]. This value will determine the density of the system.
- number of particles, an integer. The larger this value, the longer it will take to run a simulation.
- spacing for S(q), the structure factor needs a spacing value. The larger this value, the longer it will take to run a simulation.
- Number of Monte Carlo cycles, this is the value to use in order to compute the observables.

## Observables

By default, this code will compute two main observables:

- The radial distribution function, also known as g(r);
- the structure factor, from the definition, also known as S(q).

These will appear as files named `gr.dat` and `sq.dat`.

## Caveats

There is no acceleration on this code. This is plain, old O(n^2) complexity.
Although, it is very robust, and works as expected. It might take a long time, but it *will*
give you the correct, physically-relevant result.

The only potential supported is the pseudo hard sphere potential.
Maybe in the future other potentials might be added.
