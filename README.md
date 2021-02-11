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

The `mcfort` executable expects a file called `input.in`