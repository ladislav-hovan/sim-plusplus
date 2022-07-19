# sim-plusplus

## Description

C++ code for doing simple molecular dynamics of particles with Lennard-Jones potential.
It supports initial position and velocity input from simple xyz-style files, outputs atom positions and overall kinetic and potential energies.
The default parameters correspond to Argon-40 atoms.

## Compilation

A compiler compatible with the C++ 17 standard is needed (due to inline variables in `constants.h`).

```
g++ *.cpp -std=c++17 -o sim++ -O3
```

## Usage

Assuming the file `input.dat` contains the input parameters, the command to launch is:

```
./sim++ input.dat
```

The input file can also be omitted, in which case the default parameters are used.