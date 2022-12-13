# Simulate the time-dependent Schrödinger equation
We wanted to see the quantum effects of a particle moving through two slits. This program simulates a 1 x 1 box with dirichlet conditions at the boundary, so that the edges work as a reflective surface. Particle is represented by a probability distribution and the simulation utilizes the Born rule to see the effects. It gives the probability of one-slit and three-slit aswell with a wall at x = 0.8.

### Content list
- `main.cpp` : Main program to run different cases.
- `common.py` : Declaration of python functions for vizualisation.
- `plot.py` : Python script to vizualise the results.
- `src/` : Folder containing the source code.
- `include/` : Folder containing the header files.

### Run instructions (OSX and Linux)
Compile (OSX):
```sh
g++
```
Run program with (OSX):
```sh
./main
```

Compile and run (Linux):
```sh
g++ main.cpp src/* -I include -std=c++11 -larmadillo -o main; ./main
```

Produce plots with:
```sh
python plot.py
```

### Dependencies
- Armadillo 11.3
- Python 3
- GNU g++ version 12

### Notes

