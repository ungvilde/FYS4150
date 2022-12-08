# Simulate the time-dependent Schrödinger equation
We want to see the effects of the quantum waves through one, two and three slites.

### Content list
- `main.cpp` : Main program to run different cases.
- `plot.py` : Python script to vizualise the results.
- `src/` : Folder containing the source code.
- `include/` : Folder containing the header files.

### Run instructions (OSX and Linux)
Compile (OSX):
```sh
g++
```

Compile (Linux):
```sh
g++ main.cpp src/simulate_quantum.cpp -I include -std=c++14 -larmadillo -o main
```

Run program with:
```sh
./main
```

Produce plots with:
```sh
python plot.py
```

### Dependencies
- Armadillo 11.3
- Python 3.11
- GNU g++ version 12

### Notes

