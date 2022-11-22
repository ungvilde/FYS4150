# Ising model in two dimensions to simulate ferromagnet enviroment
We use the Ising model to simulate the temperature-dependent behavior in ferromagnets. We want to estimate the critical temperature where the system undergoes a phase transition, where the ferromagnets go from magnetized to a non-magnetized state.

### Content list
- `main.cpp` : main program used to run the different cases
- `plot.py` : python script to plot the results
- `common.py` : script for reading data files and calculating $|m|$, $\epsilon$, $C_V$ and $\chi$. 
- `src/` : folder containing the source code
- `include/` : folder containing the header files
- `figs/` : folder containing the figures
- `datasets/` : folder containing the data files

### Run instructions (OSX and Linux)

Compile (OSX): 
```sh
g++-12 -O3 main.cpp -std=c++11 -larmadillo src/* -I include -fopenmp -o main
```

Compile (Linux): 
```sh
g++ -O3 main.cpp -std=c++11 -larmadillo src/* -I include -fopenmp -o main
```

Run with: 
```sh
./main
```

Produce plots with: 
```sh
python plot.py
```

### Dependencies
- Armadillo
- Python 3.7
- GNU g++ version 12 (Install with `brew install gcc@12` on OSX)


### Notes

## Further information to describe the code structure
