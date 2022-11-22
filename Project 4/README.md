# Ising model in two dimensions to simulate ferromagnet enviroment
We use the Ising model to simulate the temperature-dependent behavior in ferromagnets. We want to estimate the critical temperature where the system undergoes a phase transition, where the ferromagnets go from magnetized to a non-magnetized state.

### Content list
- `main.cpp` : main program used to run the different cases
- `plot.py` : python script to plot the results
- `common.py` : script for reading data files and calculating $C_V$, $|m|$, $\epsilon$ and $\chi$. 
- `src/` : folder containing the source code
- `include/` : folder containing the header files
- `figs/` : folder containing the figures
- `datasets/` : folder containing the data files

### Run instructions

Vildes run stuff

// no parallell, don't delete..
// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ main.cpp -std=c++11 -larmadillo src/* -I include -o main

// try with optimisation!! much faster.
// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ -O2 main.cpp -std=c++11 -larmadillo src/* -I include -o main

// with parallell
// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++-12 -O2 main.cpp -std=c++11 -larmadillo src/* -I include -fopenmp -o main


### Notes

## Further information to describe the code structure
