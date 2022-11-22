# Ising model in two dimensions to simulate ferromagnet enviroment
We use the Ising model to simulate the temperature-dependent behavior in ferromagnets. We want to estimate the critical temperature when the system undergoes a phase transition, in which the ferromagnets go from a magnetized to a non-magnetized state. 
Our code is made up of one main class which runs the Markov chain Monte Carlo simulation in parallel and extracts the desired results in the `datasets/`-folder. 

### Content list
- `main.cpp` : main program used to run the different cases
- `plot.py` : python script to plot the results
- `common.py` : script for reading data files and calculating $|m|$, $\epsilon$, $C_V$ and $\chi$. 
- `analytic_numeric.py` : script for comparing the analytic and numeric results
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

For benchmarking, run
```sh
python analytic_numeric.py
```

### Dependencies
- Armadillo
- Python 3.7
- GNU g++ version 12 (Install with `brew install gcc@12` on OSX)


### Notes

When running `main.cpp` the program will produce a number of data files in the `datasets/` folder, and print to the terminal the numerical results for a $2 \times 2$-lattice. 

The `plot.py` script will produce the figures in the `figs/` folder.    
