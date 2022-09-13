# FYS4150: Project 1 by Aline, Even Tobias, and Vilde


To compile and link the C++ code, use
> g++ main.cpp -std=c++11 src/*.cpp -I include -o main.exe   

Run with 
> ./main

To generate plots and tables, run
> python plot.py

## File content overview

- **main.cpp**: runs the main program, to generate exact values, approximated values as well as timing the special and general algorithms
- **algorithms.cpp + .hpp**: includes functions to compute respective algorithms and time them
- **utils.cpp +.hpp**: includes basic reusable functions, like writing to file and printing numbers in scientific notation
- **plot.py**: script for generating plots of exact and numerical solutions, error plots, error table, and timing table from data computed in main.cpp
