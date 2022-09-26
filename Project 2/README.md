# Project 2 - Scaling equations, eigenvalue problems, code testing

## Project description
For this project, we have implemented the Jacobi rotation algorithm for solving a one-dimensional buckling beam problem. The Jacobi rotation algorithm is an iterative method that computes a numerical approximation of the eigenvalues and eigenvectors of real symmetrical matrices. 

## Instructions
To compile and link the C++ code, use
> g++ main.cpp -std=c++11 -larmadillo src/*.cpp -I include -o main.exe

And run with
> ./main.exe

To generate plots, use
> python plot.py

## Overview
- **main.cpp** is the main program where we test our functions, and compute numerical and analytical solutions
- **algorithm.cpp + .hpp** contains the functions relevant for executing the jacobi rotation algorithm, i.e. finding the largest off-diagonal element, doing the rotation, and iterating the rotation until the threshold requirement is met.
- **utils.cpp + .hpp** contains reusable functions, such as creating a tridiagonal matrix, finding analytical solutions, and comparing the numerical and analytical solutions
- **plot.py** is a script for generating the plots
