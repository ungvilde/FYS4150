#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>
#include <time.h>
#include <armadillo>

#include "omp.h" 
#include "IsingModel.hpp"


void run_ensemble(arma::vec temperatures, int n_cycles, int n0, int L);
void compute_2x2_analytic_solution();
void compute_2x2_numeric_solution();

int main()
{
    // ---------
    // Problem 4
    // ---------
    
    // compute_2x2_analytic_solution();
    // compute_2x2_numeric_solution();

    // ----------
    // Problem 5:
    // ----------
    // arma::arma_rng::set_seed(2022);
    // int L = 20;
    // IsingModel ising_model_T1(L, 1.0);
    // arma::mat data0 = ising_model_T1.run_n_MC_cycles(50000, 0, "random");
    // data0.save("datasets/data_L20_initRand_MC50000_T1.txt", arma::raw_ascii);
    // arma::mat data2 = ising_model_T1.run_n_MC_cycles(50000, 0, "ordered");
    // data2.save("datasets/data_L20_initOrdered_MC50000_T1.txt", arma::raw_ascii);

    // IsingModel ising_model_T24(L, 2.4);
    // arma::mat data1 = ising_model_T24.run_n_MC_cycles(50000, 0, "random");
    // data1.save("datasets/data_L20_initRand_MC50000_T2.4.txt", arma::raw_ascii);
    // arma::mat data3 = ising_model_T24.run_n_MC_cycles(50000, 0, "ordered");
    // data3.save("datasets/data_L20_initOrdered_MC50000_T2.4.txt", arma::raw_ascii);

    // ----------
    // Problem 6:
    // ----------

    // IsingModel ising_model(L, T=1.0);
    // arma::mat data = ising_model.run_n_MC_cycles(100000, 10000, "ordered");
    // data.save("datasets/hist_data_L20_initOrdered_MC100000_T1.txt", arma::raw_ascii);

    // IsingModel ising_model1(20, 2.4);
    // arma::mat data1 = ising_model1.run_n_MC_cycles(100000, 10000, "random");
    // data1.save("datasets/hist_data_L20_initRand_MC100000_T2.4.txt", arma::raw_ascii);

    // IsingModel ising_model2(20, 10.0);
    // arma::mat data2 = ising_model2.run_n_MC_cycles(100000, 10000, "random");
    // data2.save("datasets/hist_data_L20_initRand_MC100000_T10.txt", arma::raw_ascii);

    // ----------
    // Problem 7:
    // ----------

    // arma::vec temperatures = arma::linspace(1, 10, 1);
    // run_ensemble(temperatures, 50000, 0, 20);
    // IsingModel ising_model(20, 2.0);
    // arma::mat data = ising_model.run_n_MC_cycles(50000, 0, "random");

    // ----------
    // Problem 7:
    // ----------
    
    std::cout << "Search in critical range" << std::endl;
    arma::vec temperatures = arma::linspace(2.1, 2.4, 40);
    // std::cout << temperatures << std::endl;
    int n_cycles = 1000000;
    int n0 = 20000;
    
    std::cout << "L = " << 40 << std::endl;
    run_ensemble(temperatures, n_cycles, n0, 40);
    std::cout << "L = " << 60 << std::endl;
    run_ensemble(temperatures, n_cycles, n0, 60);
    std::cout << "L = " << 80 << std::endl;
    run_ensemble(temperatures, n_cycles, n0, 80);
    std::cout << "L = " << 100 << std::endl;
    run_ensemble(temperatures, n_cycles, n0, 100);

    return 0;
}

void run_ensemble(arma::vec temperatures, int n_cycles, int n0, int L)
{
    int K = temperatures.size(); //num. temperatures to run

    // here we do the parallelisation
    #pragma omp parallel for
    for(int k=0; k < K; k++)
    {
        double T = temperatures(k);

        IsingModel ising_model(L, T);

        arma::mat data = ising_model.run_n_MC_cycles(n_cycles, n0, "random");
        std::ostringstream oss;
        oss << "_L" << L << "initRand_MC" << n_cycles << "T" << T;
        std::string info = oss.str();
        data.save("datasets/data_parallel" + info + ".txt", arma::raw_ascii);
        
    }  
}

void compute_2x2_numeric_solution()
{
    arma::mat values = arma::mat(3, 4, arma::fill::zeros);
    int L = 2;
    double T = 1.0;
    int N = L*L;
    IsingModel ising_model(L, T);
    int n = 500;

    n = 5000;
    std::cout << " n = " << n <<std::endl;
    arma::mat data = ising_model.run_n_MC_cycles(n, 0, "random");
    std::cout << "Numeric expected energy = " << mean(data.col(0))/N << std::endl;
    std::cout << "Numeric expected heat capacity = " << 1.0/ N * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0))) << std::endl;
    std::cout << "Numeric expected abs. magnetisation = " << mean(abs(data.col(1)))/N << std::endl;
    std::cout << "Numeric expected susceptibility = " << 1.0 / N * (mean( data.col(1) % data.col(1)) - mean(abs(data.col(1)))*mean(abs(data.col(1)))) << std::endl;
    values(0,0) =  mean(data.col(0))/N;
    values(0,1) =  1.0/ N * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0)));
    values(0,2) = mean(abs(data.col(1)))/N;
    values(0,3) = 1.0 / N * (mean( data.col(1) % data.col(1)) - mean(abs(data.col(1)))*mean(abs(data.col(1))));

    n = 50000;
    std::cout << " n = " << n <<std::endl;
    data = ising_model.run_n_MC_cycles(n, 0, "random");
    std::cout << "Numeric expected energy = " << mean(data.col(0))/N << std::endl;
    std::cout << "Numeric expected heat capacity = " << 1.0/ N * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0))) << std::endl;
    std::cout << "Numeric expected abs. magnetisation = " << mean(abs(data.col(1)))/N << std::endl;
    std::cout << "Numeric expected susceptibility = " << 1.0 / N * (mean( data.col(1) % data.col(1)) - mean(abs(data.col(1)))*mean(abs(data.col(1)))) << std::endl;
    values(1,0) = mean(data.col(0))/N;
    values(1,1) = 1.0/ N * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0)));
    values(1,2) = mean(abs(data.col(1)))/N;
    values(1,3) = 1.0 / N * (mean( data.col(1) % data.col(1)) - mean(abs(data.col(1)))*mean(abs(data.col(1))));

    n = 500000;
    std::cout << " n = " << n <<std::endl;
    data = ising_model.run_n_MC_cycles(n, 0, "random");
    std::cout << "Numeric expected energy = " << mean(data.col(0))/N << std::endl;
    std::cout << "Numeric expected heat capacity = " << 1.0/ N * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0))) << std::endl;
    std::cout << "Numeric expected abs. magnetisation = " << mean(abs(data.col(1)))/N << std::endl;
    std::cout << "Numeric expected susceptibility = " << 1.0 / N * (mean( data.col(1) % data.col(1)) - mean(abs(data.col(1)))*mean(abs(data.col(1)))) << std::endl;
    values(2,0) = mean(data.col(0))/N;
    values(2,1) = 1.0/ N * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0)));
    values(2,2) = mean(abs(data.col(1)))/N;
    values(2,3) = 1.0 / N * (mean( data.col(1) % data.col(1)) - mean(abs(data.col(1)))*mean(abs(data.col(1))));

    values.save("datasets/numeric_sol.txt", arma::raw_ascii);
}