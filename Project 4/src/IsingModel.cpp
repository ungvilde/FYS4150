# include "IsingModel.hpp"

IsingModel::IsingModel(int lattice_size, double temperature)
{
    L = lattice_size;
    N = L*L;
    T = temperature;
    beta = 1.0 / (T*kb);
}

// for initialising the lattice
void IsingModel::initialize_lattice(std::string initialisation)
{
    if(initialisation == "random")
    {
        lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0, 1));
        lattice = 2.0 * lattice - 1;
    } else if(initialisation == "ordered")
    {
        lattice = arma::mat(L, L, arma::fill::ones);
    } else
    {
        throw std::invalid_argument( "Invalid initialisation." );
    }
}

// for running through a single monte carlo cycle
void IsingModel::run_MC_cycle()
{
    for(int k=0; k<N; k++)
    {
        int i = arma::randi(arma::distr_param(0, L-1)); // choose random indeces
        int j = arma::randi(arma::distr_param(0, L-1));
        
        // compute energy at chosen spin
        double energy_noflip = - J * lattice(i, j) * (
            lattice(i, ((j - 1) % L + L) % L) +
            lattice(i, ((j + 1) % L + L) % L) +
            lattice(((i - 1) % L + L) % L, j) +
            lattice(((i + 1) % L + L) % L, j)
        );
    
        // compute energy if spin is flipped
        double energy_flip = -1 * energy_noflip;
        double u = arma::randu(arma::distr_param(0, 1));
        double dE = energy_flip - energy_noflip;

        // this is the Metropolis algorithm
        if(dE <= 0)
        {
            lattice(i, j) *= -1; // do the flip
            energy += dE;
            magnetisation += 2*lattice(i, j);
        }
        else if(u <= exp(-beta * dE)) // could hardcode this exp-value, there are only two cases to check
        {
            lattice(i, j) *= -1; // do the flip
            energy += dE;
            magnetisation += 2*lattice(i, j);
        }
    }
}

// for running through n monte carlo cycles, storing the E, M, etc for each cycle
arma::mat IsingModel::run_n_MC_cycles(int n_cycles, int n0, std::string initialisation)
{
    int k = 0;
    int i = 0;
    initialize_lattice(initialisation);

    // compute energy and magnetisation of initial lattice
    // these variables are updated whenever we do a flip
    energy = compute_energy();
    magnetisation = compute_magnetisation();

    // for storing the values    
    arma::mat data = arma::mat(n_cycles - n0, 2);

    while(k < n_cycles)
    {
        run_MC_cycle();
        
        if(k >= n0) //burn-in
        {
            // here we store the values
            data(i, 0) = energy;
            data(i, 1) = magnetisation;

            i += 1;
        }

        k += 1;
    }
    
    return data;
}

// for computing the energy of a lattice
double IsingModel::compute_energy()
{
    energy = 0;
    for(int i = 0; i<L; i++)
    {
        for(int j = 0; j < L; j++)
        {
            energy += -J * lattice(i, j)*(
                lattice(i, ((j - 1) % L + L) % L) +
                lattice(i, ((j + 1) % L + L) % L) +
                lattice(((i - 1) % L + L) % L, j) +
                lattice(((i + 1) % L + L) % L, j)
            );
        }
    }

    energy = energy / 2.0; // to account for double counting 
    
    return energy;
}

// for computing magnetisation of a lattice
double IsingModel::compute_magnetisation()
{
    double magnetization = 0;
    for(int i = 0; i < L; i++)
    {
        for(int j = 0; j < L; j++)
        {
            magnetization += lattice(i, j);
        }
    }

    return magnetization;
}