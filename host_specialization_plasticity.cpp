// C++ source code for a Gillespie simulation
// for a 2 species model with virulence evolution due to genetics,
// horizontal transmission and plasticity
// Written by Bram Kuijper, October/November 2018
//
// Used in the paper
// How do bacterial pathogens emerge in novel hosts? 
// Camille Bonneaud, Lucy Weinert & Bram Kuijper


#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "auxiliary.h"

// whether we want to output the whole distribution of all trait loci
// We do so if we are making plots of the branching events
//#define DISTRIBUTION

gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_r; // gnu scientific rng 

using namespace std;

const int Npop = 20000;

double lambda = 0.0; // immigration of susceptible hosts
double d[2] = {0.0,0.0}; // susceptible host death rate
double gamma_par[2] = {0.0,0.0}; // susceptible host clearance
double p = 0.0; // proportion one host versus second
double rho = 0.0; // strenght of trade-off of virulence in host vs another
double tau2 = 0.0;
double phi[2][2] = {{0.0,0.0},{0.0,0.0}};

double mu_a = 0; // mutation rate of the elevation
double mu_b = 0; // mutation rate of reaction norm slope
double mu_h = 0; // mutation rate of reaction norm slope
double init_v = 0;
double sdmu_a = 0; // mutational variance elevation
double sdmu_b = 0; // mutational variance plasticity
double sdmu_h = 0; // mutational variance plasticity

int init_s[2] = {0,0};
int init_i[2] = {0,0};
double kappa = 0; // strength with which transmissibility is reduced in envt 2

double U = 2; // increase in environmental signal in new host
double sigma_e = 2; // random environmental noise in virulence



// stats
int Ns[2] = {0,0}; // number of susceptible of species 1
int Ni[2] = {0,0}; // number of infected of species 2

long int t_max = 0; // maximum number of timesteps in this simulation
long int t = 0; 

double seed = 0; // the random seed

int sample_transmit = 0;
int sample_burst = 0;

long int skip = 0;
long int skip_distribution = 0;

// enum with the different types
enum State { susceptible, infected };
enum Species { native_host, new_host };
enum TraitType { genetic, plasticity, horizontal };

// the individual struct
struct Individual
{
    double a; // genetic adaptation of virulence
    double b; // reaction norm of virulence
    double h; // horizontal transmission rate
    double envt; // current environment experienced by pathogen
    double vphen; // virulence phenotype

    State host_type;

    Species host_species;
};

// the numbers of the individuals
// that are sampled to calculate
// the transmission rate
// we need these guys in case we 
// are going to do some actual
// new nfections

typedef Individual Population[Npop];
Population Susceptible[2];
Population Infected[2];

string filename("sim_sir_evo");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

#ifdef DISTRIBUTION
string filename_dist(filename_new + "_dist");
ofstream DataFileDist(filename_dist.c_str());
#endif 

void init_arguments(int argc, char *argv[])
{
    lambda = atof(argv[1]);
    d[native_host] = atof(argv[2]);
    d[new_host] = atof(argv[3]);
    gamma_par[native_host] = atof(argv[4]);
    gamma_par[new_host] = atof(argv[5]);
    p = atof(argv[6]);
    rho = atof(argv[7]);
    
    phi[native_host][native_host] = atof(argv[8]);
    phi[native_host][new_host] = atof(argv[9]);
    phi[new_host][native_host] = atof(argv[10]);
    phi[new_host][new_host] = atof(argv[11]);

    mu_a = atof(argv[12]);
    mu_b = atof(argv[13]);
    mu_h = atof(argv[14]);
    sdmu_a = atof(argv[15]);
    sdmu_b = atof(argv[16]);
    sdmu_h = atof(argv[17]);
    init_v= atof(argv[18]);

    init_s[native_host] = atoi(argv[19]);
    init_s[new_host] = atoi(argv[20]);
    init_i[native_host] = atoi(argv[21]);
    init_i[new_host] = atoi(argv[22]);

    t_max = atoi(argv[23]);
    
    kappa = atof(argv[24]);
    U = atof(argv[25]);
    sigma_e = atof(argv[26]);
    tau2 = atof(argv[27]);

    skip = t_max / 5000;

    skip_distribution = t_max / 5000;
}

void write_parameters()
{
	DataFile << endl << endl
        << "type;sir_evolutionary_gandon" << endl
        << "seed;" << setprecision(10) << seed << endl
        << "lambda;" <<  lambda << endl
        << "d1;" <<  d[native_host] << endl
        << "d2;" <<  d[new_host] << endl
        << "gamma_par1;" <<  gamma_par[native_host] << endl
        << "gamma_par2;" <<  gamma_par[new_host] << endl
        << "rho;" <<  rho << endl
        << "p;" <<  p << endl
        << "tau2;" <<  tau2 << endl
        << "phi11;" <<  phi[native_host][native_host] << endl
        << "phi12;" <<  phi[native_host][new_host] << endl
        << "phi21;" <<  phi[new_host][native_host] << endl
        << "phi22;" <<  phi[new_host][new_host] << endl
        << "mu_a;" <<  mu_a << endl
        << "mu_b;" <<  mu_b << endl
        << "mu_h;" <<  mu_h << endl
        << "sdmu_a;" << sdmu_a << endl
        << "sdmu_b;" << sdmu_b << endl
        << "sdmu_h;" << sdmu_h << endl
        << "init_v;" << init_v << endl
        << "init_s1;" << init_s[0] << endl
        << "init_s2;" << init_s[1] << endl
        << "init_i1;" << init_i[0] << endl
        << "init_i2;" << init_i[1] << endl
        << "kappa;" << kappa << endl
        << "U;" << U << endl
        << "sigma_e;" << sigma_e << endl;
}

// calculate probability of transmission
double tau(double const v, 
        Species const host_species)
{
    double eps = host_species == new_host ? (1.0 - rho) * v : v;

    return(eps / (tau2 + eps));
}

// initialize all the phenotypes
void init_population()
{
    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();
    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng_r = gsl_rng_alloc(T);
    gsl_rng_set(rng_r, seed);

    Ns[native_host] = init_s[native_host];
    Ns[new_host] = init_s[new_host];
    Ni[native_host] = init_i[native_host];
    Ni[new_host] = init_i[new_host];

    for (int i = 0; i < Ns[native_host]; ++i)
    {
        Susceptible[native_host][i].a = 0;
        Susceptible[native_host][i].b = 0;
        Susceptible[native_host][i].h = 0;
        Susceptible[native_host][i].envt = 0;
        Susceptible[native_host][i].vphen = 0;
        Susceptible[native_host][i].host_type = susceptible;
    }
    
    for (int i = 0; i < Ns[new_host]; ++i)
    {
        Susceptible[new_host][i].a = 0;
        Susceptible[new_host][i].b = 0;
        Susceptible[new_host][i].h = 0;
        Susceptible[new_host][i].envt = 0;
        Susceptible[new_host][i].vphen = 0;
        Susceptible[new_host][i].host_type = susceptible;
    }
    
    for (int i = 0; i < Ni[native_host]; ++i)
    {
        Infected[native_host][i].a = init_v;
        Infected[native_host][i].b = 0;
        Infected[native_host][i].h = 0;
        Infected[native_host][i].envt = gsl_ran_gaussian(rng_r, sigma_e);
        Infected[native_host][i].vphen = init_v;
        Infected[native_host][i].host_type = infected;
    }
    
    for (int i = 0; i < Ni[new_host]; ++i)
    {
        Infected[new_host][i].a = init_v;
        Infected[new_host][i].b = 0;
        Infected[new_host][i].h = 0;
        Infected[new_host][i].envt = U + gsl_ran_gaussian(rng_r, sigma_e);
        Infected[new_host][i].vphen = init_v;
        Infected[new_host][i].host_type = infected;
    }
}

// create a newborn individual
// this individual is always born as a susceptible
void newborn(Species host_species)
{
    // if already at carrying capacity 
    // forget producing a newborn
    if (Ns[native_host] + Ni[native_host] 
            + Ns[new_host] + Ni[new_host] < Npop-4)
    {
        // set all trait values to zero
        Individual Kid;
        Kid.a = 0;
        Kid.b = 0;
        Kid.h = 0;
        Kid.envt = 0;
        Kid.vphen = 0;
        Kid.host_type = susceptible;

        // if parent was type 1, host is type 1
        Susceptible[host_species][Ns[host_species]++] = Kid;
    }
}

// kill random individual
void death(State const host_state, Species const host_species)
{
    if (host_state == susceptible)
    {
        int n = Ns[host_species];

        assert(n > 0);
        assert(n < Npop);

        // delete susceptible
        Susceptible[host_species][gsl_rng_uniform_int(rng_r, n)] 
            = Susceptible[host_species][n - 1];

        --Ns[host_species];
    }
    else
    {
        int n = Ni[host_species];

        assert(n > 0);
        assert(n < Npop);

        // delete infected
        int rand_inf = gsl_rng_uniform_int(rng_r, n);
        Infected[host_species][rand_inf] = Infected[host_species][n - 1];

        --Ni[host_species];
    }
}

// infection of susceptible
void infect(
        int const host_id, // index number of original infected who infects susceptible
        Species const host_species_from, // the host species who infects
        Species const host_species_to) // the host species who infects)
{
    assert(Ns[host_species_to] > 0);
    assert(Ni[host_species_from] > 0);
    assert(host_id >= 0);
    assert(host_id < Ni[host_species_from]);

    int to_be_infected = gsl_rng_uniform_int(rng_r, Ns[host_species_to]);

    // remove the susceptible individual
    Susceptible[host_species_to][to_be_infected] = 
        Susceptible[host_species_to][Ns[host_species_to] - 1];

    --Ns[host_species_to];
    
    // add the newly infected to the infected
    // individuals
    Individual newInfected;

    // do some error checking
    assert(Infected[host_species_from][host_id].host_type == infected);

    // assign values to the infected individual from the host
    newInfected.a = Infected[host_species_from][host_id].a;
    newInfected.b = Infected[host_species_from][host_id].b;
    newInfected.h = Infected[host_species_from][host_id].h;

    // determine the environment (either 1 or 2)
    // background variation is given by sigma_envt around 0
    newInfected.envt = gsl_ran_gaussian(rng_r, sigma_e);

    // steep change in the environment when it is in the new host
    if (host_species_to == new_host)
    {
        newInfected.envt += U;
    }

    // make the virulence phenotype
    newInfected.vphen = newInfected.a + newInfected.b * newInfected.envt;

    if (newInfected.vphen < 0)
    {
        newInfected.vphen = 0;
    }

    // set the individual's state to infected
    newInfected.host_type = infected;
    newInfected.host_species = host_species_to;

    // add newly infected to stack
    Infected[host_species_to][Ni[host_species_to]++] = newInfected;

    assert(Ni[native_host] + Ns[native_host] + 
            Ni[new_host] + Ns[new_host] <= Npop);
}

// clearance of infection
void clearance(Species const host_species)
{
    assert(Ni[host_species]>0);

    Infected[host_species][gsl_rng_uniform_int(rng_r, Ni[host_species])] = 
        Infected[host_species][Ni[host_species]-1];

    --Ni[host_species];

    Individual newSusceptible;
    newSusceptible.host_type = susceptible;
    Susceptible[host_species][++Ns[host_species]] = newSusceptible;
    
    assert(Ni[native_host] + Ns[native_host] + 
            Ni[new_host] + Ns[new_host] <= Npop);
}

// death of an infected 'host_id', due to
// virulence of the parasite
void death_due_to_virulence(int const host_id, Species const host_species)
{
    assert(host_id >= 0 && host_id < Ni[host_species]);

    Infected[host_species][host_id] = Infected[host_species][Ni[host_species]-1];
    Ni[host_species]--;
}

// horizontally transmit a and b within a host species
void horizontal_transmission(int const host_id, Species const host_species)
{
    // if there is only one infected host, don't bother
    if (Ni[host_species] <= 1)
    {
        return;
    }

    int other_individual;

    // randomly sample another host
    do 
    {
        other_individual = gsl_rng_uniform_int(rng_r, Ni[host_species]);
    }
    while(other_individual == host_id);

    // horizontal transmission of a or b
    if (gsl_rng_uniform(rng_r) < 0.5)
    {
        Infected[host_species][host_id].a = Infected[host_species][other_individual].a;
    }
    else
    {
        Infected[host_species][host_id].b = Infected[host_species][other_individual].b;
    }
}


// mutation of the parasite
void mutate(
        TraitType const trait_type, 
        Species const host_species)
{
    assert(Ni[host_species]>0);

    int parasite = gsl_rng_uniform_int(rng_r, Ni[host_species]);

    assert(Infected[host_species][parasite].host_type == infected);

    double dev;

    // mutate the elevation
    if (trait_type == genetic)
    {
        dev = sdmu_a * gsl_rng_uniform(rng_r);

        Infected[host_species][parasite].a += 
            gsl_rng_uniform(rng_r) < 0.5 ? -dev : dev;
    }
    else if (trait_type == plasticity)
    {
        dev = sdmu_b * gsl_rng_uniform(rng_r);

        Infected[host_species][parasite].b += 
            gsl_rng_uniform(rng_r) < 0.5 ? -dev : dev;
    }
    else
    {
        dev = sdmu_h * gsl_rng_uniform(rng_r);

        Infected[host_species][parasite].h += 
            gsl_rng_uniform(rng_r) < 0.5 ? -dev : dev;

        // check boundaries
        if (Infected[host_species][parasite].h < 0)
        {
            Infected[host_species][parasite].h = 0.0;
        }
    }
}

void write_data()
{
    double meanv_overall = 0;
    double ssv_overall = 0;

    double meanv[2] = {0,0};
    double ssv[2] = {0,0};

    double meana_overall = 0;
    double ssa_overall = 0;

    double meana[2] = {0,0};
    double ssa[2] = {0,0};
    
    double meanb_overall = 0;
    double ssb_overall = 0;

    double meanb[2] = {0,0};
    double ssb[2] = {0,0};
    
    double meanh_overall = 0;
    double ssh_overall = 0;

    double meanh[2] = {0,0};
    double ssh[2] = {0,0};


    double varv_overall = 0;
    double vara_overall = 0;
    double varb_overall = 0;
    double varh_overall = 0;

    double varv[2] = {0,0};
    double vara[2] = {0,0};
    double varb[2] = {0,0};
    double varh[2] = {0,0};


    double v,b,a,h;

    for (int i = 0; i < Ni[native_host]; ++i)
    {
        // stats for v
        v = Infected[native_host][i].vphen;
        
        meanv_overall += v;
        ssv_overall += v*v;

        assert(v >= 0);
        meanv[native_host] += v;
        ssv[native_host]+= v * v;
       
        // stats for a 
        a = Infected[native_host][i].a;

        meana_overall += a;
        ssa_overall += a*a;

        meana[native_host] += a;
        ssa[native_host] += a*a;

        // stats for b
        b = Infected[native_host][i].b;
        meanb_overall += b;
        ssb_overall += b*b;
        
        meanb[native_host] += b;
        ssb[native_host]  += b*b;

        // stats for h
        h = Infected[native_host][i].h;
        meanh_overall += h;
        ssh_overall += h*h;
        
        meanh[native_host] += h;
        ssh[native_host]  += h*h;
    }
    
    for (int i = 0; i < Ni[new_host]; ++i)
    {
        // stats for v
        v = Infected[new_host][i].vphen;
        
        meanv_overall += v;
        ssv_overall += v*v;

        assert(v >= 0);
        meanv[new_host] += v;
        ssv[new_host]+= v * v;
       
        // stats for a 
        a = Infected[new_host][i].a;

        meana_overall += a;
        ssa_overall += a*a;

        meana[new_host] += a;
        ssa[new_host] += a*a;

        // stats for b
        b = Infected[new_host][i].b;
        meanb_overall += b;
        ssb_overall += b*b;
        
        meanb[new_host] += b;
        ssb[new_host]  += b*b;
        
        // stats for h
        h = Infected[new_host][i].h;
        meanh_overall += h;
        ssh_overall += h*h;
        
        meanh[new_host] += h;
        ssh[new_host]  += h*h;
    }

    int Ni_overall = 0;

    if (Ni[native_host] > 0)
    {
        Ni_overall += Ni[native_host];

        meanv[native_host]/= Ni[native_host];
        meana[native_host] /= Ni[native_host];
        meanb[native_host] /= Ni[native_host];
        meanh[native_host] /= Ni[native_host];

        varv[native_host] = ssv[native_host]/Ni[native_host] 
            - meanv[native_host]*meanv[native_host];

        vara[native_host] = ssa[native_host]/Ni[native_host] 
            - meana[native_host]*meana[native_host];

        varb[native_host] = ssb[native_host]/Ni[native_host] 
            - meanb[native_host]*meanb[native_host];
        
        varh[native_host] = ssh[native_host]/Ni[native_host] 
            - meanh[native_host]*meanh[native_host];
    }
    
    if (Ni[new_host] > 0)
    {
        Ni_overall += Ni[new_host];

        meanv[new_host] /= Ni[new_host];
        meana[new_host] /= Ni[new_host];
        meanb[new_host] /= Ni[new_host];
        meanh[new_host] /= Ni[new_host];

        varv[new_host] = ssv[new_host]/Ni[new_host] 
            - meanv[new_host]*meanv[new_host];

        vara[new_host] = ssa[new_host]/Ni[new_host] 
            - meana[new_host]*meana[new_host];

        varb[new_host] = ssb[new_host]/Ni[new_host] 
            - meanb[new_host]*meanb[new_host];
        
        varh[new_host] = ssh[new_host]/Ni[new_host] 
            - meanh[new_host]*meanh[new_host];
    }

    if (Ni_overall > 0)
    {
        meanv_overall /= Ni_overall;
        meana_overall /= Ni_overall;
        meanb_overall /= Ni_overall;
        meanh_overall /= Ni_overall;
        varv_overall = ssv_overall / Ni_overall - meanv_overall * meanv_overall;
        vara_overall = ssa_overall / Ni_overall - meana_overall * meana_overall;
        varb_overall = ssb_overall / Ni_overall - meanb_overall * meanb_overall;
        varh_overall = ssh_overall / Ni_overall - meanh_overall * meanh_overall;
    }


    DataFile << t << ";" 
        << meana_overall << ";"
        << meanb_overall << ";"
        << meanv_overall << ";"
        << meanh_overall << ";"
        << vara_overall << ";"
        << varb_overall << ";"
        << varv_overall << ";"
        << varh_overall << ";"
        << Ns[native_host] << ";" 
        << Ni[native_host] << ";" 
        << meana[native_host] << ";" 
        << meanb[native_host] << ";" 
        << meanv[native_host] << ";" 
        << meanh[native_host] << ";" 
        << vara[native_host] << ";" 
        << varb[native_host] << ";" 
        << varv[native_host] << ";" 
        << varh[native_host] << ";" 
        << Ns[new_host] << ";" 
        << Ni[new_host] << ";" 
        << meana[new_host] << ";" 
        << meanb[new_host] << ";" 
        << meanv[new_host] << ";" 
        << meanh[new_host] << ";" 
        << vara[new_host] << ";" 
        << varb[new_host] << ";" 
        << varv[new_host] << ";" 
        << varh[new_host] << ";" 
        << endl;
}

void write_distribution_data(long int const t)
{
#ifdef DISTRIBUTION
    for (int i = 0; i < Ni[native_host]; ++i)
    {
        DataFileDist << t << ";"
            << Infected[native_host][i].a << ";"
            << Infected[native_host][i].b << ";"
            << Infected[native_host][i].h << ";"
            << Infected[native_host][i].vphen << ";" << endl;
    }
    
    for (int i = 0; i < Ni[new_host]; ++i)
    {
        DataFileDist << t << ";"
            << Infected[new_host][i].a << ";"
            << Infected[new_host][i].b << ";"
            << Infected[new_host][i].h << ";"
            << Infected[new_host][i].vphen << ";" << endl;
    }
#endif
}

void write_data_headers()
{
    DataFile << 
        "t;meana;meanb;meanv;meanh;vara;varb;varv;varh;" << 
        "Ns1;Ni1;meana1;meanb1;meanv1;meanh1;vara1;varb1;varv1;varh1;" << 
        "Ns2;Ni2;meana2;meanb2;meanv2;meanh2;vara2;varb2;varv2;varh2;" << endl;
}

void write_distribution_data_headers()
{
#ifdef DISTRIBUTION
    DataFileDist << "t;a;b;h;v;" << endl;
#endif
}


int main(int argc, char ** argv)
{
	init_arguments(argc, argv);
	write_data_headers();
	init_population();

#ifdef DISTRIBUTION
    write_distribution_data_headers();
#endif

    double sumprob;

    int mod, id, nprobs;

	for (t = 0; t < t_max; ++t)
	{
        // total number of probabilities:
        // 1. immigration species 1
        // 2. immigration species 2
        // 3. death of susceptible species 1
        // 4. death of susceptible species 2
        // 5. clearance of infected species 1
        // 6. clearance of infected species 2
        // 7. death of infected species 1
        // 8. death of infected species 2
        // 9. mutation of elevation infected species 1
        // 10. mutation of elevation infected species 2
        // 11. mutation of slope infected species 1
        // 12. mutation of slope infected species 2
        // 13. mutation of horizontal transmission species 1
        // 14. mutation of horizontal transmission species 2
        nprobs = 14 + 
            4*Ni[native_host] 
            + 4*Ni[new_host];

        double probs[nprobs];

        assert(Ns[native_host] + Ni[native_host] + 
                Ns[new_host] + Ni[new_host] < Npop);

        assert(Ns[native_host] >= 0);
        assert(Ns[new_host] >= 0);
        assert(Ni[native_host] >= 0);
        assert(Ni[new_host] >= 0);
        
        sumprob = probs[0] = lambda * (1.0 - p); // immigration of susceptible hosts 1
        sumprob = probs[1] = sumprob + lambda * p; // immigration of susceptible hosts 2

        sumprob = probs[2] = sumprob + d[native_host] * Ns[native_host]; // death of susceptible 1
        sumprob = probs[3] = sumprob + d[new_host] * Ns[new_host]; // death of susceptible 2

        sumprob = probs[4] = sumprob + 
            gamma_par[native_host] * Ni[native_host]; // clearance of random infected 1

        sumprob = probs[5] = sumprob + 
            gamma_par[new_host] * Ni[new_host]; // clearance of random infected 2

        sumprob = probs[6] = sumprob + d[native_host] * Ni[native_host]; // normal death of infected 1
        sumprob = probs[7] = sumprob + d[new_host] * Ni[new_host]; // normal death of infected 2

        sumprob = probs[8] = sumprob + mu_a * Ni[native_host]; // mutation of a 1
        sumprob = probs[9] = sumprob + mu_a * Ni[new_host]; // mutation of a 2

        sumprob = probs[10] = sumprob + mu_b * Ni[native_host]; // mutation of b 1
        sumprob = probs[11] = sumprob + mu_b * Ni[new_host]; // mutation of b 2
        
        sumprob = probs[12] = sumprob + mu_h * Ni[native_host]; // mutation of h 1
        sumprob = probs[13] = sumprob + mu_h * Ni[new_host]; // mutation of h 2

        int prob_counter = 14;

        // now let's deal with infection, 
        // virulence related deaths 
        // and mutation of species 1
        //
        // these are functions of the individual virulence levels
        // in the population
        for (int i = 0; i < Ni[native_host]; ++i)
        {
            // infection of susceptible species 1
            sumprob = probs[prob_counter + i*4] = sumprob
                + phi[native_host][native_host]
                    * tau(Infected[native_host][i].vphen, 
                            native_host) * Ns[native_host];
           
            // infection of susceptible species 2
            sumprob = probs[prob_counter + i*4 + 1] = sumprob
                + phi[native_host][new_host]
                    * tau(Infected[native_host][i].vphen, 
                            native_host) * Ns[new_host];

            // virulence related death 
            sumprob = probs[prob_counter + i*4 + 2] = sumprob
                + Infected[native_host][i].vphen;

            // horizontal transmission within a native host
            sumprob = probs[prob_counter + i*4 + 3] = sumprob
                + Infected[native_host][i].h; 
        }

        prob_counter = prob_counter + Ni[native_host]*4;
        
        // same for species 2
        for (int i = 0; i < Ni[new_host]; ++i)
        {
            // infection of susceptible species 1
            sumprob = probs[prob_counter + i*4] = sumprob
                + phi[new_host][native_host] *
                    tau(Infected[new_host][i].vphen, 
                            new_host) * Ns[native_host];
           
            // infection of susceptible species 2
            sumprob = probs[prob_counter + i*4 + 1] = sumprob
                + phi[new_host][new_host] *
                    tau(Infected[new_host][i].vphen, 
                            new_host) * Ns[new_host];

            // virulence related death 
            sumprob = probs[prob_counter + i*4 + 2] = sumprob
                + (1-rho) * Infected[new_host][i].vphen;
            
            // horizontal transmission within a new host
            sumprob = probs[prob_counter + i*4 + 3] = sumprob
                + Infected[new_host][i].h; 
        }

        assert(prob_counter + Ni[new_host] * 4 == nprobs);

        double prob = gsl_rng_uniform(rng_r) * sumprob;

        int choice;

        for (choice = 0; choice < nprobs; ++choice)
        {
            if (prob <= probs[choice])
            {
                break;
            }
        }

        // ok, we get into the domain of individual
        // values of virulence
        if (choice >= 14 && choice < 14 + Ni[native_host]*4)
        {
            // virulence concerning native host

            // find out which of the 4 events has been chosen
            // - infection susceptible species 1,
            // - infection susceptible species 2,
            // - death due to virulence occurs here
            // - horizontal transmission within a native host
            mod = (choice - 14) % 4;

            // find out the position of the individual in the stack
            // which is the current number minus the first 14 events
            // minus the nth of the three events (mod)
            // if this number is 0 we have the first individual
            // however is this number is >0 we have to divide 
            //id = choice - 14 - mod == 0 ? 0 : (choice - 14 - mod)/4;
            id = floor((choice - 14)/4);

            assert(id >= 0 && id < Ni[native_host]);

            // ok, infection of susceptible species 1
            if (mod == 0)
            {
                infect(id, native_host, native_host);
            }
            else if (mod == 1)  // ok infection of susceptible species 2
            {
                infect(id, native_host, new_host);
            }
            else if (mod == 2) // virulence related death
            {
                death_due_to_virulence(id, native_host);
            }
            else 
            {
                horizontal_transmission(id, native_host);
            }
        }
        else if (choice >= 14 + Ni[native_host]*4)
        {
            // virulence concerning new host

            // find out which of the three events has been chosen
            // - infection susceptible species 1,
            // - infection susceptible species 2,
            // - death due to virulence 
            // - horizontal transmission
            mod = (choice - Ni[native_host]*4 - 14) % 4;

            // find out the position of the individual in the stack
            // which is the current number minus the first 12 events
            // minus the nth of the three events (mod)
            // if this number is 0 we have the first individual
            // however is this number is >0 we have to divide 
            //id = choice - 12 - mod == 0 ? 0 : (choice - 12 - mod)/3;
            id = floor((choice - Ni[native_host]*4 - 14)/4);

            assert(id >= 0 && id < Ni[new_host]);

            // infection of susceptible species 1
            if (mod == 0)
            {
                infect(id, new_host, native_host);
            }
            else if (mod == 1)  // infection of susceptible species 2
            {
                infect(id, new_host, new_host);
            }
            else if (mod == 2)
            {
                death_due_to_virulence(id, new_host);
            }
            else
            {
                horizontal_transmission(id, new_host);
            }
        }
        else 
        {
            assert(choice < 14);

            switch(choice)
            {
                case 0:
                    newborn(native_host);
                    break;
                case 1:
                    newborn(new_host);
                    break;
                case 2:
                    death(susceptible, native_host);
                    break;
                case 3:
                    death(susceptible, new_host);
                    break;
                case 4:
                    clearance(native_host);
                    break;
                case 5:
                    clearance(new_host);
                    break;
                case 6:
                    death(infected, native_host);
                    break;
                case 7:
                    death(infected, new_host);
                    break;
                case 8:
                    mutate(genetic, native_host);
                    break;
                case 9:
                    mutate(genetic, new_host);
                    break;
                case 10:
                    mutate(plasticity, native_host);
                    break;
                case 11:
                    mutate(plasticity, new_host);
                    break;
                case 12:
                    mutate(horizontal, native_host);
                    break;
                case 13:
                    mutate(horizontal, new_host);
                    break;
                default:
                    cout << "something is going quite wrong here." << endl;
                    exit(1);
                    break;
            }
        }

        if (t % skip == 0)
        {
            write_data();
        }

#ifdef DISTRIBUTION
        if (t % skip_distribution == 0)
        {
            write_distribution_data(t);
        }
#endif
    }

    write_data();
	write_parameters();

    return(0);
}
