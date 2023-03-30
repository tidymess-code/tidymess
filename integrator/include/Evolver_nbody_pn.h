#include "Evolver_nbody_base.h"

#ifndef __Evolver_nbody_pn_h
#define __Evolver_nbody_pn_h

class Evolver_nbody_pn : public Evolver_nbody_base {
    public:

    // Evolve function
    void leapfrog_step(vector<Body> &bodies, double dt, Force *force);
    void mclachlan_step(vector<Body> &bodies, double dt, Force *force);
 
    void initialize(vector<Body> &bodies, Force *force);
    void evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep);
    void evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep);
};

#endif


