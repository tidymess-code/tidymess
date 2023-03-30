#include "Evolver_shape_base.h"

#ifndef __Evolver_linear_pn_mb_h
#define __Evolver_linear_pn_mb_h

class Evolver_linear_pn_mb : public Evolver_shape_base {
    public:
    
    // Evolve function
    void leapfrog_step(vector<Body> &bodies, double dt, Force *force);
    void sync_step(vector<Body> &bodies, Force *force, double dt);
    
    void mclachlan_step(vector<Body> &bodies, double dt, Force *force);
    
    void initialize(vector<Body> &bodies, Force *force);
    void evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep);
    void evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep);    
};

#endif


