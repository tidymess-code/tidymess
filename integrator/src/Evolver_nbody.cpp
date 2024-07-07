#include "Evolver_nbody.h"

// Evolve function
void Evolver_nbody::leapfrog_step(vector<Body> &bodies, double dt, Force *force) {
    double dth = 0.5*dt;
    orbit.drift_r(bodies, dth);
    force->update_acceleration(bodies);
    orbit.kick_v(bodies, dt);
    orbit.drift_r(bodies, dth);
}
void Evolver_nbody::mclachlan_step(vector<Body> &bodies, double dt, Force *force) {
    leapfrog_step(bodies, w1*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w3*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w1*dt, force);
}

void Evolver_nbody::initialize(vector<Body> &bodies, Force *force) {
    force->update_acceleration(bodies);
}
void Evolver_nbody::evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
    bool toContinue = true;

    num_integration_step = 0;

    toContinue = stopping_condition_time(t, t_end, timestep->dt_sgn);
        
    while(toContinue) {
        dt_step = timestep->get_timestep(bodies);

        if(dt_step > 0 && t+dt_step > t_end) dt_step = t_end-t;
        else if(dt_step < 0 && t+dt_step < t_end) dt_step = t_end-t;    

        //----------------------------------------------

        mclachlan_step(bodies, dt_step, force);

        //----------------------------------------------

        t += dt_step;
        num_integration_step++;

        //----------------------------------------------

        toContinue = stopping_condition_time(t, t_end, timestep->dt_sgn);
    }
}
void Evolver_nbody::evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
    num_integration_step = 0;
        
    for(int i=0; i<N_step; i++) {
        dt_step = timestep->get_timestep(bodies);

        //----------------------------------------------

        mclachlan_step(bodies, dt_step, force);

        //----------------------------------------------

        t += dt_step;
        num_integration_step++;
    }
}


