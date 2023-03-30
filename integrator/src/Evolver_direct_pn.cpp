#include "Evolver_direct_pn.h"

// Evolve function

void Evolver_direct_pn::leapfrog_step(vector<Body> &bodies, double dt, Force *force) {
    double dth = 0.5*dt;

    orbit.drift_r(bodies, dth);

    orbit.kick_vv(bodies, dth);
    spin.kick_K(bodies, dth);
    shape.kick_direct_J(bodies, dth);

    shape.calculate_J_inv_and_w_with_K(bodies);

    force->update_acceleration_and_torque_pn_and_deformation_tensor(bodies, true, true);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    force->update_deformation_tensor(bodies, true);

    orbit.kick_v(bodies, dt);
    spin.kick_L(bodies, dt);
    shape.kick_direct(bodies, dt);

    orbit.drift_r(bodies, dth);      

    shape.calculate_I_inv_and_w(bodies);

    force->update_acceleration_and_torque_pn_and_deformation_tensor(bodies, false, false);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    force->update_deformation_tensor(bodies, false);

    orbit.kick_vv(bodies, dth);
    spin.kick_K(bodies, dth);
    shape.kick_direct_J(bodies, dth);
}

void Evolver_direct_pn::mclachlan_step(vector<Body> &bodies, double dt, Force *force) {   
    leapfrog_step(bodies, w1*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w3*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w1*dt, force);
}

void Evolver_direct_pn::initialize(vector<Body> &bodies, Force *force) {
    shape.calculate_I_inv_and_w(bodies);

    force->update_acceleration_and_torque_pn_and_deformation_tensor(bodies, false, false);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    force->update_deformation_tensor(bodies, false);
}
void Evolver_direct_pn::evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
    bool toContinue = true;

    num_integration_step = 0;
        
    while(toContinue) {
        reset_aux(bodies);

        //----------------------------------------------

        dt_step = timestep->get_timestep(bodies);

        //----------------------------------------------

        mclachlan_step(bodies, dt_step, force);

        //----------------------------------------------

        t += dt_step;
        num_integration_step++;

        //----------------------------------------------

        toContinue = stopping_condition_time(t, t_end, timestep->dt_sgn);
    }    
}
void Evolver_direct_pn::evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
    num_integration_step = 0;
        
    for(int i=0; i<N_step; i++) {
        reset_aux(bodies);

        //----------------------------------------------

        dt_step = timestep->get_timestep(bodies);

        //----------------------------------------------

        mclachlan_step(bodies, dt_step, force);

        //----------------------------------------------

        t += dt_step;
        num_integration_step++;
    }    
}



