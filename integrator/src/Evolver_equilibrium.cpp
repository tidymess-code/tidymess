#include "Evolver_equilibrium.h"

// Evolve function
void Evolver_equilibrium::leapfrog_step(vector<Body> &bodies, double dt, Force *force) {
    double dth = 0.5*dt;

    orbit.drift_r(bodies, dth);

    force->update_tidal_deformation_save(bodies);
    force->update_equilibrium_tensor(bodies);

    shape.set_to_equilibrium_shape(bodies);

    shape.calculate_I_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.set_to_equilibrium_shape(bodies);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);
        
        shape.set_to_equilibrium_shape(bodies);
    }

    force->update_torque_load(bodies, false);

    spin.kick_K(bodies, dth);

    shape.calculate_I_inv_and_w_with_K(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.set_to_equilibrium_shape(bodies);

    for(int n=0; n<n_iter; n++) {
        shape.calculate_I_inv_and_w_with_K(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);
        
        shape.set_to_equilibrium_shape(bodies);
    }
    
    force->update_acceleration_and_torque_load(bodies, false);

    orbit.kick_v(bodies, dt);

    spin.kick_L(bodies, dt);

    shape.calculate_I_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.set_to_equilibrium_shape(bodies);

    for(int n=0; n<n_iter; n++) {
        shape.calculate_I_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);
        
        shape.set_to_equilibrium_shape(bodies);
    }

    force->update_torque_load(bodies, false);

    spin.kick_K(bodies, dth);

    orbit.drift_r(bodies, dth);    
}

void Evolver_equilibrium::sync_step(vector<Body> &bodies, Force *force) {
    force->update_tidal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.set_to_equilibrium_shape(bodies);

    shape.calculate_I_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.set_to_equilibrium_shape(bodies);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);
        
        shape.set_to_equilibrium_shape(bodies);
    }
    
    shape.calculate_I_inv_and_w(bodies);
}

void Evolver_equilibrium::mclachlan_step(vector<Body> &bodies, double dt, Force *force) {
    leapfrog_step(bodies, w1*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w3*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w1*dt, force);        
}

void Evolver_equilibrium::initialize(vector<Body> &bodies, Force *force) {
    shape.calculate_I_inv_and_w(bodies);

    force->update_tidal_deformation(bodies);
    force->update_centrifugal_deformation(bodies);

    force->update_acceleration_and_torque(bodies, false);
}
void Evolver_equilibrium::evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
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

    sync_step(bodies, force);              
}
void Evolver_equilibrium::evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
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

    sync_step(bodies, force);              
}



