#include "Evolver_linear_mb.h"

// Evolve function
void Evolver_linear_mb::leapfrog_step(vector<Body> &bodies, double dt, Force *force) {
    double dth = 0.5*dt;

    orbit.drift_r(bodies, dth);

    force->update_tidal_deformation_save(bodies);
    force->update_equilibrium_tensor_and_memorize(bodies);   

    shape.kick_linear_discrete(bodies, dth);

    shape.calculate_I_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.kick_linear_discrete(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);
        
        shape.kick_linear_discrete(bodies, dth);
    }

    force->update_torque_load(bodies, false);

    spin.kick_K_mb(bodies, dth);
    spin.kick_K(bodies, dth);

    shape.calculate_I_inv_and_w_with_K(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.kick_linear_discrete(bodies, dth);

    for(int n=0; n<n_iter; n++) {
        shape.calculate_I_inv_and_w_with_K(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);
        
        shape.kick_linear_discrete(bodies, dth);
    }
    
    spin.memorize_w(bodies);
    shape.backup_I_e(bodies);
    
    force->update_acceleration_and_torque_load(bodies, false);
    
    orbit.kick_v(bodies, dt);

    spin.kick_L(bodies, dt);
    spin.kick_L_mb(bodies, dt);

    orbit.drift_r(bodies, dth);      

    sync_step(bodies, force, dt_step);              
}

void Evolver_linear_mb::sync_step(vector<Body> &bodies, Force *force, double dt) {
    double dth = 0.5*dt;
    
    force->update_tidal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    shape.set_I_e_prev_to_backup(bodies);
    
    shape.kick_linear_discrete(bodies, dth);

    shape.calculate_I_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    spin.reset_w(bodies);
    
    shape.kick_linear_discrete(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);

        spin.reset_w(bodies);
        
        shape.kick_linear_discrete(bodies, dth);
    }
    
    shape.calculate_I_inv_and_w(bodies);

    reset_aux(bodies);
}

void Evolver_linear_mb::mclachlan_step(vector<Body> &bodies, double dt, Force *force) {
    leapfrog_step(bodies, w1*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w3*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w1*dt, force);   
}

void Evolver_linear_mb::initialize(vector<Body> &bodies, Force *force) {
    shape.calculate_I_inv_and_w(bodies);

    force->update_tidal_deformation(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    force->update_acceleration_and_torque(bodies, false);        
}
void Evolver_linear_mb::evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
    bool toContinue = true;

    num_integration_step = 0;
        
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
void Evolver_linear_mb::evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
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


