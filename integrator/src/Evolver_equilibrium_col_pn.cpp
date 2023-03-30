#include "Evolver_equilibrium_col_pn.h"

// Evolve function
void Evolver_equilibrium_col_pn::leapfrog_step(vector<Body> &bodies, double dt, Force *force) {
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

    force->update_acceleration_and_torque_pn_load(bodies, false, false);

    orbit.kick_vv(bodies, dth);
    
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
    
    force->update_acceleration_and_torque_pn_load(bodies, true, false);

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

    force->update_acceleration_and_torque_pn_load(bodies, false, false);

    orbit.kick_vv(bodies, dth);
    
    spin.kick_K(bodies, dth);

    orbit.drift_r(bodies, dth);    
}
void Evolver_equilibrium_col_pn::leapfrog_step_with_collisions(vector<Body> &bodies, double dt, Force *force) {
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

    force->update_acceleration_and_torque_pn_load(bodies, false, false);

    orbit.kick_vv(bodies, dth);
    
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
    
    force->update_acceleration_and_torque_pn_col_load(bodies, true, false);

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

    force->update_acceleration_and_torque_pn_load(bodies, false, false);

    orbit.kick_vv(bodies, dth);
    
    spin.kick_K(bodies, dth);

    orbit.drift_r(bodies, dth);    
}

void Evolver_equilibrium_col_pn::sync_step(vector<Body> &bodies, Force *force) {
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

void Evolver_equilibrium_col_pn::mclachlan_step(vector<Body> &bodies, double dt, Force *force) {
    leapfrog_step(bodies, w1*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step(bodies, w3*dt, force);
    leapfrog_step(bodies, w2*dt, force);
    leapfrog_step_with_collisions(bodies, w1*dt, force);
}

void Evolver_equilibrium_col_pn::initialize(vector<Body> &bodies, Force *force) {
    shape.calculate_I_inv_and_w(bodies);

    force->update_tidal_deformation(bodies);
    force->update_centrifugal_deformation(bodies);
    
    force->update_acceleration_and_torque_pn(bodies, false, false);    
}
void Evolver_equilibrium_col_pn::evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
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

        //----------------------------------------------

        if(force->collision_detected) {
            collision_flag = true;
        }
        if(force->roche_detected) {
            roche_flag = true;
        }        
        
        if(force->collision_detected) {        
            if(collision_mode >= 2) {            
                break;
            }
        }
        if(force->roche_detected) {
            if(roche_mode >= 2) {
                break;
            }
        }
    }
    
    sync_step(bodies, force);          
}
void Evolver_equilibrium_col_pn::evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) {
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

        //----------------------------------------------

        if(force->collision_detected) {
            collision_flag = true;
        }
        if(force->roche_detected) {
            roche_flag = true;
        }        
        
        if(force->collision_detected) {        
            if(collision_mode >= 2) {            
                break;
            }
        }
        if(force->roche_detected) {
            if(roche_mode >= 2) {
                break;
            }
        }
    }
    
    sync_step(bodies, force);          
}




