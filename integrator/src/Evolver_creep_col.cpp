#include "Evolver_creep_col.h"

// Evolve function
void Evolver_creep_col::leapfrog_step(vector<Body> &bodies, double dt, Force *force) {
    double dth = 0.5*dt;

    orbit.drift_r(bodies, dth);

    force->update_tidal_deformation_save(bodies);
    force->update_equilibrium_tensor_and_memorize(bodies);   

    shape.memorize_J(bodies);
    shape.deform_and_rotate_J_half(bodies, dth);

    shape.calculate_J_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    shape.reset_J(bodies);
    shape.deform_and_rotate_J_half(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_J_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);

        shape.reset_J(bodies);
        shape.deform_and_rotate_J_half(bodies, dth);
    }

    force->update_torque_load(bodies, true);

    spin.kick_K(bodies, dth);

    shape.calculate_J_inv_and_w_with_K(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    shape.memorize_I(bodies);
    shape.deform_and_rotate_half(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w_with_K(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);

        shape.reset_I(bodies);
        shape.deform_and_rotate_half(bodies, dth);
    }

    spin.memorize_w(bodies);
    shape.backup_I_e(bodies);

    force->update_acceleration_and_torque_load(bodies, false);

    orbit.kick_v(bodies, dt);

    spin.kick_L(bodies, dt);

    orbit.drift_r(bodies, dth);

    sync_step(bodies, force, dth);
}
void Evolver_creep_col::leapfrog_step_with_collisions(vector<Body> &bodies, double dt, Force *force) {
    double dth = 0.5*dt;

    orbit.drift_r(bodies, dth);

    force->update_tidal_deformation_save(bodies);
    force->update_equilibrium_tensor_and_memorize(bodies);   

    shape.memorize_J(bodies);
    shape.deform_and_rotate_J_half(bodies, dth);

    shape.calculate_J_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    shape.reset_J(bodies);
    shape.deform_and_rotate_J_half(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_J_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);

        shape.reset_J(bodies);
        shape.deform_and_rotate_J_half(bodies, dth);
    }

    force->update_torque_load(bodies, true);

    spin.kick_K(bodies, dth);

    shape.calculate_J_inv_and_w_with_K(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    shape.memorize_I(bodies);
    shape.deform_and_rotate_half(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w_with_K(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);

        shape.reset_I(bodies);
        shape.deform_and_rotate_half(bodies, dth);
    }

    spin.memorize_w(bodies);
    shape.backup_I_e(bodies);

    force->update_acceleration_and_torque_col_load(bodies, false);

    orbit.kick_v(bodies, dt);

    spin.kick_L(bodies, dt);

    orbit.drift_r(bodies, dth);

    sync_step(bodies, force, dth);
}

void Evolver_creep_col::sync_step(vector<Body> &bodies, Force *force, double dth) {
    force->update_tidal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);   

    shape.set_I_e_prev_to_backup(bodies);

    shape.memorize_I(bodies);
    
    shape.deform_and_rotate_half(bodies, dth);

    shape.calculate_I_inv_and_w(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);

    spin.reset_w(bodies);
    shape.reset_I(bodies);
    
    shape.deform_and_rotate_half(bodies, dth);

    for(int n=0; n<n_iter; n++) {    
        shape.calculate_I_inv_and_w(bodies);
        force->update_centrifugal_deformation(bodies);
        force->update_equilibrium_tensor(bodies);

        spin.reset_w(bodies);
        shape.reset_I(bodies);
    
        shape.deform_and_rotate_half(bodies, dth);
    }

    shape.calculate_I_inv_and_w(bodies);    

    reset_aux(bodies);    
}

void Evolver_creep_col::mclachlan_step(vector<Body> &bodies, double dt, Force *force) {        
    shape.calc_expos(bodies, w1*dt);
    shape.copy_expos_to_1(bodies);

    leapfrog_step(bodies, w1*dt, force);

    shape.calc_expos(bodies, w2*dt);
    shape.copy_expos_to_2(bodies);

    leapfrog_step(bodies, w2*dt, force);

    shape.calc_expos(bodies, w3*dt);

    leapfrog_step(bodies, w3*dt, force);

    shape.revert_expos_2(bodies);

    leapfrog_step(bodies, w2*dt, force);

    shape.revert_expos_1(bodies);

    leapfrog_step_with_collisions(bodies, w1*dt, force);
}

void Evolver_creep_col::initialize(vector<Body> &bodies, Force *force) {
    shape.calculate_I_inv_and_w(bodies);

    force->update_tidal_deformation(bodies);
    force->update_centrifugal_deformation(bodies);
    force->update_equilibrium_tensor(bodies);
    
    force->update_acceleration_and_torque(bodies, false);    
}
void Evolver_creep_col::evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) {  
    bool toContinue = true;

    num_integration_step = 0;
        
    while(toContinue) {
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
}
void Evolver_creep_col::evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) {  
    num_integration_step = 0;
        
    for(int i=0; i<N_step; i++) {
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
}



