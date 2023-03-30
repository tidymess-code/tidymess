#include "Evolver_base.h"

// Initializers
Evolver_base::Evolver_base() {
    dt_step = 1.0;
    
    w1 = 0.28;
    w2 = 0.62546642846767004501;
    w3 = 1 - 2*w1 - 2*w2;
    
    n_iter = 1;
    num_integration_step = 0;
    
    collision_mode = 0;
    roche_mode = 0;
    
    collision_flag = false;
    roche_flag = false;
}
    
// Get diagnostics
array<double, 3> Evolver_base::get_center_of_mass(vector<Body> &bodies) {
    return orbit.get_center_of_mass(bodies);
}
array<double, 3> Evolver_base::get_center_of_mass_velocity(vector<Body> &bodies) {
    return orbit.get_center_of_mass_velocity(bodies);
}
array<double, 3> Evolver_base::get_orbital_angular_momentum(vector<Body> &bodies) {
    return orbit.get_angular_momentum(bodies);
}
double Evolver_base::get_orbital_kinetic_energy(vector<Body> &bodies) {
    return orbit.get_kinetic_energy(bodies);
}
    
void Evolver_base::assign_vectors(vector<Body> &bodies, Force *force) {    
    int N = bodies.size();
    int Np = 0.5*N*(N-1); // this only works for N < few k in terms of memory
    force->dr2_vec.assign(Np, 0);
    force->dr_3_vec.assign(Np, 0);
}    

void Evolver_base::set_n_iter(int n_iter) {
    this->n_iter = n_iter; 
}
int Evolver_base::get_n_iter() {
    return n_iter;
}

void Evolver_base::reset_aux(vector<Body> &bodies) {
    orbit.init_vv(bodies);
    spin.init_K(bodies);
    shape.init_J(bodies);
}
       
void Evolver_base::set_collision_mode(int collision_mode) {
    this->collision_mode = collision_mode;
}
void Evolver_base::set_roche_mode(int roche_mode) {
    this->roche_mode = roche_mode;
}

bool Evolver_base::stopping_condition_time(double t0, double t1, int dt_sgn) {
    bool toContinue = true;
    if(dt_sgn == 1) {
        if(t0 >= t1) toContinue = false;
    }
    else { 
        if(t0 <= t1) toContinue = false;
    }
    return toContinue;
}


    
