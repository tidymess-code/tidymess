#include "Timestep_base.h"

// Initializer
Timestep_base::Timestep_base() {
    set_default_values();
}

void Timestep_base::set_default_values() {
    this->dt_const = 0.015625;
    this->eta = 0.0625;
    this->eta2 = this->eta*this->eta;
    
    this->h = 0;   
    this->h2 = 0; 
    this->dh = 0;
    this->h_spin = 0;
    this->h_shape = 0;

    this->dt_prev = 0;
    this->dt_next = 0;

    this->fp_bound = 4;
    this->fm_bound = 0.25;
    this->fp2_bound = this->fp_bound*this->fp_bound;
    this->fm2_bound = this->fm_bound*this->fm_bound;
    this->alpha = 10;
    this->dt_sgn = 1;
}
    
// Setters and Getters
void Timestep_base::set_dt_const(double dt_const) {
    this->dt_const = dt_const;
}
double Timestep_base::get_dt_const() {
    return dt_const;
}
    
void Timestep_base::set_eta(double eta) {
    this->eta = eta;
    this->eta2 = this->eta*this->eta;
}
double Timestep_base::get_eta() {
    return eta;
}

void Timestep_base::commit() {
    dt_const = abs(dt_const);

    if(dt_sgn > 0) {
        if(dt_prev < 0) dt_prev = abs(dt_prev);
    }
    else {
        if(dt_prev > 0) dt_prev = -dt_prev;
    }
}

void Timestep_base::set_dt_sgn(double t0, double t1) {
    dt_sgn = (t1 >= t0) ? 1 : -1;
    commit();
}
void Timestep_base::set_dt_sgn(int dt_sgn) {
    this->dt_sgn = dt_sgn;
    commit();
}

void Timestep_base::set_dt_prev(double dt_prev) {
    this->dt_prev = abs(dt_prev);

    if(this->dt_sgn < 0) {
        this->dt_prev *= -1;
    }
}


    

