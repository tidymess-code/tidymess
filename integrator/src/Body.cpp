#include "Body.h"

// Initializers
Body::Body() {
    this->reset();
}

Body::Body(double m, double R, double xi, double kf, double tau, double a_mb, double wx, double wy, double wz, double x, double y, double z, double vx, double vy, double vz) {
    this->reset();
    
    // Internal properties
    this->m = m;
    this->R = R;
    this->xi = xi;

    this->kf = kf;
    this->tau = tau;

    this->a_mb = a_mb;
            
    // Orbital properties
    this->r[0] = x;
    this->r[1] = y;
    this->r[2] = z;

    this->v[0] = vx;
    this->v[1] = vy;
    this->v[2] = vz;

    // Spin properties
    this->w[0] = wx;
    this->w[1] = wy;
    this->w[2] = wz;
    
    // Aux properties
    this->update_aux_properties(); 
}
Body::Body(double m, double R, double xi, double kf, double tau, double a_mb, array<double, 3> w, array<double, 3> r, array<double, 3> v) {
    this->reset();

    // Internal properties
    this->m = m;
    this->R = R;
    this->xi = xi;

    this->kf = kf;
    this->tau = tau;
        
    this->a_mb = a_mb;
        
    // Orbital properties
    for(int k=0; k<3; k++) {
        this->r[k] = r[k];
        this->v[k] = v[k];
    }

    // Spin properties
    for(int k=0; k<3; k++) {
        this->w[k] = w[k];
    }

    // Aux properties
    this->update_aux_properties(); 
}
Body::Body(array<double, 15> &d) {
    this->reset();

    // Internal properties
    this->m = d[0];
    this->R = d[1];
    this->xi = d[2];
    
    this->kf = d[3];
    this->tau = d[4];

    this->a_mb = d[5];
    
    // Spin properties
    this->w[0] = d[6];
    this->w[1] = d[7];
    this->w[2] = d[8];
        
    // Orbital properties
    this->r[0] = d[9];
    this->r[1] = d[10];
    this->r[2] = d[11];

    this->v[0] = d[12];
    this->v[1] = d[13];
    this->v[2] = d[14];

    // Aux properties
    this->update_aux_properties();  
}

void Body::setup(array<double, 15> &d) {
    this->reset();

    // Internal properties
    this->m = d[0];
    this->R = d[1];
    this->xi = d[2];
    
    this->kf = d[3];
    this->tau = d[4];
    
    this->a_mb = d[5];

    // Spin properties
    this->w[0] = d[6];
    this->w[1] = d[7];
    this->w[2] = d[8];
        
    // Orbital properties
    this->r[0] = d[9];
    this->r[1] = d[10];
    this->r[2] = d[11];

    this->v[0] = d[12];
    this->v[1] = d[13];
    this->v[2] = d[14];
    
    // Aux properties
    this->update_aux_properties(); 
}
void Body::set_id(int id) {
    this->id = id;
}
void Body::set_name(string name) {
    this->name = name;
}
    
void Body::reset() {
    // Tag
    this->id = 0;
    this->name = "0";

    this->particle_type = 3;

    // Internal properties
    this->m = 0.;
    this->R = 0.;
    this->xi = 0.;
    this->kf = 0.;
    this->tau = 0.;    
    this->a_mb = 0.;

    this->R3 = 0.;
    this->R5 = 0.;
    this->R5_3 = 0.;
    this->kf_R5 = 0.;
    this->kf_R5_3 = 0.;
    this->tau_inv = 0.;
    this->rho = 0.;
    
    this->roche_factor = 2.44; // fluid body
    this->roche_factor3 = this->roche_factor*this->roche_factor*this->roche_factor;
    
    // Orbital properties
    for(int i=0; i<3; i++) {
        this->r[i] = 0.;
        this->v[i] = 0.;
        this->a[i] = 0.;
        
        this->r_prev[i] = 0.;
        this->vv[i] = 0.;
    }

    // Spin properties
    for(int i=0; i<3; i++) {
        this->L[i] = 0.;
        this->w[i] = 0.;
        this->T[i] = 0.;
        
        this->L_prev[i] = 0.;
        this->w_prev[i] = 0.;
        this->K[i] = 0.;
    }

    // Shape properties
    for(int i=0; i<6; i++) {
        this->I_p[i] = 0;
        this->I_n[i] = 0;
        this->I[i] = 0;
        this->I_inv[i] = 0;
        this->I_e_r[i] = 0;
        this->dI_e_r[i] = 0;
        this->I_e_w[i] = 0;
        this->I_e[i] = 0;
        this->dI_e[i] = 0;
        this->dI_n[i] = 0;

        this->I_e_rh[i] = 0;        
        this->I_e_prev[i] = 0;        
        this->I_n_prev[i] = 0;        
        this->I_e_prev_bu[i] = 0;        

        this->J_n[i] = 0;
        this->J[i] = 0;        
        this->J_inv[i] = 0;        
    }         
    
    // Integration variables
    this->isLinear = false;
    
    this->expo = 0;
    this->expoh = 0;
    this->dt_rot = 0;
    this->dth_rot = 0;

    this->isLinear_1 = false;
    this->expo_1 = 0;
    this->expoh_1 = 0;
    this->dt_rot_1 = 0;
    this->dth_rot_1 = 0;

    this->isLinear_2 = false;
    this->expo_2 = 0;
    this->expoh_2 = 0;
    this->dt_rot_2 = 0;
    this->dth_rot_2 = 0;     
}

void Body::update_aux_properties() {
    this->R3 = pow(this->R, 3.);
    this->R5 = pow(this->R, 5.);
    this->R5_3 = this->R5 / 3.;    
    this->kf_R5 = this->kf*this->R5;
    this->kf_R5_3 = this->kf*this->R5_3;

    if(this->tau > 0) this->tau_inv = 1./this->tau;
    
    double V = 4./3*M_PI * this->R3;
    if(V > 0) this->rho = this->m / V;
    
    for(int k=0; k<3; k++) {
        this->vv[k] = v[k];
    }    
    
    if(this->m <= 0) {
        this->particle_type = 0;        
    }
    else {
        if(this->R <= 0 || this->xi <= 0) {
            this->particle_type = 1;           
        }
        else {
            if(this->kf <= 0) {
                this->particle_type = 2;               
            }
            else {
                this->particle_type = 3;                               
            }
        }
    }
}
    

