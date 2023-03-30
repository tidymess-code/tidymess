#include "Spin.h"

// Initializers
Spin::Spin() {
    ;
}

// Calculators
void Spin::calculate_w(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {
            b->w[0] = b->I_inv[0]*b->L[0] + b->I_inv[1]*b->L[1] + b->I_inv[2]*b->L[2];
            b->w[1] = b->I_inv[1]*b->L[0] + b->I_inv[3]*b->L[1] + b->I_inv[4]*b->L[2];
            b->w[2] = b->I_inv[2]*b->L[0] + b->I_inv[4]*b->L[1] + b->I_inv[5]*b->L[2];
        }
        else {
            b->w[0] = 0;
            b->w[1] = 0;
            b->w[2] = 0;
        }
    }
}
void Spin::calculate_w_with_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {
            b->w[0] = b->J_inv[0]*b->L[0] + b->J_inv[1]*b->L[1] + b->J_inv[2]*b->L[2];
            b->w[1] = b->J_inv[1]*b->L[0] + b->J_inv[3]*b->L[1] + b->J_inv[4]*b->L[2];
            b->w[2] = b->J_inv[2]*b->L[0] + b->J_inv[4]*b->L[1] + b->J_inv[5]*b->L[2];
        }
        else {
            b->w[0] = 0;
            b->w[1] = 0;
            b->w[2] = 0;        
        }
    }
}
void Spin::calculate_w_with_K(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {
            b->w[0] = b->I_inv[0]*b->K[0] + b->I_inv[1]*b->K[1] + b->I_inv[2]*b->K[2];
            b->w[1] = b->I_inv[1]*b->K[0] + b->I_inv[3]*b->K[1] + b->I_inv[4]*b->K[2];
            b->w[2] = b->I_inv[2]*b->K[0] + b->I_inv[4]*b->K[1] + b->I_inv[5]*b->K[2];
        }
        else {
            b->w[0] = 0;
            b->w[1] = 0;
            b->w[2] = 0;                
        }
    }
}

void Spin::calculate_L(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {
            b->L[0] = b->I[0]*b->w[0]+b->I[1]*b->w[1]+b->I[2]*b->w[2];
            b->L[1] = b->I[1]*b->w[0]+b->I[3]*b->w[1]+b->I[4]*b->w[2];
            b->L[2] = b->I[2]*b->w[0]+b->I[4]*b->w[1]+b->I[5]*b->w[2];
        }
        else {
            b->L[0] = 0;
            b->L[1] = 0;
            b->L[2] = 0;                        
        }
    }
}
void Spin::init_K(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->K[k] = b->L[k];
        }
    }
}

void Spin::memorize_w(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->w_prev[k] = b->w[k];
        }
    }
}
void Spin::reset_w(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->w[k] = b->w_prev[k];
        }
    }
}

void Spin::kick_L(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {  
            for(int k=0; k<3; k++) {
                b->L[k] += b->T[k]*dt;
            }
        }
    }
}
void Spin::kick_K(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {  
            for(int k=0; k<3; k++) {
                b->K[k] += b->T[k]*dt;
            }
        }
    }
}

void Spin::kick_L_mb(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {  
            if(b->a_mb > 0) {
                double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);

                //double expo = exp(-b->a_mb * w2 * dt);
                //for(int k=0; k<3; k++) { 
                //    b->L[k] = b->L[k] * expo;
                //}

                double x = b->a_mb * w2 * dt;
                for(int k=0; k<3; k++) { 
                    b->L[k] = b->L[k] * (1-x);
                }            
            } 
        }
    }
}
void Spin::kick_K_mb(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {  
            if(b->a_mb > 0) {
                double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);

                //double expo = exp(-b->a_mb * w2 * dt);
                //for(int k=0; k<3; k++) { 
                //    b->K[k] = b->K[k] * expo;
                //}

                double x = b->a_mb * w2 * dt;
                for(int k=0; k<3; k++) { 
                    b->K[k] = b->K[k] * (1-x);
                }            
            }
        }
    }
}

// Diagnostics
double Spin::get_kinetic_energy(vector<Body> &bodies) {
    double E = 0;
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {  
            double dwxIdwx = b->w[0]*(b->w[0]*b->I[0] + b->w[1]*b->I[1] + b->w[2]*b->I[2]);
            double dwyIdwy = b->w[1]*(b->w[0]*b->I[1] + b->w[1]*b->I[3] + b->w[2]*b->I[4]);
            double dwzIdwz = b->w[2]*(b->w[0]*b->I[2] + b->w[1]*b->I[4] + b->w[2]*b->I[5]);
            double myE = 0.5 * (dwxIdwx + dwyIdwy + dwzIdwz);
            E += myE;
        }
    }
    return E;
}

array<double, 3> Spin::get_angular_momentum(vector<Body> &bodies) {
    array<double, 3> L = {};
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type >= 2) {  
            double Lx = b->I[0]*b->w[0] + b->I[1]*b->w[1] + b->I[2]*b->w[2];
            double Ly = b->I[1]*b->w[0] + b->I[3]*b->w[1] + b->I[4]*b->w[2];
            double Lz = b->I[2]*b->w[0] + b->I[4]*b->w[1] + b->I[5]*b->w[2];
            L[0] += Lx;
            L[1] += Ly;
            L[2] += Lz;
        }
    }
    return L;
}


