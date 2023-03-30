#include "Shape.h"

// Initializers
Shape::Shape() {
    f_lim = 7; // transition ratio h/tau from full equation to linear approximation
}

void Shape::set_permanent_shape(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        double I0 = b->xi*b->m*b->R*b->R; 
        b->I_p.fill(0);
        b->I_p[0] = I0;
        b->I_p[3] = I0;
        b->I_p[5] = I0;
    }
}
void Shape::set_to_spherical_shape(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->I_n.fill(0);
        for(int k=0; k<6; k++) {
            b->I[k] = b->I_p[k] + b->I_n[k];
        }
    }
}
void Shape::set_to_equilibrium_shape(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<5; k++) {
            b->I_n[k] = b->I_e[k];
            b->I[k] = b->I_p[k] + b->I_n[k];
        }
        b->I_n[5] = -(b->I_n[0]+b->I_n[3]);
        b->I[5] = b->I_p[5] + b->I_n[5];
    }
}
void Shape::set_to_equilibrium_shape_with_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<5; k++) {
            b->J_n[k] = b->I_e[k];
            b->J[k] = b->I_p[k] + b->J_n[k];
        }
        b->J_n[5] = -(b->J_n[0]+b->J_n[3]);
        b->J[5] = b->I_p[5] + b->J_n[5];
    }
}
void Shape::set_to_linear_shape(vector<Body> &bodies) {
    kick_linear(bodies);
}
void Shape::init_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->J_n[k] = b->I_n[k];
            b->J[k]   = b->I[k];
        }
    }
}

void Shape::calc_expos(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3 && b->tau > 0) {
            double dt_abs       = abs(dt);
        
            double f_ratio      = dt*b->tau_inv; 
            double f_ratio_abs  = dt_abs*b->tau_inv; 

            double fh_ratio     = 0.5*f_ratio; 
            double fh_ratio_abs = 0.5*f_ratio_abs;

            double expoh_abs = exp(-fh_ratio_abs);
            double expo_abs  = expoh_abs*expoh_abs;

            if(f_ratio_abs < f_lim) {
                double expoh = (dt > 0) ? expoh_abs : 1./expoh_abs;
                        
                double expo      = expoh*expoh;
                
                int dt_sgn       = dt/dt_abs;

                double dt_step   = b->tau*(1-expo_abs);            
                dt_step         *= dt_sgn;

                double dth_step  = b->tau*(1-expoh_abs);            
                dth_step        *= dt_sgn;

                b->expo    = expo;
                b->expoh   = expoh;
                b->dt_rot  = dt_step;
                b->dth_rot = dth_step;
                b->isLinear = false;
            }
            else {
                double dth_step  = b->tau*(1-expoh_abs);             
                double dt_step   = b->tau*(1-expo_abs);                            
 
                b->dt_rot  = dt_step;
                b->dth_rot = dth_step;
                b->isLinear = true;
            }
        }
    }
}    
void Shape::copy_expos_to_1(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->expo_1     = b->expo;
        b->expoh_1    = b->expoh;
        b->dt_rot_1   = b->dt_rot;
        b->dth_rot_1  = b->dth_rot;
        b->isLinear_1 = b->isLinear;
    }            
}
void Shape::copy_expos_to_2(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->expo_2     = b->expo;
        b->expoh_2    = b->expoh;
        b->dt_rot_2   = b->dt_rot;
        b->dth_rot_2  = b->dth_rot;
        b->isLinear_2 = b->isLinear;
    }            
}
void Shape::revert_expos_1(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->expo     = b->expo_1;
        b->expoh    = b->expoh_1;
        b->dt_rot   = b->dt_rot_1;
        b->dth_rot  = b->dth_rot_1;
        b->isLinear = b->isLinear_1;
    }            
}
void Shape::revert_expos_2(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->expo     = b->expo_2;
        b->expoh    = b->expoh_2;
        b->dt_rot   = b->dt_rot_2;
        b->dth_rot  = b->dth_rot_2;
        b->isLinear = b->isLinear_2;
    }            
}
            
// Calculators
void Shape::calculate_I_inv(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double I11 = b->I[1]*b->I[1];
            double I22 = b->I[2]*b->I[2];
            double I44 = b->I[4]*b->I[4];

            double I24 = b->I[2]*b->I[4];
            double I35 = b->I[3]*b->I[5];

            double D = 1./(b->I[0]*I35 + 2*b->I[1]*I24 - b->I[0]*I44 - b->I[3]*I22 - b->I[5]*I11);

            b->I_inv[0] = D * (I35-I44);
            b->I_inv[1] = D * (I24-b->I[1]*b->I[5]);
            b->I_inv[2] = D * (b->I[1]*b->I[4]-b->I[3]*b->I[2]);
            b->I_inv[3] = D * (b->I[0]*b->I[5]-I22);
            b->I_inv[4] = D * (b->I[1]*b->I[2]-b->I[0]*b->I[4]);
            b->I_inv[5] = D * (b->I[0]*b->I[3]-I11);
        }
        else if(b->particle_type == 2) {
            b->I_inv[0] = 1./b->I[0];
            b->I_inv[1] = 0;
            b->I_inv[2] = 0;
            b->I_inv[3] = 1./b->I[3];
            b->I_inv[4] = 0;
            b->I_inv[5] = 1./b->I[5];
        }
    }
}
void Shape::calculate_I_inv_and_w(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double I11 = b->I[1]*b->I[1];
            double I22 = b->I[2]*b->I[2];
            double I44 = b->I[4]*b->I[4];

            double I24 = b->I[2]*b->I[4];
            double I35 = b->I[3]*b->I[5];

            double D = 1./(b->I[0]*I35 + 2*b->I[1]*I24 - b->I[0]*I44 - b->I[3]*I22 - b->I[5]*I11);

            b->I_inv[0] = D * (I35-I44);
            b->I_inv[1] = D * (I24-b->I[1]*b->I[5]);
            b->I_inv[2] = D * (b->I[1]*b->I[4]-b->I[3]*b->I[2]);
            b->I_inv[3] = D * (b->I[0]*b->I[5]-I22);
            b->I_inv[4] = D * (b->I[1]*b->I[2]-b->I[0]*b->I[4]);
            b->I_inv[5] = D * (b->I[0]*b->I[3]-I11);
        
            b->w[0] = b->I_inv[0]*b->L[0] + b->I_inv[1]*b->L[1] + b->I_inv[2]*b->L[2];
            b->w[1] = b->I_inv[1]*b->L[0] + b->I_inv[3]*b->L[1] + b->I_inv[4]*b->L[2];
            b->w[2] = b->I_inv[2]*b->L[0] + b->I_inv[4]*b->L[1] + b->I_inv[5]*b->L[2];        
        }
        else if(b->particle_type == 2) {
            b->I_inv[0] = 1./b->I[0];
            b->I_inv[1] = 0;
            b->I_inv[2] = 0;
            b->I_inv[3] = 1./b->I[3];
            b->I_inv[4] = 0;
            b->I_inv[5] = 1./b->I[5];

            b->w[0] = b->I_inv[0]*b->L[0] + b->I_inv[1]*b->L[1] + b->I_inv[2]*b->L[2];
            b->w[1] = b->I_inv[1]*b->L[0] + b->I_inv[3]*b->L[1] + b->I_inv[4]*b->L[2];
            b->w[2] = b->I_inv[2]*b->L[0] + b->I_inv[4]*b->L[1] + b->I_inv[5]*b->L[2];                 
        }
    }
}

void Shape::calculate_J_inv_and_w(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double I11 = b->J[1]*b->J[1];
            double I22 = b->J[2]*b->J[2];
            double I44 = b->J[4]*b->J[4];

            double I24 = b->J[2]*b->J[4];
            double I35 = b->J[3]*b->J[5];

            double D = 1./(b->J[0]*I35 + 2*b->J[1]*I24 - b->J[0]*I44 - b->J[3]*I22 - b->J[5]*I11);

            b->J_inv[0] = D * (I35-I44);
            b->J_inv[1] = D * (I24-b->J[1]*b->J[5]);
            b->J_inv[2] = D * (b->J[1]*b->J[4]-b->J[3]*b->J[2]);
            b->J_inv[3] = D * (b->J[0]*b->J[5]-I22);
            b->J_inv[4] = D * (b->J[1]*b->J[2]-b->J[0]*b->J[4]);
            b->J_inv[5] = D * (b->J[0]*b->J[3]-I11);
        
            b->w[0] = b->J_inv[0]*b->L[0] + b->J_inv[1]*b->L[1] + b->J_inv[2]*b->L[2];
            b->w[1] = b->J_inv[1]*b->L[0] + b->J_inv[3]*b->L[1] + b->J_inv[4]*b->L[2];
            b->w[2] = b->J_inv[2]*b->L[0] + b->J_inv[4]*b->L[1] + b->J_inv[5]*b->L[2];        
        }
        else if(b->particle_type == 2) {
            b->J_inv[0] = 1./b->J[0];
            b->J_inv[1] = 0;
            b->J_inv[2] = 0;
            b->J_inv[3] = 1./b->J[3];
            b->J_inv[4] = 0;
            b->J_inv[5] = 1./b->J[5];

            b->w[0] = b->J_inv[0]*b->L[0] + b->J_inv[1]*b->L[1] + b->J_inv[2]*b->L[2];
            b->w[1] = b->J_inv[1]*b->L[0] + b->J_inv[3]*b->L[1] + b->J_inv[4]*b->L[2];
            b->w[2] = b->J_inv[2]*b->L[0] + b->J_inv[4]*b->L[1] + b->J_inv[5]*b->L[2];        
        }        
    }
}

void Shape::calculate_I_inv_and_w_with_K(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double I11 = b->I[1]*b->I[1];
            double I22 = b->I[2]*b->I[2];
            double I44 = b->I[4]*b->I[4];

            double I24 = b->I[2]*b->I[4];
            double I35 = b->I[3]*b->I[5];

            double D = 1./(b->I[0]*I35 + 2*b->I[1]*I24 - b->I[0]*I44 - b->I[3]*I22 - b->I[5]*I11);

            b->I_inv[0] = D * (I35-I44);
            b->I_inv[1] = D * (I24-b->I[1]*b->I[5]);
            b->I_inv[2] = D * (b->I[1]*b->I[4]-b->I[3]*b->I[2]);
            b->I_inv[3] = D * (b->I[0]*b->I[5]-I22);
            b->I_inv[4] = D * (b->I[1]*b->I[2]-b->I[0]*b->I[4]);
            b->I_inv[5] = D * (b->I[0]*b->I[3]-I11);
        
            b->w[0] = b->I_inv[0]*b->K[0] + b->I_inv[1]*b->K[1] + b->I_inv[2]*b->K[2];
            b->w[1] = b->I_inv[1]*b->K[0] + b->I_inv[3]*b->K[1] + b->I_inv[4]*b->K[2];
            b->w[2] = b->I_inv[2]*b->K[0] + b->I_inv[4]*b->K[1] + b->I_inv[5]*b->K[2];        
        }
        else if(b->particle_type == 2) {
            b->I_inv[0] = 1./b->I[0];
            b->I_inv[1] = 0;
            b->I_inv[2] = 0;
            b->I_inv[3] = 1./b->I[3];
            b->I_inv[4] = 0;
            b->I_inv[5] = 1./b->I[5];

            b->w[0] = b->I_inv[0]*b->K[0] + b->I_inv[1]*b->K[1] + b->I_inv[2]*b->K[2];
            b->w[1] = b->I_inv[1]*b->K[0] + b->I_inv[3]*b->K[1] + b->I_inv[4]*b->K[2];
            b->w[2] = b->I_inv[2]*b->K[0] + b->I_inv[4]*b->K[1] + b->I_inv[5]*b->K[2];        
        }        
    }
}
void Shape::calculate_J_inv_and_w_with_K(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double I11 = b->J[1]*b->J[1];
            double I22 = b->J[2]*b->J[2];
            double I44 = b->J[4]*b->J[4];

            double I24 = b->J[2]*b->J[4];
            double I35 = b->J[3]*b->J[5];

            double D = 1./(b->J[0]*I35 + 2*b->J[1]*I24 - b->J[0]*I44 - b->J[3]*I22 - b->J[5]*I11);

            b->J_inv[0] = D * (I35-I44);
            b->J_inv[1] = D * (I24-b->J[1]*b->J[5]);
            b->J_inv[2] = D * (b->J[1]*b->J[4]-b->J[3]*b->J[2]);
            b->J_inv[3] = D * (b->J[0]*b->J[5]-I22);
            b->J_inv[4] = D * (b->J[1]*b->J[2]-b->J[0]*b->J[4]);
            b->J_inv[5] = D * (b->J[0]*b->J[3]-I11);
        
            b->w[0] = b->J_inv[0]*b->K[0] + b->J_inv[1]*b->K[1] + b->J_inv[2]*b->K[2];
            b->w[1] = b->J_inv[1]*b->K[0] + b->J_inv[3]*b->K[1] + b->J_inv[4]*b->K[2];
            b->w[2] = b->J_inv[2]*b->K[0] + b->J_inv[4]*b->K[1] + b->J_inv[5]*b->K[2];        
        }
        else if(b->particle_type == 2) {
            b->J_inv[0] = 1./b->J[0];
            b->J_inv[1] = 0;
            b->J_inv[2] = 0;
            b->J_inv[3] = 1./b->J[3];
            b->J_inv[4] = 0;
            b->J_inv[5] = 1./b->J[5];

            b->w[0] = b->J_inv[0]*b->K[0] + b->J_inv[1]*b->K[1] + b->J_inv[2]*b->K[2];
            b->w[1] = b->J_inv[1]*b->K[0] + b->J_inv[3]*b->K[1] + b->J_inv[4]*b->K[2];
            b->w[2] = b->J_inv[2]*b->K[0] + b->J_inv[4]*b->K[1] + b->J_inv[5]*b->K[2];        
        }                
    }
}

// Evolvers
void Shape::kick_linear(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            b->I_n[0] = b->I_e[0] + b->tau*( 2*b->I_e[2]*b->w[1] - 2*b->I_e[1]*b->w[2] - b->dI_e[0] );
            b->I_n[3] = b->I_e[3] + b->tau*(-2*b->I_e[4]*b->w[0] + 2*b->I_e[1]*b->w[2] - b->dI_e[3] );
            b->I_n[1] = b->I_e[1] + b->tau*(  -b->I_e[2]*b->w[0] + b->I_e[4]*b->w[1] + (b->I_e[0]-b->I_e[3])*b->w[2] - b->dI_e[1]);
            b->I_n[2] = b->I_e[2] + b->tau*(    b->I_e[1]*b->w[0] - (2*b->I_e[0]+b->I_e[3])*b->w[1] - b->I_e[4]*b->w[2] - b->dI_e[2] );
            b->I_n[4] = b->I_e[4] + b->tau*(  (b->I_e[0]+2*b->I_e[3])*b->w[0] - b->I_e[1]*b->w[1] + b->I_e[2]*b->w[2] - b->dI_e[4] );

            b->I_n[5] = -(b->I_n[0]+b->I_n[3]);

            for(int i=0; i<6; i++) {
                b->I[i] = b->I_p[i] + b->I_n[i];
            }
        }
    }
}
void Shape::kick_linear_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            b->J_n[0] = b->I_e[0] + b->tau*( 2*b->I_e[2]*b->w[1] - 2*b->I_e[1]*b->w[2] - b->dI_e[0] );
            b->J_n[3] = b->I_e[3] + b->tau*(-2*b->I_e[4]*b->w[0] + 2*b->I_e[1]*b->w[2] - b->dI_e[3] );
            b->J_n[1] = b->I_e[1] + b->tau*(  -b->I_e[2]*b->w[0] + b->I_e[4]*b->w[1] + (b->I_e[0]-b->I_e[3])*b->w[2] - b->dI_e[1]);
            b->J_n[2] = b->I_e[2] + b->tau*(    b->I_e[1]*b->w[0] - (2*b->I_e[0]+b->I_e[3])*b->w[1] - b->I_e[4]*b->w[2] - b->dI_e[2] );
            b->J_n[4] = b->I_e[4] + b->tau*(  (b->I_e[0]+2*b->I_e[3])*b->w[0] - b->I_e[1]*b->w[1] + b->I_e[2]*b->w[2] - b->dI_e[4] );

            b->J_n[5] = -(b->J_n[0]+b->J_n[3]);

            for(int i=0; i<6; i++) {
                b->J[i] = b->I_p[i] + b->J_n[i];
            }
        }
    }
}

void Shape::kick_linear_discrete(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            for(int k=0; k<5; k++) {
                b->dI_e[k] = (b->I_e[k]-b->I_e_prev[k]) / dt;
            }
    
            b->I_n[0] = b->I_e[0] + b->tau*( 2*b->I_e[2]*b->w[1] - 2*b->I_e[1]*b->w[2] - b->dI_e[0] );
            b->I_n[3] = b->I_e[3] + b->tau*(-2*b->I_e[4]*b->w[0] + 2*b->I_e[1]*b->w[2] - b->dI_e[3] );
            b->I_n[1] = b->I_e[1] + b->tau*(  -b->I_e[2]*b->w[0] + b->I_e[4]*b->w[1] + (b->I_e[0]-b->I_e[3])*b->w[2] - b->dI_e[1]);
            b->I_n[2] = b->I_e[2] + b->tau*(    b->I_e[1]*b->w[0] - (2*b->I_e[0]+b->I_e[3])*b->w[1] - b->I_e[4]*b->w[2] - b->dI_e[2] );
            b->I_n[4] = b->I_e[4] + b->tau*(  (b->I_e[0]+2*b->I_e[3])*b->w[0] - b->I_e[1]*b->w[1] + b->I_e[2]*b->w[2] - b->dI_e[4] );

            b->I_n[5] = -(b->I_n[0]+b->I_n[3]);

            for(int i=0; i<6; i++) {
                b->I[i] = b->I_p[i] + b->I_n[i];
            }
        }
    }
}
void Shape::kick_linear_discrete_J(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            for(int k=0; k<5; k++) {
                b->dI_e[k] = (b->I_e[k]-b->I_e_prev[k]) / dt;
            }

            b->J_n[0] = b->I_e[0] + b->tau*( 2*b->I_e[2]*b->w[1] - 2*b->I_e[1]*b->w[2] - b->dI_e[0] );
            b->J_n[3] = b->I_e[3] + b->tau*(-2*b->I_e[4]*b->w[0] + 2*b->I_e[1]*b->w[2] - b->dI_e[3] );
            b->J_n[1] = b->I_e[1] + b->tau*(  -b->I_e[2]*b->w[0] + b->I_e[4]*b->w[1] + (b->I_e[0]-b->I_e[3])*b->w[2] - b->dI_e[1]);
            b->J_n[2] = b->I_e[2] + b->tau*(    b->I_e[1]*b->w[0] - (2*b->I_e[0]+b->I_e[3])*b->w[1] - b->I_e[4]*b->w[2] - b->dI_e[2] );
            b->J_n[4] = b->I_e[4] + b->tau*(  (b->I_e[0]+2*b->I_e[3])*b->w[0] - b->I_e[1]*b->w[1] + b->I_e[2]*b->w[2] - b->dI_e[4] );

            b->J_n[5] = -(b->J_n[0]+b->J_n[3]);

            for(int i=0; i<6; i++) {
                b->J[i] = b->I_p[i] + b->J_n[i];
            }
        }
    }
}

void Shape::kick_direct(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            for(int k=0; k<5; k++) {
                b->I_n[k] += b->dI_n[k]*dt;
            }
            b->I_n[5] = -(b->I_n[0]+b->I_n[3]);

            for(int k=0; k<6; k++) {
                b->I[k] = b->I_p[k] + b->I_n[k];
            }
        }
    }
}
void Shape::kick_direct_J(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            for(int k=0; k<5; k++) {
                b->J_n[k] += b->dI_n[k]*dt;
            }
            b->J_n[5] = -(b->J_n[0]+b->J_n[3]);

            for(int k=0; k<6; k++) {
                b->J[k] = b->I_p[k] + b->J_n[k];
            }
        }
    }
}

void Shape::deform_and_rotate(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            if(b->tau > 0) {
                double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);
                double w  = sqrt(w2);
            
                double dt_step = b->dt_rot; 

                double wdt = 0.5*w*dt_step;
                double cosa = cos(wdt);
                double sina = sin(wdt);
            
                array<double, 4> q0 = {};
                q0[0] = 1.;

                if(w > 0) {
                    double w_1 = 1/w;

                    q0[0] = cosa;
                    q0[1] = sina*b->w[0]*w_1;
                    q0[2] = sina*b->w[1]*w_1;
                    q0[3] = sina*b->w[2]*w_1;
                }

                if(b->isLinear) {
                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->I_n[k] = b->I_e[k] - b->tau*dIe;
                    }
                    b->I_n[5] = -(b->I_n[0]+b->I_n[3]);

                    array<double, 5> I2 = {};
                    rotate_tensor(q0, b->I_n, I2); 
                    for(int k=0; k<5; k++) {
                        b->I_n[k] = I2[k]; 
                    }      
                    b->I_n[5] = -(b->I_n[0]+b->I_n[3]);                      
                }
                else {
                    double expo = b->expo;

                    if(dt < 0) {
                        array<double, 5> I1 = {};
                        rotate_tensor(q0, b->I_n, I1); 
                        for(int k=0; k<5; k++) {
                            b->I_n[k] = I1[k]; 
                        }     
                        b->I_n[5] = -(b->I_n[0]+b->I_n[3]);                          
                    }

                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->I_n[k] = b->I_e[k] + (b->I_n[k] - b->I_e_prev[k])*expo - b->tau*dIe*(1-expo);
                    }
                    b->I_n[5] = -(b->I_n[0]+b->I_n[3]);            

                    if(dt > 0) {
                        array<double, 5> I2 = {};
                        rotate_tensor(q0, b->I_n, I2); 
                        for(int k=0; k<5; k++) {
                            b->I_n[k] = I2[k]; 
                        }     
                        b->I_n[5] = -(b->I_n[0]+b->I_n[3]);                          
                    }                
                }            
            }
            else {
                for(int k=0; k<5; k++) {
                    b->I_n[k] = b->I_e[k];
                }
                b->I_n[5] = -(b->I_n[0]+b->I_n[3]);
            }

            for(int k=0; k<6; k++) {
                b->I[k] = b->I_p[k] + b->I_n[k];
            }
        }
    }
}
void Shape::deform_and_rotate_half(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            if(b->tau > 0) {
                double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);
                double w  = sqrt(w2);
            
                double dt_step = b->dth_rot; 

                double wdt = 0.5*w*dt_step;
                double cosa = cos(wdt);
                double sina = sin(wdt);
            
                array<double, 4> q0 = {};
                q0[0] = 1.;

                if(w > 0) {
                    double w_1 = 1/w;

                    q0[0] = cosa;
                    q0[1] = sina*b->w[0]*w_1;
                    q0[2] = sina*b->w[1]*w_1;
                    q0[3] = sina*b->w[2]*w_1;
                }

                if(b->isLinear) {
                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->I_n[k] = b->I_e[k] - b->tau*dIe;
                    }
                    b->I_n[5] = -(b->I_n[0]+b->I_n[3]);

                    array<double, 5> I2 = {};
                    rotate_tensor(q0, b->I_n, I2); 
                    for(int k=0; k<5; k++) {
                        b->I_n[k] = I2[k]; 
                    }     
                    b->I_n[5] = -(b->I_n[0]+b->I_n[3]);                      
                }
                else {
                    double expo = b->expoh;

                    if(dt < 0) {
                        array<double, 5> I1 = {};
                        rotate_tensor(q0, b->I_n, I1); 
                        for(int k=0; k<5; k++) {
                            b->I_n[k] = I1[k]; 
                        }     
                        b->I_n[5] = -(b->I_n[0]+b->I_n[3]);                          
                    }

                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->I_n[k] = b->I_e[k] + (b->I_n[k] - b->I_e_prev[k])*expo - b->tau*dIe*(1-expo);
                    }
                    b->I_n[5] = -(b->I_n[0]+b->I_n[3]);            

                    if(dt > 0) {
                        array<double, 5> I2 = {};
                        rotate_tensor(q0, b->I_n, I2); 
                        for(int k=0; k<5; k++) {
                            b->I_n[k] = I2[k]; 
                        }     
                        b->I_n[5] = -(b->I_n[0]+b->I_n[3]);                          
                    }                
                }            
            }
            else {
                for(int k=0; k<5; k++) {
                    b->I_n[k] = b->I_e[k];
                }
                b->I_n[5] = -(b->I_n[0]+b->I_n[3]);
            }

            for(int k=0; k<6; k++) {
                b->I[k] = b->I_p[k] + b->I_n[k];
            }
        }
    }
}

void Shape::deform_and_rotate_J(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            if(b->tau > 0) {
                double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);
                double w  = sqrt(w2);

                double dt_step = b->dt_rot; 

                double wdt = 0.5*w*dt_step;
                double cosa = cos(wdt);
                double sina = sin(wdt);
            
                array<double, 4> q0 = {};
                q0[0] = 1.;

                if(w > 0) {
                    double w_1 = 1/w;

                    q0[0] = cosa;
                    q0[1] = sina*b->w[0]*w_1;
                    q0[2] = sina*b->w[1]*w_1;
                    q0[3] = sina*b->w[2]*w_1;
                }

                if(b->isLinear) {
                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->J_n[k] = b->I_e[k] - b->tau*dIe;
                    }
                    b->J_n[5] = -(b->J_n[0]+b->J_n[3]);

                    array<double, 5> I2 = {};
                    rotate_tensor(q0, b->J_n, I2); 
                    for(int k=0; k<5; k++) {
                        b->J_n[k] = I2[k]; 
                    }     
                    b->J_n[5] = -(b->J_n[0]+b->J_n[3]);                      
                }
                else {
                    double expo = b->expo;

                    if(dt < 0) {
                        array<double, 5> I1 = {};
                        rotate_tensor(q0, b->J_n, I1); 
                        for(int k=0; k<5; k++) {
                            b->J_n[k] = I1[k]; 
                        }     
                        b->J_n[5] = -(b->J_n[0]+b->J_n[3]);                          
                    }

                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->J_n[k] = b->I_e[k] + (b->J_n[k] - b->I_e_prev[k])*expo - b->tau*dIe*(1-expo);
                    }
                    b->J_n[5] = -(b->J_n[0]+b->J_n[3]);            

                    if(dt > 0) {
                        array<double, 5> I2 = {};
                        rotate_tensor(q0, b->J_n, I2); 
                        for(int k=0; k<5; k++) {
                            b->J_n[k] = I2[k]; 
                        }     
                        b->J_n[5] = -(b->J_n[0]+b->J_n[3]);                          
                    }                
                }            
            }
            else {
                for(int k=0; k<5; k++) {
                    b->J_n[k] = b->I_e[k];
                }
                b->J_n[5] = -(b->J_n[0]+b->J_n[3]);
            }

            for(int k=0; k<6; k++) {
                b->J[k] = b->I_p[k] + b->J_n[k];
            }
        }
    }
}
void Shape::deform_and_rotate_J_half(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            if(b->tau > 0) {
                double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);
                double w  = sqrt(w2);

                double dt_step = b->dth_rot; 

                double wdt = 0.5*w*dt_step;
                double cosa = cos(wdt);
                double sina = sin(wdt);
            
                array<double, 4> q0 = {};
                q0[0] = 1.;

                if(w > 0) {
                    double w_1 = 1/w;

                    q0[0] = cosa;
                    q0[1] = sina*b->w[0]*w_1;
                    q0[2] = sina*b->w[1]*w_1;
                    q0[3] = sina*b->w[2]*w_1;
                }

                if(b->isLinear) {
                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->J_n[k] = b->I_e[k] - b->tau*dIe;
                    }
                    b->J_n[5] = -(b->J_n[0]+b->J_n[3]);

                    array<double, 5> I2 = {};
                    rotate_tensor(q0, b->J_n, I2); 
                    for(int k=0; k<5; k++) {
                        b->J_n[k] = I2[k]; 
                    }     
                    b->J_n[5] = -(b->J_n[0]+b->J_n[3]);                      
                }
                else {
                    double expo = b->expoh;

                    if(dt < 0) {
                        array<double, 5> I1 = {};
                        rotate_tensor(q0, b->J_n, I1); 
                        for(int k=0; k<5; k++) {
                            b->J_n[k] = I1[k]; 
                        }     
                        b->J_n[5] = -(b->J_n[0]+b->J_n[3]);                          
                    }

                    for(int k=0; k<5; k++) {
                        double dIe = (b->I_e[k]-b->I_e_prev[k])/dt;   
                        b->J_n[k] = b->I_e[k] + (b->J_n[k] - b->I_e_prev[k])*expo - b->tau*dIe*(1-expo);
                    }
                    b->J_n[5] = -(b->J_n[0]+b->J_n[3]);            

                    if(dt > 0) {
                        array<double, 5> I2 = {};
                        rotate_tensor(q0, b->J_n, I2); 
                        for(int k=0; k<5; k++) {
                            b->J_n[k] = I2[k]; 
                        }     
                        b->J_n[5] = -(b->J_n[0]+b->J_n[3]);                          
                    }                
                }            
            }
            else {
                for(int k=0; k<5; k++) {
                    b->J_n[k] = b->I_e[k];
                }
                b->J_n[5] = -(b->J_n[0]+b->J_n[3]);
            }

            for(int k=0; k<6; k++) {
                b->J[k] = b->I_p[k] + b->J_n[k];
            }
        }
    }
}

void Shape::rotate_tensor(array<double, 4> &q, array<double, 6> &A, array<double, 5> &B) {
    // Quaternion to Rotation matrix
    double q00 = q[0] * q[0];
    double q11 = q[1] * q[1];
    double q22 = q[2] * q[2];
    double q33 = q[3] * q[3];
    double q01 = q[0] * q[1];
    double q02 = q[0] * q[2];
    double q03 = q[0] * q[3];
    double q12 = q[1] * q[2];
    double q13 = q[1] * q[3];
    double q23 = q[2] * q[3];

    double r00 = 2 * (q00 + q11) - 1;
    double r01 = 2 * (q12 - q03);
    double r02 = 2 * (q13 + q02);
     
    double r10 = 2 * (q12 + q03);
    double r11 = 2 * (q00 + q22) - 1;
    double r12 = 2 * (q23 - q01);
     
    double r20 = 2 * (q13 - q02);
    double r21 = 2 * (q23 + q01);
    double r22 = 2 * (q00 + q33) - 1;
        
    // RA
    double x00 = r00*A[0]+r01*A[1]+r02*A[2];
    double x01 = r00*A[1]+r01*A[3]+r02*A[4];
    double x02 = r00*A[2]+r01*A[4]+r02*A[5];

    double x10 = r10*A[0]+r11*A[1]+r12*A[2];
    double x11 = r10*A[1]+r11*A[3]+r12*A[4];
    double x12 = r10*A[2]+r11*A[4]+r12*A[5];
    
    // RAR^T
    B[0] = x00*r00+x01*r01+x02*r02;
    B[1] = x00*r10+x01*r11+x02*r12;
    B[2] = x00*r20+x01*r21+x02*r22;
    B[3] = x10*r10+x11*r11+x12*r12;
    B[4] = x10*r20+x11*r21+x12*r22;
}

void Shape::copy_I_e_r_to_I_e_rh(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->I_e_rh[k] = b->I_e_r[k];
        }    
    }
}

void Shape::copy_I_to_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->J_n[k] = b->I_n[k];
            b->J[k] = b->I[k];
        }    
    }
}
        
void Shape::memorize_I(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->I_n_prev[k] = b->I_n[k];
        }    
    }
}
void Shape::reset_I(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->I_n[k] = b->I_n_prev[k];
            b->I[k] = b->I_p[k] + b->I_n[k];
        }    
    }
}

void Shape::memorize_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->I_n_prev[k] = b->J_n[k];
        }    
    }
}
void Shape::reset_J(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->J_n[k] = b->I_n_prev[k];
            b->J[k] = b->I_p[k] + b->J_n[k];
        }    
    }
}

void Shape::backup_I_e(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->I_e_prev_bu[k] = b->I_e[k];
        }    
    }
}
void Shape::set_I_e_prev_to_backup(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<6; k++) {
            b->I_e_prev[k] = b->I_e_prev_bu[k];
        }    
    }
}




        
        
        
