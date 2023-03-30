#include "Timestep_direct.h"

// Calculators
void Timestep_direct::calculate_shared_adaptive_minimum_timestep_orbital(vector<Body> &bodies) {
    this->h = 0;

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            if(bi->particle_type > 0 || bj->particle_type > 0) {

                double mu = bi->m + bj->m;

                array<double, 3> dr, dv;

                for(int k=0; k<3; k++) {
                    dr[k] = bi->r[k] - bj->r[k];             
                    dv[k] = bi->v[k] - bj->v[k];             
                }

                double dr2 = inner_product(dr.begin(), dr.end(), dr.begin(), 0.);
                double dv2 = inner_product(dv.begin(), dv.end(), dv.begin(), 0.);
        
                double dr4 = dr2*dr2;
        
                double dt1_4 = mu*mu/(dr4*dr2);
                double dt2_4 = dv2*dv2/dr4;
        
                double P_4 = max(dt1_4, dt2_4);
        
                if(P_4 > this->h) this->h = P_4;        
            }
        }
    }

    this->h = 1./sqrt(this->h);
}

void Timestep_direct::calculate_minimum_timestep_spin(vector<Body> &bodies) {
    h_spin = 0.;

    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);
            if(w2 > h_spin) {
                h_spin = w2;
            }
        }
    }
    
    if(h_spin == 0) {
        h_spin = 1e100;
    }
    else {
        h_spin = 1/h_spin;
    }
}
void Timestep_direct::calculate_minimum_timestep_shape(vector<Body> &bodies) {
    h_shape = 1e100;

    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            if(b->tau > 0) {
                if(b->tau < h_shape) {
                    h_shape = b->tau;
                }
            }
        }
    }  
    
    h_shape = h_shape*h_shape;  
}

void Timestep_direct::initialize(vector<Body> &bodies) {
    ;
}
double Timestep_direct::get_timestep(vector<Body> &bodies) {
    calculate_shared_adaptive_minimum_timestep_orbital(bodies);
    calculate_minimum_timestep_spin(bodies);
    calculate_minimum_timestep_shape(bodies);
    
    h = min(min(h_spin, h_shape), h);
    h = sqrt(h);
    
    h *= eta;
    h *= dt_sgn;

    return h;
}

    
    


