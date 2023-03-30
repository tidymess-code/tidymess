#include "Timestep_adapt.h"

// Calculators
void Timestep_adapt::calculate_shared_adaptive_minimum_timestep_orbital(vector<Body> &bodies) {
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

    this->h = 1./sqrt(sqrt(this->h));
}

void Timestep_adapt::initialize(vector<Body> &bodies) {
    ;
}
double Timestep_adapt::get_timestep(vector<Body> &bodies) {
    calculate_shared_adaptive_minimum_timestep_orbital(bodies);
    h *= eta;
    h *= dt_sgn;
    return h;
}

    
    


