#include "Timestep_adapt_weight.h"

// Calculators
void Timestep_adapt_weight::calculate_shared_adaptive_weighted_timestep_orbital(vector<Body> &bodies) {
    double wsum = 0;
    double wh2sum = 0;
      
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
            
                //---------------------------------------------------------------
            
                double P_4  = dt1_4 + dt2_4;
                double P_2  = sqrt(P_4);
                double P_8  = P_4*P_4;
                double P_10 = P_8*P_2;

                //---------------------------------------------------------------

                double wh2 = P_8;
                double w   = P_10;
            
                //---------------------------------------------------------------

                wsum   += w;
                wh2sum += wh2;   
            }
        }
    }
    
    this->h2 = eta2 * wh2sum/wsum;
}
void Timestep_adapt_weight::calculate_shared_adaptive_weighted_timestep_orbital_and_derivative(vector<Body> &bodies) {
    double h_8  = 0; 
    double h_10 = 0; 
    double h_9_dh  = 0; 
    double h_11_dh = 0; 

    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            if(bi->particle_type > 0 || bj->particle_type > 0) {
            
                double mu = bi->m + bj->m;

                array<double, 3> dr, dv, da;

                for(int k=0; k<3; k++) {
                    dr[k] = bj->r[k] - bi->r[k];             
                    dv[k] = bj->v[k] - bi->v[k];
                    da[k] = bj->a[k] - bi->a[k];             
                }

                double dr2 = inner_product(dr.begin(), dr.end(), dr.begin(), 0.);
                double dv2 = inner_product(dv.begin(), dv.end(), dv.begin(), 0.);

                double rdotv = inner_product(dr.begin(), dr.end(), dv.begin(), 0.);
                double vdota = inner_product(dv.begin(), dv.end(), da.begin(), 0.);
                        
                //------------------------------------------

                double h_4_ra = mu*mu/(dr2*dr2*dr2);

                double h_1_dh_ra = 1.5*rdotv/dr2;
                double h_5_dh_ra  = h_4_ra * h_1_dh_ra;

                //------------------------------------------

                double h_2_rv  = dv2/dr2;
                double h_4_rv  = h_2_rv*h_2_rv;

                double h_3_dh_rv = (rdotv * h_2_rv - vdota) / dr2;
                double h_5_dh_rv  = h_2_rv * h_3_dh_rv;

                //------------------------------------------
      
                double T_4 = h_4_ra + h_4_rv;

                double T_8 = T_4*T_4;
                double T_2 = sqrt(T_4);
                double T_10 = T_8*T_2;
            
                double T_5_dT = h_5_dh_ra + h_5_dh_rv;      
            
                double T_9_dT  = T_4 * T_5_dT;      
                double T_11_dT = T_2 * T_9_dT;      
            
                //------------------------------------------
      
                h_8  += T_8;
                h_10 += T_10;            

                h_9_dh  += T_9_dT;
                h_11_dh += T_11_dT;            
            }
        }
    }    
  
    double T  = sqrt(h_8 / h_10);  
    double dT = 0.5 / T * ( (2-alpha)*h_9_dh/h_10 - h_8/(h_10*h_10)*-alpha*h_11_dh );
  
    T *= eta;
    dT *= eta;
    
    this->h = T;
    this->dh = dT;
}

void Timestep_adapt_weight::initialize(vector<Body> &bodies) {
    calculate_shared_adaptive_weighted_timestep_orbital_and_derivative(bodies); 

    if(dt_sgn > 0) {
        this->dt_prev = this->h / (1+0.5*this->dh);
    }
    else {
        this->dt_prev = this->h / (1-0.5*this->dh);    
    }

    this->dt_prev *= dt_sgn;
}
double Timestep_adapt_weight::get_timestep(vector<Body> &bodies) {
    calculate_shared_adaptive_weighted_timestep_orbital(bodies);

    this->dt_next = this->h2 / this->dt_prev;

    double f2 = this->dt_next*this->dt_next / this->h2;
    if(f2 > this->fp2_bound) this->dt_next = dt_sgn * sqrt(this->h2);
    else if(f2 < this->fm2_bound) this->dt_next = dt_sgn * sqrt(this->h2);

    this->dt_prev = this->dt_next;

    return this->dt_next;
}


