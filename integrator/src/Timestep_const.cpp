#include "Timestep_const.h"

// Calculators
void Timestep_const::initialize(vector<Body> &bodies) {
    ; 
}
double Timestep_const::get_timestep(vector<Body> &bodies) {
    return dt_sgn*dt_const;
}

    
    


