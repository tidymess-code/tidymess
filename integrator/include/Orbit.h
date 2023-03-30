#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric>

#include "Body.h"

#ifndef __Orbit_h
#define __Orbit_h

class Orbit {
    public:

    // Initializers
    Orbit();

    void drift_r(vector<Body> &bodies, double dt);
    
    void kick_v(vector<Body> &bodies, double dt);
    void kick_vv(vector<Body> &bodies, double dt);

    void memorize_r(vector<Body> &bodies);
    void swap_r_and_r_prev(vector<Body> &bodies);

    void init_vv(vector<Body> &bodies);
    void swap_v_and_vv(vector<Body> &bodies);
        
    // Diagnostics
    double get_kinetic_energy(vector<Body> &bodies);
    double get_potential_energy(vector<Body> &bodies);
    double get_potential_energy_nbody(vector<Body> &bodies);
    
    array<double, 3> get_center_of_mass(vector<Body> &bodies);
    array<double, 3> get_center_of_mass_velocity(vector<Body> &bodies);
    array<double, 3> get_angular_momentum(vector<Body> &bodies);
};

#endif


