#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric>

#include "Body.h"

#ifndef __Spin_h
#define __Spin_h

class Spin {
    public:

    // Initializers
    Spin();

    // Calculators
    void calculate_w(vector<Body> &bodies);
    void calculate_w_with_J(vector<Body> &bodies);
    void calculate_w_with_K(vector<Body> &bodies);

    void calculate_L(vector<Body> &bodies);
    void init_K(vector<Body> &bodies);    

    void memorize_w(vector<Body> &bodies);
    void reset_w(vector<Body> &bodies);
    
    void kick_L(vector<Body> &bodies, double dt);
    void kick_K(vector<Body> &bodies, double dt);

    void kick_L_mb(vector<Body> &bodies, double dt);
    void kick_K_mb(vector<Body> &bodies, double dt);

    // Diagnostics
    double get_kinetic_energy(vector<Body> &bodies);
    array<double, 3> get_angular_momentum(vector<Body> &bodies);
};

#endif


