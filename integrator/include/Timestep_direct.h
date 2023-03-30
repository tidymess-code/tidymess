#include <iostream>
using namespace std;

#include <array>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include "Body.h"
#include "Timestep_base.h"

#ifndef __Timestep_direct_h
#define __Timestep_direct_h

class Timestep_direct : public Timestep_base {
    public:

    // Calculators
    void calculate_shared_adaptive_minimum_timestep_orbital(vector<Body> &bodies);
    
    void calculate_minimum_timestep_spin(vector<Body> &bodies);
    void calculate_minimum_timestep_shape(vector<Body> &bodies);
    
    void initialize(vector<Body> &bodies);
    double get_timestep(vector<Body> &bodies);
};

#endif



