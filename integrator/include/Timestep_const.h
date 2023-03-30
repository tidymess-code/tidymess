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

#ifndef __Timestep_const_h
#define __Timestep_const_h

class Timestep_const : public Timestep_base {
    public:

    // Calculators
    void initialize(vector<Body> &bodies);
    double get_timestep(vector<Body> &bodies);
};

#endif



