#include <iostream>
using namespace std;

#include <array>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include "Body.h"

#ifndef __Timestep_base_h
#define __Timestep_base_h

class Timestep_base {
    public:

    // Variables
    double dt_const;
    double eta, eta2;

    double h, h2, dh, h_spin, h_shape;    
    double dt_prev, dt_next;
    double fp_bound, fm_bound, fp2_bound, fm2_bound;
    double alpha;
    int dt_sgn;
   
    // Initializer
    Timestep_base();

    void set_default_values();

    // Setters and Getters
    void set_dt_const(double dt_const);
    double get_dt_const();
    
    void set_eta(double eta);
    double get_eta();

    void commit();
    void set_dt_sgn(double t0, double t1);
    void set_dt_sgn(int dt_sgn);

    void set_dt_prev(double dt_prev);

    // Calculators
    virtual void initialize(vector<Body> &bodies) = 0;
    virtual double get_timestep(vector<Body> &bodies) = 0;
};

#endif



