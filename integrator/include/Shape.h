#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <cstdlib>
#include <numeric>

#include "Body.h"

#ifndef __Shape_h
#define __Shape_h

class Shape {
    double f_lim;

    public:

    // Initializers
    Shape();

    void set_permanent_shape(vector<Body> &bodies);
    void set_to_spherical_shape(vector<Body> &bodies);
    void set_to_equilibrium_shape(vector<Body> &bodies);
    void set_to_equilibrium_shape_with_J(vector<Body> &bodies);
    void set_to_linear_shape(vector<Body> &bodies);
    void init_J(vector<Body> &bodies);
    
    void calc_expos(vector<Body> &bodies, double dt);    
    void copy_expos_to_1(vector<Body> &bodies);
    void copy_expos_to_2(vector<Body> &bodies);
    void revert_expos_1(vector<Body> &bodies);
    void revert_expos_2(vector<Body> &bodies);
        
    // Calculators
    void calculate_I_inv(vector<Body> &bodies);
    void calculate_I_inv_and_w(vector<Body> &bodies);

    void calculate_J_inv_and_w(vector<Body> &bodies);
    void calculate_w_with_J(vector<Body> &bodies);

    void calculate_I_inv_and_w_with_K(vector<Body> &bodies);
    void calculate_J_inv_and_w_with_K(vector<Body> &bodies);

    // Evolvers
    void kick_linear(vector<Body> &bodies);
    void kick_linear_J(vector<Body> &bodies);
    void kick_linear_discrete(vector<Body> &bodies, double dt);
    void kick_linear_discrete_J(vector<Body> &bodies, double dt);

    void kick_direct(vector<Body> &bodies, double dt);
    void kick_direct_J(vector<Body> &bodies, double dt);

    void deform_and_rotate(vector<Body> &bodies, double dt);
    void deform_and_rotate_half(vector<Body> &bodies, double dt);

    void deform_and_rotate_J(vector<Body> &bodies, double dt);
    void deform_and_rotate_J_half(vector<Body> &bodies, double dt);

    void rotate_tensor(array<double, 4> &q, array<double, 6> &A, array<double, 5> &B);
                
    void copy_I_e_r_to_I_e_rh(vector<Body> &bodies);
    void copy_I_to_J(vector<Body> &bodies);
    
    void memorize_I(vector<Body> &bodies);
    void reset_I(vector<Body> &bodies);

    void memorize_J(vector<Body> &bodies);
    void reset_J(vector<Body> &bodies);

    void backup_I_e(vector<Body> &bodies);
    void set_I_e_prev_to_backup(vector<Body> &bodies);
};

#endif


