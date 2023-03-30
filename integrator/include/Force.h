#include <iostream>
using namespace std;

#include <vector>
#include <array>
#include <cmath>
#include <numeric>

#include "Body.h"

#ifndef __Force_h
#define __Force_h

class Force {
    public:

    // Variables
    int pn_order;
    double c, c_2, c_4, c_5;

    bool collision_detected;
    vector< array<int, 2> > index_collisions;

    bool roche_detected;
    vector< array<int, 2> > index_roche;

    vector<double> dr2_vec, dr_3_vec;

    // Initializers
    Force();
    Force(int pn_order, double c);

    void set_default_settings();

    // Setters
    void set_pn_order(int pn_order);
    void set_speed_of_light(double c);

    // Collision checker
    void check_for_collisions(Body &bi, Body &bj, double &dr_1);    
    void check_for_collisions(Body &bi, Body &bj, double &dr_3, double &dr2);

    void check_for_roche(Body &bi, Body &bj, double &dr_1);    
    void check_for_roche(Body &bi, Body &bj, double &dr_3, double &dr2);

    void reset_collisions();
    void reset_roche();    
        
    // Precalculation functions

    void precalc_sep_vec(Body &bi, Body &bj, array<double, 3> &dr);
    void precalc_sep_vec_norm(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &drn, double &dr_4, double &dr_3);
    void precalc_sep_vec_norm_load(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &drn, double &dr_4, double &dr2, double &dr_3);
    void precalc_sep_vec_norm_save(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &drn, double &dr_4, double &dr2, double &dr_3);
    void precalc_sep_mag(Body &bi, Body &bj, array<double, 3> &dr, double &dr2, double &dr_3);

    void precalc_pn(Body &bi, Body &bj, bool &use_vv, array<double, 3> &drn, double &dr_1, double &vivi, double &vjvj, double &vivj, double &drnvi, double &drnvj, double &dv2, array<double, 3> &dv, double &drnvidrnvi, double &drnvjdrnvj, double &bimdr_1, double &bjmdr_1);

    void precalc_tidal(double &dr2, double &dr_3, double &dr_1, double &dr_2, double &dr_4, double &dr_5, array<double, 3> &dr, array<double, 3> &drn, double &dyndyn, double &dzndzn, double &dxndyn, double &dxndzn, double &dyndzn);

    // Force functions

    void calculate_nbody_force(Body &bi, Body &bj, array<double, 3> &dr, double &dr2, double &dr_3, array<double, 3> &ai, array<double, 3> &aj);
    void calculate_nbody_force_load(Body &bi, Body &bj, array<double, 3> &dr, double &dr2, double &dr_3, array<double, 3> &ai, array<double, 3> &aj);

    void calculate_pn_force(Body &bi, Body &bj, double &dr2, double &dr_1, double &dr_2, double &dr_3, double &dr_4, array<double, 3> &dr, array<double, 3> &drn, double &dv2, array<double, 3> &dv, array<double, 3> &ai_pn, array<double, 3> &aj_pn, bool use_vv);

    void calculate_tidal_force(Body &bi, Body &bj, double &dr2, double &dr_3, double &dr_1, double &dr_2, double &dr_4, double &dr_5, array<double, 3> &dr, array<double, 3> &drn, array<double, 3> &ai_tidal, array<double, 3> &aj_tidal, array<double, 3> &hij, array<double, 3> &hji, bool use_J);

    void calculate_tidal_h(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji, bool use_J);
    void calculate_tidal_h_load(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji, bool use_J, double &dr2, double &dr_3);
    void calculate_tidal_h_save(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji, bool use_J, double &dr2, double &dr_3);
    
    // Torque functions
    
    void calculate_torque(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji); 

    // Inertia tensor functions

    void calculate_deformation_tensor(Body &bi, Body &bj, array<double, 3> &drn, double &dr_3);

    // Acceleration updaters
    
    void update_acceleration(vector<Body> &bodies);
    void update_acceleration_pn(vector<Body> &bodies, bool use_vv);
    void update_acceleration_tidal(vector<Body> &bodies, bool use_J);
    void update_acceleration_col(vector<Body> &bodies);
    void update_acceleration_pn_col(vector<Body> &bodies, bool use_vv);
    void update_acceleration_tidal_col(vector<Body> &bodies, bool use_J);
    void update_acceleration_tidal_pn(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_tidal_pn_col(vector<Body> &bodies, bool use_vv, bool use_J);
            
    void update_acceleration_save(vector<Body> &bodies);
    void update_acceleration_pn_save(vector<Body> &bodies, bool use_vv);
    void update_acceleration_tidal_save(vector<Body> &bodies, bool use_J);
    void update_acceleration_col_save(vector<Body> &bodies);
    void update_acceleration_pn_col_save(vector<Body> &bodies, bool use_vv);
    void update_acceleration_tidal_col_save(vector<Body> &bodies, bool use_J);
    void update_acceleration_tidal_pn_save(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_tidal_pn_col_save(vector<Body> &bodies, bool use_vv, bool use_J);
    
    void update_acceleration_load(vector<Body> &bodies);
    void update_acceleration_pn_load(vector<Body> &bodies, bool use_vv);
    void update_acceleration_tidal_load(vector<Body> &bodies, bool use_J);
    void update_acceleration_col_load(vector<Body> &bodies);
    void update_acceleration_pn_col_load(vector<Body> &bodies, bool use_vv);
    void update_acceleration_tidal_col_load(vector<Body> &bodies, bool use_J);
    void update_acceleration_tidal_pn_load(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_tidal_pn_col_load(vector<Body> &bodies, bool use_vv, bool use_J);

    // Acceleration and torque updaters

    void update_acceleration_and_torque(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_and_torque_col(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn_col(vector<Body> &bodies, bool use_vv, bool use_J);

    void update_acceleration_and_torque_save(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn_save(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_and_torque_col_save(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn_col_save(vector<Body> &bodies, bool use_vv, bool use_J);

    void update_acceleration_and_torque_load(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn_load(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_and_torque_col_load(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn_col_load(vector<Body> &bodies, bool use_vv, bool use_J);

    void update_acceleration_and_torque_and_deformation_tensor(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_and_deformation_tensor_col(vector<Body> &bodies, bool use_J);
    void update_acceleration_and_torque_pn_and_deformation_tensor(vector<Body> &bodies, bool use_vv, bool use_J);
    void update_acceleration_and_torque_pn_and_deformation_tensor_col(vector<Body> &bodies, bool use_vv, bool use_J);    
    
    // Torque updaters
    
    void update_torque(vector<Body> &bodies, bool use_J);
    void update_torque_save(vector<Body> &bodies, bool use_J);
    void update_torque_load(vector<Body> &bodies, bool use_J);

    // Deformation updaters

    void update_tidal_deformation(vector<Body> &bodies);
    void update_tidal_deformation_save(vector<Body> &bodies);
    void update_tidal_deformation_load(vector<Body> &bodies);

    void update_centrifugal_deformation(vector<Body> &bodies);    
    
    void update_equilibrium_tensor(vector<Body> &bodies);
    void update_equilibrium_tensor_and_memorize(vector<Body> &bodies);

    void update_deformation_tensor(vector<Body> &bodies, bool use_J);    
};

#endif

