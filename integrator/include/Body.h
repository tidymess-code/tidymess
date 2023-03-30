#include <cmath>
#include <array>
#include <string>
using namespace std;

#ifndef __Body_h
#define __Body_h

class Body {
    public:

    // Initializers
    Body();

    Body(double m, double R, double xi, double kf, double tau, double a_mb, double wx, double wy, double wz, double x, double y, double z, double vx, double vy, double vz);
    Body(double m, double R, double xi, double kf, double tau, double a_mb, array<double, 3> w, array<double, 3> r, array<double, 3> v);
    Body(array<double, 15> &d);

    void setup(array<double, 15> &d);
    void set_id(int id);
    void set_name(string name);

    void reset();
    void update_aux_properties();

    // Tag
    int id;
    string name;

    int particle_type; // 0 = test-particle, 1 = point-particle, 2 = rigid sphere, 3 = deformable ellipsoid

    // Internal properties
    double m, R, xi;
    double kf, tau;
    double a_mb;

    double R3, R5, R5_3, kf_R5, kf_R5_3;
    double rho, tau_inv;
    
    double roche_factor, roche_factor3;
    
    // Orbital properties
    array<double, 3> r, v;
    array<double, 3> a;
    
    array<double, 3> r_prev;
    array<double, 3> vv;
    
    // Spin properties
    array<double, 3> L, w;
    array<double, 3> T;

    array<double, 3> L_prev, w_prev;
    array<double, 3> K;
        
    // Shape properties
    array<double, 6> I_p, I_n, I, I_inv;
    array<double, 6> I_e_r, dI_e_r, I_e_w, I_e, dI_e, dI_n;

    array<double, 6> I_e_rh, I_e_prev, I_n_prev, I_e_prev_bu;
    array<double, 6> J_n, J, J_inv;

    // Integration variables
    bool isLinear;
    double expo, expoh, dt_rot, dth_rot;

    bool isLinear_1;
    double expo_1, expoh_1, dt_rot_1, dth_rot_1;

    bool isLinear_2;
    double expo_2, expoh_2, dt_rot_2, dth_rot_2;
};

#endif


