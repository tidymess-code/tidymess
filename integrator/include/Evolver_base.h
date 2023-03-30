#include "Force.h"

#include "Timestep_base.h"

#include "Timestep_const.h"
#include "Timestep_adapt.h"
#include "Timestep_direct.h"
#include "Timestep_adapt_weight.h"
#include "Timestep_direct_weight.h"

#include "Orbit.h"
#include "Spin.h"
#include "Shape.h"

#ifndef __Evolver_base_h
#define __Evolver_base_h

class Evolver_base {
    public:

    // Variables and Objects    
    Orbit orbit;
    Spin spin;
    Shape shape;

    double dt_step;
    double w1, w2, w3;
    int n_iter;
    int num_integration_step;
  
    int collision_mode, roche_mode;
    bool collision_flag, roche_flag;

    // Initializers
    Evolver_base();

    void assign_vectors(vector<Body> &bodies, Force *force);
    
    void set_n_iter(int n_iter);
    int get_n_iter();

    void reset_aux(vector<Body> &bodies);

    void set_collision_mode(int collision_mode);
    void set_roche_mode(int roche_mode);
 
    // Get diagnostics
    array<double, 3> get_center_of_mass(vector<Body> &bodies);
    array<double, 3> get_center_of_mass_velocity(vector<Body> &bodies);
    array<double, 3> get_orbital_angular_momentum(vector<Body> &bodies);
    double get_orbital_kinetic_energy(vector<Body> &bodies);

    virtual array<double, 3> get_spin_angular_momentum(vector<Body> &bodies) = 0;
    virtual double get_spin_kinetic_energy(vector<Body> &bodies) = 0;
    virtual double get_potential_energy(vector<Body> &bodies) = 0;
    virtual array<double, 3> get_angular_momentum(vector<Body> &bodies) = 0;
    virtual double get_energy(vector<Body> &bodies) = 0;
    
    // Evolve function
    virtual void initialize(vector<Body> &bodies, Force *force) = 0;

    virtual void evolve_model(double &t, double &t_end, vector<Body> &bodies, Force *force, Timestep_base *timestep) = 0;
    virtual void evolve_model(double &t, int &N_step, vector<Body> &bodies, Force *force, Timestep_base *timestep) = 0;

    bool stopping_condition_time(double t0, double t1, int dt_sgn);
};

#endif


