#include "Force.h"

#include "Orbit.h"
#include "Spin.h"
#include "Shape.h"

#include "Timestep_base.h"

#include "Timestep_const.h"
#include "Timestep_adapt.h"
#include "Timestep_direct.h"
#include "Timestep_adapt_weight.h"
#include "Timestep_direct_weight.h"

#include "Evolver_base.h"

#include "Evolver_nbody.h"
#include "Evolver_nbody_col.h"
#include "Evolver_nbody_pn.h"
#include "Evolver_nbody_col_pn.h"

#include "Evolver_equilibrium.h"
#include "Evolver_equilibrium_col.h"
#include "Evolver_equilibrium_pn.h"
#include "Evolver_equilibrium_mb.h"
#include "Evolver_equilibrium_col_pn.h"
#include "Evolver_equilibrium_col_mb.h"
#include "Evolver_equilibrium_pn_mb.h"
#include "Evolver_equilibrium_col_pn_mb.h"

#include "Evolver_linear.h"
#include "Evolver_linear_col.h"
#include "Evolver_linear_pn.h"
#include "Evolver_linear_mb.h"
#include "Evolver_linear_col_pn.h"
#include "Evolver_linear_col_mb.h"
#include "Evolver_linear_pn_mb.h"
#include "Evolver_linear_col_pn_mb.h"

#include "Evolver_direct.h"
#include "Evolver_direct_col.h"
#include "Evolver_direct_pn.h"
#include "Evolver_direct_mb.h"
#include "Evolver_direct_col_pn.h"
#include "Evolver_direct_col_mb.h"
#include "Evolver_direct_pn_mb.h"
#include "Evolver_direct_col_pn_mb.h"

#include "Evolver_creep.h"
#include "Evolver_creep_col.h"
#include "Evolver_creep_pn.h"
#include "Evolver_creep_mb.h"
#include "Evolver_creep_col_pn.h"
#include "Evolver_creep_col_mb.h"
#include "Evolver_creep_pn_mb.h"
#include "Evolver_creep_col_pn_mb.h"

#ifndef __Tidy_h
#define __Tidy_h

class Tidy {
    public:

    // User parameters
    double t;
    vector<Body> bodies;

    int tidal_model;

    int collision_mode, roche_mode, encounter_mode;

    int pn_order;
    double speed_of_light;

    int magnetic_braking;

    int dt_mode;
    double dt_const, eta, dt_prev;
    int num_integration_step;

    int n_iter;
    int dt_sgn;

    // Setters and Getters
    void set_model_time(double t);
    double get_model_time();

    void set_particles(vector<Body> &bodies);
    void set_particles(vector< array<double, 15> > &d);
    void set_particles(vector<string> &name, vector<int> &id, vector< array<double, 15> > &d);
    vector<Body> get_particles();
    
    void set_tidal_model(int tidal_model);
    int get_tidal_model();
   
    void set_collision_mode(int collision_mode);
    int get_collision_mode();

    void set_roche_mode(int roche_mode);
    int get_roche_mode();

    void set_encounter_mode();
    int get_encounter_mode();
    
    void set_pn_order(int pn_order);
    int get_pn_order();

    void set_speed_of_light(double c);
    double get_speed_of_light();
    
    void set_magnetic_braking(int b);
    int get_magnetic_braking();
    
    void set_dt_mode(int dt_mode);
    double get_dt_mode();

    void set_dt_const(double dt_const);
    double get_dt_const();
    
    void set_eta(double eta);
    double get_eta();

    void set_dt_prev(double dt_prev);
    double get_dt_prev();

    void set_n_iter(int n_iter);
    int get_n_iter();

    void set_num_integration_step(int num_integration_step);
    int get_num_integration_step();

    void set_dt_sgn(int dt_sgn);
            
    // Constructor/Destructor
    Tidy();
    ~Tidy();

    // Timesteppers and integrators
    Timestep_const timestep_const;
    Timestep_adapt timestep_adapt;
    Timestep_direct timestep_direct;
    Timestep_adapt_weight timestep_adapt_weight;
    Timestep_direct_weight timestep_direct_weight;
    
    Evolver_nbody evolver_nbody;
    Evolver_nbody_col evolver_nbody_col;
    Evolver_nbody_pn evolver_nbody_pn;
    Evolver_nbody_col_pn evolver_nbody_col_pn;

    Evolver_equilibrium evolver_equilibrium;
    Evolver_equilibrium_col evolver_equilibrium_col;
    Evolver_equilibrium_pn evolver_equilibrium_pn;
    Evolver_equilibrium_mb evolver_equilibrium_mb;
    Evolver_equilibrium_col_pn evolver_equilibrium_col_pn;
    Evolver_equilibrium_col_mb evolver_equilibrium_col_mb;
    Evolver_equilibrium_pn_mb evolver_equilibrium_pn_mb;
    Evolver_equilibrium_col_pn_mb evolver_equilibrium_col_pn_mb;

    Evolver_linear evolver_linear;
    Evolver_linear_col evolver_linear_col;
    Evolver_linear_pn evolver_linear_pn;
    Evolver_linear_mb evolver_linear_mb;
    Evolver_linear_col_pn evolver_linear_col_pn;
    Evolver_linear_col_mb evolver_linear_col_mb;
    Evolver_linear_pn_mb evolver_linear_pn_mb;
    Evolver_linear_col_pn_mb evolver_linear_col_pn_mb;

    Evolver_direct evolver_direct;
    Evolver_direct_col evolver_direct_col;
    Evolver_direct_pn evolver_direct_pn;
    Evolver_direct_mb evolver_direct_mb;
    Evolver_direct_col_pn evolver_direct_col_pn;
    Evolver_direct_col_mb evolver_direct_col_mb;
    Evolver_direct_pn_mb evolver_direct_pn_mb;
    Evolver_direct_col_pn_mb evolver_direct_col_pn_mb;

    Evolver_creep evolver_creep;
    Evolver_creep_col evolver_creep_col;
    Evolver_creep_pn evolver_creep_pn;
    Evolver_creep_mb evolver_creep_mb;
    Evolver_creep_col_pn evolver_creep_col_pn;
    Evolver_creep_col_mb evolver_creep_col_mb;
    Evolver_creep_pn_mb evolver_creep_pn_mb;
    Evolver_creep_col_pn_mb evolver_creep_col_pn_mb;

    // Pointers to objects to be used
    Force *force_ptr;
    Timestep_base *timestep_ptr;
    Evolver_base *evolver_ptr;
     
    // Committers
    void commit_parameters();
    void commit_particles();
    void set_pointers();
    void upload_parameters();    
    void initialize();      
             
    // Collision handling
    bool is_collision_detected();
    vector< array<int, 2> > get_collision_indices();
    bool get_collision_flag();

    bool is_roche_detected();
    vector< array<int, 2> > get_roche_indices();
    bool get_roche_flag();

    // Shape handling
    Force force;
    Spin spin;
    Shape shape;

    void set_to_spherical_shape();
    void set_to_equilibrium_shape();
    void update_angular_momentum();

    // Get diagnostics
    array<double, 3> get_center_of_mass();
    array<double, 3> get_center_of_mass_velocity();
    array<double, 3> get_orbital_angular_momentum();
    array<double, 3> get_spin_angular_momentum();
    double get_orbital_kinetic_energy();
    double get_spin_kinetic_energy();
    double get_potential_energy();
    array<double, 3> get_angular_momentum();
    double get_energy();

    // Top level evolve function
    void evolve_model(double t_end);    
    void evolve_model(int N_step);    
};

#endif


