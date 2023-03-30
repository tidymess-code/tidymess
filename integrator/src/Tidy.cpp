#include "Tidy.h"

// Setters and Getters
void Tidy::set_model_time(double t) {
    this->t = t;
}
double Tidy::get_model_time() {
    return t;
}

void Tidy::set_particles(vector<Body> &bodies) {
    this->bodies.clear();
    this->bodies = bodies;
}
void Tidy::set_particles(vector< array<double, 15> > &d) {
    this->bodies.clear();
    
    int N = d.size();
    this->bodies.resize(N);
    
    for(int i=0; i<N; i++) {
        bodies[i].setup(d[i]);
    }
}
void Tidy::set_particles(vector<string> &name, vector<int> &id, vector< array<double, 15> > &d) {
    this->bodies.clear();
    
    int N = d.size();
    this->bodies.resize(N);
    
    for(int i=0; i<N; i++) {
        bodies[i].setup(d[i]);
    }
    for(int i=0; i<N; i++) {
        bodies[i].name = name[i];
        bodies[i].id = id[i];
    }
}
vector<Body> Tidy::get_particles() {
    return this->bodies;
}
    
void Tidy::set_tidal_model(int tidal_model) {
    this->tidal_model = tidal_model;
}
int Tidy::get_tidal_model() {
    return tidal_model;
}

void Tidy::set_collision_mode(int collision_mode) {
    this->collision_mode = collision_mode;
}
int Tidy::get_collision_mode() {
    return collision_mode;
}

void Tidy::set_roche_mode(int roche_mode) {
    this->roche_mode = roche_mode;
}
int Tidy::get_roche_mode() {
    return roche_mode;
}
    
void Tidy::set_encounter_mode() {
    encounter_mode = 0;
    if(collision_mode > 0 || roche_mode > 0) encounter_mode = 1;
}    
int Tidy::get_encounter_mode() {
    return encounter_mode;
}
   
void Tidy::set_pn_order(int pn_order) {
    this->pn_order = pn_order;
}
int Tidy::get_pn_order() {
    return pn_order;
}

void Tidy::set_speed_of_light(double c) {
    this->speed_of_light = c;   
}
double Tidy::get_speed_of_light() {
    return speed_of_light;
}
    
void Tidy::set_magnetic_braking(int b) {
    this->magnetic_braking = b;
}
int Tidy::get_magnetic_braking() {
    return magnetic_braking;
}

void Tidy::set_dt_mode(int dt_mode) {
    this->dt_mode = dt_mode;
}
double Tidy::get_dt_mode() {
    return dt_mode;
}

void Tidy::set_dt_const(double dt_const) {
    this->dt_const = dt_const;
}
double Tidy::get_dt_const() {
    return dt_const;
}
    
void Tidy::set_eta(double eta) {
    this->eta = eta;
}
double Tidy::get_eta() {
    return eta;
}

void Tidy::set_dt_prev(double dt_prev) {
    timestep_ptr->set_dt_prev(dt_prev);
    this->dt_prev = timestep_ptr->dt_prev;
}
double Tidy::get_dt_prev() {
    this->dt_prev = timestep_ptr->dt_prev;
    return dt_prev;
}

void Tidy::set_n_iter(int n_iter) {
    this->n_iter = n_iter;
}
int Tidy::get_n_iter() {
    this->n_iter = evolver_ptr->n_iter;
    return this->n_iter;
}

void Tidy::set_num_integration_step(int num_integration_step) {
    this->num_integration_step = num_integration_step;
}
int Tidy::get_num_integration_step() {
    return num_integration_step;
}
   
void Tidy::set_dt_sgn(int dt_sgn) {
    this->dt_sgn = dt_sgn;
    timestep_ptr->set_dt_sgn(this->dt_sgn);    
}
       
// Initializers
Tidy::Tidy() {
    t = 0.;
    bodies.clear();
    
    tidal_model = 0;
    
    collision_mode = 0;
    roche_mode = 0;
    encounter_mode = 0;

    pn_order = 0;
    speed_of_light = 1e100;
        
    magnetic_braking = 0;
    
    dt_mode = 1;
    dt_const = 0.015625;
    eta = 0.0625;
    num_integration_step = 0;

    n_iter = 1;
}
Tidy::~Tidy() {
    ;
}  

// Committers
void Tidy::commit_parameters() {
    set_encounter_mode();
    set_pointers();
    upload_parameters();
    initialize();
    this->dt_prev = timestep_ptr->dt_prev;
}
void Tidy::commit_particles() {
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        bodies[i].update_aux_properties();
    }

    if(evolver_ptr && force_ptr) {
        evolver_ptr->initialize(bodies, force_ptr);
    }
    else {
        commit_parameters();
    }
}
void Tidy::set_pointers() {
    force_ptr = &force;
        
    switch(dt_mode) {
        case 0:
            timestep_ptr = &timestep_const;            
            break;
        case 1:
            if(this->tidal_model == 3) {
                timestep_ptr = &timestep_direct;                        
            }
            else {
                timestep_ptr = &timestep_adapt;            
            }
            break;
        default:
            if(this->tidal_model == 3) {
                timestep_ptr = &timestep_direct_weight;                        
            }
            else {
                timestep_ptr = &timestep_adapt_weight;            
            }
            break;
    }
    
    switch(tidal_model) {
        case 0:
           switch(encounter_mode) {
               case 0:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_nbody;
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_nbody;
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_nbody_pn;                               
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_nbody_pn;
                                   break;
                           }
                           break;
                   }
                   break;               
               default:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_nbody_col;                               
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_nbody_col;
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_nbody_col_pn;                               
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_nbody_col_pn;
                                   break;
                           }
                           break;
                   }
                   break;
           } 
           break;        
        case 1:
           switch(encounter_mode) {
               case 0:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_equilibrium;                               
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_equilibrium_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_equilibrium_pn;                                
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_equilibrium_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;               
               default:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_equilibrium_col; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_equilibrium_col_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_equilibrium_col_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_equilibrium_col_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;
           } 
           break;        
        case 2:
           switch(encounter_mode) {
               case 0:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_linear; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_linear_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_linear_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_linear_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;               
               default:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_linear_col; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_linear_col_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_linear_col_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_linear_col_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;
           } 
           break;        
        case 3:
           switch(encounter_mode) {
               case 0:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_direct; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_direct_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_direct_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_direct_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;               
               default:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_direct_col; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_direct_col_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_direct_col_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_direct_col_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;
           } 
           break;  
        default:
           switch(encounter_mode) {
               case 0:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_creep; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_creep_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_creep_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_creep_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;               
               default:
                   switch(pn_order) {
                       case 0:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_creep_col; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_creep_col_mb; 
                                   break;
                           }
                           break;                       
                       default:
                           switch(magnetic_braking) {
                               case 0:
                                   evolver_ptr = &evolver_creep_col_pn; 
                                   break;                               
                               default:
                                   evolver_ptr = &evolver_creep_col_pn_mb; 
                                   break;
                           }
                           break;
                   }
                   break;
           } 
           break;
    }
}

void Tidy::upload_parameters() {
    force_ptr->set_pn_order(pn_order);
    force_ptr->set_speed_of_light(speed_of_light);
    
    timestep_ptr->set_dt_const(dt_const);
    timestep_ptr->set_eta(eta);   
    
    evolver_ptr->set_n_iter(n_iter); 

    evolver_ptr->set_collision_mode(collision_mode);
    evolver_ptr->set_roche_mode(roche_mode);
}  
void Tidy::initialize() {
    evolver_ptr->initialize(bodies, force_ptr);
    timestep_ptr->initialize(bodies);
}      
       
// Collision handling
bool Tidy::is_collision_detected() {
    return force_ptr->collision_detected;
}
vector< array<int, 2> > Tidy::get_collision_indices() {
    return force_ptr->index_collisions;
}
bool Tidy::get_collision_flag() {
    return evolver_ptr->collision_flag;
}

bool Tidy::is_roche_detected() {
    return force_ptr->roche_detected;
}
vector< array<int, 2> > Tidy::get_roche_indices() {
    return force_ptr->index_roche;
}
bool Tidy::get_roche_flag() {
    return evolver_ptr->roche_flag;
}
    
// Shape handling
void Tidy::set_to_spherical_shape() {
    shape.set_permanent_shape(bodies);
    shape.set_to_spherical_shape(bodies);
    shape.calculate_I_inv(bodies);
}
void Tidy::set_to_equilibrium_shape() {
    force.update_tidal_deformation(bodies);
    force.update_centrifugal_deformation(bodies);
    force.update_equilibrium_tensor(bodies);

    shape.set_permanent_shape(bodies);
    shape.set_to_equilibrium_shape(bodies);
    shape.calculate_I_inv(bodies);
}
void Tidy::update_angular_momentum() {
    spin.calculate_L(bodies);
    spin.init_K(bodies);
    shape.init_J(bodies);
}
    
// Get diagnostics
array<double, 3> Tidy::get_center_of_mass() {
    return evolver_ptr->get_center_of_mass(bodies);
}
array<double, 3> Tidy::get_center_of_mass_velocity() {
    return evolver_ptr->get_center_of_mass_velocity(bodies);
}
array<double, 3> Tidy::get_orbital_angular_momentum() {
    return evolver_ptr->get_orbital_angular_momentum(bodies);
}
array<double, 3> Tidy::get_spin_angular_momentum() {
    return evolver_ptr->get_spin_angular_momentum(bodies);
}
double Tidy::get_orbital_kinetic_energy() {
    return evolver_ptr->get_orbital_kinetic_energy(bodies);
}
double Tidy::get_spin_kinetic_energy() {
    return evolver_ptr->get_spin_kinetic_energy(bodies);
}
double Tidy::get_potential_energy() {
    return evolver_ptr->get_potential_energy(bodies);
}
array<double, 3> Tidy::get_angular_momentum() {
    return evolver_ptr->get_angular_momentum(bodies);
}
double Tidy::get_energy() {
    return evolver_ptr->get_energy(bodies);
}
    
// Top level evolve functions
void Tidy::evolve_model(double t_end) {
    evolver_ptr->assign_vectors(bodies, force_ptr);
    timestep_ptr->set_dt_sgn(dt_sgn);

    evolver_ptr->evolve_model(t, t_end, bodies, force_ptr, timestep_ptr);
    
    num_integration_step += evolver_ptr->num_integration_step;
}
void Tidy::evolve_model(int N_step) {
    evolver_ptr->assign_vectors(bodies, force_ptr);
    timestep_ptr->set_dt_sgn(dt_sgn);

    evolver_ptr->evolve_model(t, N_step, bodies, force_ptr, timestep_ptr);
    
    num_integration_step += evolver_ptr->num_integration_step;
}

    



