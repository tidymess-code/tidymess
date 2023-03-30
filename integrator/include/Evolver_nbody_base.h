#include "Evolver_base.h"

#ifndef __Evolver_nbody_base_h
#define __Evolver_nbody_base_h

class Evolver_nbody_base : public Evolver_base {
    public:

    // Get diagnostics
    array<double, 3> get_spin_angular_momentum(vector<Body> &bodies);
    double get_spin_kinetic_energy(vector<Body> &bodies);
    double get_potential_energy(vector<Body> &bodies);
    array<double, 3> get_angular_momentum(vector<Body> &bodies);
    double get_energy(vector<Body> &bodies);
};

#endif


