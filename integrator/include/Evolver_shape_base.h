#include "Force.h"

#include "Timestep_base.h"

#include "Timestep_const.h"
#include "Timestep_adapt.h"
#include "Timestep_adapt_weight.h"
#include "Timestep_direct_weight.h"

#include "Evolver_base.h"

#include "Orbit.h"
#include "Spin.h"
#include "Shape.h"

#ifndef __Evolver_shape_base_h
#define __Evolver_shape_base_h

class Evolver_shape_base : public Evolver_base {
    public:

    // Get diagnostics
    array<double, 3> get_spin_angular_momentum(vector<Body> &bodies);
    double get_spin_kinetic_energy(vector<Body> &bodies);
    double get_potential_energy(vector<Body> &bodies);
    array<double, 3> get_angular_momentum(vector<Body> &bodies);
    double get_energy(vector<Body> &bodies);
};

#endif


