#include "Evolver_nbody_base.h"

// Get diagnostics
array<double, 3> Evolver_nbody_base::get_spin_angular_momentum(vector<Body> &bodies) {
    array<double, 3> L = {};
    return L;
}
double Evolver_nbody_base::get_spin_kinetic_energy(vector<Body> &bodies) {
    return 0;
}
double Evolver_nbody_base::get_potential_energy(vector<Body> &bodies) {
    return orbit.get_potential_energy_nbody(bodies);
}
array<double, 3> Evolver_nbody_base::get_angular_momentum(vector<Body> &bodies) {
    return orbit.get_angular_momentum(bodies);
}
double Evolver_nbody_base::get_energy(vector<Body> &bodies) {
    double EK = get_orbital_kinetic_energy(bodies);
    double EP = get_potential_energy(bodies);
    return EK + EP;
}

