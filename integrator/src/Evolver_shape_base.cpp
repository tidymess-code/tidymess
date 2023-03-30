#include "Evolver_shape_base.h"

// Get diagnostics
array<double, 3> Evolver_shape_base::get_spin_angular_momentum(vector<Body> &bodies) {
    return spin.get_angular_momentum(bodies);
}
double Evolver_shape_base::get_spin_kinetic_energy(vector<Body> &bodies) {
    return spin.get_kinetic_energy(bodies);
}
double Evolver_shape_base::get_potential_energy(vector<Body> &bodies) {
    return orbit.get_potential_energy(bodies);
}
array<double, 3> Evolver_shape_base::get_angular_momentum(vector<Body> &bodies) {
    array<double, 3> Lorbit = get_orbital_angular_momentum(bodies);
    array<double, 3> Lspin = get_spin_angular_momentum(bodies);
    array<double, 3> Ltot = {};
    for(int k=0; k<3; k++) {
        Ltot[k] = Lorbit[k] + Lspin[k];
    }
    return Ltot;
}
double Evolver_shape_base::get_energy(vector<Body> &bodies) {
    double EK_orbit = get_orbital_kinetic_energy(bodies);
    double EK_spin = get_spin_kinetic_energy(bodies);
    double EP = get_potential_energy(bodies);
    double E = EK_orbit + EK_spin + EP;
    return E;
}

