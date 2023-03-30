#include "Orbit.h"

// Initializers
Orbit::Orbit() {
    ;
}

void Orbit::drift_r(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->r[k] += b->v[k]*dt;
        }
    }
}

void Orbit::kick_v(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->v[k] += b->a[k]*dt;
        }
    }
}
void Orbit::kick_vv(vector<Body> &bodies, double dt) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->vv[k] += b->a[k]*dt;
        }
    }
}

void Orbit::memorize_r(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->r_prev[k] = b->r[k];
        }
    }
}
void Orbit::swap_r_and_r_prev(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            double d = b->r_prev[k];
            b->r_prev[k] = b->r[k];
            b->r[k] = d;
        }
    }
}

void Orbit::init_vv(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            b->vv[k] = b->v[k];
        }
    }
}
void Orbit::swap_v_and_vv(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        for(int k=0; k<3; k++) {
            double dummy = b->v[k];
            b->v[k] = b->vv[k];
            b->vv[k] = dummy;
        }
    }
}
    
// Diagnostics
double Orbit::get_kinetic_energy(vector<Body> &bodies) {
    double EK = 0;
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type > 0) {
            double v2 = inner_product(b->v.begin(), b->v.end(), b->v.begin(), 0.);
            EK += 0.5*b->m*v2;
        }
    }
    return EK;
}
double Orbit::get_potential_energy(vector<Body> &bodies) {
    double EP = 0;

    array<double, 3> dr, drn;
    array<double, 6> J;

    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
 
        if(bi->particle_type > 0) {
 
            for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

                if(bj->particle_type > 0) {

                    for(int k=0; k<3; k++) {
                        dr[k] = bi->r[k]-bj->r[k];
                    }

                    double dr2 = inner_product(dr.begin(), dr.end(), dr.begin(), 0.);
                    double dr_3 = 1./sqrt(dr2*dr2*dr2);
                    double dr_1 = dr_3*dr2;

                    double EP1 = -bi->m*bj->m*dr_1;

                    EP += EP1;

                    if(bi->particle_type >= 2 || bj->particle_type >= 2) {

                        for(int k=0; k<3; k++) {
                            drn[k] = dr[k]*dr_1;
                        }

                        for(int k=0; k<6; k++) {
                            J[k] = bi->m*bj->I[k] + bj->m*bi->I[k];
                        }

                        double dxJdx = drn[0]*(drn[0]*J[0] + drn[1]*J[1] + drn[2]*J[2]);
                        double dyJdy = drn[1]*(drn[0]*J[1] + drn[1]*J[3] + drn[2]*J[4]);
                        double dzJdz = drn[2]*(drn[0]*J[2] + drn[1]*J[4] + drn[2]*J[5]);

                        double EP2 = 1.5*dr_3*(dxJdx + dyJdy + dzJdz - (J[0]+J[3]+J[5])/3.);

                        EP += EP2;
                    }
                
                }
            }
        
        }
    }

    return EP;
}
double Orbit::get_potential_energy_nbody(vector<Body> &bodies) {
    double EP = 0;

    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {

        if(bi->particle_type > 0) {

            for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

                if(bj->particle_type > 0) {

                    array<double, 3> dr = {};
                    for(int k=0; k<3; k++) {
                        dr[k] = bi->r[k]-bj->r[k];
                    }

                    double dr2 = inner_product(dr.begin(), dr.end(), dr.begin(), 0.);
                    double dr_1 = 1./sqrt(dr2);

                    double EP1 = -bi->m*bj->m*dr_1;

                    EP += EP1;
                }

            }

        }

    }

    return EP;
}

array<double, 3> Orbit::get_center_of_mass(vector<Body> &bodies) {
    array<double, 3> r = {};
    int N = bodies.size();
    double M = 0.;
    for(int i=0; i<N; i++) {
        if(bodies[i].particle_type > 0) { 
            M += bodies[i].m;
            for(int k=0; k<3; k++) {
                r[k] += bodies[i].m*bodies[i].r[k];
            }
        }
    }
    for(int k=0; k<3; k++) {
        r[k] /= M;
    }
    return r;
}
array<double, 3> Orbit::get_center_of_mass_velocity(vector<Body> &bodies) {
    array<double, 3> v = {};
    int N = bodies.size();
    double M = 0.;
    for(int i=0; i<N; i++) {
        if(bodies[i].particle_type > 0) { 
            M += bodies[i].m;
            for(int k=0; k<3; k++) {
                v[k] += bodies[i].m*bodies[i].v[k];
            }
        }
    }
    for(int k=0; k<3; k++) {
        v[k] /= M;
    }
    return v;
}
array<double, 3> Orbit::get_angular_momentum(vector<Body> &bodies) {
    array<double, 3> L = {};
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        if(bodies[i].particle_type > 0) { 
            double Lx = bodies[i].r[1]*bodies[i].v[2] - bodies[i].r[2]*bodies[i].v[1];
            double Ly = bodies[i].r[2]*bodies[i].v[0] - bodies[i].r[0]*bodies[i].v[2];
            double Lz = bodies[i].r[0]*bodies[i].v[1] - bodies[i].r[1]*bodies[i].v[0];
            L[0] += bodies[i].m*Lx;
            L[1] += bodies[i].m*Ly;
            L[2] += bodies[i].m*Lz;
        }
    }
    return L;
}


