#include "Force.h"

// Initializers

Force::Force() {
    this->set_default_settings();
}

Force::Force(int pn_order, double c) {
    this->set_default_settings();

    this->pn_order = pn_order;
    this->c = c;
    this->set_speed_of_light(this->c);
}

void Force::set_default_settings() {
    this->pn_order = 0;

    this->c = 1e100;
    this->set_speed_of_light(this->c);

    collision_detected = false;
    index_collisions.clear();

    roche_detected = false;
    index_roche.clear();
}

// Setters

void Force::set_pn_order(int pn_order) {
    this->pn_order = pn_order;
}
void Force::set_speed_of_light(double c) {
    this->c = c;
    this->c_2 = 1./(this->c*this->c);
    this->c_4 = this->c_2*this->c_2;
    this->c_5 = this->c_4/this->c;
}

// Collision checker

void Force::check_for_collisions(Body &bi, Body &bj, double &dr_1) { 
    double Rsum = bi.R+bj.R;
    double fcol_R = Rsum*dr_1;

    if(fcol_R > 1) {
        this->collision_detected = true;

        array<int, 2> couple;
        couple[0] = bi.id;
        couple[1] = bj.id;

        this->index_collisions.push_back(couple);
    }
}
void Force::check_for_collisions(Body &bi, Body &bj, double &dr_3, double &dr2) { 
    double dr_1 = dr_3*dr2;
    this->check_for_collisions(bi, bj, dr_1);            
}

void Force::check_for_roche(Body &bi, Body &bj, double &dr_1) { 
    bool is_roche = false;

    if(bi.rho > 0 && bj.m > 0) {
        double C = bi.roche_factor;
        double m = bi.m;
        double r = bi.R;
        double M = bj.m;
        
        double d = C * r * pow(M/m, 1./3);
        double f = d*dr_1;
                                         
        if(f > 1) {
            this->roche_detected = true;  

            array<int, 2> couple;
            couple[0] = bi.id;
            couple[1] = bj.id;

            this->index_roche.push_back(couple);        

            is_roche = true;
        }
    }
    
    if(!is_roche) {
        if(bj.rho > 0 && bi.m > 0) {
            double C = bj.roche_factor;
            double m = bj.m;
            double r = bj.R;
            double M = bi.m;
        
            double d = C * r * pow(M/m, 1./3);
            double f = d*dr_1;
            
            if(f > 1) {
                this->roche_detected = true;  

                array<int, 2> couple;
                couple[0] = bi.id;
                couple[1] = bj.id;

                this->index_roche.push_back(couple);        

                is_roche = true;
            }
        }
    }    
}
void Force::check_for_roche(Body &bi, Body &bj, double &dr_3, double &dr2) { 
    double dr_1 = dr_3*dr2;
    this->check_for_roche(bi, bj, dr_1);            
}

void Force::reset_collisions() { 
    this->collision_detected = false;
    this->index_collisions.clear();
}
void Force::reset_roche() { 
    this->roche_detected = false;
    this->index_roche.clear();
}

// Force functions

void Force::precalc_sep_vec(Body &bi, Body &bj, array<double, 3> &dr) { 
    for(int k=0; k<3; k++) {
        dr[k] = bi.r[k]-bj.r[k];
    }    
}        
void Force::precalc_sep_mag(Body &bi, Body &bj, array<double, 3> &dr, double &dr2, double &dr_3) { 
    precalc_sep_vec(bi, bj, dr);
    dr2 = inner_product(dr.begin(), dr.end(), dr.begin(), 0.);
    dr_3 = 1./sqrt(dr2*dr2*dr2);
}    
void Force::precalc_sep_vec_norm(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &drn, double &dr_4, double &dr_3) {
    double dr2;
    
    precalc_sep_mag(bi, bj, dr, dr2, dr_3);

    double dr_1 = dr2*dr_3;
    dr_4 = dr_3*dr_1;

    for(int k=0; k<3; k++) {
        drn[k] = dr[k]*dr_1;
    }
}        
void Force::precalc_sep_vec_norm_load(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &drn, double &dr_4, double &dr2, double &dr_3) {
    precalc_sep_vec(bi, bj, dr);

    double dr_1 = dr2*dr_3;
    dr_4 = dr_3*dr_1;

    for(int k=0; k<3; k++) {
        drn[k] = dr[k]*dr_1;
    }
}        
void Force::precalc_sep_vec_norm_save(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &drn, double &dr_4, double &dr2, double &dr_3) {
    precalc_sep_mag(bi, bj, dr, dr2, dr_3);

    double dr_1 = dr2*dr_3;
    dr_4 = dr_3*dr_1;

    for(int k=0; k<3; k++) {
        drn[k] = dr[k]*dr_1;
    }
}        

void Force::precalc_pn(Body &bi, Body &bj, bool &use_vv, array<double, 3> &drn, double &dr_1, double &vivi, double &vjvj, double &vivj, double &drnvi, double &drnvj, double &dv2, array<double, 3> &dv, double &drnvidrnvi, double &drnvjdrnvj, double &bimdr_1, double &bjmdr_1) { 
    if(use_vv) {
        for(int k=0; k<3; k++) dv[k] = bi.vv[k]-bj.vv[k];  
        vivi  = inner_product(bi.vv.begin(), bi.vv.end(), bi.vv.begin(), 0.);
        vjvj  = inner_product(bj.vv.begin(), bj.vv.end(), bj.vv.begin(), 0.);
        vivj  = inner_product(bi.vv.begin(), bi.vv.end(), bj.vv.begin(), 0.);
        drnvi = inner_product(drn.begin(), drn.end(), bi.vv.begin(), 0.);
        drnvj = inner_product(drn.begin(), drn.end(), bj.vv.begin(), 0.);
    }
    else {
        for(int k=0; k<3; k++) dv[k] = bi.v[k]-bj.v[k];  
        vivi  = inner_product(bi.v.begin(), bi.v.end(), bi.v.begin(), 0.);
        vjvj  = inner_product(bj.v.begin(), bj.v.end(), bj.v.begin(), 0.);
        vivj  = inner_product(bi.v.begin(), bi.v.end(), bj.v.begin(), 0.);
        drnvi = inner_product(drn.begin(), drn.end(), bi.v.begin(), 0.);
        drnvj = inner_product(drn.begin(), drn.end(), bj.v.begin(), 0.);
    }

    dv2 = inner_product(dv.begin(), dv.end(), dv.begin(), 0.);

    drnvidrnvi = drnvi*drnvi;	
    drnvjdrnvj = drnvj*drnvj;	

    bimdr_1 = bi.m*dr_1;
    bjmdr_1 = bj.m*dr_1;
}
void Force::precalc_tidal(double &dr2, double &dr_3, double &dr_1, double &dr_2, double &dr_4, double &dr_5, array<double, 3> &dr, array<double, 3> &drn, double &dyndyn, double &dzndzn, double &dxndyn, double &dxndzn, double &dyndzn) {
    dr_1 = dr2*dr_3;
    dr_2 = dr_1*dr_1;
    dr_4 = dr_2*dr_2;
    dr_5 = dr_3*dr_2;

    for(int k=0; k<3; k++) {
        drn[k] = dr[k]*dr_1;
    }

    dyndyn = 0.5*(drn[1]*drn[1] - 0.2);
    dzndzn = 0.5*(drn[2]*drn[2] - 0.2);
    dxndyn = drn[0]*drn[1];
    dxndzn = drn[0]*drn[2];
    dyndzn = drn[1]*drn[2];
}   
   
void Force::calculate_nbody_force(Body &bi, Body &bj, array<double, 3> &dr, double &dr2, double &dr_3, array<double, 3> &ai, array<double, 3> &aj) {
    bool calc_done = false;

    if(bj.particle_type > 0) {    
        precalc_sep_mag(bi, bj, dr, dr2, dr_3);
        calc_done = true;
    
        double amag = bj.m*dr_3;

        for(int k=0; k<3; k++) {
            ai[k] = -amag*dr[k];
        }  
    }
    else {
        for(int k=0; k<3; k++) {
            ai[k] = 0;
        }      
    }
    
    if(bi.particle_type > 0) {
        if(calc_done == false) {
            precalc_sep_mag(bi, bj, dr, dr2, dr_3);
            calc_done = true;
        }
    
        double amag = bi.m*dr_3;

        for(int k=0; k<3; k++) {
            aj[k] =  amag*dr[k];
        }  
    }    
    else {
        for(int k=0; k<3; k++) {
            aj[k] = 0;
        }      
    }
}
void Force::calculate_nbody_force_load(Body &bi, Body &bj, array<double, 3> &dr, double &dr2, double &dr_3, array<double, 3> &ai, array<double, 3> &aj) { 
    bool calc_done = false;

    if(bj.particle_type > 0) {    
        precalc_sep_vec(bi, bj, dr);
        calc_done = true;
    
        double amag = bj.m*dr_3;

        for(int k=0; k<3; k++) {
            ai[k] = -amag*dr[k];
        }  
    }
    else {
        for(int k=0; k<3; k++) {
            ai[k] = 0;
        }      
    }
    
    if(bi.particle_type > 0) {
        if(calc_done == false) {
            precalc_sep_vec(bi, bj, dr);
            calc_done = true;
        }
    
        double amag = bi.m*dr_3;

        for(int k=0; k<3; k++) {
            aj[k] =  amag*dr[k];
        }  
    }    
    else {
        for(int k=0; k<3; k++) {
            aj[k] = 0;
        }      
    }
}

void Force::calculate_pn_force(Body &bi, Body &bj, double &dr2, double &dr_1, double &dr_2, double &dr_3, double &dr_4, array<double, 3> &dr, array<double, 3> &drn, double &dv2, array<double, 3> &dv, array<double, 3> &ai_pn, array<double, 3> &aj_pn, bool use_vv) { 
    bool calc_done = false;

    double vivi, vjvj, vivj, drnvi, drnvj, drnvidrnvi, drnvjdrnvj, bimdr_1, bjmdr_1;

    if(bj.particle_type > 0) {    
        precalc_pn(bi, bj, use_vv, drn, dr_1, vivi, vjvj, vivj, drnvi, drnvj, dv2, dv, drnvidrnvi, drnvjdrnvj, bimdr_1, bjmdr_1);
        calc_done = true;
    
        // 1 PN
        array<double, 3> ai_1pn = {};

        double A = bj.m*dr_2;
        double B = -vivi - 2*vjvj + 4*vivj + 1.5*drnvjdrnvj + 5*bimdr_1 + 4*bjmdr_1;
        double C = 4*drnvi - 3*drnvj;

        for(int k=0; k<3; k++) {
            ai_1pn[k] = A*(B*drn[k] + C*dv[k]);
        }

        // 2 PN
        array<double, 3> ai_2pn = {};

        if(pn_order > 1) {
            A = bj.m*dr_2;
            B = vivi*drnvj + 4*vjvj*drnvi - 5*vjvj*drnvj - 4*vivj*drnvi + 4*vivj*drnvj - 6*drnvi*drnvjdrnvj + 4.5*drnvjdrnvj*drnvj + 
                bimdr_1*(-15.75*drnvi + 13.75*drnvj) + bjmdr_1*(-2*drnvi - 2*drnvj);
            C = -2*vjvj*vjvj + 4*vjvj*vivj - 2*vivj*vivj + 1.5*vivi*drnvjdrnvj + 4.5*vjvj*drnvjdrnvj - 6*vivj*drnvjdrnvj - 1.875*drnvjdrnvj*drnvjdrnvj + 
                bjmdr_1*(4*vjvj - 8*vivj + 2*drnvidrnvi - 4*drnvi*drnvj - 6*drnvjdrnvj) + 
                bimdr_1*(-3.75*vivi + 1.25*vjvj - 2.5*vivj + 19.5*drnvidrnvi - 39.*drnvi*drnvj + 8.5*drnvjdrnvj);
            double D = bj.m*dr_4*(-14.25*bi.m*bi.m - 9*bj.m*bj.m - 34.5*bi.m*bj.m);

            for(int k=0; k<3; k++) {
                ai_2pn[k] = A*( B*dv[k] + C*drn[k]) + D*drn[k];
            }
        }

        // 2.5 PN
        array<double, 3> ai_25pn = {};

        if(pn_order > 2) {
            A = 0.8*bi.m*bj.m*dr_3;
            B = -dv2 + 2*bimdr_1 - 8*bjmdr_1;
            C = (3*dv2 - 6*bimdr_1 + (17.+1./3)*bjmdr_1)*(drnvi - drnvj);
            for(int k=0; k<3; k++) {
                ai_25pn[k] = A*( B*dv[k] + C*drn[k] );
            }
        }

        for(int k=0; k<3; k++) {
            ai_pn[k] = ai_1pn[k]*c_2 + ai_2pn[k]*c_4 + ai_25pn[k]*c_5;
        }            
    }
    else {
        for(int k=0; k<3; k++) {
            ai_pn[k] = 0;
        }      
    }
    
    if(bi.particle_type > 0) {
        if(calc_done == false) {
            precalc_pn(bi, bj, use_vv, drn, dr_1, vivi, vjvj, vivj, drnvi, drnvj, dv2, dv, drnvidrnvi, drnvjdrnvj, bimdr_1, bjmdr_1);
            calc_done = true;
        }
    
        // 1 PN
        array<double, 3> aj_1pn = {};

        double A = bi.m*dr_2;
        double B = -vjvj - 2*vivi + 4*vivj + 1.5*drnvidrnvi + 5*bjmdr_1 + 4*bimdr_1;
        double C = 4*-drnvj - 3*-drnvi;

        for(int k=0; k<3; k++) {
            aj_1pn[k] = A*(B*-drn[k] + C*-dv[k] );
        }

        // 2 PN
        array<double, 3> aj_2pn = {};

        if(pn_order > 1) {
            A = bi.m*dr_2;
            B = vjvj*-drnvi + 4*vivi*-drnvj - 5*vivi*-drnvi - 4*vivj*-drnvj + 4*vivj*-drnvi - 6*-drnvj*drnvidrnvi + 4.5*drnvidrnvi*-drnvi + 
                bjmdr_1*(-15.75*-drnvj + 13.75*-drnvi) + bimdr_1*(-2*-drnvj - 2*-drnvi);
            C = -2*vivi*vivi + 4*vivi*vivj - 2*vivj*vivj + 1.5*vjvj*drnvidrnvi + 4.5*vivi*drnvidrnvi - 6*vivj*drnvidrnvi - 1.875*drnvidrnvi*drnvidrnvi + 
                bimdr_1*(4*vivi - 8*vivj + 2*drnvjdrnvj - 4*-drnvj*-drnvi - 6*drnvidrnvi) + 
                bjmdr_1*(-3.75*vjvj + 1.25*vivi - 2.5*vivj + 19.5*drnvjdrnvj - 39.*-drnvj*-drnvi + 8.5*drnvidrnvi);
            double D = bi.m*dr_4*(-14.25*bj.m*bj.m - 9*bi.m*bi.m - 34.5*bj.m*bi.m);

            for(int k=0; k<3; k++) {
                aj_2pn[k] = A*( B*-dv[k] + C*-drn[k]) + D*-drn[k];
            }
        }

        // 2.5 PN
        array<double, 3> aj_25pn = {};

        if(pn_order > 2) {
            A = 0.8*bi.m*bj.m*dr_3;
            B = -dv2 + 2*bjmdr_1 - 8*bimdr_1;
            C = (3*dv2 - 6*bjmdr_1 + (17.+1./3)*bimdr_1)*(-drnvj - -drnvi);
            for(int k=0; k<3; k++) {
                aj_25pn[k] = A*( B*-dv[k] + C*-drn[k] );
            }
        }

        for(int k=0; k<3; k++) {
            aj_pn[k] = aj_1pn[k]*c_2 + aj_2pn[k]*c_4 + aj_25pn[k]*c_5;
        }            
    }    
    else {
        for(int k=0; k<3; k++) {
            aj_pn[k] = 0;
        }      
    }  
}     

void Force::calculate_tidal_force(Body &bi, Body &bj, double &dr2, double &dr_3, double &dr_1, double &dr_2, double &dr_4, double &dr_5, array<double, 3> &dr, array<double, 3> &drn, array<double, 3> &ai_tidal, array<double, 3> &aj_tidal, array<double, 3> &hij, array<double, 3> &hji, bool use_J) {
    bool calc_done = false;

    double dyndyn, dzndzn, dxndyn, dxndzn, dyndzn;
    array<double, 6> bi_I, bj_I;
    array<double, 3> gij = {};
    array<double, 3> gji = {};
    double gmag, A;
    
    if(bj.particle_type == 3) {    
        precalc_tidal(dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, dyndyn, dzndzn, dxndyn, dxndzn, dyndzn);
        A = -3*dr_4;
        calc_done = true;

        if(use_J) {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.I[k];
            }    
        }

        gmag = 15*dr_5 * ( (bj_I[3]-bj_I[0])*dyndyn + (bj_I[5]-bj_I[0])*dzndzn + bj_I[1]*dxndyn + bj_I[2]*dxndzn + bj_I[4]*dyndzn );

        for(int k=0; k<3; k++) {
            gij[k] = gmag*dr[k];
        }
        
        hij[0] = A * ( bj_I[2]*drn[2] + bj_I[1]*drn[1] );
        hij[1] = A * ( (bj_I[3]-bj_I[0])*drn[1] + bj_I[1]*drn[0] + bj_I[4]*drn[2] );
        hij[2] = A * ( (bj_I[5]-bj_I[0])*drn[2] + bj_I[2]*drn[0] + bj_I[4]*drn[1] );
    }
    else {
        for(int k=0; k<3; k++) {
            gij[k] = 0;
            hij[k] = 0;
        }         
    }
    
    if(bi.particle_type == 3) {
        if(calc_done == false) {
            precalc_tidal(dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, dyndyn, dzndzn, dxndyn, dxndzn, dyndzn);
            A = -3*dr_4;
            calc_done = true;
        }

        if(use_J) {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.I[k];
            }    
        }

        gmag = 15*dr_5 * ( (bi_I[3]-bi_I[0])*dyndyn + (bi_I[5]-bi_I[0])*dzndzn + bi_I[1]*dxndyn + bi_I[2]*dxndzn + bi_I[4]*dyndzn );

        for(int k=0; k<3; k++) {
            gji[k] = gmag*dr[k];
        }
            
        hji[0] = A * ( bi_I[2]*drn[2] + bi_I[1]*drn[1] );
        hji[1] = A * ( (bi_I[3]-bi_I[0])*drn[1] + bi_I[1]*drn[0] + bi_I[4]*drn[2] );
        hji[2] = A * ( (bi_I[5]-bi_I[0])*drn[2] + bi_I[2]*drn[0] + bi_I[4]*drn[1] );
    }    
    else {
        for(int k=0; k<3; k++) {
            gji[k] = 0;
            hji[k] = 0;
        }          
    }

    if(bi.particle_type > 0 && bj.particle_type > 0) {
        for(int k=0; k<3; k++) {
            double F = bi.m*gij[k] + bi.m*hij[k] + bj.m*gji[k] + bj.m*hji[k];
            ai_tidal[k] =  F / bi.m;
            aj_tidal[k] = -F / bj.m;
        }
    }
    else {
        if(bi.particle_type > 0) {
            for(int k=0; k<3; k++) {
                double acc = gji[k] + hji[k];
                ai_tidal[k] = 0;
                aj_tidal[k] = -acc;
            }
        }
        else if(bj.particle_type > 0) {
            for(int k=0; k<3; k++) {
                double acc = gij[k] + hij[k];
                ai_tidal[k] = acc;
                aj_tidal[k] = 0;
            }
        }
        else {
            for(int k=0; k<3; k++) {
                ai_tidal[k] = 0;
                aj_tidal[k] = 0;
            }                
        }
    }
    
    for(int k=0; k<3; k++) {
        hij[k] *= bi.m;
        hji[k] *= bj.m;
    }
}

void Force::calculate_tidal_h(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji, bool use_J) {
    bool calc_done = false;

    array<double, 3> drn;
    array<double, 6> bi_I, bj_I;
    double A, dr_4, dr_3;
    
    if(bj.particle_type == 3) {    
        precalc_sep_vec_norm(bi, bj, dr, drn, dr_4, dr_3);
        A = -3*dr_4;
        calc_done = true;

        if(use_J) {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.I[k];
            }    
        }
        
        hij[0] = A * ( bj_I[2]*drn[2] + bj_I[1]*drn[1] );
        hij[1] = A * ( (bj_I[3]-bj_I[0])*drn[1] + bj_I[1]*drn[0] + bj_I[4]*drn[2] );
        hij[2] = A * ( (bj_I[5]-bj_I[0])*drn[2] + bj_I[2]*drn[0] + bj_I[4]*drn[1] );
    }
    else {
        for(int k=0; k<3; k++) {
            hij[k] = 0;
        }         
    }
    
    if(bi.particle_type == 3) {
        if(calc_done == false) {
            precalc_sep_vec_norm(bi, bj, dr, drn, dr_4, dr_3);
            A = -3*dr_4;
            calc_done = true;
        }

        if(use_J) {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.I[k];
            }    
        }
            
        hji[0] = A * ( bi_I[2]*drn[2] + bi_I[1]*drn[1] );
        hji[1] = A * ( (bi_I[3]-bi_I[0])*drn[1] + bi_I[1]*drn[0] + bi_I[4]*drn[2] );
        hji[2] = A * ( (bi_I[5]-bi_I[0])*drn[2] + bi_I[2]*drn[0] + bi_I[4]*drn[1] );
    }    
    else {
        for(int k=0; k<3; k++) {
            hji[k] = 0;
        }          
    }
    
    for(int k=0; k<3; k++) {
        hij[k] *= bi.m;
        hji[k] *= bj.m;
    }
}
void Force::calculate_tidal_h_load(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji, bool use_J, double &dr2, double &dr_3) {
    bool calc_done = false;

    array<double, 3> drn;
    array<double, 6> bi_I, bj_I;
    double A, dr_4;
    
    if(bj.particle_type == 3) {    
        precalc_sep_vec_norm_load(bi, bj, dr, drn, dr_4, dr2, dr_3);
        A = -3*dr_4;
        calc_done = true;

        if(use_J) {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.I[k];
            }    
        }
        
        hij[0] = A * ( bj_I[2]*drn[2] + bj_I[1]*drn[1] );
        hij[1] = A * ( (bj_I[3]-bj_I[0])*drn[1] + bj_I[1]*drn[0] + bj_I[4]*drn[2] );
        hij[2] = A * ( (bj_I[5]-bj_I[0])*drn[2] + bj_I[2]*drn[0] + bj_I[4]*drn[1] );
    }
    else {
        for(int k=0; k<3; k++) {
            hij[k] = 0;
        }         
    }
    
    if(bi.particle_type == 3) {
        if(calc_done == false) {
            precalc_sep_vec_norm_load(bi, bj, dr, drn, dr_4, dr2, dr_3);
            A = -3*dr_4;
            calc_done = true;
        }

        if(use_J) {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.I[k];
            }    
        }
            
        hji[0] = A * ( bi_I[2]*drn[2] + bi_I[1]*drn[1] );
        hji[1] = A * ( (bi_I[3]-bi_I[0])*drn[1] + bi_I[1]*drn[0] + bi_I[4]*drn[2] );
        hji[2] = A * ( (bi_I[5]-bi_I[0])*drn[2] + bi_I[2]*drn[0] + bi_I[4]*drn[1] );
    }    
    else {
        for(int k=0; k<3; k++) {
            hji[k] = 0;
        }          
    }
    
    for(int k=0; k<3; k++) {
        hij[k] *= bi.m;
        hji[k] *= bj.m;
    }    
}
void Force::calculate_tidal_h_save(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji, bool use_J, double &dr2, double &dr_3) {
    bool calc_done = false;

    array<double, 3> drn;
    array<double, 6> bi_I, bj_I;
    double A, dr_4;
    
    if(bi.particle_type > 0 || bj.particle_type > 0) {
        precalc_sep_vec_norm_save(bi, bj, dr, drn, dr_4, dr2, dr_3);    
        A = -3*dr_4;
        calc_done = true;
    }
    
    if(bj.particle_type == 3) {    
        if(calc_done == false) {
            precalc_sep_vec_norm_save(bi, bj, dr, drn, dr_4, dr2, dr_3);
            A = -3*dr_4;
            calc_done = true;
        }

        if(use_J) {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bj_I[k] = bj.I[k];
            }    
        }
        
        hij[0] = A * ( bj_I[2]*drn[2] + bj_I[1]*drn[1] );
        hij[1] = A * ( (bj_I[3]-bj_I[0])*drn[1] + bj_I[1]*drn[0] + bj_I[4]*drn[2] );
        hij[2] = A * ( (bj_I[5]-bj_I[0])*drn[2] + bj_I[2]*drn[0] + bj_I[4]*drn[1] );
    }
    else {
        for(int k=0; k<3; k++) {
            hij[k] = 0;
        }         
    }
    
    if(bi.particle_type == 3) {
        if(calc_done == false) {
            precalc_sep_vec_norm_save(bi, bj, dr, drn, dr_4, dr2, dr_3);
            A = -3*dr_4;
            calc_done = true;
        }

        if(use_J) {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.J[k];
            }
        }
        else {
            for(int k=0; k<6; k++) {
                bi_I[k] = bi.I[k];
            }    
        }
            
        hji[0] = A * ( bi_I[2]*drn[2] + bi_I[1]*drn[1] );
        hji[1] = A * ( (bi_I[3]-bi_I[0])*drn[1] + bi_I[1]*drn[0] + bi_I[4]*drn[2] );
        hji[2] = A * ( (bi_I[5]-bi_I[0])*drn[2] + bi_I[2]*drn[0] + bi_I[4]*drn[1] );
    }    
    else {
        for(int k=0; k<3; k++) {
            hji[k] = 0;
        }          
    }
    
    for(int k=0; k<3; k++) {
        hij[k] *= bi.m;
        hji[k] *= bj.m;
    }    
}

// Torque functions
    
void Force::calculate_torque(Body &bi, Body &bj, array<double, 3> &dr, array<double, 3> &hij, array<double, 3> &hji) {    
    if(bi.particle_type == 3) {
        bi.T[0] += -dr[1]*hji[2] + dr[2]*hji[1];
        bi.T[1] += -dr[2]*hji[0] + dr[0]*hji[2];
        bi.T[2] += -dr[0]*hji[1] + dr[1]*hji[0];
    }
    if(bj.particle_type == 3) {    
        bj.T[0] += -dr[1]*hij[2] + dr[2]*hij[1];
        bj.T[1] += -dr[2]*hij[0] + dr[0]*hij[2];
        bj.T[2] += -dr[0]*hij[1] + dr[1]*hij[0];
    }
}  

// Inertia tensor functions

void Force::calculate_deformation_tensor(Body &bi, Body &bj, array<double, 3> &drn, double &dr_3) {    
    bool calc_done = false; 
 
    double third, Xi, Xj, dr_4;
    array<double, 5> Ie_tidal = {};
    array<double, 3> dr = {};

    if(bi.particle_type == 3 && bj.particle_type > 0) {
        third = 1./3;
                             
        Ie_tidal[0] = (drn[0]*drn[0]-third);
        Ie_tidal[3] = (drn[1]*drn[1]-third);
        Ie_tidal[1] = drn[0]*drn[1];
        Ie_tidal[2] = drn[0]*drn[2];
        Ie_tidal[4] = drn[1]*drn[2];

        calc_done = true;

        Xi  = -bj.m*bi.kf_R5*dr_3;

        for(int k=0; k<5; k++) {
            bi.I_e_r[k] += Xi*Ie_tidal[k];
        }    
    }

    if(bj.particle_type == 3 && bi.particle_type > 0) {
        if(calc_done == false) {
            third = 1./3;
                        
            Ie_tidal[0] = (drn[0]*drn[0]-third);
            Ie_tidal[3] = (drn[1]*drn[1]-third);
            Ie_tidal[1] = drn[0]*drn[1];
            Ie_tidal[2] = drn[0]*drn[2];
            Ie_tidal[4] = drn[1]*drn[2];        

            calc_done = true;
        }

        Xj  = -bi.m*bj.kf_R5*dr_3;

        for(int k=0; k<5; k++) {
            bj.I_e_r[k] += Xj*Ie_tidal[k];
        }    
    }
}

// Acceleration updaters

void Force::update_acceleration(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k];
                bj->a[k] += aj[k];
            }
        }
    }
}

void Force::update_acceleration_pn(vector<Body> &bodies, bool use_vv) { 
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_pn
            double dr_1 = dr2*dr_3;
            double dr_2 = dr_1*dr_1;
            double dr_4 = dr_2*dr_2;
            array<double, 3> drn = {};
            for(int k=0; k<3; k++) {
                drn[k] = dr[k]*dr_1;
            }
            
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);
            
            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_pn[k];
            }
        }
    }
}

void Force::update_acceleration_tidal(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
        }
    }
}

void Force::update_acceleration_col(vector<Body> &bodies) { 
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k];
                bj->a[k] += aj[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}


void Force::update_acceleration_pn_col(vector<Body> &bodies, bool use_vv) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_pn
            double dr_1 = dr2*dr_3;
            double dr_2 = dr_1*dr_1;
            double dr_4 = dr_2*dr_2;
            array<double, 3> drn = {};
            for(int k=0; k<3; k++) {
                drn[k] = dr[k]*dr_1;
            }
            
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_pn[k];   
                bj->a[k] += aj[k] + aj_pn[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_tidal_col(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_tidal_pn(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }
        }
    }
}

void Force::update_acceleration_tidal_pn_col(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_acceleration_and_torque_pn(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_acceleration_and_torque_col(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
            
            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_pn_col(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_and_deformation_tensor(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
        b->I_e_r.fill(0);        
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Tidal deformation
            calculate_deformation_tensor(*bi, *bj, drn, dr_3);    
        }
    }
}

void Force::update_acceleration_and_torque_and_deformation_tensor_col(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
        b->I_e_r.fill(0);        
    }

    reset_collisions();
    reset_roche();

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Tidal deformation
            calculate_deformation_tensor(*bi, *bj, drn, dr_3);    

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_pn_and_deformation_tensor(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
        b->I_e_r.fill(0);        
    }
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Tidal deformation
            calculate_deformation_tensor(*bi, *bj, drn, dr_3);    
        }
    }
}
void Force::update_acceleration_and_torque_pn_and_deformation_tensor_col(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
        b->I_e_r.fill(0);        
    }

    reset_collisions();
    reset_roche();
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Tidal deformation
            calculate_deformation_tensor(*bi, *bj, drn, dr_3);    

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_load(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k];
                bj->a[k] += aj[k];
            }
        }
    }
}

void Force::update_acceleration_pn_load(vector<Body> &bodies, bool use_vv) { 
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_pn
            double dr_1 = dr2*dr_3;
            double dr_2 = dr_1*dr_1;
            double dr_4 = dr_2*dr_2;
            array<double, 3> drn = {};
            for(int k=0; k<3; k++) {
                drn[k] = dr[k]*dr_1;
            }
            
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);
            
            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_pn[k];   
                bj->a[k] += aj[k] + aj_pn[k];
            }
        }
    }
}

void Force::update_acceleration_tidal_load(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
        }
    }
}

void Force::update_acceleration_col_load(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k];
                bj->a[k] += aj[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_pn_col_load(vector<Body> &bodies, bool use_vv) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_pn
            double dr_1 = dr2*dr_3;
            double dr_2 = dr_1*dr_1;
            double dr_4 = dr_2*dr_2;
            array<double, 3> drn = {};
            for(int k=0; k<3; k++) {
                drn[k] = dr[k]*dr_1;
            }
            
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_pn[k];   
                bj->a[k] += aj[k] + aj_pn[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_tidal_col_load(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_tidal_pn_load(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }
        }
    }
}

void Force::update_acceleration_tidal_pn_col_load(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_load(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }    
}

void Force::update_acceleration_and_torque_pn_load(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_acceleration_and_torque_col_load(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
            
            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_pn_col_load(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // Load
            double dr2 = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];            
            cnt++;

            // f
            array<double, 3> dr = {};
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force_load(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_save(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;
                                    
            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k];
                bj->a[k] += aj[k];
            }
        }
    }
}

void Force::update_acceleration_pn_save(vector<Body> &bodies, bool use_vv) { 
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;
            
            // f_pn
            double dr_1 = dr2*dr_3;
            double dr_2 = dr_1*dr_1;
            double dr_4 = dr_2*dr_2;
            array<double, 3> drn = {};
            for(int k=0; k<3; k++) {
                drn[k] = dr[k]*dr_1;
            }
            
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);
            
            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_pn[k];   
                bj->a[k] += aj[k] + aj_pn[k];
            }
        }
    }
}

void Force::update_acceleration_tidal_save(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
        }
    }
}

void Force::update_acceleration_col_save(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }
    
    reset_collisions();
    reset_roche();

    int cnt = 0;

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k];
                bj->a[k] += aj[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_pn_col_save(vector<Body> &bodies, bool use_vv) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {
            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_pn
            double dr_1 = dr2*dr_3;
            double dr_2 = dr_1*dr_1;
            double dr_4 = dr_2*dr_2;
            array<double, 3> drn = {};
            for(int k=0; k<3; k++) {
                drn[k] = dr[k]*dr_1;
            }
            
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_pn[k];   
                bj->a[k] += aj[k] + aj_pn[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_tidal_col_save(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_tidal_pn_save(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }
        }
    }
}

void Force::update_acceleration_tidal_pn_col_save(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_save(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_acceleration_and_torque_pn_save(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_acceleration_and_torque_col_save(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    reset_collisions();
    reset_roche();
    
    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};

            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k];
                bj->a[k] += aj[k] + aj_tidal[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
            
            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_acceleration_and_torque_pn_col_save(vector<Body> &bodies, bool use_vv, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->a.fill(0);
        b->T.fill(0);
    }

    reset_collisions();
    reset_roche();

    int cnt = 0;
        
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // f
            array<double, 3> dr = {};
            double dr2, dr_3;
            array<double, 3> ai = {};
            array<double, 3> aj = {};
 
            calculate_nbody_force(*bi, *bj, dr, dr2, dr_3, ai, aj);

            // Save
            dr2_vec[cnt] = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;

            // f_tidal
            array<double, 3> drn = {};
            double dr_1, dr_2, dr_4, dr_5;
            array<double, 3> ai_tidal = {};
            array<double, 3> aj_tidal = {};
            array<double, 3> hij = {};
            array<double, 3> hji = {};
      
            calculate_tidal_force(*bi, *bj, dr2, dr_3, dr_1, dr_2, dr_4, dr_5, dr, drn, ai_tidal, aj_tidal, hij, hji, use_J);

            // f_pn
            double dv2;
            array<double, 3> dv = {};
            array<double, 3> ai_pn = {};
            array<double, 3> aj_pn = {};

            calculate_pn_force(*bi, *bj, dr2, dr_1, dr_2, dr_3, dr_4, dr, drn, dv2, dv, ai_pn, aj_pn, use_vv);

            // Acceleration
            for(int k=0; k<3; k++) {
                bi->a[k] += ai[k] + ai_tidal[k] + ai_pn[k];
                bj->a[k] += aj[k] + aj_tidal[k] + aj_pn[k];
            }

            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 

            // Check for collisions 
            check_for_collisions(*bi, *bj, dr_3, dr2);

            // Check for Roche limit 
            check_for_roche(*bi, *bj, dr_3, dr2);
        }
    }
}

void Force::update_torque(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->T.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // Tidal h function
            array<double, 3> dr;
            array<double, 3> hij, hji;    
                
            calculate_tidal_h(*bi, *bj, dr, hij, hji, use_J);
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_torque_load(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->T.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            double dr2  = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];
            cnt++;

            // Tidal h function
            array<double, 3> dr;
            array<double, 3> hij, hji;    
                
            calculate_tidal_h_load(*bi, *bj, dr, hij, hji, use_J, dr2, dr_3);
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_torque_save(vector<Body> &bodies, bool use_J) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->T.fill(0);
    }

    int cnt = 0;
    
    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            // Tidal h function
            array<double, 3> dr;
            array<double, 3> hij, hji;    
            double dr2, dr_3;
                
            calculate_tidal_h_save(*bi, *bj, dr, hij, hji, use_J, dr2, dr_3);

            // Save
            dr2_vec[cnt]  = dr2;
            dr_3_vec[cnt] = dr_3;
            cnt++;
            
            // Torques
            calculate_torque(*bi, *bj, dr, hij, hji); 
        }
    }
}

void Force::update_tidal_deformation(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->I_e_r.fill(0);
    }

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            bool calc_done = false; 
 
            double third, Xi, Xj, dr_4, dr_3;
            array<double, 5> Ie_tidal = {};
            array<double, 3> dr = {};
            array<double, 3> drn = {};

            if(bi->particle_type == 3 && bj->particle_type > 0) {
                third = 1./3;
         
                precalc_sep_vec_norm(*bi, *bj, dr, drn, dr_4, dr_3);
                    
                Ie_tidal[0] = (drn[0]*drn[0]-third);
                Ie_tidal[3] = (drn[1]*drn[1]-third);
                Ie_tidal[1] = drn[0]*drn[1];
                Ie_tidal[2] = drn[0]*drn[2];
                Ie_tidal[4] = drn[1]*drn[2];

                calc_done = true;

                Xi  = -bj->m*bi->kf_R5*dr_3;

                for(int k=0; k<5; k++) {
                    bi->I_e_r[k] += Xi*Ie_tidal[k];
                }    
            }

            if(bj->particle_type == 3 && bi->particle_type > 0) {
                if(calc_done == false) {
                    third = 1./3;

                    precalc_sep_vec_norm(*bi, *bj, dr, drn, dr_4, dr_3);
                                        
                    Ie_tidal[0] = (drn[0]*drn[0]-third);
                    Ie_tidal[3] = (drn[1]*drn[1]-third);
                    Ie_tidal[1] = drn[0]*drn[1];
                    Ie_tidal[2] = drn[0]*drn[2];
                    Ie_tidal[4] = drn[1]*drn[2];        

                    calc_done = true;
                }

                Xj  = -bi->m*bj->kf_R5*dr_3;

                for(int k=0; k<5; k++) {
                    bj->I_e_r[k] += Xj*Ie_tidal[k];
                }    
            }
        }
    }
}    
    
void Force::update_centrifugal_deformation(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            double w00 = b->w[0]*b->w[0];
            double w11 = b->w[1]*b->w[1];
            double w22 = b->w[2]*b->w[2];

            double w2 = w00 + w11 + w22;
            double w23 = w2/3;

            b->I_e_w[0] = b->kf_R5_3*(w00 - w23);
            b->I_e_w[3] = b->kf_R5_3*(w11 - w23);
            b->I_e_w[1] = b->kf_R5_3*b->w[0]*b->w[1];
            b->I_e_w[2] = b->kf_R5_3*b->w[0]*b->w[2];
            b->I_e_w[4] = b->kf_R5_3*b->w[1]*b->w[2];
        }
    }
}

void Force::update_equilibrium_tensor(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            for(int i=0; i<5; i++) {
                b->I_e[i] = b->I_e_r[i] + b->I_e_w[i];
            }
        }
    }
}

void Force::update_equilibrium_tensor_and_memorize(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        if(b->particle_type == 3) {
            for(int i=0; i<5; i++) {
                b->I_e_prev[i] = b->I_e[i];            
                b->I_e[i] = b->I_e_r[i] + b->I_e_w[i];
            }
        }
    }
}

void Force::update_deformation_tensor(vector<Body> &bodies, bool use_J) {
    if(use_J) {
        for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
            if(b->particle_type == 3) {
                b->dI_n[0] =  2*b->J_n[2]*b->w[1] - 2*b->J_n[1]*b->w[2] + (b->I_e[0]-b->J_n[0])*b->tau_inv;
                b->dI_n[3] = -2*b->J_n[4]*b->w[0] + 2*b->J_n[1]*b->w[2] + (b->I_e[3]-b->J_n[3])*b->tau_inv;
                b->dI_n[1] =   -b->J_n[2]*b->w[0] + b->J_n[4]*b->w[1] + (b->J_n[0]-b->J_n[3])*b->w[2] + (b->I_e[1]-b->J_n[1])*b->tau_inv;
                b->dI_n[2] =    b->J_n[1]*b->w[0] - (2*b->J_n[0]+b->J_n[3])*b->w[1] - b->J_n[4]*b->w[2] + (b->I_e[2]-b->J_n[2])*b->tau_inv;
                b->dI_n[4] =   (b->J_n[0]+2*b->J_n[3])*b->w[0] - b->J_n[1]*b->w[1] + b->J_n[2]*b->w[2] + (b->I_e[4]-b->J_n[4])*b->tau_inv;
            }
        }    
    }
    else {
        for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
            if(b->particle_type == 3) {
                b->dI_n[0] =  2*b->I_n[2]*b->w[1] - 2*b->I_n[1]*b->w[2] + (b->I_e[0]-b->I_n[0])*b->tau_inv;
                b->dI_n[3] = -2*b->I_n[4]*b->w[0] + 2*b->I_n[1]*b->w[2] + (b->I_e[3]-b->I_n[3])*b->tau_inv;
                b->dI_n[1] =   -b->I_n[2]*b->w[0] + b->I_n[4]*b->w[1] + (b->I_n[0]-b->I_n[3])*b->w[2] + (b->I_e[1]-b->I_n[1])*b->tau_inv;
                b->dI_n[2] =    b->I_n[1]*b->w[0] - (2*b->I_n[0]+b->I_n[3])*b->w[1] - b->I_n[4]*b->w[2] + (b->I_e[2]-b->I_n[2])*b->tau_inv;
                b->dI_n[4] =   (b->I_n[0]+2*b->I_n[3])*b->w[0] - b->I_n[1]*b->w[1] + b->I_n[2]*b->w[2] + (b->I_e[4]-b->I_n[4])*b->tau_inv;
            }
        }
    }
}

void Force::update_tidal_deformation_load(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->I_e_r.fill(0);
    }

    int cnt = 0;

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            bool calc_done = false; 
 
            double third, Xi, Xj, dr_4;
            array<double, 5> Ie_tidal = {};
            array<double, 3> dr = {};
            array<double, 3> drn = {};

            double dr2  = dr2_vec[cnt];
            double dr_3 = dr_3_vec[cnt];
            cnt++;

            if(bi->particle_type == 3 && bj->particle_type > 0) {
                third = 1./3;
         
                precalc_sep_vec_norm_load(*bi, *bj, dr, drn, dr_4, dr2, dr_3);
                    
                Ie_tidal[0] = (drn[0]*drn[0]-third);
                Ie_tidal[3] = (drn[1]*drn[1]-third);
                Ie_tidal[1] = drn[0]*drn[1];
                Ie_tidal[2] = drn[0]*drn[2];
                Ie_tidal[4] = drn[1]*drn[2];

                calc_done = true;

                Xi  = -bj->m*bi->kf_R5*dr_3;

                for(int k=0; k<5; k++) {
                    bi->I_e_r[k] += Xi*Ie_tidal[k];
                }    
            }

            if(bj->particle_type == 3 && bi->particle_type > 0) {
                if(calc_done == false) {
                    third = 1./3;

                    precalc_sep_vec_norm_load(*bi, *bj, dr, drn, dr_4, dr2, dr_3);
                                        
                    Ie_tidal[0] = (drn[0]*drn[0]-third);
                    Ie_tidal[3] = (drn[1]*drn[1]-third);
                    Ie_tidal[1] = drn[0]*drn[1];
                    Ie_tidal[2] = drn[0]*drn[2];
                    Ie_tidal[4] = drn[1]*drn[2];        

                    calc_done = true;
                }

                Xj  = -bi->m*bj->kf_R5*dr_3;

                for(int k=0; k<5; k++) {
                    bj->I_e_r[k] += Xj*Ie_tidal[k];
                }    
            }
        }
    }
}

void Force::update_tidal_deformation_save(vector<Body> &bodies) {
    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end(); ++b) {
        b->I_e_r.fill(0);
    }

    int cnt = 0;

    // O(N^2) operations
    for(vector<Body>::iterator bi = bodies.begin(); bi != bodies.end()-1; ++bi) {
        for(vector<Body>::iterator bj = bi+1; bj != bodies.end(); ++bj) {

            bool calc_done = false; 
 
            double third, Xi, Xj, dr_4, dr2, dr_3;
            array<double, 5> Ie_tidal = {};
            array<double, 3> dr = {};
            array<double, 3> drn = {};

            if(bi->particle_type == 3 && bj->particle_type > 0) {
                third = 1./3;
         
                precalc_sep_vec_norm_save(*bi, *bj, dr, drn, dr_4, dr2, dr_3);
                    
                dr2_vec[cnt] = dr2;    
                dr_3_vec[cnt] = dr_3;    
                    
                Ie_tidal[0] = (drn[0]*drn[0]-third);
                Ie_tidal[3] = (drn[1]*drn[1]-third);
                Ie_tidal[1] = drn[0]*drn[1];
                Ie_tidal[2] = drn[0]*drn[2];
                Ie_tidal[4] = drn[1]*drn[2];

                calc_done = true;

                Xi  = -bj->m*bi->kf_R5*dr_3;

                for(int k=0; k<5; k++) {
                    bi->I_e_r[k] += Xi*Ie_tidal[k];
                }    
            }

            if(bj->particle_type == 3 && bi->particle_type > 0) {
                if(calc_done == false) {
                    third = 1./3;

                    precalc_sep_vec_norm_save(*bi, *bj, dr, drn, dr_4, dr2, dr_3);

                    dr2_vec[cnt] = dr2;    
                    dr_3_vec[cnt] = dr_3;    
                                        
                    Ie_tidal[0] = (drn[0]*drn[0]-third);
                    Ie_tidal[3] = (drn[1]*drn[1]-third);
                    Ie_tidal[1] = drn[0]*drn[1];
                    Ie_tidal[2] = drn[0]*drn[2];
                    Ie_tidal[4] = drn[1]*drn[2];        

                    calc_done = true;
                }

                Xj  = -bi->m*bj->kf_R5*dr_3;

                for(int k=0; k<5; k++) {
                    bj->I_e_r[k] += Xj*Ie_tidal[k];
                }    
            }
            
            if(calc_done == false) {
                precalc_sep_vec_norm_save(*bi, *bj, dr, drn, dr_4, dr2, dr_3);

                dr2_vec[cnt] = dr2;    
                dr_3_vec[cnt] = dr_3;                
            }
            
            cnt++;
        }
    }
}


