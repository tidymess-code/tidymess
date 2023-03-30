#include "Collision.h"

// Initializers
Collision::Collision() {
    this->collision_mode = 0;
    this->roche_mode = 0;
    setup();
}

// Getters and Setters
void Collision::set_collision_mode(int collision_mode) {
    switch(collision_mode) {
        case 0:
            this->collision_mode = 0;
            break;
        case 1:
            this->collision_mode = 1;
            break;
        case 2:
            this->collision_mode = 2;
            break;
        case 3:
            this->collision_mode = 3;
            break;
        default:
            cerr << "Invalid collision mode: " << collision_mode << endl;
            exit(1);
            break;
    }
}
void Collision::set_roche_mode(int roche_mode) {
    switch(roche_mode) {
        case 0:
            this->roche_mode = 0;
            break;
        case 1:
            this->roche_mode = 1;
            break;
        case 2:
            this->roche_mode = 2;
            break;
        default:
            cerr << "Invalid roche mode: " << roche_mode << endl;
            exit(1);
            break;
    }
}

void Collision::setup() {
    if(collision_mode > 0) {
        to_detect_collision = true;
    }
    else {
        to_detect_collision = false;
    }

    if(roche_mode > 0) {
        to_detect_roche = true;
    }
    else {
        to_detect_roche = false;
    }
}
    
// Handlers
void Collision::process_collision_chains(vector< array<int, 2> > &cindex) {
    int Nc = cindex.size();

    chain_id.clear();

    vector<int> couple(2);
    couple[0] = cindex[0][0];
    couple[1] = cindex[0][1];

    chain_id.push_back(couple);

    for(int i=1; i<Nc; i++) {

        bool isNewChain = false;

        for(int j=0; j<chain_id.size(); j++) {
            bool isPresent0 = false, isPresent1 = false;

            if( find(chain_id[j].begin(), chain_id[j].end(), cindex[i][0]) != chain_id[j].end() ) {
                isPresent0 = true; 
            }
            if( find(chain_id[j].begin(), chain_id[j].end(), cindex[i][1]) != chain_id[j].end() ) {
                isPresent1 = true; 
            }

            if(isPresent0 == true && isPresent1 == true) {
                isNewChain = false;
                break;
            }
            else if(isPresent0 == true && isPresent1 == false) {
                chain_id[j].push_back(cindex[i][1]);
                isNewChain = false;
                break;
            }
            else if(isPresent0 == false && isPresent1 == true) {
                chain_id[j].push_back(cindex[i][0]);
                isNewChain = false;
                break;
            }
            else {
                isNewChain = true;
            }
   
        }

        if(isNewChain == true) {
            couple[0] = cindex[i][0];
            couple[1] = cindex[i][1];
            chain_id.push_back(couple);
        }
    }
}    
void Collision::convert_id_to_index(vector<Body> &bodies) {
    this->Ng = this->chain_id.size();
    this->chain_index.clear();

    // Get chains and ids of collision partners
    int N = bodies.size();
    for(int i=0; i<Ng; i++) {
        int Ns = chain_id[i].size();
        vector<int> mychain(Ns, 0);
        
        for(int j=0; j<Ns; j++) {
            for(int k=0; k<N; k++) {
                if(bodies[k].id == chain_id[i][j]) {
                    mychain[j] = k;
                    break;
                }
            }
        }
        
        this->chain_index.push_back(mychain);
    }        
}

void Collision::replace_test_body(vector<Body> &bodies, int i, int Ns, vector<int> &index_sec) {
    // Get center of mass quantities

    double Mcm = 0;
    array<double, 3> rcm = {};
    array<double, 3> vcm = {};

    for(int j=0; j<Ns; j++) {
        for(int k=0; k<3; k++) {
            rcm[k] += bodies[chain_index[i][j]].r[k];
            vcm[k] += bodies[chain_index[i][j]].v[k];
        }
    }

    for(int k=0; k<3; k++) {
        rcm[k] /= Ns;
        vcm[k] /= Ns;
    }

    // Identify largest collision partner

    int index_max = 0;
    double Rmax = -1;

    for(int j=0; j<Ns; j++) {
        if(bodies[chain_index[i][j]].R > Rmax) {
            Rmax = bodies[chain_index[i][j]].R;
            index_max = j; 
        }
    }

    // Derive new internal properties based on assumptions

    // New radius is largest radius
    double Rcm = Rmax;

    // Update primary collision partner
    bodies[chain_index[i][index_max]].m = Mcm;
    bodies[chain_index[i][index_max]].R = Rcm;
       
    for(int k=0; k<3; k++) {
        bodies[chain_index[i][index_max]].r[k]  = rcm[k];
        bodies[chain_index[i][index_max]].v[k]  = vcm[k];
    }        
         
    // Set auxiliary quantities
    bodies[chain_index[i][index_max]].update_aux_properties(); 
         
    // Book keep indices of secondary collision partners  
    for(int j=0; j<Ns; j++) {
        if(j != index_max) {   
            index_sec.push_back(chain_index[i][j]);
        }
    }  
}
void Collision::replace_point_body(vector<Body> &bodies, int i, int Ns, vector<int> &index_sec) {
    // Get center of mass quantities

    double Mcm = 0;
    array<double, 3> rcm = {};
    array<double, 3> vcm = {};

    for(int j=0; j<Ns; j++) {
        Mcm += bodies[chain_index[i][j]].m;
        for(int k=0; k<3; k++) {
            rcm[k] += bodies[chain_index[i][j]].m*bodies[chain_index[i][j]].r[k];
            vcm[k] += bodies[chain_index[i][j]].m*bodies[chain_index[i][j]].v[k];
        }
    }

    for(int k=0; k<3; k++) {
        rcm[k] /= Mcm;
        vcm[k] /= Mcm;
    }

    // Identify most massive collision partner

    int index_max = 0;
    double Mmax = 0;

    for(int j=0; j<Ns; j++) {
        if(bodies[chain_index[i][j]].m > Mmax) {
            Mmax = bodies[chain_index[i][j]].m;
            index_max = j; 
        }
    }

    // Derive new internal properties based on assumptions

    // New radius based on constancy of density
    double Rcm = 0;

    if(bodies[chain_index[i][index_max]].R > 0) {
        double rho = bodies[chain_index[i][index_max]].rho;
        if(rho <= 0) {
            bodies[chain_index[i][index_max]].update_aux_properties();                
            rho = bodies[chain_index[i][index_max]].rho;
        }
        Rcm = pow(3*Mcm / rho / (4*M_PI), 1./3);                    
    }

    // Update primary collision partner
    bodies[chain_index[i][index_max]].m = Mcm;
    bodies[chain_index[i][index_max]].R = Rcm;
       
    for(int k=0; k<3; k++) {
        bodies[chain_index[i][index_max]].r[k]  = rcm[k];
        bodies[chain_index[i][index_max]].v[k]  = vcm[k];
    }        
         
    // Set auxiliary quantities
    bodies[chain_index[i][index_max]].update_aux_properties(); 
         
    // Book keep indices of secondary collision partners  
    for(int j=0; j<Ns; j++) {
        if(j != index_max) {   
            index_sec.push_back(chain_index[i][j]);
        }
    }  
}
void Collision::replace_tidal_body(vector<Body> &bodies, int i, int Ns, vector<int> &index_sec) {
    // Get center of mass quantities

    double Mcm = 0;
    array<double, 3> rcm = {};
    array<double, 3> vcm = {};

    for(int j=0; j<Ns; j++) {
        Mcm += bodies[chain_index[i][j]].m;
        for(int k=0; k<3; k++) {
            rcm[k] += bodies[chain_index[i][j]].m*bodies[chain_index[i][j]].r[k];
            vcm[k] += bodies[chain_index[i][j]].m*bodies[chain_index[i][j]].v[k];
        }
    }

    for(int k=0; k<3; k++) {
        rcm[k] /= Mcm;
        vcm[k] /= Mcm;
    }

    // Got total L involved

    array<double, 3> Lorb = {};
    array<double, 3> Lspin = {};

    for(int j=0; j<Ns; j++) {
        int q = chain_index[i][j];

        if(bodies[q].particle_type > 0) {
            double Lx_orb = bodies[q].m*(bodies[q].r[1]*bodies[q].v[2] - bodies[q].r[2]*bodies[q].v[1]);
            double Ly_orb = bodies[q].m*(bodies[q].r[2]*bodies[q].v[0] - bodies[q].r[0]*bodies[q].v[2]);
            double Lz_orb = bodies[q].m*(bodies[q].r[0]*bodies[q].v[1] - bodies[q].r[1]*bodies[q].v[0]);

            Lorb[0] += Lx_orb;
            Lorb[1] += Ly_orb;
            Lorb[2] += Lz_orb;
        }
        if(bodies[q].particle_type >= 2) {        
            double Lx_spin = 0;
            double Ly_spin = 0;
            double Lz_spin = 0;

            Lx_spin = bodies[q].I[0]*bodies[q].w[0] + bodies[q].I[1]*bodies[q].w[1] + bodies[q].I[2]*bodies[q].w[2];
            Ly_spin = bodies[q].I[1]*bodies[q].w[0] + bodies[q].I[3]*bodies[q].w[1] + bodies[q].I[4]*bodies[q].w[2];
            Lz_spin = bodies[q].I[2]*bodies[q].w[0] + bodies[q].I[4]*bodies[q].w[1] + bodies[q].I[5]*bodies[q].w[2];

            Lspin[0] += Lx_spin;
            Lspin[1] += Ly_spin;
            Lspin[2] += Lz_spin;
        }
    }    
        
    array<double, 3> Ltot = {};
  
    for(int k=0; k<3; k++) {
        Ltot[k] += Lorb[k] + Lspin[k];
    }

    // Get new spin angular momentum

    array<double, 3> Lorb_new = {};

    Lorb_new[0] = Mcm*(rcm[1]*vcm[2] - rcm[2]*vcm[1]);
    Lorb_new[1] = Mcm*(rcm[2]*vcm[0] - rcm[0]*vcm[2]);
    Lorb_new[2] = Mcm*(rcm[0]*vcm[1] - rcm[1]*vcm[0]);
           
    array<double, 3> Lspin_new = {};

    for(int k=0; k<3; k++) {
        Lspin_new[k] = Ltot[k]-Lorb_new[k];
    }

    // Identify most massive collision partner

    int index_max = 0;
    double Mmax = 0;

    for(int j=0; j<Ns; j++) {
        if(bodies[chain_index[i][j]].m > Mmax) {
            Mmax = bodies[chain_index[i][j]].m;
            index_max = j; 
        }
    }

    // Derive new internal properties based on assumptions

    // New radius
    double Rcm = 0;

    if(bodies[chain_index[i][index_max]].R > 0) {
        double rho = bodies[chain_index[i][index_max]].rho;
        if(rho <= 0) {
            bodies[chain_index[i][index_max]].update_aux_properties();                
            rho = bodies[chain_index[i][index_max]].rho;
        }
        Rcm = pow(3*Mcm / rho / (4*M_PI), 1./3);                    
    }

    // Assume xi remains unchanged
    // Assume kf remains unchanged
    // Assume tau remains unchanged
    // Assume a_mb remains unchanged

    // Shape: assume sphere after collision has relaxed
    double Imag = bodies[chain_index[i][index_max]].xi * Mcm * Rcm * Rcm;

    array<double, 6> Icm_p = {};
    Icm_p.fill(0);
    Icm_p[0] = Imag;
    Icm_p[3] = Imag;
    Icm_p[5] = Imag;

    array<double, 6> Icm_n = {};
    Icm_n.fill(0);

    array<double, 6> Icm = {};
    for(int k=0; k<6; k++) {
        Icm[k] = Icm_p[k] + Icm_n[k];
    }

    array<double, 6> Icm_inv = {};
    if(Imag > 0) {
        Icm_inv.fill(0);
        Icm_inv[0] = 1./Imag;
        Icm_inv[3] = 1./Imag;
        Icm_inv[5] = 1./Imag;
    }

    // Spin
    array<double, 3> wcm = {};

    wcm[0] = Icm_inv[0]*Lspin_new[0] + Icm_inv[1]*Lspin_new[1] + Icm_inv[2]*Lspin_new[2];
    wcm[1] = Icm_inv[1]*Lspin_new[0] + Icm_inv[3]*Lspin_new[1] + Icm_inv[4]*Lspin_new[2];
    wcm[2] = Icm_inv[2]*Lspin_new[0] + Icm_inv[4]*Lspin_new[1] + Icm_inv[5]*Lspin_new[2];

    // Update primary collision partner
    bodies[chain_index[i][index_max]].m = Mcm;
    bodies[chain_index[i][index_max]].R = Rcm;

    for(int k=0; k<6; k++) {
        bodies[chain_index[i][index_max]].I_p[k] = Icm_p[k];
        bodies[chain_index[i][index_max]].I_n[k] = Icm_n[k];
        bodies[chain_index[i][index_max]].I[k] = Icm[k];
        bodies[chain_index[i][index_max]].I_inv[k] = Icm_inv[k];
    }
        
    for(int k=0; k<3; k++) {
        bodies[chain_index[i][index_max]].w[k] = wcm[k];
        bodies[chain_index[i][index_max]].L[k] = Lspin_new[k];
    }
       
    for(int k=0; k<3; k++) {
        bodies[chain_index[i][index_max]].r[k]  = rcm[k];
        bodies[chain_index[i][index_max]].v[k]  = vcm[k];
    }        
         
    // Set auxiliary quantities
    bodies[chain_index[i][index_max]].update_aux_properties(); 
         
    // Book keep indices of secondary collision partners  
    for(int j=0; j<Ns; j++) {
        if(j != index_max) {   
            index_sec.push_back(chain_index[i][j]);
        }
    }  
}
void Collision::replace(vector<Body> &bodies, vector< array<int, 2> > &cindex) {
    process_collision_chains(cindex);
    convert_id_to_index(bodies);

    vector<int> index_sec;

    for(int i=0; i<Ng; i++) {
        int Ns = chain_index[i].size();
        
        int collision_type = 0;

        for(int j=0; j<Ns; j++) {
            if( bodies[chain_index[i][j]].particle_type >= 2 ) {
                collision_type = 2;
                break;
            }
            else if( bodies[chain_index[i][j]].particle_type == 1 ) {
                collision_type = 1;
            }
        }
        
        if(collision_type == 0) {
            replace_test_body(bodies, i, Ns, index_sec);
        }
        else if(collision_type == 1) {
            replace_point_body(bodies, i, Ns, index_sec);
        }
        else {
            replace_tidal_body(bodies, i, Ns, index_sec);
        }
    }

    // Remove merger components
    std::sort(index_sec.begin(), index_sec.end(), std::greater<int>());

    for(int i=0; i<index_sec.size(); i++) {
        bodies.erase(bodies.begin()+index_sec[i]);
    }
}

// Printer    
void Collision::print_indices(vector< array<int, 2> > &cindex) {
    int Nc = cindex.size();
    for(int i=0; i<Nc; i++) {
        cout << cindex[i][0] << " " << cindex[i][1] << endl;
    }
    cout << endl;
}



