#include "Breakup.h"

// Initializers
Breakup::Breakup() {
    this->mode = 0;
    setup();

    breakup_id.clear();
    N_breakup = 0;
}

// Getters and Setters
void Breakup::setup() {
    to_detect = false;
    
    if(mode > 0) {
        to_detect = true;
    }
}

void Breakup::set_breakup_mode(int breakup_mode) {
    switch(breakup_mode) {
        case 0:
            this->mode = 0;
            break;
        case 1:
            this->mode = 1;
            break;
        case 2:
            this->mode = 2;
            break;
        default:
            cerr << "Invalid breakup mode: " << breakup_mode << endl;
            exit(1);
            break;
    }
}

int Breakup::get_N_breakup() {
    return N_breakup;
}
vector<int> Breakup::get_breakup_indices() {
    return breakup_id;
}

// Handlers
bool Breakup::detect_breakup(vector<Body> &bodies) {
    N_breakup = 0;
    breakup_id.clear();

    for(vector<Body>::iterator b = bodies.begin(); b != bodies.end()-1; ++b) {
        if(b->particle_type >= 2) {
            double w2 = inner_product(b->w.begin(), b->w.end(), b->w.begin(), 0.);
            double Re = 1.5 * b->R;
            double w2_breakup = b->m / (Re*Re*Re);
            if(w2 > w2_breakup) {
                int id = b-bodies.begin();
                breakup_id.push_back(id);
                N_breakup++;
            }
        }    
    }

    bool detection = false;
    if(N_breakup > 0) detection = true;
    
    return detection;
}    

// Printer    
void Breakup::print_id() {
    int Nc = breakup_id.size();
    for(int i=0; i<Nc; i++) {
        cout << breakup_id[i] << " ";
    }
    cout << endl;
}




