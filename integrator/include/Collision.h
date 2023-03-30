#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "Body.h"

#ifndef __Collision_h
#define __Collision_h

class Collision {
    public:

    // Variables and structures
    int collision_mode, roche_mode;
    bool to_detect_collision, to_detect_roche;

    vector< vector<int> > chain_id, chain_index;
    int Ng;

    // Initializers
    Collision();

    // Getters and Setters
    void set_collision_mode(int collision_mode);
    void set_roche_mode(int roche_mode);

    void setup();

    // Handlers
    void process_collision_chains(vector< array<int, 2> > &cindex);    
    void convert_id_to_index(vector<Body> &bodies);

    void replace_test_body(vector<Body> &bodies, int i, int Ns, vector<int> &index_sec);
    void replace_point_body(vector<Body> &bodies, int i, int Ns, vector<int> &index_sec);
    void replace_tidal_body(vector<Body> &bodies, int i, int Ns, vector<int> &index_sec);
    
    void replace(vector<Body> &bodies, vector< array<int, 2> > &cindex);

    // Printer    
    void print_indices(vector< array<int, 2> > &cindex);
};

#endif


