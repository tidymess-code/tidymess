#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>

#include "Body.h"

#ifndef __Breakup_h
#define __Breakup_h

class Breakup {
    public:

    // Variables and structures
    bool to_detect;
    int mode;

    vector<int> breakup_id;
    int N_breakup;

    // Initializers
    Breakup();

    // Getters and Setters
    void set_breakup_mode(int breakup_mode);
    void setup();

    int get_N_breakup();
    vector<int> get_breakup_indices();

    // Handlers
    bool detect_breakup(vector<Body> &bodies);    

    // Printer    
    void print_id();
};

#endif


