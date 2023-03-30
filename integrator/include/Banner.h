#include <iostream>
using namespace std;

#ifndef __Banner_h
#define __Banner_h

class Banner {
    public:
    
    void print_banner();
    
    void print_intro(bool to_continue);
    void print_outro(string output_dir);

    void print_reference();
    void print_stars();

    void print_header();
};

#endif

