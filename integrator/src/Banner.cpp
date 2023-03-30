#include "Banner.h"

void Banner::print_banner() {
    cout << "  _____  ___  ____ __   __ __  __  _____  ____  ____   " << endl;
    cout << " |_   _||_ _||  _ \\\\ \\ / /|  \\/  || ____|/ ___|/ ___|   " << endl;
    cout << "   | |   | | | | | |\\ V / | |\\/| ||  _|  \\___ \\\\___ \\  " << endl;
    cout << "   | |   | | | |_| | | |  | |  | || |___  ___) |___) | " << endl;
    cout << "   |_|  |___||____/  |_|  |_|  |_||_____||____/|____/  " << endl;
}

void Banner::print_intro(bool to_continue) {
    cout << endl;
    if(to_continue) {
        cout << "Continuing simulation ..." << endl;    
    }
    else {
        cout << "Starting new simulation ..." << endl;    
    }
    cout << " " << endl;
}

void Banner::print_outro(string output_dir) {
    cout << endl;
    cout << "Simulation finished!" << endl;    
    cout << "Output files available in directory: " << output_dir << endl;
}

void Banner::print_reference() {
    cout << endl;
    cout << "Citation: Boekholt, T. C. N. and Correia, A. C. M.            " << endl;
    cout << "          Monthly Notices of the Royal Astronomical Society   " << endl;
    cout << "          Volume 503, Issue 1, pp.xxxx-yyyy, November 2022   " << endl;
}

void Banner::print_stars() {
    cout << endl;
    cout << "**************************************************************" << endl;
}

void Banner::print_header() {
    cout << " " << endl;
    cout << "***************************************************************************" << endl;
    cout << "***  TIDYMESS 1.0 - TIdal DYnamics of Multi-body ExtraSolar Systems     ***" << endl;
    cout << "***  Authors      : Tjarda C. N. Boekholt and Alexandre C. M. Correia   ***" << endl;
    cout << "***  Contact      : tjardaboekholt@gmail.com, acm.correia@gmail.com     ***" << endl;
    cout << "***  Citation     : Boekholt, T. C. N. and Correia, A. C. M.            ***" << endl;
    cout << "***                 Monthly Notices of the Royal Astronomical Society   ***" << endl;
    cout << "***                 Volume 503, Issue 1, pp.xxxx-yyyy, November 2021    ***" << endl;
    cout << "***************************************************************************" << endl;
    cout << " " << endl;
}

                                                     

