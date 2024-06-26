#include <iostream>
using namespace std;

#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iomanip> 
#include <cmath>
#include <array>

#include "Units.h"
#include "Banner.h"

#ifndef __Initializer_h
#define __Initializer_h

class Initializer {
    public:

    // Simulation parameters

    bool to_continue;		// 0=new simulation, 1=continue simulation
    double max_cpu_time;	// Maximum CPU running time in seconds. If 0 (default), then cpu time has no limit.

    // Physical model parameters

    int tidal_model;		// 0=none, 1=conservative, 2=linear, 3=creep direct, 4=creep tidymess (default)
    int pn_order;		// Post-Newtonian order: 0=none, 1=1pn, 2=1+2pn, 25=1+2+2.5pn    
    int magnetic_braking;	// Magnetic braking. 0=off, 1=on

    int collisions;		// 0=off, 1=flag, 2=exception, 3=replace
    int roche_limit;		// 0=off, 1=flag, 2=exception
    int breakup_speed;		// Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception

    // Unit system

    string mass_unit;		// Unit of mass in output: []=Nbody unit, [g], [kg], [Mearth], [Mjupiter], [Msun]
    string length_unit;	// Unit of length in output: []=Nbody unit, [m], [km], [Rsun], [AU], [pc]
    string time_unit;		// Unit of time for 1) t_begin, t_end and 2) unit of time in output: []=Nbody unit, [s], [hr], [day], [yr], [Myr], [Gyr]  
    string speed_unit;		// Unit of speed in output: []=Nbody unit, [m/s], [km/s], [km/hr], [km/hour], [AU/day], [AU/yr], [AU/year], [pc/Myr]

    double speed_of_light;	// Speed of light in N-body units. Only used in conjunction with N-body units and pn_order>0, otherwise equal to c.

    // Initial condition parameters

    string file_ic;		// initial condition file

    int orbit_coor;		// 0=cartesian inertial, 1=elliptical astrocentric, 2=elliptical jacobian
    int spin_coor;		// 0=absolute in inertial frame, 1=relative to its orbit; body 0 in the inertial frame, 2=relative to its orbit; body 0 relative to innermost orbit

    int initial_shape;		// 0=sphere, 1=equilibrium (future, 2=specified in initial condition)
    int num_body;		// 0=all, num_body+1=number of bodies to include  

    // Output parameters

    int snapshot_mode;		// 0=linear interval (default), 1=logarithmic interval, 2=variable interval 
    int n_snapshot;		// Total number of snapshots between t_begin and t_end (linear or in log10), or output a snapshot every fixed number (n_snapshot) of integration steps (variable) 

    string output_dir;		// Output directory; default is 'data/'. If left blank or set to '/', then file_ic will be adopted without the extension.  
    int overwrite;		// overwrite existing files: 0=no, 1=yes

    int output_format;		// 0=file per body, 1=file per snapshot, 2=single file
    int output_info;		// 0=time-varying quantities, 1=all quantities
    int output_coor;		// 0=cartesian inertial

    int output_diag;		// 0=no (default), 1=yes: output diagnostics, such as E and L, are written to a separate diagnostics file with extension '.diag'
    int output_terminal;	// Display progress of simulation in terminal window. 0=no, 1=yes

    // Integration parameters

    double t_begin;		// begin time in units given by time_unit 
    double t_end;		// final time in units given by time_unit

    int dt_mode;		// 0=constant dt, 1=adaptive dt, 2=adaptive, weighted dt
    double dt_const;		// constant time step in units given by time_unit (only used if dt_mode=0)
    double eta;		// accuracy parameter; timestep multiplication factor, default=0.0625 (only used if dt_mode>0)

    int n_iter;		// Number of iterations to improve reversibility (default=1)

    // Program parameters

    string file_par;		// Name of setup file with parameters. Set to "none" if no setup file is required.

    bool toSetOutputDir;	// Boolean for copying output dir from initial condition file
    string file_out;		// Initial condition file without the extension

    // Program variables

    Units units;
    bool physical_units;	// boolean for using physical units or N-body units    

    Banner banner;

    // Snapshot variables

    int N;
    vector< array<double, 15> > data_ic;    
    vector<int> id;
    vector<string> name;
    double t_end_phys;

    // Constructor
    Initializer();
    void set_to_default_values();
    
    double get_Cm();
    double get_Cr();
    double get_Cv();
    double get_Ct();
    
    // Readme
    void display_version();
    void display_help();
    void display_readme();
    void display_parameters();
    void display_units();
    void display_ic_args();
    
    // Argument parser
    void parse_arguments(int argc, char* argv[]);
    bool validate_command_line_args(int argc, vector<string> &args);
    bool validate_setup_file_args(int argc, vector<string> &args);
    bool check_valid_argument(string &str);
    bool check_valid_value(string &str, string &str_val);
    void process_args(int argc, vector<string> &args);

    void print_arguments_as_setup_file();
    void print_arguments_as_setup_file(ofstream &fo);
    void print_initial_condition_args();

    int check_var(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine, string x);

    bool validate_initial_condition_cartesian(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine);
    bool validate_initial_condition_elliptical(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine);
    bool validate_initial_condition_spin(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine);

    void determine_physical_units(vector<string> &var, vector<string> &unit);
    bool validate_units_internal_properties(vector<string> &var, vector<string> &unit);
    bool validate_units_orbital_cartesian(vector<string> &var, vector<string> &unit);
    bool validate_units_orbital_elliptical(vector<string> &var, vector<string> &unit);
    bool validate_units_spin(vector<string> &var, vector<string> &unit); 
            
    // Aux functions    
    inline bool check_file_exists(const std::string& path);
    bool check_dir_exists(string &dir);
    vector<string> get_words(string line);
    vector< vector<string> > read_file(string file_in);
    
    // Initial condition functions
    vector<double> rotZrotX(double anglez, double anglex, vector<double> vin);

    void convert_spin_vectors_to_inertial();
    void convert_spin_vectors_from_elliptical_body0abs();
    void convert_spin_vectors_from_elliptical_body0rel();
    
    double mean_to_eccentric_anomaly(double MA, double e);
    double eccentric_to_true_anomaly(double EA, double e);
    double convert_mean_to_true_anomaly(double MA, double ecc);  

    vector<double> get_relative_posvel_from_orbital_elements(double m1, double m2, double a, double ecc, double inc, double O, double w, double TA, double G);

    void move_to_center();
    void move_to_center(int Nenc);	

    void convert_astrocentric_elements_to_cartesian_coordinates();
    void convert_jacobian_elements_to_cartesian_coordinates();     
};

#endif

