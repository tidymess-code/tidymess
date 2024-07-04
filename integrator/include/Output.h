#include <iostream>
using namespace std;

#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include "Body.h"
#include "Initializer.h"
#include "Units.h"

#ifndef __Output_h
#define __Output_h

class Output {
    public:

    int format;
    int info;
    int coor;
    string dir;
    bool overwrite;

    Units units;
    bool physical_units;
    string time_unit_output;
    string mass_unit_output;
    string length_unit_output;
    string speed_unit_output;
    string frequency_unit_output;
    string angular_momentum_unit_output;
    string energy_unit_output;  
    string inertia_unit_output;
           
    int tidal_model;
    
    string file_out;
    string file_parbu;
        
    int output_class;
    
    bool firstWrite;
    bool firstWriteDiag;
    int file_counter;
    int bin_number;

    int n_width;

    // Constructor
    Output();

    // Functions
    void set_format(int output_format);
    void set_info(int output_info);
    void set_coor(int output_coor);
    void set_dir(string output_dir);
    void set_overwrite(int overwrite);

    void set_units(bool physical_units, string mass_unit, string length_unit, string time_unit, string speed_unit);
    void set_conversion_factors(double Cm, double Cr, double Cv, double Ct);
    void convert_time_to_frequency_unit(string time_unit);
    void construct_angular_momentum_unit(string mass_unit, string length_unit, string speed_unit);
    void construct_energy_unit(string mass_unit, string speed_unit);
    void construct_inertia_unit(string mass_unit, string length_unit);
 
    void set_tidal_model(int tidal_model);

    void set_file_out(string fo);

    void determine_output_mode();        
    
    void write_snapshot(double t, double tcpu, vector<Body> &bodies);
    void create_output_file(string outputFile, bool isFirstWrite);

    // Writers
    void write_snapshot_to_file(double t, double tcpu, vector<Body> &bodies);
    void write_snapshot_to_file_compact(double t, double tcpu, vector<Body> &bodies);
    void write_new_snapshot_file(double t, double tcpu, vector<Body> &bodies);
    void write_new_snapshot_file_compact(double t, double tcpu, vector<Body> &bodies);
    void write_snapshot_per_body(double t, double tcpu, vector<Body> &bodies);
    void write_snapshot_per_body_compact(double t, double tcpu, vector<Body> &bodies);

    void write_snapshot_to_file_nbody(double t, double tcpu, vector<Body> &bodies);
    void write_snapshot_to_file_compact_nbody(double t, double tcpu, vector<Body> &bodies);
    void write_new_snapshot_file_nbody(double t, double tcpu, vector<Body> &bodies);
    void write_new_snapshot_file_compact_nbody(double t, double tcpu, vector<Body> &bodies);
    void write_snapshot_per_body_nbody(double t, double tcpu, vector<Body> &bodies);
    void write_snapshot_per_body_compact_nbody(double t, double tcpu, vector<Body> &bodies);

    void save_to_binary(double &t, int &N, double &tcpu, double &dt_prev, double &t_end, int &num_integration_steps, int &collision_flag, int &roche_flag, int &breakup_flag, vector<Body> &bodies, double &dt_snapshot, double &dt0_log, double &fmul_log, int &num_snapshot);
    void read_binary_backup(double &t, int &N, double &tcpu, double &dt_prev, double &t_end, int &num_integration_steps, vector<Body> &bodies, int &collision_flag, int &roche_flag, int &breakup_flag, double &dt_snapshot, double &dt0_log, double &fmul_log, int &num_snapshot);
    void read_from_binary(double &t, int &N, double &tcpu, double &dt_prev, double &t_end, int &num_integration_steps, vector<Body> &bodies, string &bin_name, bool &myfirstWrite, bool &myfirstWriteDiag, int &myfile_counter, int &mybin_number, int &collision_flag, int &roche_flag, int &breakup_flag, double &dt_snapshot, double &dt0_log, double &fmul_log, int &num_snapshot);
            
    // Diagnostics log
    void write_log(int argc, char* argv[], double t, double tcpu, int &num_integration_steps, int N0, int N1, double dr, double dv, double dLx, double dLy, double dLz, double dLx_rel, double dLy_rel, double dLz_rel, double dL_abs, double dL_rel, double dE_abs, double dE_rel, double x0, int outcome_type, int collision_flag, int roche_flag, int breakup_flag);
    void print_log(int argc, char* argv[], double t, double tcpu, int &num_integration_steps, int N0, int N1, double dr, double dv, double dLx, double dLy, double dLz, double dLx_rel, double dLy_rel, double dLz_rel, double dL_abs, double dL_rel, double dE_abs, double dE_rel, double x0, int outcome_type, int collision_flag, int roche_flag, int breakup_flag);

    void write_diag(double &t, double &t_cpu, int &num_integration_steps, int &N, array<double, 3> &r0, array<double, 3> &v0, array<double, 3> &L0_orb, array<double, 3> &L0_spin, double &Ekin0_orb, double Epot0, double &Ekin0_spin); 

    void backup_file(string f1);
    void copy_file(string f1, string f2);
};

#endif

