#include <iostream>
using namespace std;

#include <string>
#include <cmath>
#include <vector>

#ifndef __Units_h
#define __Units_h

class Units {
    bool physical_units;

    public:

    // Constants
    double c_m_s, c_km_s, c_AU_yr, c_two_pi_AU_yr;
    double G_standard;

    // Time: s, hr, day, yr, Myr, Gyr
    double s_to_yr, min_to_yr, hr_to_yr, day_to_yr, Myr_to_yr, Gyr_to_yr;
    double yr_to_s, yr_to_min, yr_to_hr, yr_to_day, yr_to_Myr, yr_to_Gyr;

    // Length: m, km, Rsun/RSun, AU, pc
    double m_to_AU, km_to_AU, Rsun_to_AU, pc_to_AU;
    double AU_to_m, AU_to_km, AU_to_Rsun, AU_to_pc;

    double km_in_m, Rsun_in_m, AU_in_m, pc_in_m; 

    // Mass: g, kg, Mearth/MEarth, Mjupiter/MJupiter/Mjup/MJup, Msun/MSun
    double g_to_Msun, kg_to_Msun, Mearth_to_Msun, Mjupiter_to_Msun;
    double Msun_to_g, Msun_to_kg, Msun_to_Mearth, Msun_to_Mjupiter;

    double GMjupiter_in_m_s, GMearth_in_m_s, GMjupiter_in_AU_yr, GMearth_in_AU_yr;
    double Mjupiter_in_Msun, Mearth_in_Msun, Msun_in_kg, Msun_in_g;

    // Velocity: m/s, km/s, km/hr, AU/day, pc/Myr 
    double m_s_to_AU_yr, km_s_to_AU_yr, km_hr_to_AU_yr, AU_day_to_AU_yr, pc_Myr_to_AU_yr;
    double AU_yr_to_m_s, AU_yr_to_km_s, AU_yr_to_km_hr, AU_yr_to_AU_day, AU_yr_to_pc_Myr;    

    // Angles: deg, rad
    double rad_to_deg, deg_to_rad;

    // G=4pi <-> G=1 conversion
    double yr_to_yr_two_pi, yr_two_pi_to_yr;
    double AU_yr_to_two_pi_AU_yr, two_pi_AU_yr_to_AU_yr; 
    
    // Conversion factors to dimensionless units
    double Cm, Cr, Cv, Ct;
    
    Units();

    void set_physical_units(bool physical_units);
    
    bool validate_mass_unit(string u);
    bool validate_length_unit(string u);
    bool validate_time_unit(string u);
    bool validate_speed_unit(string u);
    bool validate_angular_unit(string u);

    double convert_mass_to_standard(double value, string unit);
    double convert_length_to_standard(double value, string unit);
    double convert_time_to_standard(double value, string unit);
    double convert_speed_to_standard(double value, string unit);
    double convert_angle_to_standard(double value, string unit);

    double convert_time_from_standard_to_code(double value);
    double convert_speed_from_standard_to_code(double value);
    double convert_frequency_from_standard_to_code(double value);

    double convert_time_from_code_to_output(double value, string unit);
    double convert_frequency_from_code_to_output(double value, string unit);
    double convert_mass_from_code_to_output(double value, string unit);
    double convert_length_from_code_to_output(double value, string unit);
    double convert_speed_from_code_to_output(double value, string unit);
    double convert_angular_momentum_from_code_to_output(double value, string unit);
    double convert_energy_from_code_to_output(double value, string unit);
    double convert_inertia_from_code_to_output(double value, string unit);

    vector<string> split_units(string &u);
    void print_units();
};

#endif


