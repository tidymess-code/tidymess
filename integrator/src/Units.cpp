#include "Units.h"

Units::Units() {
    // Constants
    c_m_s   = 299792458;
    
    c_km_s  = c_m_s / 1.e3;
    c_AU_yr = c_m_s * m_s_to_AU_yr;
    c_two_pi_AU_yr = c_AU_yr * AU_yr_to_two_pi_AU_yr;

    G_standard = 39.47841760435743; //39.48524924651484; // Gravitational constant in units of MSun, AU and yr
    
    // Time: s, min, hr, day, yr, Myr, Gyr
    s_to_yr   = 1./60/60/24/365.25;
    min_to_yr = 1./60/24/365.25;
    hr_to_yr  = 1./24/365.25;
    day_to_yr = 1./365.25;
    Myr_to_yr = 1.e6;
    Gyr_to_yr = 1.e9;

    yr_to_s   = 1./s_to_yr;
    yr_to_min = 1./min_to_yr;
    yr_to_hr  = 1./hr_to_yr;
    yr_to_day = 1./day_to_yr;
    yr_to_Myr = 1./Myr_to_yr;
    yr_to_Gyr = 1./Gyr_to_yr;

    // Length: m, km, Rsun/RSun, AU, pc
    km_in_m   = 1.e3;
    Rsun_in_m = 6.957e8;
    AU_in_m   = 149597870700;
    pc_in_m   = 180*60*60*AU_in_m / M_PI;

    m_to_AU    = 1/AU_in_m;
    km_to_AU   = 1e3/AU_in_m;
    Rsun_to_AU = Rsun_in_m/AU_in_m;
    pc_to_AU   = pc_in_m/AU_in_m;

    AU_to_m    = 1./m_to_AU;
    AU_to_km   = 1./km_to_AU;
    AU_to_Rsun = 1./Rsun_to_AU;
    AU_to_pc   = 1./pc_to_AU;

    // Mass: g, kg, Mearth/MEarth, Mjupiter/MJupiter/Mjup/MJup, Msun/MSun

    GMjupiter_in_m_s   = 1.2668653e17; // m^3 s^-2        
    GMearth_in_m_s     = 3.986004e14; // m^3 s^-2

    GMjupiter_in_AU_yr = GMjupiter_in_m_s * pow(AU_in_m, -3) * pow(s_to_yr, -2); // AU^3 yr^-2
    GMearth_in_AU_yr   = GMearth_in_m_s * pow(AU_in_m, -3) * pow(s_to_yr, -2); // AU^3 yr^-2

    Mjupiter_in_Msun   = GMjupiter_in_AU_yr / (4*M_PI*M_PI);
    Mearth_in_Msun     = GMearth_in_AU_yr / (4*M_PI*M_PI);

    Msun_in_kg         = 1.98847e30;
    Msun_in_g          = Msun_in_kg*1.e3;

    g_to_Msun          = 1./Msun_in_g;
    kg_to_Msun         = 1./Msun_in_kg;
    Mearth_to_Msun     = Mearth_in_Msun;
    Mjupiter_to_Msun   = Mjupiter_in_Msun;

    Msun_to_g          = 1./g_to_Msun;
    Msun_to_kg         = 1./kg_to_Msun;
    Msun_to_Mearth     = 1./Mearth_to_Msun;
    Msun_to_Mjupiter   = 1./Mjupiter_to_Msun;

    // Velocity: m/s, km/s, km/hr, AU/day, pc/Myr
    m_s_to_AU_yr    = m_to_AU / s_to_yr;
    km_s_to_AU_yr   = km_to_AU / s_to_yr;
    km_hr_to_AU_yr  = km_to_AU / hr_to_yr;
    AU_day_to_AU_yr = 1. / day_to_yr;
    pc_Myr_to_AU_yr = pc_to_AU / Myr_to_yr;

    AU_yr_to_m_s    = 1./m_s_to_AU_yr;
    AU_yr_to_km_s   = 1./km_s_to_AU_yr;
    AU_yr_to_km_hr  = 1./km_hr_to_AU_yr;
    AU_yr_to_AU_day = 1./AU_day_to_AU_yr;
    AU_yr_to_pc_Myr = 1./pc_Myr_to_AU_yr;    

    // Angles: deg, rad
    rad_to_deg = 1. / M_PI * 180;
    deg_to_rad = 1. / 180 * M_PI;
    
    // G=4pi <-> G=1 conversion
    yr_to_yr_two_pi = 2*M_PI;
    yr_two_pi_to_yr = 1./yr_to_yr_two_pi;

    AU_yr_to_two_pi_AU_yr = 1./(2*M_PI);
    two_pi_AU_yr_to_AU_yr = 1./AU_yr_to_two_pi_AU_yr; 

    // Conversion factors to dimensionless units
    Cm = 1.;
    Cr = 1.;
    Cv = 1.;
    Ct = 1.;
    
    physical_units = true;
}

void Units::set_physical_units(bool physical_units) {
    this->physical_units = physical_units;
}

bool Units::validate_mass_unit(string u) {
    if(u == "" || u == "[]") return true;
    else if(u == "[g]") return true;
    else if(u == "[kg]") return true;
    else if(u == "[Mearth]" || u == "[MEarth]") return true;
    else if(u == "[Mjupiter]" || u == "[MJupiter]") return true;
    else if(u == "[Mjup]" || u == "[MJup]") return true;
    else if(u == "[Msun]" || u == "[MSun]") return true;
    return false;
}
bool Units::validate_length_unit(string u) {
    if(u == "" || u == "[]") return true;
    else if(u == "[m]") return true;
    else if(u == "[km]") return true;
    else if(u == "[Rsun]" || u == "[RSun]") return true;
    else if(u == "[AU]" || u == "[au]") return true;
    else if(u == "[parsec]" || u == "[pc]") return true;
    return false;
}
bool Units::validate_time_unit(string u) {
    if(u == "" || u == "[]") return true;
    else if(u == "[s]") return true;
    else if(u == "[hour]" || u == "[hr]") return true;
    else if(u == "[day]") return true;
    else if(u == "[year]" || u == "[yr]") return true;
    else if(u == "[Myr]") return true;
    else if(u == "[Gyr]") return true;
    return false;
}
bool Units::validate_speed_unit(string u) {
    if(u == "" || u == "[]") return true;
    else if(u == "[m/s]") return true;
    else if(u == "[km/s]") return true;
    else if(u == "[km/hr]" || u == "[km/hour]") return true;
    else if(u == "[AU/day]" || u == "[au/day]") return true;
    else if(u == "[AU/yr]" || u == "[AU/year]") return true;
    else if(u == "[au/yr]" || u == "[au/year]") return true;
    else if(u == "[pc/Myr]") return true;
    return false;
}
bool Units::validate_angular_unit(string u) {
    if(u == "" || u == "[]") return true;
    else if(u == "[deg]") return true;
    else if(u == "[rad]") return true;
    return false;
}

double Units::convert_mass_to_standard(double value, string unit) {
    if(physical_units) {
        if(unit == "[g]") return value * g_to_Msun;
        else if(unit == "[kg]") return value * kg_to_Msun;
        else if(unit == "[Mearth]" || unit == "[MEarth]") return value * Mearth_to_Msun;
        else if(unit == "[Mjupiter]" || unit == "[MJupiter]") return value * Mjupiter_to_Msun;
        else if(unit == "[Mjup]" || unit == "[MJup]") return value * Mjupiter_to_Msun;
        else if(unit == "[Msun]" || unit == "[MSun]") return value;
    }
    return value;
}
double Units::convert_length_to_standard(double value, string unit) {
    if(physical_units) {
        if(unit == "[m]") return value * m_to_AU;
        else if(unit == "[km]") return value * km_to_AU;
        else if(unit == "[Rsun]" || unit == "[RSun]") return value * Rsun_to_AU;
        else if(unit == "[AU]" || unit == "[au]") return value;
        else if(unit == "[parsec]" || unit == "[pc]") return value * pc_to_AU;
    }
    return value;
}
double Units::convert_time_to_standard(double value, string unit) {
    if(physical_units) {
        if(unit == "[s]") return value * s_to_yr;
        else if(unit == "[hour]" || unit == "[hr]") return value * hr_to_yr;
        else if(unit == "[day]") return value * day_to_yr;
        else if(unit == "[year]" || unit == "[yr]") return value;
        else if(unit == "[Myr]") return value * Myr_to_yr;
        else if(unit == "[Gyr]") return value * Gyr_to_yr;
    }
    return value;
}
double Units::convert_speed_to_standard(double value, string unit) {
    if(physical_units) {
        if(unit == "[m/s]") return value * m_s_to_AU_yr;
        else if(unit == "[km/s]") return value * km_s_to_AU_yr;
        else if(unit == "[km/hr]" || unit == "[km/hour]") return value * km_hr_to_AU_yr;
        else if(unit == "[AU/day]" || unit == "[au/day]") return value * AU_day_to_AU_yr;
        else if(unit == "[AU/yr]" || unit == "[AU/year]") return value;
        else if(unit == "[au/yr]" || unit == "[au/year]") return value;
        else if(unit == "[pc/Myr]") return value * pc_Myr_to_AU_yr;
    }
    return value;
}
double Units::convert_angle_to_standard(double value, string unit) {
    if(unit == "" || unit == "[]") return value * deg_to_rad;
    else if(unit == "[deg]") return value * deg_to_rad;
    else if(unit == "[rad]") return value;
    return value;
}

double Units::convert_time_from_standard_to_code(double value) {
    if(physical_units) {
        return value * yr_to_yr_two_pi;
    }
    return value;
}
double Units::convert_speed_from_standard_to_code(double value) {
    if(physical_units) {
        return value * AU_yr_to_two_pi_AU_yr;
    }
    return value;
}
double Units::convert_frequency_from_standard_to_code(double value) {
    if(physical_units) {
        return value / yr_to_yr_two_pi;
    }
    return value;
}

double Units::convert_time_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Ct * val_G1_nodim;
        double val_standard = yr_two_pi_to_yr * val_G1_dim;
        if(unit == "[s]") return val_standard * yr_to_s;
        else if(unit == "[hour]" || unit == "[hr]") return val_standard * yr_to_hr;
        else if(unit == "[day]") return val_standard * yr_to_day;
        else if(unit == "[year]" || unit == "[yr]") return val_standard;
        else if(unit == "[Myr]") return val_standard * yr_to_Myr;
        else if(unit == "[Gyr]") return val_standard * yr_to_Gyr;
    }
    return value;
}
double Units::convert_frequency_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = val_G1_nodim / Ct;
        double val_standard = val_G1_dim / yr_two_pi_to_yr;
        if(unit == "[1/s]") return val_standard / yr_to_s;
        else if(unit == "[1/hour]" || unit == "[1/hr]") return val_standard / yr_to_hr;
        else if(unit == "[1/day]") return val_standard / yr_to_day;
        else if(unit == "[1/year]" || unit == "[1/yr]") return val_standard;
        else if(unit == "[1/Myr]") return val_standard / yr_to_Myr;
        else if(unit == "[1/Gyr]") return val_standard / yr_to_Gyr;
    }
    return value;
}
double Units::convert_mass_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Cm * val_G1_nodim;
        double val_standard = val_G1_dim;
        if(unit == "[g]") return val_standard * Msun_to_g;
        else if(unit == "[kg]") return val_standard * Msun_to_kg;
        else if(unit == "[Mearth]" || unit == "[MEarth]") return val_standard * Msun_to_Mearth;
        else if(unit == "[Mjupiter]" || unit == "[MJupiter]") return val_standard * Msun_to_Mjupiter;
        else if(unit == "[Mjup]" || unit == "[MJup]") return val_standard * Msun_to_Mjupiter;
        else if(unit == "[Msun]" || unit == "[MSun]") return val_standard;
    }
    return value;
}
double Units::convert_length_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Cr * val_G1_nodim;
        double val_standard = val_G1_dim;
        if(unit == "[m]") return val_standard * AU_to_m;
        else if(unit == "[km]") return val_standard * AU_to_km;
        else if(unit == "[Rsun]" || unit == "[RSun]") return val_standard * AU_to_Rsun;
        else if(unit == "[AU]" || unit == "[au]") return val_standard;
        else if(unit == "[parsec]" || unit == "[pc]") return val_standard * AU_to_pc;
    }
    return value;
}
double Units::convert_speed_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Cv * val_G1_nodim;
        double val_standard = two_pi_AU_yr_to_AU_yr * val_G1_dim;
        if(unit == "[m/s]") return val_standard * AU_yr_to_m_s;
        else if(unit == "[km/s]") return val_standard * AU_yr_to_km_s;
        else if(unit == "[km/hr]" || unit == "[km/hour]") return val_standard * AU_yr_to_km_hr;
        else if(unit == "[AU/day]" || unit == "[au/day]") return val_standard * AU_yr_to_AU_day;
        else if(unit == "[AU/yr]" || unit == "[AU/year]") return val_standard;
        else if(unit == "[au/yr]" || unit == "[au/year]") return val_standard;
        else if(unit == "[pc/Myr]") return val_standard * AU_yr_to_pc_Myr;
    }
    return value;
}
double Units::convert_angular_momentum_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Cm*Cr*Cv * val_G1_nodim;
        double val_standard = two_pi_AU_yr_to_AU_yr * val_G1_dim;

        vector<string> str = split_units(unit);

        for(int i=0; i<str.size(); i++) {
            if(str[i] == "g") val_standard *= Msun_to_g;
            else if(str[i] == "kg") val_standard *= Msun_to_kg;
            else if(str[i] == "Mearth" || str[i] == "MEarth") val_standard *= Msun_to_Mearth;
            else if(str[i] == "Mjupiter" || str[i] == "MJupiter") val_standard *= Msun_to_Mjupiter;
            else if(str[i] == "Mjup" || str[i] == "MJup") val_standard *= Msun_to_Mjupiter;
            else if(str[i] == "Msun" || str[i] == "MSun") val_standard *= 1;

            if(str[i] == "m") val_standard *= AU_to_m;
            else if(str[i] == "km") val_standard *= AU_to_km;
            else if(str[i] == "Rsun" || str[i] == "RSun") val_standard *= AU_to_Rsun;
            else if(str[i] == "AU" || str[i] == "au") val_standard *= 1;
            else if(str[i] == "parsec" || str[i] == "pc") val_standard *= AU_to_pc;

            if(str[i] == "m/s") val_standard *= AU_yr_to_m_s;
            else if(str[i] == "km/s") val_standard *= AU_yr_to_km_s;
            else if(str[i] == "km/hr" || str[i] == "km/hour") val_standard *= AU_yr_to_km_hr;
            else if(str[i] == "AU/day" || str[i] == "au/day") val_standard *= AU_yr_to_AU_day;
            else if(str[i] == "AU/yr" || str[i] == "AU/year") val_standard *= 1;
            else if(str[i] == "au/yr" || str[i] == "au/year") val_standard *= 1;
            else if(str[i] == "pc/Myr") val_standard *= AU_yr_to_pc_Myr;
        }
        return val_standard;
    }
    return value;
}
double Units::convert_energy_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Cm*Cv*Cv * val_G1_nodim;
        double val_standard = two_pi_AU_yr_to_AU_yr*two_pi_AU_yr_to_AU_yr * val_G1_dim;

        vector<string> str = split_units(unit);
        
        for(int i=0; i<str.size(); i++) {
            if(str[i] == "g") val_standard *= Msun_to_g;
            else if(str[i] == "kg") val_standard *= Msun_to_kg;
            else if(str[i] == "Mearth" || str[i] == "MEarth") val_standard *= Msun_to_Mearth;
            else if(str[i] == "Mjupiter" || str[i] == "MJupiter") val_standard *= Msun_to_Mjupiter;
            else if(str[i] == "Mjup" || str[i] == "MJup") val_standard *= Msun_to_Mjupiter;
            else if(str[i] == "Msun" || str[i] == "MSun") val_standard *= 1;

            if(str[i] == "m/s") val_standard *= AU_yr_to_m_s;
            else if(str[i] == "km/s") val_standard *= AU_yr_to_km_s;
            else if(str[i] == "km/hr" || str[i] == "km/hour") val_standard *= AU_yr_to_km_hr;
            else if(str[i] == "AU/day" || str[i] == "au/day") val_standard *= AU_yr_to_AU_day;
            else if(str[i] == "AU/yr" || str[i] == "AU/year") val_standard *= 1;
            else if(str[i] == "au/yr" || str[i] == "au/year") val_standard *= 1;
            else if(str[i] == "pc/Myr") val_standard *= AU_yr_to_pc_Myr;
        }
        return val_standard;
    }
    return value;
}

double Units::convert_inertia_from_code_to_output(double value, string unit) {
    if(physical_units) {
        double val_G1_nodim = value;
        double val_G1_dim = Cm*Cr*Cr * val_G1_nodim;
        double val_standard = val_G1_dim;

        vector<string> str = split_units(unit);

        for(int i=0; i<str.size(); i++) {
            if(str[i] == "g") val_standard *= Msun_to_g;
            else if(str[i] == "kg") val_standard *= Msun_to_kg;
            else if(str[i] == "Mearth" || str[i] == "MEarth") val_standard *= Msun_to_Mearth;
            else if(str[i] == "Mjupiter" || str[i] == "MJupiter") val_standard *= Msun_to_Mjupiter;
            else if(str[i] == "Mjup" || str[i] == "MJup") val_standard *= Msun_to_Mjupiter;
            else if(str[i] == "Msun" || str[i] == "MSun") val_standard *= 1;

            if(str[i] == "m") val_standard *= AU_to_m;
            else if(str[i] == "km") val_standard *= AU_to_km;
            else if(str[i] == "Rsun" || str[i] == "RSun") val_standard *= AU_to_Rsun;
            else if(str[i] == "AU" || str[i] == "au") val_standard *= 1;
            else if(str[i] == "parsec" || str[i] == "pc") val_standard *= AU_to_pc;
        }
        return val_standard;
    }
    return value;
}

vector<string> Units::split_units(string &u) {
    vector<string> v;
    
    vector<int> ind;
    ind.push_back(-1);
        
    int N_str = u.length();
    for(int k=0; k<N_str; k++) {
        if( isspace( u[k]) ) {
            ind.push_back(k);
        }
    }
    
    ind.push_back(N_str);
        
    int Np = ind.size();    
    for(int k=0; k<Np-1; k++) {
        int i0 = ind[k]+1;
        int di = ind[k+1]-(ind[k]+1);
        string s = u.substr(i0, di);            

        v.push_back(s);
    }

    Np = v.size();

    v[0] = v[0].substr(1, v[0].size()-1);
    v[Np-1] = v[Np-1].substr(0, v[Np-1].size()-1);
        
    return v;
}    
void Units::print_units() {
    cout << "Mass units:" << endl;
  
    cout << "[]         = N-body unit" << endl;
    cout << "[g]        = gram" << endl;
    cout << "[kg]       = kilogram" << endl;
    cout << "[Mearth]   = Earth mass" << endl;
    cout << "[MEarth]   = Earth mass" << endl;
    cout << "[Mjupiter] = Jupiter mass" << endl;
    cout << "[MJupiter] = Jupiter mass" << endl;
    cout << "[Mjup]     = Jupiter mass" << endl;
    cout << "[MJup]     = Jupiter mass" << endl;
    cout << "[Msun]     = Solar mass" << endl;
    cout << "[MSun]     = Solar mass" << endl;

    cout << " " << endl;

    cout << "Length units:" << endl;

    cout << "[]         = N-body unit" << endl;
    cout << "[m]        = meter" << endl;
    cout << "[km]       = kilometer" << endl;
    cout << "[Rsun]     = Solar radius" << endl;
    cout << "[RSun]     = Solar radius" << endl;
    cout << "[AU]       = Astronomical Unit" << endl;
    cout << "[au]       = Astronomical Unit" << endl;
    cout << "[parsec]   = parsec" << endl;
    cout << "[pc]       = parsec" << endl;

    cout << " " << endl;

    cout << "Time units:" << endl;

    cout << "[]         = N-body unit" << endl;
    cout << "[s]        = second" << endl;
    cout << "[hour]     = hour" << endl;
    cout << "[hr]       = hour" << endl;
    cout << "[day]      = day" << endl;
    cout << "[year]     = year" << endl;
    cout << "[yr]       = year" << endl;
    cout << "[Myr]      = Million years" << endl;
    cout << "[Gyr]      = Billion years" << endl;

    cout << " " << endl;

    cout << "Speed units:" << endl;

    cout << "[]         = N-body unit" << endl;
    cout << "[m/s]      = meters per second" << endl;
    cout << "[km/s]     = kilometers per second" << endl;
    cout << "[km/hr]    = kilometers per hour" << endl;
    cout << "[km/hour]  = kilometers per hour" << endl;
    cout << "[AU/day]   = Astronomical Units per day" << endl;
    cout << "[au/day]   = Astronomical Units per day" << endl;
    cout << "[AU/yr]    = Astronomical Units per year" << endl;
    cout << "[au/yr]    = Astronomical Units per year" << endl;
    cout << "[AU/year]  = Astronomical Units per year" << endl;
    cout << "[au/year]  = Astronomical Units per year" << endl;
    cout << "[pc/Myr]   = parsecs per million years" << endl;

    cout << " " << endl;
    
    cout << "Angular units:" << endl;
    
    cout << "[]         = degrees" << endl;
    cout << "[deg]      = degrees" << endl;
    cout << "[rad]      = radians" << endl;
}





