#include "Initializer.h"

// Constructor
Initializer::Initializer() {
    set_to_default_values();
}

void Initializer::set_to_default_values() {
    // Simulation parameters

    to_continue	= false;		// 0=new simulation, 1=continue simulation
    max_cpu_time	= 0;			// Maximum CPU running time in seconds. If 0 (default), then cpu time has no limit.

    // Physical model parameters

    tidal_model	= 4;			// 0=none, 1=conservative, 2=linear, 3=creep direct, 4=creep tidymess (default)
    pn_order		= 0;			// Post-Newtonian order: 0=none, 1=1pn, 2=1+2pn, 25=1+2+2.5pn    
    magnetic_braking	= 0;			// Magnetic braking. 0=off, 1=on

    collisions		= 0;			// 0=off, 1=flag, 2=exception, 3=replace
    roche_limit	= 0;			// 0=off, 1=flag, 2=exception
    breakup_speed	= 0;			// Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception
    
    // Unit system

    mass_unit		= "[Msun]";		// Unit of mass in output: []=Nbody unit, [g], [kg], [Mearth], [Mjupiter], [Msun]
    length_unit	= "[au]";		// Unit of length in output: []=Nbody unit, [m], [km], [Rsun], [au], [pc]
    time_unit		= "[yr]";		// Unit of time for 1) t_begin, t_end and 2) unit of time in output: []=Nbody unit, [s], [hr], [day], [yr], [Myr], [Gyr]  
    speed_unit		= "[au/yr]";		// Unit of speed in output: []=Nbody unit, [m/s], [km/s], [km/hr], [km/hour], [au/day], [au/yr], [au/year], [pc/Myr]

    speed_of_light	= 0;			// Speed of light in N-body units. Only used in conjunction with N-body units and pn_order>0, otherwise equal to c.

    // Initial condition parameters

    file_ic		= "tidymess.ic";	// initial condition file

    orbit_coor		= 1;			// 0=cartesian inertial, 1=elliptical astrocentric, 2=elliptical jacobian
    spin_coor		= 2;			// 0=absolute in inertial frame, 1=relative to its orbit; body 0 in the inertial frame, 2=relative to its orbit; body 0 relative to innermost orbit

    initial_shape	= 0;			// 0=sphere, 1=equilibrium (future, 2=specified in initial condition)
    num_body		= 0;			// 0=all, num_body+1=number of bodies to include  

    // Output parameters

    snapshot_mode      = 0;			// 0=linear interval (default), 1=logarithmic interval, 2=variable interval 
    n_snapshot         = 100;			// Total number of snapshots between t_begin and t_end (linear or in log10), or output a snapshot every fixed number (n_snapshot) of integration steps (variable) 

    output_dir		= "data/";		// Output directory; default is 'data/'. If left blank or set to '/', then file_ic will be adopted without the extension.  
    overwrite		= 1;			// overwrite existing files: 0=no, 1=yes

    output_format	= 0;			// 0=file per body, 1=file per snapshot, 2=single file
    output_info	= 1;			// 0=time-varying quantities, 1=all quantities
    output_coor	= 0;			// 0=cartesian inertial

    output_diag        = 0;			// 0=no (default), 1=yes: output diagnostics, such as E and L, are written to a separate diagnostics file with extension '.diag'
    output_terminal	= 1;			// Display progress of simulation in terminal window. 0=no, 1=yes

    // Integration parameters

    t_begin		= 0;			// begin time in units given by time_unit 
    t_end		= 1e2;			// final time in units given by time_unit

    dt_mode		= 2;			// 0=constant dt, 1=adaptive dt, 2=adaptive, weighted dt
    dt_const		= 0.015625;		// constant time step in units given by time_unit (only used if dt_mode=0)
    eta		= 0.0625;		// accuracy parameter; timestep multiplication factor, default=0.125 (only used if dt_mode>0)

    n_iter		= 1;			// Number of iterations to improve reversibility (default=1)

    // Program parameters

    file_par		= "tidymess.par";	// Name of setup file with parameters. Set to "none" if no setup file is required.

    toSetOutputDir	= false;		// Boolean for copying output dir from initial condition file

    Units units;
    physical_units	= true;		// boolean for using physical units or N-body units    

    Banner banner;

    // Snapshot variables

    N = 0;
    data_ic.clear();    
    id.clear();
    name.clear();
    
    t1_code = 1;    
}

double Initializer::get_Cm() {
    return units.Cm;
}
double Initializer::get_Cr() {
    return units.Cr;
}
double Initializer::get_Cv() {
    return units.Cv;
}
double Initializer::get_Ct() {
    return units.Ct;
}

// Readme
void Initializer::display_version() {
    banner.print_header();
}
void Initializer::display_help() {
    cout << " " << endl;
    display_readme();
    cout << " " << endl;
}
void Initializer::display_readme() {
    display_version();

    cout << "To run TIDYMESS successfully, please make sure to have " << endl;
    cout << "  1) a proper initial condition file (e.g. tidymess.ic), and " << endl;
    cout << "  2) a valid setup file (e.g. tidymess.par) and/or command line arguments." << endl;
    cout << " " << endl;
    cout << "Then run ./tidymess.exe, which by default will look for tidymess.par. " << endl;
    cout << "User specified parameter files can be used by ./tidymess.exe -f name_of_par_file." << endl;
    cout << "The initial condition file is specified in the setup file, or on the command line." << endl;
    cout << " " << endl;
    cout << "The default output directory is 'data/', but this can be changed by setting the output_dir parameter." << endl;
    cout << "When set to '/' or left blank, then a folder will be created with the name of the initial condition file." << endl;
    cout << " " << endl;
    cout << "For more info, run ./tidymess.exe with the following flags:" << endl;
    cout << "--parameters or -p to see the list of simulation parameters." << endl;
    cout << "--ic to see an overview of properties to be assigned to the bodies in the initial condition file." << endl;
    cout << "--units or -u to see the list of recognized units in TIDYMESS." << endl;
    cout << " " << endl;
    cout << "Please see the examples folder for TIDYMESS tutorials." << endl;
    cout << endl;
}
void Initializer::display_parameters() {
    cout << " " << endl;
    cout << "Parameters to be specified in the setup file (e.g. tidymess.par)." << endl;
    cout << "They can also be set manually in the command line. In that case, add '--' in front and remove '=' (e.g. --t_begin 0)" << endl;
    cout << "If a value is defined in both the setup file and command line, then TIDYMESS adopts the command line value." << endl;
    cout << "The values below are the default values. " << endl;
    cout << " " << endl;
    print_arguments_as_setup_file();
    cout << " " << endl;
}
void Initializer::display_units() {
    cout << " " << endl;
    cout << "The following units are recognized by TIDYMESS, and should be used in the" << endl;
    cout << "parameter setup file (e.g. tidymess.par) and initial condition file (e.g. tidymess.ic)." << endl;
    cout << " " << endl;
    units.print_units();
    cout << " " << endl;
}
void Initializer::display_ic_args() {
    cout << " " << endl;
    print_initial_condition_args();
    cout << " " << endl;
    cout << "Use the -u or --units flag to see the list of recognized units in TIDYMESS." << endl;
    cout << " " << endl;
}

// Argument parser
void Initializer::parse_arguments(int argc, char* argv[]) {
    // Convert command line arguments to strings
    vector<string> args(argv+1, argv+argc);
    argc = args.size();

    // Check for support flags
    if(argc > 0) {
        for(int i=0; i<argc; i++) {
            if(args[i] == "--help" || args[i] == "-h") {
                display_readme();
                exit(1);
            }
        }
        for(int i=0; i<argc; i++) {
            if(args[i] == "--parameters" || args[i] == "-p") {
                display_parameters();
                exit(1);
            }
        }
        for(int i=0; i<argc; i++) {
            if(args[i] == "--ic") {
                display_ic_args();
                exit(1);
            }
        }
        for(int i=0; i<argc; i++) {
            if(args[i] == "--units" || args[i] == "-u") {
                display_units();
                exit(1);
            }
        }
        for(int i=0; i<argc; i++) {
            if(args[i] == "--version" || args[i] == "-v") {
                display_version();
                exit(1);
            }
        }
    }

    // Validate command line arguments
    bool valid_command_line_args = validate_command_line_args(argc, args);
    if(!valid_command_line_args) {
        exit(1);
    }

    // Determine setup file usage
    for(int i=0; i<argc; i++) {
        if(args[i] == "--file_par" || args[i] == "-f") {
            string f = args[i+1];
            if(f == "none" || f == "None") {
                file_par = "none";        
            }
            else {
                file_par = f;
            }            
            break;
        }
    }    

    // Read in, validate and process setup file
    if(file_par != "none") {

        // Read setup file
        vector< vector<string> > setup_line = read_file(file_par);

        vector<string> args_f;
        int numLine = setup_line.size();
        for(int i=0; i<numLine; i++) {
            int numWord = setup_line[i].size();
            if(numWord >= 3 && setup_line[i][1] == "=") {
                args_f.push_back("--"+setup_line[i][0]);
                args_f.push_back(setup_line[i][2]);
            }
        }
        int argc_f = args_f.size();
        
        bool valid_setup_file_args = validate_setup_file_args(argc_f, args_f);
        if(!valid_setup_file_args) {
            exit(1);
        }

        // Process setup file arguments
        process_args(argc_f, args_f);
    }

    // Override setup file args with command line args
    process_args(argc, args);    

    // Extract initial condition file name without extension and paths
    string f1 = file_ic;
    string f2 = "";

    int N_str = f1.length();
    for(int k=0; k<N_str; k++) {
        f2 = f1.substr(N_str-1-k, 1+k);            
        if(f1.at(N_str-1-k) == '/') {
            f2 = f1.substr(N_str-1-k+1, 1+k-1);            
            break;
        }
    }

    string f3 = "";

    N_str = f2.length();
    for(int k=0; k<N_str; k++) {
        f3 = f2.substr(0, k+1);            
        if(f2.at(k) == '.') {
            f3 = f2.substr(0, k);            
            break;
        }
    }
    
    file_out = f3;

    // Set output directory and try creating it
    if(toSetOutputDir) {
        output_dir = file_out;    
    }

    bool dir_exist = check_dir_exists(output_dir);
    if(!dir_exist) {
        // try make dir
        int status = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if(status != 0) {
            cerr << "Cannot create output directory: " << output_dir << endl;
            exit(1);
        }
    }

    // Read in initial condition file
    vector< vector<string> > data = read_file(file_ic);
    
    vector<string> var, unit;
    int numLine = data.size();
    for(int i=0; i<numLine; i++) {
        string str = data[i][0];
        int N_str = str.length();

        string myvar, myunit;
        if(str.at(N_str-1) == ']') {
            for(int k=0; k<N_str; k++) {
                if(str.at(k) == '[') {
                    myvar = str.substr(0, k);
                    myunit = str.substr(k, N_str-k);             
                    break;
                }
            }
        }
        else {
            myvar = str;
            myunit = "";
        }
        
        var.push_back(myvar);
        unit.push_back(myunit);
    }

    // Determine unit system and number of bodies

    int num_mass = check_var(data, var, unit, numLine, "mass");
    if(num_mass == 0) {
        exit(1);
    } 

    N = num_mass;
    if(N == 1) {
        cerr << "Only a single mass value detected. TIDYMESS requires at least two bodies." << endl;
        exit(1);
    }

    // Determine physical units or N-body units
    determine_physical_units(var, unit);
    units.set_physical_units(physical_units);
    
    // Check consistency output units
    if(mass_unit == "" || mass_unit == "[]") {
        if(physical_units) {
            cerr << endl;
            cerr << "Incompatible units: mass_unit = " << mass_unit << endl;
            cerr << endl;
            exit(1);
        }
    }
    else {
        if(!physical_units) {
            cerr << endl;
            cerr << "Incompatible units: mass_unit = " << mass_unit << endl;
            cerr << endl;
            exit(1);
        }    
    }

    if(length_unit == "" || length_unit == "[]") {
        if(physical_units) {
            cerr << endl;
            cerr << "Incompatible units: length_unit = " << length_unit << endl;
            cerr << endl;
            exit(1);
        }
    }
    else {
        if(!physical_units) {
            cerr << endl;
            cerr << "Incompatible units: length_unit = " << length_unit << endl;
            cerr << endl;
            exit(1);
        }    
    }

    if(time_unit == "" || time_unit == "[]") {
        if(physical_units) {
            cerr << endl;
            cerr << "Incompatible units: time_unit = " << time_unit << endl;
            cerr << endl;
            exit(1);
        }
    }
    else {
        if(!physical_units) {
            cerr << endl;
            cerr << "Incompatible units: time_unit = " << time_unit << endl;
            cerr << endl;
            exit(1);
        }    
    }

    if(speed_unit == "" || speed_unit == "[]") {
        if(physical_units) {
            cerr << endl;
            cerr << "Incompatible units: speed_unit = " << speed_unit << endl;
            cerr << endl;
            exit(1);
        }
    }
    else {
        if(!physical_units) {
            cerr << endl;
            cerr << "Incompatible units: speed_unit = " << speed_unit << endl;
            cerr << endl;
            exit(1);
        }    
    }

    // Validate initial condition table
    
    if(orbit_coor == 0) { // Cartesian coordinates
        bool valid_cartesian = validate_initial_condition_cartesian(data, var, unit, numLine);
        if(!valid_cartesian) {
            exit(1);          
        }
    }
    else { // Elliptical coordinates
        bool valid_elliptical = validate_initial_condition_elliptical(data, var, unit, numLine);
        if(!valid_elliptical) {
            exit(1);          
        }
    }

    if(tidal_model >= 1) {
        int num_R = check_var(data, var, unit, numLine, "R");
        if(num_R == 0) {
            exit(1);
        }
        else if(num_R != N) {
            cerr << "Missing data in the initial condition table: R" << endl;
            exit(1);
        }
    
        int num_xi = check_var(data, var, unit, numLine, "xi");
        if(num_xi == 0) {
            exit(1);
        }
        else if(num_xi != N) {
            cerr << "Missing data in the initial condition table: xi" << endl; 
            exit(1);
        }

        int num_kf = check_var(data, var, unit, numLine, "kf");
        if(num_kf == 0) {
            exit(1);
        }
        else if(num_kf != N) {
            cerr << "Missing data in the initial condition table: kf" << endl; 
            exit(1);
        }
    }
    if(tidal_model >= 2) {
        int num_tau = check_var(data, var, unit, numLine, "tau");
        if(num_tau == 0) {
            exit(1);
        }
        else if(num_tau != N) {
            cerr << "Missing data in the initial condition table: tau" << endl; 
            exit(1);
        }
    }
    
    if(collisions > 0 || roche_limit > 0 || breakup_speed > 0) {
        int num_R = check_var(data, var, unit, numLine, "R");
        if(num_R == 0) {
            exit(1);
        }
        else if(num_R != N) {
            cerr << "Missing data in the initial condition table: R" << endl; 
            exit(1);
        }
    }

    if(magnetic_braking == 1) {
        int num_a_mb = check_var(data, var, unit, numLine, "a_mb");
        if(num_a_mb == 0) {
            exit(1);
        }
        else if(num_a_mb != N) {
            cerr << "Missing data in the initial condition table: a_mb" << endl; 
            exit(1);
        }
    }

    // Validate Post-Newtonian parameters
    if(pn_order > 0) {
        if(physical_units) {
            speed_of_light = 299792.458; // c = 299792.458 km/s
            speed_of_light = units.convert_speed_to_standard(speed_of_light, "[km/s]");
        }
        else {
            if(speed_of_light == 0) {
                speed_of_light = 1e100;
            }
        }
    }

    // Validate units
    bool valid_units = validate_units_internal_properties(var, unit); 
    if(!valid_units) {
        exit(1);
    }

    if(orbit_coor == 0) { // Cartesian coordinates
        valid_units = validate_units_orbital_cartesian(var, unit); 
        if(!valid_units) {
            exit(1);
        }
    }
    else { // Elliptical coordinates
        valid_units = validate_units_orbital_elliptical(var, unit); 
        if(!valid_units) {
            exit(1);
        }
    }    

    // If spin properties are specified in the initial conditions, validate them. 
    // If no spin properties are specified, TIDYMESS assumes no initial spin. 
    bool valid_spin = validate_initial_condition_spin(data, var, unit, numLine);
    if(!valid_spin) {
        exit(1);          
    }

    // Make sure simulation time is not zero
    double t_sim = t_end - t_begin;
    if(t_sim == 0) {
        cerr << "t_begin = " << t_begin << endl;
        cerr << "t_end   = " << t_end << endl;
        cerr << "Please set the end time to a different value than the begin time." << endl;
        exit(1);
    }

    // Set units of time integration parameters
    t_begin     = units.convert_time_to_standard(t_begin, time_unit);    
    t_end       = units.convert_time_to_standard(t_end, time_unit);    
    dt_const    = units.convert_time_to_standard(dt_const, time_unit);    
    
    // Determine first logarithmic step
    if(physical_units) {
        if(t_begin == 0) {
            if(t_end > 1) {
                t1_code = 1;
            }
            else {
                if(t_end * 365.25 > 1) {
                    t1_code = 1./365.25;
                }
                else if(t_end * 365.25 * 24 > 1) {
                    t1_code = 1./(365.25 * 24);
                }
                else if(t_end * 365.25 * 24 * 60 * 60 > 1) {
                    t1_code = 1./(365.25 * 24 * 60 * 60);
                }
                else {
                    t1_code = 1e-6 * t_end;
                }
            }
        }
    }
    else {
        if(t_begin == 0) {
            if(t_end > 1) {
                t1_code = 1;
            }
            else {
                t1_code = 1e-6 * t_end;
            }
        }
    }    
        
    // Process initial conditions and store in data_ic container
    data_ic.clear();
    for(int i=0; i<N; i++) {
        array<double, 15> dd = {};
        data_ic.push_back(dd);
    }

    // Upload internal properties in standard units
    for(int i=0; i<numLine; i++) {
        if(var[i] == "mass") {
            for(int j=0; j<N; j++) {
                data_ic[j][0] = units.convert_mass_to_standard( stod(data[i][1+j]), unit[i] );
            } 
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "R") {
            for(int j=0; j<N; j++) {
                data_ic[j][1] = units.convert_length_to_standard( stod(data[i][1+j]), unit[i] );
            } 
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "xi") {
            for(int j=0; j<N; j++) {
                data_ic[j][2] = stod(data[i][1+j]);
            } 
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "kf") {
            for(int j=0; j<N; j++) {
                data_ic[j][3] = stod(data[i][1+j]);
            } 
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "tau") {
            for(int j=0; j<N; j++) {
                data_ic[j][4] = units.convert_time_to_standard( stod(data[i][1+j]), unit[i] );
            } 
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "a_mb") {
            for(int j=0; j<N; j++) {
                data_ic[j][5] = units.convert_time_to_standard( stod(data[i][1+j]), unit[i] );
            } 
            break;
        }
    }

    // Upload orbital coordinates in standard units
    if(orbit_coor == 0) {
        for(int i=0; i<numLine; i++) {
            if(var[i] == "x") {
                for(int j=0; j<N; j++) {
                    data_ic[j][9] = units.convert_length_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "y") {
                for(int j=0; j<N; j++) {
                    data_ic[j][10] = units.convert_length_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "z") {
                for(int j=0; j<N; j++) {
                    data_ic[j][11] = units.convert_length_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "vx") {
                for(int j=0; j<N; j++) {
                    data_ic[j][12] = units.convert_speed_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "vy") {
                for(int j=0; j<N; j++) {
                    data_ic[j][13] = units.convert_speed_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "vz") {
                for(int j=0; j<N; j++) {
                    data_ic[j][14] = units.convert_speed_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
    }
    else {
        for(int i=0; i<numLine; i++) {
            if(var[i] == "a") {
                for(int j=0; j<N-1; j++) {
                    data_ic[j+1][9] = units.convert_length_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "e") {
                for(int j=0; j<N-1; j++) {
                    data_ic[j+1][10] = stod(data[i][1+j]);
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "I") {
                for(int j=0; j<N-1; j++) {
                    data_ic[j+1][11] = units.convert_angle_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "O") {
                for(int j=0; j<N-1; j++) {
                    data_ic[j+1][12] = units.convert_angle_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "w") {
                for(int j=0; j<N-1; j++) {
                    data_ic[j+1][13] = units.convert_angle_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
        for(int i=0; i<numLine; i++) {
            if(var[i] == "M") {
                for(int j=0; j<N-1; j++) {
                    data_ic[j+1][14] = units.convert_angle_to_standard( stod(data[i][1+j]), unit[i] );
                }  
                break;
            }
        }
    }
        
    // Upload spins in standard units
    for(int i=0; i<numLine; i++) {
        if(var[i] == "lod") {
            for(int j=0; j<N; j++) {
                data_ic[j][6] = units.convert_time_to_standard( stod(data[i][1+j]), unit[i] );
            }  
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "obl") {
            for(int j=0; j<N; j++) {
                data_ic[j][7] = units.convert_angle_to_standard( stod(data[i][1+j]), unit[i] );
            }
            break;
        }
    }
    for(int i=0; i<numLine; i++) {
        if(var[i] == "psi") {
            for(int j=0; j<N; j++) {
                data_ic[j][8] = units.convert_angle_to_standard( stod(data[i][1+j]), unit[i] );
            }
            break;
        }
    }

    // Convert orbital and spin coordinates to Cartesian inertial

    if(orbit_coor == 0) { // Cartesian inertial
        if(spin_coor == 0) { // absolute in inertial frame
            convert_spin_vectors_to_inertial();
        }
        else {
            cerr << " " << endl;
            cerr << "Incompatibility detected: orbit_coor = " << orbit_coor << ", spin_coor = " << spin_coor << endl;        
            cerr << " " << endl;
            exit(1);
        }
    }
    else if(orbit_coor == 1) { // Elliptical astrocentric
        if(spin_coor == 0) { // absolute in inertial frame
            convert_spin_vectors_to_inertial();                      
        }
        else if(spin_coor == 1) { // relative to orbit and body 0 absolute in inertial frame
            convert_spin_vectors_from_elliptical_body0abs();
        }
        else if(spin_coor == 2) { // relative to orbit and body 0 relative to inner orbit
            convert_spin_vectors_from_elliptical_body0rel();
        }    

        convert_astrocentric_elements_to_cartesian_coordinates();          
    }
    else if(orbit_coor == 2) { // Elliptical Jacobian
        if(spin_coor == 0) { // absolute in inertial frame
            convert_spin_vectors_to_inertial();                              
        }
        else if(spin_coor == 1) { // relative to orbit and body 0 absolute in inertial frame
            convert_spin_vectors_from_elliptical_body0abs();        
        }
        else if(spin_coor == 2) { // relative to orbit and body 0 relative to inner orbit
            convert_spin_vectors_from_elliptical_body0rel();        
        }    
        
        convert_jacobian_elements_to_cartesian_coordinates();                  
    }

    // Convert units to G=1 system
    if(physical_units) {    
        t_begin     = units.convert_time_from_standard_to_code(t_begin);    
        t_end       = units.convert_time_from_standard_to_code(t_end);    

        dt_const    = units.convert_time_from_standard_to_code(dt_const);   
        
        speed_of_light = units.convert_speed_from_standard_to_code(speed_of_light);
        
        for(int i=0; i<N; i++) {
            data_ic[i][4] = units.convert_time_from_standard_to_code(data_ic[i][4]);
            data_ic[i][5] = units.convert_time_from_standard_to_code(data_ic[i][5]);

            data_ic[i][6] = units.convert_frequency_from_standard_to_code(data_ic[i][6]);
            data_ic[i][7] = units.convert_frequency_from_standard_to_code(data_ic[i][7]);
            data_ic[i][8] = units.convert_frequency_from_standard_to_code(data_ic[i][8]);

            data_ic[i][12] = units.convert_speed_from_standard_to_code(data_ic[i][12]);
            data_ic[i][13] = units.convert_speed_from_standard_to_code(data_ic[i][13]);
            data_ic[i][14] = units.convert_speed_from_standard_to_code(data_ic[i][14]);
        }        

        // Logarithmic first step
        t1_code = units.convert_time_from_standard_to_code(t1_code);   
    }    

    // Convert units to dimensionless units and store scale factors
    if(physical_units) {
        units.Cm = 0;
        for(int i=0; i<N; i++) {
            units.Cm += data_ic[i][0];
        }

        if(units.Cm == 0) {
            cerr << "The cumulative mass of the N-body system is zero! Please include a non-zero mass body." << endl;
            exit(1);
        }

        units.Cr = 0;
        int cnt = 0;
        for(int i=0; i<N-1; i++) {
            for(int j=i+1; j<N; j++) {
                double dx = data_ic[i][9]-data_ic[j][9];
                double dy = data_ic[i][10]-data_ic[j][10];
                double dz = data_ic[i][11]-data_ic[j][11];
                double dr2 = dx*dx + dy*dy + dz*dz;
                double dr1 = sqrt(dr2);
                units.Cr += dr1;
                cnt++;
            }
        }
        units.Cr /= cnt;

        units.Cv = sqrt(units.Cm/units.Cr);
        units.Ct = units.Cr/units.Cv;
    }

    t_begin /= units.Ct;
    t_end /= units.Ct;

    dt_const /= units.Ct;

    speed_of_light /= units.Cv;

    // Logarithmic first step
    t1_code /= units.Ct;

    for(int i=0; i<N; i++) {
        data_ic[i][0] /= units.Cm;
        data_ic[i][1] /= units.Cr;
        data_ic[i][4] /= units.Ct;
        data_ic[i][5] /= units.Ct;
        data_ic[i][6] *= units.Ct;
        data_ic[i][7] *= units.Ct;
        data_ic[i][8] *= units.Ct;
        data_ic[i][9] /= units.Cr;
        data_ic[i][10] /= units.Cr;
        data_ic[i][11] /= units.Cr;
        data_ic[i][12] /= units.Cv;
        data_ic[i][13] /= units.Cv;
        data_ic[i][14] /= units.Cv;    
    }

    // Name and id tags        
    for(int i=0; i<N; i++) {
        id.push_back(i);
        string str = "body" + to_string(i);
        name.push_back(str);
    }

    for(int i=0; i<numLine; i++) {
        if(var[i] == "name" || var[i] == "Name") {
            for(int j=0; j<data[i].size()-1; j++) {
                name[j] = data[i][1+j];  
            }  
            break;
        }
    }    

    // Only include 1+num_body bodies
    if(num_body > 0) {
        for(int i=0; i<N; i++) {
            if(N-1-i >= 1+num_body) {
                data_ic.erase (data_ic.begin()+N-1-i);
            }
        }
        N = data_ic.size();

        // move to center
        move_to_center();
    }    
}
bool Initializer::validate_command_line_args(int argc, vector<string> &args) {
    vector<bool> arg_validated(argc, false);

    for(int i=0; i<argc; i++) {
        string str = args[i];
        if(str.length() > 1) {
            if(str.at(0) == '-' && str.at(1) == '-') {
                str.erase(str.begin()+0, str.begin()+2);
                bool isValid = check_valid_argument(str);
                if(!isValid) {
                    cerr << " " << endl;
                    cerr << "Unrecognized command line argument: " << args[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
                arg_validated[i] = true;
                                
                if(i+1 < argc) {
                    string str_val = args[i+1];
                    isValid = check_valid_value(str, str_val);
                    if(!isValid) {
                        cerr << " " << endl;
                        cerr << "Invalid command line value: " << args[i] << " " << args[i+1] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                    arg_validated[i+1] = true;
                }
                else {                
                    if(str == "output_dir") {
                        ;
                    }
                    else {
                        cerr << " " << endl;
                        cerr << "Please specify a value after: " << args[i] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                }
            }
            else if(str == "-f") {
                str.erase(str.begin()+0, str.begin()+1);
                bool isValid = check_valid_argument(str);
                if(!isValid) {
                    cerr << " " << endl;
                    cerr << "Unrecognized command line argument: " << args[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
                arg_validated[i] = true;
                
                if(i+1 < argc) {
                    string str_val = args[i+1];
                    isValid = check_valid_value(str, str_val);
                    if(!isValid) {
                        cerr << " " << endl;
                        cerr << "Invalid command line value: " << args[i] << " " << args[i+1] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                    arg_validated[i+1] = true;
                }
                else {
                    cerr << " " << endl;
                    cerr << "Please specify a value after: " << args[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }
        }
    }

    for(int i=0; i<argc; i++) {
        if(!arg_validated[i]) {
            cerr << " " << endl;
            cerr << "Unrecognized command line entry: " << args[i] << endl;
            cerr << " " << endl;
            return false;
        }
    }

    return true;
}
bool Initializer::validate_setup_file_args(int argc, vector<string> &args) {
    vector<bool> arg_validated(argc, false);

    bool ignoring = false;

    for(int i=0; i<argc; i++) {
        string str = args[i];
        if(str.length() > 1) {
            if(str.at(0) == '-' && str.at(1) == '-') {
                str.erase(str.begin()+0, str.begin()+2);
                bool isValid = check_valid_argument(str);
                if(!isValid) {
                    cerr << " " << endl;
                    cerr << "Warning: ignoring line in parameter setup file starting with: " << str << endl;
                    ignoring = true;
                }
                else {
                    if(i+1 < argc) {
                        string str_val = args[i+1];
                        isValid = check_valid_value(str, str_val);
                        if(!isValid) {
                            cerr << " " << endl;
                            cerr << "Invalid value in parameter setup file: " << args[i] << " " << args[i+1] << endl;
                            cerr << " " << endl;
                            return false;
                        }
                        arg_validated[i+1] = true;
                    }
                    else {
                        if(str == "output_dir") {
                            ;
                        }
                        else {
                            cerr << " " << endl;
                            cerr << "Please specify a value after: " << args[i] << endl;
                            cerr << " " << endl;
                            return false;
                        }
                    }
                
                }
                arg_validated[i] = true;
            }
        }
    }

    if(ignoring) cerr << " " << endl;

    return true;
}
bool Initializer::check_valid_argument(string &str) {
    if(str == "to_continue") {
        return true;
    }
    else if(str == "max_cpu_time") {
        return true;
    }
    else if(str == "file_par" || str == "f") {
        return true;
    }
    else if(str == "file_ic") {
        return true;
    }
    else if(str == "orbit_coor") {
        return true;
    }
    else if(str == "num_body") {
        return true;
    }
    else if(str == "output_format") {
        return true;
    }
    else if(str == "output_info") {
        return true;
    }
    else if(str == "output_coor") {
        return true;
    }
    else if(str == "output_dir") {
        return true;
    }
    else if(str == "overwrite") {
        return true;
    }
    else if(str == "output_diag") {
        return true;
    }
    else if(str == "output_terminal") {
        return true;
    }    
    else if(str == "mass_unit") {
        return true;
    }
    else if(str == "length_unit") {
        return true;
    }
    else if(str == "time_unit") {
        return true;
    }
    else if(str == "speed_unit") {
        return true;
    }
    else if(str == "t_begin") {
        return true;
    }
    else if(str == "t_end") {
        return true;
    }
    else if(str == "dt_mode") {
        return true;
    }
    else if(str == "dt_const") {
        return true;
    }
    else if(str == "eta") {
        return true;
    }
    else if(str == "snapshot_mode") {
        return true;
    }
    else if(str == "n_snapshot") {
        return true;
    }
    else if(str == "tidal_model") {
        return true;
    }
    else if(str == "initial_shape") {
        return true;
    }
    else if(str == "spin_coor") {
        return true;
    }
    else if(str == "pn_order") {
        return true;
    }
    else if(str == "speed_of_light") {
        return true;
    }
    else if(str == "B_braking") {
        return true;
    }
    else if(str == "collisions") {
        return true;
    }
    else if(str == "roche_limit") {
        return true;
    }
    else if(str == "breakup_speed") {
        return true;
    }
    else if(str == "n_iter") {
        return true;
    }
    return false;
}
bool Initializer::check_valid_value(string &str, string &str_val) {
    if(str == "to_continue") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "max_cpu_time") {
        try {
            size_t offset = 0;
            double t = stod(str_val, &offset);
            if(offset == str_val.length()) {
                if(t >= 0) return true;
                else {
                    cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "file_par" || str == "f") {
        //string f = str_val;
        //if(f == "none" || f == "None") return true;
        //else {
        //    bool f_exist = check_file_exists(f);
        //    if(f_exist) return true;
        //    else {
                //cerr << "Cannot open " << f << "!" << endl;
                //return false;                      
        //    }
        //}
        return true;
    }
    else if(str == "file_ic") {
        string f = str_val;
        //bool f_exist = check_file_exists(f);
        //if(f_exist) return true;
        //else {
            //cerr << "Cannot open " << f << "!" << endl;
            //return false;                      
        //}
        return true;
    }
    else if(str == "orbit_coor") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "num_body") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                if(num >= 0) return true;
                else {
                    cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "output_format") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "output_info") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "output_coor") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "output_dir") {
        return true;
    }
    else if(str == "overwrite") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "output_diag") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "output_terminal") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "mass_unit") {
        string u = str_val;
        bool valid_u = units.validate_mass_unit(u);
        if(valid_u) return true;
        else {
            cerr << "Please specify a valid value after: " << str << endl;        
        }
    }
    else if(str == "length_unit") {
        string u = str_val;
        bool valid_u = units.validate_length_unit(u);
        if(valid_u) return true;
        else {
            cerr << "Please specify a valid value after: " << str << endl;        
        }
    }
    else if(str == "time_unit") {
        string u = str_val;
        bool valid_u = units.validate_time_unit(u);
        if(valid_u) return true;
        else {
            cerr << "Please specify a valid value after: " << str << endl;        
        }
    }
    else if(str == "speed_unit") {
        string u = str_val;
        bool valid_u = units.validate_speed_unit(u);
        if(valid_u) return true;
        else {
            cerr << "Please specify a valid value after: " << str << endl;        
        }
    }
    else if(str == "t_begin") {
        try {
            size_t offset = 0;
            double t = stod(str_val, &offset);
            if(offset == str_val.length()) return true;
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "t_end") {
        try {
            size_t offset = 0;
            double t = stod(str_val, &offset);
            if(offset == str_val.length()) return true;
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "dt_mode") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }    
    else if(str == "dt_const") {
        try {
            size_t offset = 0;
            double t = stod(str_val, &offset);
            if(offset == str_val.length()) {
                if(t > 0) return true;
                else {
                    if(dt_mode > 0) return true;
                    else cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "eta") {
        try {
            size_t offset = 0;
            double t = stod(str_val, &offset);
            if(offset == str_val.length()) {
                if(t > 0) return true;
                else {
                    cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "snapshot_mode") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "n_snapshot") {
        try {
            size_t offset = 0;

            double numd = stod(str_val, &offset);
            int num = (int)numd;
                        
            if(offset == str_val.length()) {
                if(num > 0) return true;
                else {
                    cerr << "Please specify a valid value after: " << str << endl;
                }           
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "tidal_model") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    case 3: 
                        return true;
                    case 4: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "initial_shape") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "spin_coor") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "pn_order") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    case 25: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "speed_of_light") {
        try {
            size_t offset = 0;
            double t = stod(str_val, &offset);
            if(offset == str_val.length()) {
                if(t >= 0) return true;
                else {
                    cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }

    }
    else if(str == "B_braking") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "collisions") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    case 3: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "roche_limit") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "breakup_speed") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                switch(num) {
                    case 0:
                        return true;
                    case 1: 
                        return true;
                    case 2: 
                        return true;
                    default:
                        cerr << "Please specify a valid value after: " << str << endl;
                }
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }
    else if(str == "n_iter") {
        try {
            size_t offset = 0;
            int num = stoi(str_val, &offset);
            if(offset == str_val.length()) {
                if(num >= 0) return true;
                else {
                    cerr << "Please specify a valid value after: " << str << endl;
                }           
            }
            else {
                cerr << "Please specify a valid value after: " << str << endl;
            }
        }
        catch(exception &e) {
            cerr << "Please specify a valid value after: " << str << endl;
        }
    }

    cerr << "Type \"./tidymess.exe -p\" or \"./tidymess.exe --parameters\" for an overview of input parameters and values" << endl;
    return false;
}
void Initializer::process_args(int argc, vector<string> &args) {
    for(int i=0; i<argc; i++) {
        string str = args[i];
        if(str.length() > 1) {
            if(str.at(0) == '-' && str.at(1) == '-') {
                str.erase(str.begin()+0, str.begin()+2);
                if(str == "to_continue") {
                    to_continue = stoi(args[i+1]);
                }
                else if(str == "max_cpu_time") {
                    max_cpu_time = stod(args[i+1]);
                }
                else if(str == "file_ic") {
                    file_ic = args[i+1];
                }
                else if(str == "orbit_coor") {
                    orbit_coor = stoi(args[i+1]);
                }
                else if(str == "num_body") {
                    num_body = stoi(args[i+1]);
                }
                else if(str == "output_format") {
                    output_format = stoi(args[i+1]);
                }
                else if(str == "output_info") {
                    output_info = stoi(args[i+1]);
                }
                else if(str == "output_coor") {
                    output_coor = stoi(args[i+1]);
                }
                else if(str == "output_dir") {
                    if(i+1 < argc) {
                        string dir = args[i+1];
        
                        if(dir == "/" || dir == "//") {
                            output_dir = "";
                            toSetOutputDir = true;
                        }
                        else if(dir.at(0) == '-') {
                            output_dir = "";
                            toSetOutputDir = true;                        
                        }
                        else {
                            output_dir = dir;
                            toSetOutputDir = false;
                        }
                    }
                    else {
                        output_dir = "";
                        toSetOutputDir = true;
                    }
                }
                else if(str == "overwrite") {
                    overwrite = stoi(args[i+1]);
                }
                else if(str == "output_diag") {
                    output_diag = stoi(args[i+1]);
                }
                else if(str == "output_terminal") {
                    output_terminal = stoi(args[i+1]);
                }
                else if(str == "mass_unit") {
                    mass_unit = args[i+1];
                }
                else if(str == "length_unit") {
                    length_unit = args[i+1];
                }
                else if(str == "time_unit") {
                    time_unit = args[i+1];
                }
                else if(str == "speed_unit") {
                    speed_unit = args[i+1];
                }
                else if(str == "t_begin") {
                    t_begin = stod(args[i+1]);
                }
                else if(str == "t_end") {
                    t_end = stod(args[i+1]);
                }
                else if(str == "dt_mode") {
                    dt_mode = stoi(args[i+1]); 
                }
                else if(str == "dt_const") {
                    dt_const = stod(args[i+1]);
                }
                else if(str == "eta") {
                    eta = stod(args[i+1]);
                }
                else if(str == "snapshot_mode") {
                    snapshot_mode = stoi(args[i+1]); 
                }
                else if(str == "n_snapshot") {
                    double numd = stod(args[i+1]);
                    int num = (int)numd;
                    n_snapshot = num; 
                }
                else if(str == "tidal_model") {
                    tidal_model = stoi(args[i+1]); 
                }
                else if(str == "initial_shape") {
                    initial_shape = stoi(args[i+1]);
                }
                else if(str == "spin_coor") {
                    spin_coor = stoi(args[i+1]);
                }
                else if(str == "pn_order") {
                    pn_order = stoi(args[i+1]);
                }
                else if(str == "speed_of_light") {
                    speed_of_light = stod(args[i+1]);
                }
                else if(str == "B_braking") {
                    magnetic_braking = stoi(args[i+1]);
                }
                else if(str == "collisions") {
                    collisions = stoi(args[i+1]);
                }
                else if(str == "roche_limit") {
                    roche_limit = stoi(args[i+1]);
                }
                else if(str == "breakup_speed") {
                    breakup_speed = stoi(args[i+1]);
                }
                else if(str == "n_iter") {
                    n_iter = stoi(args[i+1]);
                }
            }
        }
    }
}

void Initializer::print_arguments_as_setup_file() {
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "// Simulation parameters" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "to_continue\t= " << to_continue << "\t\t// 0=new simulation, 1=continue simulation" << endl;
    cout << "max_cpu_time\t= " << max_cpu_time << "\t\t// Maximum CPU running time in seconds. If 0 (default), then cpu time has no limit." << endl;
    cout << "" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "// Physical model parameters" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "tidal_model\t= " << tidal_model << "\t\t// 0=none, 1=conservative, 2=linear, 3=creep direct, 4=creep tidymess (default)" << endl;
    cout << "pn_order\t= " << pn_order << "\t\t// Post-Newtonian order: 0=none, 1=1pn, 2=1+2pn, 25=1+2+2.5pn" << endl;
    cout << "B_braking\t= " << magnetic_braking << "\t\t// Magnetic braking. 0=off, 1=on" << endl;
    cout << "" << endl;
    cout << "collisions\t= " << collisions << "\t\t// 0=off, 1=flag, 2=exception, 3=replace" << endl;
    cout << "roche_limit\t= " << roche_limit << "\t\t// 0=off, 1=flag, 2=exception" << endl;
    cout << "breakup_speed\t= " << breakup_speed << "\t\t// Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception" << endl;
    cout << "" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "// Unit system" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "mass_unit\t= " << mass_unit << "\t// Unit of mass in output: []=Nbody unit, [g], [kg], [Mearth], [Mjupiter], [Msun]" << endl;
    cout << "length_unit\t= " << length_unit << "\t\t// Unit of length in output: []=Nbody unit, [m], [km], [Rsun], [au], [pc]" << endl;
    cout << "time_unit\t= " << time_unit << "\t\t// Unit of time for 1) t_begin, t_end and 2) unit of time in output: []=Nbody unit, [s], [hr], [day], [yr], [Myr], [Gyr]" << endl;  
    cout << "speed_unit\t= " << speed_unit << "\t// Unit of speed in output: []=Nbody unit, [m/s], [km/s], [km/hr], [km/hour], [au/day], [au/yr], [au/year], [pc/Myr]" << endl;
    cout << "" << endl;
    cout << "speed_of_light\t= " << speed_of_light << "\t\t// Speed of light in N-body units. Only used in conjunction with N-body units and pn_order>0, otherwise equal to c." << endl;
    cout << "" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "// Initial condition parameters" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "file_ic\t	= " << file_ic << "\t// initial condition file" << endl;
    cout << "" << endl;
    cout << "orbit_coor\t= " << orbit_coor << "\t\t// 0=cartesian inertial, 1=elliptical astrocentric, 2=elliptical jacobian" << endl;
    cout << "spin_coor\t= " << spin_coor << "\t\t// 0=absolute in inertial frame, 1=relative to its orbit; body 0 in the inertial frame, 2=relative to its orbit; body 0 relative to innermost orbit" << endl;
    cout << "" << endl;
    cout << "initial_shape\t= " << initial_shape << "\t\t// 0=sphere, 1=equilibrium" << endl;
    cout << "num_body\t= " << num_body << "\t\t// 0=all, num_body+1=number of bodies to include" << endl;
    cout << "" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "// Output parameters" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "snapshot_mode\t = " <<  snapshot_mode << "\t\t// 0=linear interval (default), 1=logarithmic interval, 2=variable interval" << endl; 
    cout << "n_snapshot\t = " << n_snapshot << "\t\t// Total number of snapshots between t_begin and t_end (linear or in log10), or output a snapshot every fixed number (n_snapshot) of integration steps (variable)" << endl; 
    cout << "" << endl;    
    cout << "output_dir\t= " << output_dir << "\t\t// Output directory; default is 'data/'. If left blank or set to '/', then file_ic will be adopted without the extension." << endl;  
    cout << "overwrite\t= " << overwrite << "\t\t// overwrite existing files: 0=no, 1=yes" << endl;
    cout << "" << endl;
    cout << "output_format\t= " << output_format << "\t\t// 0=file per body, 1=file per snapshot, 2=single file" << endl;
    cout << "output_info\t= " << output_info << "\t\t// 0=time-varying quantities, 1=all quantities" << endl;
    cout << "output_coor\t= " << output_coor << "\t\t// 0=cartesian inertial" << endl;
    cout << "" << endl;
    cout << "output_diag\t= " << output_diag << "\t\t// 0=no (default), 1=yes: output diagnostics, such as E and L, are written to a separate diagnostics file with extension '.diag'" << endl;
    cout << "output_terminal\t= " << output_terminal << "\t\t// Display progress of simulation in terminal window. 0=no, 1=yes" << endl;
    cout << "" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "// Integration parameters" << endl;
    cout << "//--------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    cout << "t_begin\t	= " << t_begin << "\t\t// begin time in units given by time_unit" << endl;
    cout << "t_end\t\t= " << t_end << "\t\t// final time in units given by time_unit" << endl;
    cout << "" << endl;
    cout << "dt_mode\t\t= " << dt_mode << "\t\t// 0=constant dt, 1=adaptive dt, 2=adaptive, weighted dt" << endl;
    cout << "dt_const\t= " << dt_const << "\t// constant time step in units given by time_unit (only used if dt_mode=0)" << endl;
    cout << "eta\t\t= " << eta << "\t\t// accuracy parameter; timestep multiplication factor, default=0.125 (only used if dt_mode>0)" << endl;
    cout << "" << endl;
    cout << "n_iter\t\t= " << n_iter << "\t\t// Number of iterations to improve reversibility (default=1)" << endl;
}
void Initializer::print_arguments_as_setup_file(ofstream &fo) {
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "// Simulation parameters" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "" << endl;
    fo << "to_continue\t= " << to_continue << "\t\t// 0=new simulation, 1=continue simulation" << endl;
    fo << "max_cpu_time\t= " << max_cpu_time << "\t\t// Maximum CPU running time in seconds. If 0 (default), then cpu time has no limit." << endl;
    fo << "" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "// Physical model parameters" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "" << endl;
    fo << "tidal_model\t= " << tidal_model << "\t\t// 0=none, 1=conservative, 2=linear, 3=creep direct, 4=creep tidymess (default)" << endl;
    fo << "pn_order\t= " << pn_order << "\t\t// Post-Newtonian order: 0=none, 1=1pn, 2=1+2pn, 25=1+2+2.5pn" << endl;
    fo << "B_braking\t= " << magnetic_braking << "\t\t// Magnetic braking. 0=off, 1=on" << endl;
    fo << "" << endl;
    fo << "collisions\t= " << collisions << "\t\t// 0=off, 1=flag, 2=exception, 3=replace" << endl;
    fo << "roche_limit\t= " << roche_limit << "\t\t// 0=off, 1=flag, 2=exception" << endl;
    fo << "breakup_speed\t= " << breakup_speed << "\t\t// Centrifugal breakup speed detection. 0=off, 1=flag, 2=exception" << endl;
    fo << "" << endl;    
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "// Unit system" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "" << endl;
    fo << "mass_unit\t= " << mass_unit << "\t// Unit of mass in output: []=Nbody unit, [g], [kg], [Mearth], [Mjupiter], [Msun]" << endl;
    fo << "length_unit\t= " << length_unit << "\t\t// Unit of length in output: []=Nbody unit, [m], [km], [Rsun], [au], [pc]" << endl;
    fo << "time_unit\t= " << time_unit << "\t\t// Unit of time for 1) t_begin, t_end and 2) unit of time in output: []=Nbody unit, [s], [hr], [day], [yr], [Myr], [Gyr]" << endl;  
    fo << "speed_unit\t= " << speed_unit << "\t// Unit of speed in output: []=Nbody unit, [m/s], [km/s], [km/hr], [km/hour], [au/day], [au/yr], [au/year], [pc/Myr]" << endl;
    fo << "" << endl;
    fo << "speed_of_light\t= " << speed_of_light << "\t\t// Speed of light in N-body units. Only used in conjunction with N-body units and pn_order>0, otherwise equal to c." << endl;
    fo << "" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "// Initial condition parameters" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "" << endl;
    fo << "file_ic\t	= " << file_ic << "\t// initial condition file" << endl;
    fo << "" << endl;
    fo << "orbit_coor\t= " << orbit_coor << "\t\t// 0=cartesian inertial, 1=elliptical astrocentric, 2=elliptical jacobian" << endl;
    fo << "spin_coor\t= " << spin_coor << "\t\t// 0=absolute in inertial frame, 1=relative to its orbit; body 0 in the inertial frame, 2=relative to its orbit; body 0 relative to innermost orbit" << endl;
    fo << "" << endl;
    fo << "initial_shape\t= " << initial_shape << "\t\t// 0=sphere, 1=equilibrium" << endl;
    fo << "num_body\t= " << num_body << "\t\t// 0=all, num_body+1=number of bodies to include" << endl;
    fo << "" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "// Output parameters" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "" << endl;
    fo << "snapshot_mode\t = " <<  snapshot_mode << "\t\t// 0=linear interval (default), 1=logarithmic interval, 2=variable interval" << endl; 
    fo << "n_snapshot\t = " << n_snapshot << "\t\t// Total number of snapshots between t_begin and t_end (linear or in log10), or output a snapshot every fixed number (n_snapshot) of integration steps (variable)" << endl; 
    fo << "" << endl;        
    fo << "output_dir\t= " << output_dir << "\t\t// Output directory; default is 'data/'. If left blank or set to '/', then file_ic will be adopted without the extension." << endl;  
    fo << "overwrite\t= " << overwrite << "\t\t// overwrite existing files: 0=no, 1=yes" << endl;
    fo << "" << endl;
    fo << "output_format\t= " << output_format << "\t\t// 0=file per body, 1=file per snapshot, 2=single file" << endl;
    fo << "output_info\t= " << output_info << "\t\t// 0=time-varying quantities, 1=all quantities" << endl;
    fo << "output_coor\t= " << output_coor << "\t\t// 0=cartesian inertial" << endl;
    fo << "" << endl;
    fo << "output_diag\t= " << output_diag << "\t\t// 0=no (default), 1=yes: output diagnostics, such as E and L, are written to a separate diagnostics file with extension '.diag'" << endl;
    fo << "output_terminal\t= " << output_terminal << "\t\t// Display progress of simulation in terminal window. 0=no, 1=yes" << endl;
    fo << "" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "// Integration parameters" << endl;
    fo << "//--------------------------------------------------------------------------------------------------" << endl;
    fo << "" << endl;
    fo << "t_begin\t	= " << t_begin << "\t\t// begin time in units given by time_unit" << endl;
    fo << "t_end\t\t= " << t_end << "\t\t// final time in units given by time_unit" << endl;
    fo << "" << endl;
    fo << "dt_mode\t\t= " << dt_mode << "\t\t// 0=constant dt, 1=adaptive dt, 2=adaptive, weighted dt" << endl;
    fo << "dt_const\t= " << dt_const << "\t// constant time step in units given by time_unit (only used if dt_mode=0)" << endl;
    fo << "eta\t\t= " << eta << "\t\t// accuracy parameter; timestep multiplication factor, default=0.125 (only used if dt_mode>0)" << endl;
    fo << "" << endl;
    fo << "n_iter\t\t= " << n_iter << "\t\t// Number of iterations to improve reversibility (default=1)" << endl;
}
void Initializer::print_initial_condition_args() {
    cout << "Overview of properties to be assigned to the bodies in the initial condition file:" << endl;
    cout << " " << endl;
    cout << "- Internal properties" << endl;
    cout << "mass = mass" << endl;
    cout << "R    = radius" << endl;
    cout << "xi   = moment of inertia factor" << endl;
    cout << " " << endl;
    cout << "- Tidal response parameters" << endl;
    cout << "kf   = Fluid Love number for potential" << endl;
    cout << "tau  = Fluid relaxation time" << endl;
    cout << " " << endl;
    cout << "- Other physics" << endl;
    cout << "a_mb = magnetic braking coefficient" << endl;
    cout << " " << endl;
    cout << "- Spin coordinates" << endl;
    cout << "lod  = length of day" << endl;
    cout << "obl  = obliquity" << endl;
    cout << "psi  = precession angle" << endl;
    cout << " " << endl;
    cout << "- Orbital coordinates: Cartesian " << endl;
    cout << "x    = position along x-axis" << endl;
    cout << "y    = position along y-axis " << endl;
    cout << "z    = position along z-axis " << endl;
    cout << "vx   = velocity along x-axis " << endl;
    cout << "vy   = velocity along y-axis " << endl;
    cout << "vz   = velocity along z-axis " << endl;
    cout << " " << endl;
    cout << "- Orbital coordinates: Elliptical " << endl;
    cout << "a    = semimajor axis" << endl;
    cout << "e    = eccentricity" << endl;
    cout << "I    = inclination" << endl;
    cout << "O    = longitude of ascending node" << endl;
    cout << "w    = argument of pericenter" << endl;
    cout << "M    = mean anomaly" << endl;
}

int Initializer::check_var(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine, string x) {
    int num_x = 0;
    bool found_x = false;
    for(int i=0; i<numLine; i++) {
        if(var[i] == x) {
            found_x = true;
            num_x = data[i].size()-1;
            break;
        }
    }
    if(!found_x || num_x == 0) {
        cerr << "Incomplete initial conditions: no " << x << " values detected." << endl;
    }
    return num_x;
}
bool Initializer::validate_initial_condition_cartesian(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine) {
    int num_x  = check_var(data, var, unit, numLine, "x");
    if(num_x == 0) {
        return false;
    }
    else if(num_x != N) {
        cerr << "Inconsistency in the initial condition table: x" << endl; 
        return false;
    }

    int num_y  = check_var(data, var, unit, numLine, "y");
    if(num_y == 0) {
        return false;
    }
    else if(num_y != N) {
        cerr << "Inconsistency in the initial condition table: y" << endl; 
        return false;
    }

    int num_z  = check_var(data, var, unit, numLine, "z");
    if(num_z == 0) {
        return false;
    }
    else if(num_z != N) {
        cerr << "Inconsistency in the initial condition table: z" << endl; 
        return false;
    }

    int num_vx = check_var(data, var, unit, numLine, "vx");
    if(num_vx == 0) {
        return false;
    }
    else if(num_vx != N) {
        cerr << "Inconsistency in the initial condition table: vx" << endl; 
        return false;
    }

    int num_vy = check_var(data, var, unit, numLine, "vy");
    if(num_vy == 0) {
        return false;
    }
    else if(num_vy != N) {
        cerr << "Inconsistency in the initial condition table: vy" << endl; 
        return false;
    }

    int num_vz = check_var(data, var, unit, numLine, "vz");
    if(num_vz == 0) {
        return false;
    }
    else if(num_vz != N) {
        cerr << "Inconsistency in the initial condition table: vz" << endl; 
        return false;
    }

    return true;
}
bool Initializer::validate_initial_condition_elliptical(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine) {
    int num_a  = check_var(data, var, unit, numLine, "a");
    if(num_a == 0) {
        return false;
    }
    else if(num_a != N-1) {
        cerr << "Inconsistency in the initial condition table: a" << endl; 
        return false;
    }

    int num_e  = check_var(data, var, unit, numLine, "e");
    if(num_e == 0) {
        return false;
    }
    else if(num_e != N-1) {
        cerr << "Inconsistency in the initial condition table: e" << endl; 
        return false;
    }

    int num_I  = check_var(data, var, unit, numLine, "I");
    if(num_I == 0) {
        return false;
    }
    else if(num_I != N-1) {
        cerr << "Inconsistency in the initial condition table: I" << endl; 
        return false;
    }

    int num_O = check_var(data, var, unit, numLine, "O");
    if(num_O == 0) {
        return false;
    }
    else if(num_O != N-1) {
        cerr << "Inconsistency in the initial condition table: O" << endl; 
        return false;
    }

    int num_w = check_var(data, var, unit, numLine, "w");
    if(num_w == 0) {
        return false;
    }
    else if(num_w != N-1) {
        cerr << "Inconsistency in the initial condition table: w" << endl; 
        return false;
    }

    int num_M = check_var(data, var, unit, numLine, "M");
    if(num_M == 0) {
        return false;
    }
    else if(num_M != N-1) {
        cerr << "Inconsistency in the initial condition table: M" << endl; 
        return false;
    }

    return true;
}
bool Initializer::validate_initial_condition_spin(vector< vector<string> > &data, vector<string> &var, vector<string> &unit, int numLine) {
    int num_lod = 0;
    bool found_lod = false;
    for(int i=0; i<numLine; i++) {
        if(var[i] == "lod") {
            found_lod = true;
            num_lod = data[i].size()-1;
            break;
        }
    }

    int num_obl = 0;
    bool found_obl = false;
    for(int i=0; i<numLine; i++) {
        if(var[i] == "obl") {
            found_obl = true;
            num_obl = data[i].size()-1;
            break;
        }
    }

    int num_psi = 0;
    bool found_psi = false;
    for(int i=0; i<numLine; i++) {
        if(var[i] == "psi") {
            found_psi = true;
            num_psi = data[i].size()-1;
            break;
        }
    }

    if(found_lod) {
        if(found_obl) {
            if(found_psi) {
                if(num_lod != N) {
                    cerr << "Missing data in the initial condition table: lod" << endl; 
                    return false;
                }
                else if(num_obl != N) {
                    cerr << "Missing data in the initial condition table: obl" << endl; 
                    return false;
                }
                else if(num_psi != N) {
                    cerr << "Missing data in the initial condition table: psi" << endl; 
                    return false;
                }
                else {
                    bool valid_units = validate_units_spin(var, unit); 
                    if(!valid_units) {
                        return false;
                    }
                }
            }
            else {
                cerr << "Spin period \"lod\" and angle \"obl\" detected, but no \"psi\" specified" << endl;
                return false;
            }
        }
        else {
            if(found_psi) {
                cerr << "Spin period \"lod\" and angle \"psi\" detected, but no \"obl\" specified" << endl;
                return false;
            }
            else {
                cerr << "Spin period \"lod\" detected, but no \"obl\" and \"psi\" specified" << endl;
                return false;
            }
        }
    }
    else {
        if(found_obl) {
            if(found_psi) {
                cerr << "Spin angles \"obl\" and \"psi\" detected, but no \"lod\" specified" << endl;
                return false;
            }
            else {
                cerr << "Spin angle \"obl\" detected, but no \"lod\" and \"psi\" specified" << endl;
                return false;
            }
        }
        else {
            if(found_psi) {
                cerr << "Spin angle \"psi\" detected, but no \"lod\" and \"obl\" specified" << endl;
                return false;
            }
        }
    }

    return true;
}

void Initializer::determine_physical_units(vector<string> &var, vector<string> &unit) { 
    int numLine = var.size();
    
    for(int i=0; i<numLine; i++) {
        if(var[i] == "mass") {
            bool valid_unit = units.validate_mass_unit(unit[i]);
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                exit(1);
            }

            if(unit[i] == "" || unit[i] == "[]") physical_units = false;
            else physical_units = true;

            break;        
        }
    }
}
    
bool Initializer::validate_units_internal_properties(vector<string> &var, vector<string> &unit) { 
    int numLine = var.size();
    
    for(int i=0; i<numLine; i++) {
        if(var[i] == "mass") {
            bool valid_unit = units.validate_mass_unit(unit[i]);
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                return false;
            }

            if(unit[i] == "" || unit[i] == "[]") {
                if(physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }    
            else {
                if(!physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }

            break;        
        }
    }

    for(int i=0; i<numLine; i++) {
        if(var[i] == "R") {
            bool valid_unit = units.validate_length_unit(unit[i]);        
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                return false;
            }
        
            if(unit[i] == "" || unit[i] == "[]") {
                if(physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }    
            else {
                if(!physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }

            break;        
        }
    }

    for(int i=0; i<numLine; i++) {
        if(var[i] == "tau") {
            bool valid_unit = units.validate_time_unit(unit[i]);        
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                return false;
            }
        
            if(unit[i] == "" || unit[i] == "[]") {
                if(physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }    
            else {
                if(!physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }

            break;        
        }
    }

    for(int i=0; i<numLine; i++) {
        if(var[i] == "a_mb") {
            bool valid_unit = units.validate_time_unit(unit[i]);        
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                return false;
            }
        
            if(unit[i] == "" || unit[i] == "[]") {
                if(physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }    
            else {
                if(!physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }

            break;        
        }
    }
  
    return true;  
} 

bool Initializer::validate_units_orbital_cartesian(vector<string> &var, vector<string> &unit) {
    int numLine = var.size();
    
    vector<string> rs;
    rs.push_back("x");
    rs.push_back("y");
    rs.push_back("z");
    
    for(int k=0; k<3; k++) {
        for(int i=0; i<numLine; i++) {
            if(var[i] == rs[k]) {

                bool valid_unit = units.validate_length_unit(unit[i]);        
                if(!valid_unit) {
                    cerr << " " << endl;
                    cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                    cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                    cerr << " " << endl;
                    return false;
                }
            
                if(unit[i] == "" || unit[i] == "[]") {
                    if(physical_units) {
                        cerr << " " << endl;
                        cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                }    
                else {
                    if(!physical_units) {
                        cerr << " " << endl;
                        cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                }

                break;        
            }
        }
    }

    vector<string> vs;
    vs.push_back("vx");
    vs.push_back("vy");
    vs.push_back("vz");
    
    for(int k=0; k<3; k++) {
        for(int i=0; i<numLine; i++) {
            if(var[i] == vs[k]) {

                bool valid_unit = units.validate_speed_unit(unit[i]);        
                if(!valid_unit) {
                    cerr << " " << endl;
                    cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                    cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                    cerr << " " << endl;
                    return false;
                }
                
                if(unit[i] == "" || unit[i] == "[]") {
                    if(physical_units) {
                        cerr << " " << endl;
                        cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                }    
                else {
                    if(!physical_units) {
                        cerr << " " << endl;
                        cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                        cerr << " " << endl;
                        return false;
                    }
                }

                break;        
            }
        }
    }
  
    return true;  
} 
bool Initializer::validate_units_orbital_elliptical(vector<string> &var, vector<string> &unit) {
    int numLine = var.size();
    
    for(int i=0; i<numLine; i++) {
        if(var[i] == "a") {

            bool valid_unit = units.validate_length_unit(unit[i]);        
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                return false;
            }

            if(unit[i] == "" || unit[i] == "[]") {
                if(physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }    
            else {
                if(!physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }

            break;        
        }
    }

    vector<string> angles;
    angles.push_back("I");
    angles.push_back("O");
    angles.push_back("w");
    angles.push_back("M");
    
    for(int k=0; k<4; k++) {
        for(int i=0; i<numLine; i++) {
            if(var[i] == angles[k]) {

                bool valid_unit = units.validate_angular_unit(unit[i]);
                if(!valid_unit) {
                    cerr << " " << endl;
                    cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                    cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                    cerr << " " << endl;
                    return false;
                }                
                
                break;        
            }
        }
    }
    
    return true; 
} 
bool Initializer::validate_units_spin(vector<string> &var, vector<string> &unit) {
    int numLine = var.size();
    
    for(int i=0; i<numLine; i++) {
        if(var[i] == "lod") {

            bool valid_unit = units.validate_time_unit(unit[i]);
            if(!valid_unit) {
                cerr << " " << endl;
                cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                cerr << " " << endl;
                return false;
            }

            if(unit[i] == "" || unit[i] == "[]") {
                if(physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }    
            else {
                if(!physical_units) {
                    cerr << " " << endl;
                    cerr << "Unit inconsistency detected: " << var[i] << unit[i] << endl;
                    cerr << " " << endl;
                    return false;
                }
            }
            break;        
        }
    }

    vector<string> angles;
    angles.push_back("obl");
    angles.push_back("psi");
    
    for(int k=0; k<2; k++) {
        for(int i=0; i<numLine; i++) {
            if(var[i] == angles[k]) {

                bool valid_unit = units.validate_angular_unit(unit[i]);
                if(!valid_unit) {
                    cerr << " " << endl;
                    cerr << "Unrecognized physical unit: " << var[i] << unit[i] << endl;
                    cerr << "For a list of recognized units, run \"./tidymess.exe -u \"" << endl;
                    cerr << " " << endl;
                    return false;
                }   

                break;        
            }
        }
    }
    
    return true;
} 
        
// Aux functions    
inline bool Initializer::check_file_exists(const std::string& path) {
    struct stat s;   
  
    if( stat(path.c_str(), &s) == 0 ) {
        if( s.st_mode & S_IFDIR ) {
            cerr << endl;
            cerr << "Attempting to read a directory: " << path << endl;
            cerr << endl;
            return false;
        }
        else if( s.st_mode & S_IFREG ) {
            return true;
        }
        else {
            cerr << endl;
            cerr << "Cannot read file: " << path << endl;
            cerr << endl;
            return false;
        }    
    }
    else {
        cerr << endl;
        cerr << "Cannot open: " << path << endl;
        cerr << endl;
        return false;
    }
    
    return true;
}
bool Initializer::check_dir_exists(string &dir) {
    struct stat info;
    if( stat( dir.c_str(), &info ) != 0 ) return false;
    else if( info.st_mode & S_IFDIR ) return true; 
    else return false;
}        
vector<string> Initializer::get_words(string line) {
    istringstream iss(line);
    vector<string> results;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(results));
    return results;
}
vector< vector<string> > Initializer::read_file(string file_in) {
    // Check if regular file
    bool file_exist = check_file_exists(file_in);
    if(!file_exist) {
        exit(1);
    }
    
    // Read file
    vector< vector<string> > data;

    ifstream str;
    str.open(file_in.c_str());
    if(!str) {
        cerr << "Cannot open " << file_in << "!" << endl;
        exit(1);
    }
    else {
        string line;
        while(!str.eof()) {
            getline(str, line);
            vector<string> words = get_words(line);
            if(words.size() > 0) {
                data.push_back(words);
            }
        }
        str.close();
    }

    return data;    
}

// Coordinate transformation functions

vector<double> Initializer::rotZrotX(double anglez, double anglex, vector<double> vin) {
    double sz = sin(anglez);
    double sx = sin(anglex);
    double cz = cos(anglez);
    double cx = cos(anglex);

    double rin = cx*vin[1]-sx*vin[2];

    vector<double> vec(3);
    vec[0] = cz*vin[0]-sz*rin;
    vec[1] = sz*vin[0]+cz*rin;
    vec[2] = sx*vin[1]+cx*vin[2];
            
    return vec;
}

void Initializer::convert_spin_vectors_to_inertial() {
    for(int i=0; i<N; i++) {
        double P   = data_ic[i][6];
        double obl = data_ic[i][7];
        double psi = data_ic[i][8];

        if(P == 0) {
            data_ic[i][6] = 0.;
            data_ic[i][7] = 0.;
            data_ic[i][8] = 0.;
        }
        else if(P < 0) {
            cerr << " " << endl;
            cerr << "Negative length of day detected: " << data_ic[i][6] << endl;
            cerr << " " << endl;
            exit(1);
        }
        else {
            double wmag = 2*M_PI/P;

            vector<double> w_vec(3);
            w_vec[0] = 0;
            w_vec[1] = 0;
            w_vec[2] = wmag;

            w_vec = rotZrotX(psi, obl, w_vec);

            data_ic[i][6] = w_vec[0];
            data_ic[i][7] = w_vec[1];
            data_ic[i][8] = w_vec[2];
        }
    }
}
void Initializer::convert_spin_vectors_from_elliptical_body0abs() {
    double P   = data_ic[0][6];
    double obl = data_ic[0][7];
    double psi = data_ic[0][8];

    if(P == 0) {
        data_ic[0][6] = 0.;
        data_ic[0][7] = 0.;
        data_ic[0][8] = 0.;
    }
    else if(P < 0) {
        cerr << " " << endl;
        cerr << "Negative length of day detected: " << data_ic[0][6] << endl;
        cerr << " " << endl;
        exit(1);
    }
    else {
        double wmag = 2*M_PI/P;
    
        vector<double> spinvec(3);
        spinvec[0] = 0;
        spinvec[1] = 0;
        spinvec[2] = wmag;

        vector<double> spinvec3 = rotZrotX(psi, obl, spinvec);

        for(int k=0; k<3; k++) {
            data_ic[0][6+k] = spinvec3[k];
        }
    }

    for(int i=1; i<N; i++) {
        double P   = data_ic[i][6];
        double obl = data_ic[i][7];
        double psi = data_ic[i][8];

        if(P == 0) {
            data_ic[i][6] = 0.;
            data_ic[i][7] = 0.;
            data_ic[i][8] = 0.;
        }
        else if(P < 0) {
            cerr << " " << endl;
            cerr << "Negative length of day detected: " << data_ic[i][6] << endl;
            cerr << " " << endl;
            exit(1);
        }
        else {
            double wmag = 2*M_PI/P;
    
            vector<double> spinvec(3);
            spinvec[0] = 0;
            spinvec[1] = 0;
            spinvec[2] = wmag;

            double inc = data_ic[i][11];
            double O   = data_ic[i][12];

            vector<double> spinvec2 = rotZrotX(psi-O, obl, spinvec);
            vector<double> spinvec3 = rotZrotX(O, inc, spinvec2);

            for(int k=0; k<3; k++) {
                data_ic[i][6+k] = spinvec3[k];
            }
        }
    }
}
void Initializer::convert_spin_vectors_from_elliptical_body0rel() {
    double P   = data_ic[0][6];
    double obl = data_ic[0][7];
    double psi = data_ic[0][8];
    
    double inc = data_ic[1][11];
    double O   = data_ic[1][12];

    if(P == 0) {
        data_ic[0][6] = 0.;
        data_ic[0][7] = 0.;
        data_ic[0][8] = 0.;
    }
    else if(P < 0) {
        cerr << " " << endl;
        cerr << "Negative length of day detected: " << data_ic[0][6] << endl;
        cerr << " " << endl;
        exit(1);
    }
    else {
        double wmag = 2*M_PI/P;
    
        vector<double> spinvec(3);
        spinvec[0] = 0;
        spinvec[1] = 0;
        spinvec[2] = wmag;

        vector<double> spinvec2 = rotZrotX(psi-O, obl, spinvec);
        vector<double> spinvec3 = rotZrotX(O, inc, spinvec2);

        for(int k=0; k<3; k++) {
            data_ic[0][6+k] = spinvec3[k];
        }
    }

    for(int i=1; i<N; i++) {
        double P   = data_ic[i][6];
        double obl = data_ic[i][7];
        double psi = data_ic[i][8];

        if(P == 0) {
            data_ic[i][6] = 0.;
            data_ic[i][7] = 0.;
            data_ic[i][8] = 0.;
        }
        else if(P < 0) {
            cerr << " " << endl;
            cerr << "Negative length of day detected: " << data_ic[i][6] << endl;
            cerr << " " << endl;
            exit(1);
        }
        else {
            double wmag = 2*M_PI/P;
    
            vector<double> spinvec(3);
            spinvec[0] = 0;
            spinvec[1] = 0;
            spinvec[2] = wmag;

            inc = data_ic[i][11];
            O   = data_ic[i][12];

            vector<double> spinvec2 = rotZrotX(psi-O, obl, spinvec);
            vector<double> spinvec3 = rotZrotX(O, inc, spinvec2);

            for(int k=0; k<3; k++) {
                data_ic[i][6+k] = spinvec3[k];
            }
        }
    }
}

double Initializer::mean_to_eccentric_anomaly(double MA, double e) {
    double delta = 1e-14;
    
    int maxIter = 10000;
    int cntIter = 0;

    double EA = MA;
    double diff = 1.;

    if(e > 0.4) {
        while(diff > delta && cntIter < maxIter) {
            double EA0 = EA;
            EA = EA0 + (MA + e*sin(EA0) - EA0) / (1. - e*cos(EA0));
            diff = abs(EA-EA0);
            cntIter += 1;
        }    
    }
    else {
        while(diff > delta && cntIter < maxIter) {
            double EA0 = EA;
            EA = MA + e*sin(EA0);
            diff = abs(EA-EA0);
            cntIter += 1;
        }    
    }

    return EA;
}
double Initializer::eccentric_to_true_anomaly(double EA, double e) {
    double TA = 2. * atan2(sqrt(1.+e)*sin(EA/2.), sqrt(1.-e)*cos(EA/2.));
    return TA;
}
double Initializer::convert_mean_to_true_anomaly(double MA, double ecc) {
    double EA = mean_to_eccentric_anomaly(MA, ecc);
    double TA = eccentric_to_true_anomaly(EA, ecc);   
    if(TA < 0) TA += 2*M_PI;
    else if(TA >= 2*M_PI) TA -= 2*M_PI; 
    return TA;
}

// This implementation is based on orbital_elements.py from the 
// AMUSE software framework ( https://amusecode.github.io/ ).
vector<double> Initializer::get_relative_posvel_from_orbital_elements(double m1, double m2, double a, double ecc, double inc, double O, double w, double TA, double G) {
    double cos_TA = cos(TA);
    double sin_TA = sin(TA);

    double cos_inc = cos(inc);
    double sin_inc = sin(inc);

    double cos_w = cos(w);
    double sin_w = sin(w);

    double cos_O = cos(O);
    double sin_O = sin(O);

    // alpha is a unit vector directed along the line of node
    vector<double> alpha(3);
    alpha[0] = ( cos_O*cos_w - sin_O*sin_w*cos_inc );
    alpha[1] = ( sin_O*cos_w + cos_O*sin_w*cos_inc );
    alpha[2] = sin_w*sin_inc;

    // beta is a unit vector perpendicular to alpha and the orbital angular momentum vector
    vector<double> beta(3);
    beta[0] = ( - cos_O*sin_w - sin_O*cos_w*cos_inc );
    beta[1] = ( - sin_O*sin_w + cos_O*cos_w*cos_inc );
    beta[2] = cos_w*sin_inc;

    // Relative position and velocity
    double r = a*(1. - ecc*ecc) / (1. + ecc*cos_TA); 

    vector<double> dr(3);
    dr[0] = r*cos_TA*alpha[0] + r*sin_TA*beta[0];
    dr[1] = r*cos_TA*alpha[1] + r*sin_TA*beta[1];
    dr[2] = r*cos_TA*alpha[2] + r*sin_TA*beta[2];

    double mu = G*(m1 + m2);

    double dv_aux = 0;
    if(a != 0) {
        dv_aux = sqrt( mu / (a*(1. - ecc*ecc)) );
    }
    
    vector<double> dv(3);
    dv[0] = -1. * dv_aux * sin_TA * alpha[0] + dv_aux*(ecc + cos_TA)*beta[0];
    dv[1] = -1. * dv_aux * sin_TA * alpha[1] + dv_aux*(ecc + cos_TA)*beta[1];
    dv[2] = -1. * dv_aux * sin_TA * alpha[2] + dv_aux*(ecc + cos_TA)*beta[2];

    vector<double> posvel(6);
    posvel[0] = dr[0];
    posvel[1] = dr[1];
    posvel[2] = dr[2];
    posvel[3] = dv[0];
    posvel[4] = dv[1];
    posvel[5] = dv[2];

    return posvel;
}
void Initializer::move_to_center() {
    double M = 0;
    vector<double> rcm(3,0);
    vector<double> vcm(3,0);
    for(int i=0; i<N; i++) {
        M += data_ic[i][0];
        for(int k=0; k<3; k++) {
            rcm[k] += data_ic[i][0]*data_ic[i][9+k];
            vcm[k] += data_ic[i][0]*data_ic[i][12+k];
        }
    }
    
    for(int k=0; k<3; k++) {
        rcm[k] /= M;
        vcm[k] /= M;
    }
    
    for(int i=0; i<N; i++) {
        for(int k=0; k<3; k++) {
            data_ic[i][9+k] -= rcm[k];
            data_ic[i][12+k] -= vcm[k];
        }
    }
}
void Initializer::move_to_center(int Nenc) {
    double M = 0;
    vector<double> rcm(3,0);
    vector<double> vcm(3,0);
    for(int i=0; i<Nenc; i++) {
        M += data_ic[i][0];
        for(int k=0; k<3; k++) {
            rcm[k] += data_ic[i][0]*data_ic[i][9+k];
            vcm[k] += data_ic[i][0]*data_ic[i][12+k];
        }
    }
    
    for(int k=0; k<3; k++) {
        rcm[k] /= M;
        vcm[k] /= M;
    }
    
    for(int i=0; i<Nenc; i++) {
        for(int k=0; k<3; k++) {
            data_ic[i][9+k] -= rcm[k];
            data_ic[i][12+k] -= vcm[k];
        }
    }
}

void Initializer::convert_astrocentric_elements_to_cartesian_coordinates() {
    double G = 1; // N-body units
    if(physical_units) G = units.G_standard; // Gravitational constant in units of MSun, AU and yr
    
    for(int k=0; k<3; k++) {
        data_ic[0][9+k] = 0;
        data_ic[0][12+k] = 0;
    }

    for(int i=1; i<N; i++) {
        double TA = convert_mean_to_true_anomaly(data_ic[i][14], data_ic[i][10]);
        data_ic[i][14] = TA;

        vector<double> posvel = get_relative_posvel_from_orbital_elements(data_ic[0][0], data_ic[i][0], data_ic[i][9], data_ic[i][10], data_ic[i][11], data_ic[i][12], data_ic[i][13], data_ic[i][14], G);
        for(int k=0; k<3; k++) {
            data_ic[i][9+k] = posvel[k];
            data_ic[i][12+k] = posvel[3+k];
        }
    }

    move_to_center();
}
void Initializer::convert_jacobian_elements_to_cartesian_coordinates() {
    double G = 1; // N-body units
    if(physical_units) G = units.G_standard; // Gravitational constant in units of MSun, AU and yr
    
    for(int k=0; k<3; k++) {
        data_ic[0][9+k] = 0;
        data_ic[0][12+k] = 0;
    }

    double TA = convert_mean_to_true_anomaly(data_ic[1][14], data_ic[1][10]);
    data_ic[1][14] = TA;

    vector<double> posvel = get_relative_posvel_from_orbital_elements(data_ic[0][0], data_ic[1][0], data_ic[1][9], data_ic[1][10], data_ic[1][11], data_ic[1][12], data_ic[1][13], data_ic[1][14], G);
    for(int k=0; k<3; k++) {
        data_ic[1][9+k] = posvel[k];
        data_ic[1][12+k] = posvel[3+k];
    }

    int Nenc = 2;
    double Menc = data_ic[0][0] + data_ic[1][0];
    move_to_center(Nenc);

    for(int i=2; i<N; i++) {
        double TA = convert_mean_to_true_anomaly(data_ic[i][14], data_ic[i][10]);
        data_ic[i][14] = TA;

        vector<double> posvel = get_relative_posvel_from_orbital_elements(Menc, data_ic[i][0], data_ic[i][9], data_ic[i][10], data_ic[i][11], data_ic[i][12], data_ic[i][13], data_ic[i][14], G);
        for(int k=0; k<3; k++) {
            data_ic[i][9+k] = posvel[k];
            data_ic[i][12+k] = posvel[3+k];
        }
        
        Nenc++;
        Menc += data_ic[i][0];
        move_to_center(Nenc);
    }

    move_to_center();
}


     

