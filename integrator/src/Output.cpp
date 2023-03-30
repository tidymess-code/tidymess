#include "Output.h"

// Constructor
Output::Output() {
    format = 0;
    info = 0;
    coor = 0;
    dir = "data/";
    overwrite = true;

    physical_units = false;
    time_unit_output = "";
    mass_unit_output = "";
    length_unit_output = "";
    speed_unit_output = "";
    frequency_unit_output = "";
    angular_momentum_unit_output = "";
    energy_unit_output = "";
    inertia_unit_output = "";

    tidal_model = 4;

    file_out = "outfile";
    file_parbu = "";

    output_class = 0;

    firstWrite = true;
    firstWriteDiag = true;
    file_counter = 0;
    bin_number = 0;
    
    n_width = 25;
}

// Functions
void Output::set_format(int output_format) {
    this->format = output_format;
}
void Output::set_info(int output_info) {
    this->info = output_info;
}
void Output::set_coor(int output_coor) {
    this->coor = output_coor;
}
void Output::set_dir(string output_dir) {
    this->dir = output_dir;
    
    int N_str = this->dir.length();
    if(this->dir.at(N_str-1) != '/') {
        string s = "/";
        this->dir = this->dir + s;    
    }
}
void Output::set_overwrite(int overwrite) {
    switch(overwrite) {
        case 0:
            this->overwrite = false;
            break;
        case 1:
            this->overwrite = true;
            break;
    }
}

void Output::set_units(bool physical_units, string mass_unit, string length_unit, string time_unit, string speed_unit) {
    this->physical_units = physical_units;

    if(physical_units) {
        time_unit_output = time_unit;
        mass_unit_output = mass_unit;
        length_unit_output = length_unit;
        speed_unit_output = speed_unit;
        convert_time_to_frequency_unit(time_unit);
        construct_angular_momentum_unit(mass_unit, length_unit, speed_unit);
        construct_energy_unit(mass_unit, speed_unit);
        construct_inertia_unit(mass_unit, length_unit);
    }
}
void Output::set_conversion_factors(double Cm, double Cr, double Cv, double Ct) {
    units.Cm = Cm;
    units.Cr = Cr;
    units.Cv = Cv;
    units.Ct = Ct;
}
void Output::convert_time_to_frequency_unit(string time_unit) {
    if(time_unit == "") frequency_unit_output = "";
    else if(time_unit == "[]") frequency_unit_output = "[]";
    else {
        string s1 = "[";
        string s2 = "1/";
        int n = time_unit.length();
        string s3 = time_unit.substr(1, n);            
        
        frequency_unit_output = s1 + s2 + s3;
    }
}
void Output::construct_angular_momentum_unit(string mass_unit, string length_unit, string speed_unit) {
    if(mass_unit == "" || mass_unit == "[]") angular_momentum_unit_output = "";
    else {
        string s1 = "[";
                
        int n = mass_unit.length();
        string sm = mass_unit.substr(1, n-2);            

        n = length_unit.length();
        string sr = length_unit.substr(1, n-2);            
    
        n = speed_unit.length();
        string sv = speed_unit.substr(1, n-1);            
    
        angular_momentum_unit_output = s1 + sm + " " + sr + " " + sv;
    }
}
void Output::construct_energy_unit(string mass_unit, string speed_unit) {
    if(mass_unit == "" || mass_unit == "[]") energy_unit_output = "";
    else {
        string s1 = "[";
                
        int n = mass_unit.length();
        string sm = mass_unit.substr(1, n-2);            
    
        n = speed_unit.length();
        string sv = speed_unit.substr(1, n-2);            
        
        string s2 = "]";
        
        energy_unit_output = s1 + sm + " " + sv + " " + sv + s2;
    }
}
void Output::construct_inertia_unit(string mass_unit, string length_unit) {
    if(mass_unit == "" || mass_unit == "[]") inertia_unit_output = "";
    else {
        string s1 = "[";
                
        int n = mass_unit.length();
        string sm = mass_unit.substr(1, n-2);            

        n = length_unit.length();
        string sr = length_unit.substr(1, n-2);            

        string s2 = "]";
        
        inertia_unit_output = s1 + sm + " " + sr + " " + sr + s2;
    }
}

void Output::set_tidal_model(int tidal_model) {
    this->tidal_model = tidal_model;
}
    
void Output::set_file_out(string fo) {        
    this->file_out = dir + fo;
}

void Output::determine_output_mode() {
    if(format == 0) {
        output_class = 0;
    }
    else if(format == 1) {
        output_class = 1;
    }
    else if(format == 2) {
        output_class = 2;
    }
    else {
        cerr << endl;
        cerr << "Invalid output_format: " << format << endl;
        cerr << "Choose from [0=file_per_body, 1=file_per_snapshot, 2=single_file]." << endl;
        cerr << endl;
        exit(1);
    }

    if(info == 0) {
        output_class = output_class*10 + 0;
    }
    else if(info == 1) {
        output_class = output_class*10 + 1;
    }
    else {
        cerr << endl;
        cerr << "Invalid output_info: " << info << endl;
        cerr << "Choose from [0=time-varying quantities, 1=all]." << endl;
        cerr << endl;
        exit(1);
    }
    
    if(coor == 0) {
        output_class = output_class*10 + 0;
    }
    else {
        cerr << endl;
        cerr << "Invalid output_coor: " << info << endl;
        cerr << "Choose from [0=Cartesian]." << endl;
        cerr << endl;
        exit(1);
    }
     
    if(tidal_model == 0) {
        output_class = output_class*10 + 0;
    }
    else {
        output_class = output_class*10 + 1;
    }
}        
           
void Output::write_snapshot(double t, double tcpu, vector<Body> &bodies) {
    switch(output_class) {
    
        //0	0=file per body + 0=time-varying quantities + 0=cartesian inertial + N-body
        //100	0=file per body + 1=all quantities + 0=cartesian inertial + N-body	
        
	//1000	1=file per snapshot + 0=time-varying quantities + 0=cartesian inertial + N-body
	//1100	1=file per snapshot + 1=all quantities + 0=cartesian inertial + N-body 

    	//2000	2=single file + 0=time-varying quantities + 0=cartesian inertial + N-body
    	//2100	2=single file + 1=all quantities + 0=cartesian inertial + N-body

        //1	0=file per body + 0=time-varying quantities + 0=cartesian inertial + Tidal
        //101	0=file per body + 1=all quantities + 0=cartesian inertial + Tidal	
        
	//1001	1=file per snapshot + 0=time-varying quantities + 0=cartesian inertial + Tidal
	//1101	1=file per snapshot + 1=all quantities + 0=cartesian inertial + Tidal 

    	//2001	2=single file + 0=time-varying quantities + 0=cartesian inertial + Tidal
    	//2101	2=single file + 1=all quantities + 0=cartesian inertial + Tidal

        case 0:
            write_snapshot_per_body_compact_nbody(t, tcpu, bodies);
            break;
        case 100:
            write_snapshot_per_body_nbody(t, tcpu, bodies);
            break;

        case 1000:
            write_new_snapshot_file_compact_nbody(t, tcpu, bodies);
            break;
        case 1100:
            write_new_snapshot_file_nbody(t, tcpu, bodies);
            break;
        
        case 2000:
            write_snapshot_to_file_compact_nbody(t, tcpu, bodies);
            break;
        case 2100:
            write_snapshot_to_file_nbody(t, tcpu, bodies);
            break;

        case 1:
            write_snapshot_per_body_compact(t, tcpu, bodies);
            break;
        case 101:
            write_snapshot_per_body(t, tcpu, bodies);
            break;

        case 1001:
            write_new_snapshot_file_compact(t, tcpu, bodies);
            break;
        case 1101:
            write_new_snapshot_file(t, tcpu, bodies);
            break;
        
        case 2001:
            write_snapshot_to_file_compact(t, tcpu, bodies);
            break;
        case 2101:
            write_snapshot_to_file(t, tcpu, bodies);
            break;
    }    
}
void Output::create_output_file(string outputFile, bool isFirstWrite) {
    // Check if file already exists
    std::ifstream reader(outputFile);
    bool exist = reader.good();
    reader.close();

    if(exist == false) {
        // Create the file
        ofstream writer;
        writer.open(outputFile, ios::out);
        writer.close();
    }
    else {
        if(!overwrite) {
            cerr << " " << endl;
            cerr << outputFile << " already exists!" << endl;
            cerr << " " << endl;
            exit(1);        
        }
        else if(isFirstWrite && overwrite) {
            ofstream writer;
            writer.open(outputFile, ios::out);
            writer.close();
        }
    }
}

// Writers
void Output::write_snapshot_to_file(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".run";

    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << endl;
        cerr << "Cannot open output file: " << outputFile << endl;
        cerr << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    writer << units.convert_time_from_code_to_output(t, time_unit_output) << " " << N << " " << tcpu << endl;
    for(int i=0; i<N; i++) {
        writer << std::left;

        writer << setw(n_width) << bodies[i].name;
        writer << setw(n_width) << bodies[i].id;

        writer << std::right;
        
        writer << setw(n_width) << units.convert_mass_from_code_to_output(bodies[i].m, mass_unit_output);
        writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].R, length_unit_output);
        writer << setw(n_width) << bodies[i].xi;

        writer << setw(n_width) << bodies[i].kf;
        writer << setw(n_width) << units.convert_time_from_code_to_output(bodies[i].tau, time_unit_output);

        writer << setw(n_width) << units.convert_time_from_code_to_output(bodies[i].a_mb, time_unit_output);

        for(int j=0; j<6; j++) writer << setw(n_width) << units.convert_inertia_from_code_to_output(bodies[i].I[j], inertia_unit_output);
        
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_frequency_from_code_to_output(bodies[i].w[j], frequency_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    if(firstWrite) firstWrite = false;
}
void Output::write_snapshot_to_file_compact(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".run";

    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << endl;
        cerr << "Cannot open output file: " << outputFile << endl;
        cerr << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    writer << units.convert_time_from_code_to_output(t, time_unit_output) << " " << N << " " << tcpu << endl;
    for(int i=0; i<N; i++) {
        writer << std::right;
        
        for(int j=0; j<5; j++) writer << setw(n_width) << units.convert_inertia_from_code_to_output(bodies[i].I[j], inertia_unit_output);
                
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_frequency_from_code_to_output(bodies[i].w[j], frequency_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    if(firstWrite) firstWrite = false;
}

void Output::write_new_snapshot_file(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".s" + to_string(file_counter);

    firstWrite = true;
    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << "Cannot open output file: " << outputFile << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    string s =  "# t";
    writer << s + time_unit_output << " = ";
    writer << units.convert_time_from_code_to_output(t, time_unit_output) << endl;

    writer << "# N = " << N << endl;
    writer << "# tcpu[s] = " << tcpu << endl;

    if(firstWrite) {
        s = "# ";
        writer << s;

        writer << std::left;
            
        s =  "name";
        writer << setw(n_width-2) << s;
        s =  "id";
        writer << setw(n_width) << s;

        writer << std::right;
 
        s = "mass";
        writer << setw(n_width) << s + mass_unit_output;            
        s = "R";
        writer << setw(n_width) << s + length_unit_output;
        s = "xi";
        writer << setw(n_width) << s;
            
        s = "kf";
        writer << setw(n_width) << s;            
        s = "tau";
        writer << setw(n_width) << s + time_unit_output;
            
        s = "a_mb";
        writer << setw(n_width) << s + time_unit_output;
            
        s = "Ixx";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Ixy";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Ixz";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Iyy";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Iyz";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Izz";
        writer << setw(n_width) << s + inertia_unit_output;
            
        s = "wx";
        writer << setw(n_width) << s + frequency_unit_output;
        s = "wy";
        writer << setw(n_width) << s + frequency_unit_output;
        s = "wz";
        writer << setw(n_width) << s + frequency_unit_output;

        s = "x";
        writer << setw(n_width) << s + length_unit_output;
        s = "y";
        writer << setw(n_width) << s + length_unit_output;
        s = "z";
        writer << setw(n_width) << s + length_unit_output;

        s = "vx";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vy";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vz";
        writer << setw(n_width) << s + speed_unit_output;
            
        writer << endl;
    }
    
    for(int i=0; i<N; i++) {
        writer << std::left;
        
        writer << setw(n_width) << bodies[i].name;
        writer << setw(n_width) << bodies[i].id;

        writer << std::right;

        writer << setw(n_width) << units.convert_mass_from_code_to_output(bodies[i].m, mass_unit_output);
        writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].R, length_unit_output);
        writer << setw(n_width) << bodies[i].xi;

        writer << setw(n_width) << bodies[i].kf;
        writer << setw(n_width) << units.convert_time_from_code_to_output(bodies[i].tau, time_unit_output);

        writer << setw(n_width) << units.convert_time_from_code_to_output(bodies[i].a_mb, time_unit_output);

        for(int j=0; j<6; j++) writer << setw(n_width) << units.convert_inertia_from_code_to_output(bodies[i].I[j], inertia_unit_output);
                
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_frequency_from_code_to_output(bodies[i].w[j], frequency_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    file_counter++;
}
void Output::write_new_snapshot_file_compact(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".s" + to_string(file_counter);

    firstWrite = true;
    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << "Cannot open output file: " << outputFile << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    string s =  "# t";
    writer << s + time_unit_output << " = ";
    writer << units.convert_time_from_code_to_output(t, time_unit_output) << endl;

    writer << "# N = " << N << endl;
    writer << "# tcpu[s] = " << tcpu << endl;

    if(firstWrite) {
        s = "#";
        writer << s;
                                    
        s = "Ixx";
        writer << setw(n_width-1) << s + inertia_unit_output;
        s = "Ixy";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Ixz";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Iyy";
        writer << setw(n_width) << s + inertia_unit_output;
        s = "Iyz";
        writer << setw(n_width) << s + inertia_unit_output;
            
        s = "wx";
        writer << setw(n_width) << s + frequency_unit_output;
        s = "wy";
        writer << setw(n_width) << s + frequency_unit_output;
        s = "wz";
        writer << setw(n_width) << s + frequency_unit_output;

        s = "x";
        writer << setw(n_width) << s + length_unit_output;
        s = "y";
        writer << setw(n_width) << s + length_unit_output;
        s = "z";
        writer << setw(n_width) << s + length_unit_output;

        s = "vx";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vy";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vz";
        writer << setw(n_width) << s + speed_unit_output;
            
        writer << endl;
    }
    
    for(int i=0; i<N; i++) {
        for(int j=0; j<5; j++) writer << setw(n_width) << units.convert_inertia_from_code_to_output(bodies[i].I[j], inertia_unit_output);
                
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_frequency_from_code_to_output(bodies[i].w[j], frequency_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    file_counter++;
}
void Output::write_snapshot_per_body(double t, double tcpu, vector<Body> &bodies) {
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        string outputFile = dir + bodies[i].name + ".out";

        create_output_file(outputFile, firstWrite);

        ofstream writer;
        writer.open(outputFile, ios::app);

        if (!writer) {
            cerr << "Cannot open output file: " << outputFile << endl;
            exit(1);
        }

        writer << setprecision(16);
        writer << scientific;

        if(firstWrite) {
            string s = "#";
            writer << s;
            
            s =  "t";
            writer << setw(n_width-1) << s + time_unit_output;

            s = "mass";
            writer << setw(n_width) << s + mass_unit_output;            
            s = "R";
            writer << setw(n_width) << s + length_unit_output;
            s = "xi";
            writer << setw(n_width) << s;
            
            s = "kf";
            writer << setw(n_width) << s;            
            s = "tau";
            writer << setw(n_width) << s + time_unit_output;
            
            s = "a_mb";
            writer << setw(n_width) << s + time_unit_output;
            
            s = "Ixx";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Ixy";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Ixz";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Iyy";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Iyz";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Izz";
            writer << setw(n_width) << s + inertia_unit_output;
            
            s = "wx";
            writer << setw(n_width) << s + frequency_unit_output;
            s = "wy";
            writer << setw(n_width) << s + frequency_unit_output;
            s = "wz";
            writer << setw(n_width) << s + frequency_unit_output;

            s = "x";
            writer << setw(n_width) << s + length_unit_output;
            s = "y";
            writer << setw(n_width) << s + length_unit_output;
            s = "z";
            writer << setw(n_width) << s + length_unit_output;

            s = "vx";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vy";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vz";
            writer << setw(n_width) << s + speed_unit_output;
            
            writer << endl;
        }

        writer << setw(n_width) << units.convert_time_from_code_to_output(t, time_unit_output);

        writer << setw(n_width) << units.convert_mass_from_code_to_output(bodies[i].m, mass_unit_output);
        writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].R, length_unit_output);
        writer << setw(n_width) << bodies[i].xi;

        writer << setw(n_width) << bodies[i].kf;
        writer << setw(n_width) << units.convert_time_from_code_to_output(bodies[i].tau, time_unit_output);

        writer << setw(n_width) << units.convert_time_from_code_to_output(bodies[i].a_mb, time_unit_output);

        for(int j=0; j<6; j++) writer << setw(n_width) << units.convert_inertia_from_code_to_output(bodies[i].I[j], inertia_unit_output);
                
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_frequency_from_code_to_output(bodies[i].w[j], frequency_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;

        writer.close();
    }

    if(firstWrite) firstWrite = false;
}
void Output::write_snapshot_per_body_compact(double t, double tcpu, vector<Body> &bodies) {
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        string outputFile = dir + bodies[i].name + ".out";

        create_output_file(outputFile, firstWrite);

        ofstream writer;
        writer.open(outputFile, ios::app);

        if (!writer) {
            cerr << "Cannot open output file: " << outputFile << endl;
            exit(1);
        }

        writer << setprecision(16);
        writer << scientific;

        if(firstWrite) {
            string s = "#";
            writer << s;
            
            s =  "t";
            writer << setw(n_width-1) << s + time_unit_output;
            
            s = "Ixx";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Ixy";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Ixz";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Iyy";
            writer << setw(n_width) << s + inertia_unit_output;
            s = "Iyz";
            writer << setw(n_width) << s + inertia_unit_output;
            
            s = "wx";
            writer << setw(n_width) << s + frequency_unit_output;
            s = "wy";
            writer << setw(n_width) << s + frequency_unit_output;
            s = "wz";
            writer << setw(n_width) << s + frequency_unit_output;

            s = "x";
            writer << setw(n_width) << s + length_unit_output;
            s = "y";
            writer << setw(n_width) << s + length_unit_output;
            s = "z";
            writer << setw(n_width) << s + length_unit_output;

            s = "vx";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vy";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vz";
            writer << setw(n_width) << s + speed_unit_output;
            
            writer << endl;
        }
        
        writer << setw(n_width) << units.convert_time_from_code_to_output(t, time_unit_output);

        for(int j=0; j<5; j++) writer << setw(n_width) << units.convert_inertia_from_code_to_output(bodies[i].I[j], inertia_unit_output);
                
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_frequency_from_code_to_output(bodies[i].w[j], frequency_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;

        writer.close();
    }

    if(firstWrite) firstWrite = false;
}

void Output::write_snapshot_to_file_nbody(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".run";

    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << endl;
        cerr << "Cannot open output file: " << outputFile << endl;
        cerr << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    writer << units.convert_time_from_code_to_output(t, time_unit_output) << " " << N << " " << tcpu << endl;
    for(int i=0; i<N; i++) {
        writer << std::left;

        writer << setw(n_width) << bodies[i].name;
        writer << setw(n_width) << bodies[i].id;

        writer << std::right;
        
        writer << setw(n_width) << units.convert_mass_from_code_to_output(bodies[i].m, mass_unit_output);
        writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].R, length_unit_output);
 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    if(firstWrite) firstWrite = false;
}
void Output::write_snapshot_to_file_compact_nbody(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".run";

    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << endl;
        cerr << "Cannot open output file: " << outputFile << endl;
        cerr << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    writer << units.convert_time_from_code_to_output(t, time_unit_output) << " " << N << " " << tcpu << endl;
    for(int i=0; i<N; i++) {
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    if(firstWrite) firstWrite = false;
}
void Output::write_new_snapshot_file_nbody(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".s" + to_string(file_counter);

    firstWrite = true;
    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << "Cannot open output file: " << outputFile << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    string s =  "# t";
    writer << s + time_unit_output << " = ";
    writer << units.convert_time_from_code_to_output(t, time_unit_output) << endl;

    writer << "# N = " << N << endl;
    writer << "# tcpu[s] = " << tcpu << endl;

    if(firstWrite) {
        s = "# ";
        writer << s;

        writer << std::left;
            
        s =  "name";
        writer << setw(n_width-2) << s;
        s =  "id";
        writer << setw(n_width) << s;

        writer << std::right;
 
        s = "mass";
        writer << setw(n_width) << s + mass_unit_output;            
        s = "R";
        writer << setw(n_width) << s + length_unit_output;

        s = "x";
        writer << setw(n_width) << s + length_unit_output;
        s = "y";
        writer << setw(n_width) << s + length_unit_output;
        s = "z";
        writer << setw(n_width) << s + length_unit_output;

        s = "vx";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vy";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vz";
        writer << setw(n_width) << s + speed_unit_output;
            
        writer << endl;
    }
    
    for(int i=0; i<N; i++) {
        writer << std::left;
        
        writer << setw(n_width) << bodies[i].name;
        writer << setw(n_width) << bodies[i].id;

        writer << std::right;

        writer << setw(n_width) << units.convert_mass_from_code_to_output(bodies[i].m, mass_unit_output);
        writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].R, length_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    file_counter++;
}
void Output::write_new_snapshot_file_compact_nbody(double t, double tcpu, vector<Body> &bodies) {
    string outputFile = file_out + ".s" + to_string(file_counter);

    firstWrite = true;
    create_output_file(outputFile, firstWrite);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << "Cannot open output file: " << outputFile << endl;
        exit(1);
    }

    int N = bodies.size();

    writer << setprecision(16);
    writer << scientific;

    string s =  "# t";
    writer << s + time_unit_output << " = ";
    writer << units.convert_time_from_code_to_output(t, time_unit_output) << endl;

    writer << "# N = " << N << endl;
    writer << "# tcpu[s] = " << tcpu << endl;

    if(firstWrite) {
        s = "# ";
        writer << s;

        s = "x";
        writer << setw(n_width-2) << s + length_unit_output;
        s = "y";
        writer << setw(n_width) << s + length_unit_output;
        s = "z";
        writer << setw(n_width) << s + length_unit_output;

        s = "vx";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vy";
        writer << setw(n_width) << s + speed_unit_output;
        s = "vz";
        writer << setw(n_width) << s + speed_unit_output;
            
        writer << endl;
    }
    
    for(int i=0; i<N; i++) {
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;
    }
    writer << endl;

    writer.close();

    file_counter++;
}
void Output::write_snapshot_per_body_nbody(double t, double tcpu, vector<Body> &bodies) {
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        string outputFile = dir + bodies[i].name + ".out";

        create_output_file(outputFile, firstWrite);

        ofstream writer;
        writer.open(outputFile, ios::app);

        if (!writer) {
            cerr << "Cannot open output file: " << outputFile << endl;
            exit(1);
        }

        writer << setprecision(16);
        writer << scientific;

        if(firstWrite) {
            string s = "#";
            writer << s;
            
            s =  "t";
            writer << setw(n_width-1) << s + time_unit_output;

            s = "mass";
            writer << setw(n_width) << s + mass_unit_output;            
            s = "R";
            writer << setw(n_width) << s + length_unit_output;

            s = "x";
            writer << setw(n_width) << s + length_unit_output;
            s = "y";
            writer << setw(n_width) << s + length_unit_output;
            s = "z";
            writer << setw(n_width) << s + length_unit_output;

            s = "vx";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vy";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vz";
            writer << setw(n_width) << s + speed_unit_output;
            
            writer << endl;
        }

        writer << setw(n_width) << units.convert_time_from_code_to_output(t, time_unit_output);

        writer << setw(n_width) << units.convert_mass_from_code_to_output(bodies[i].m, mass_unit_output);
        writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].R, length_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;

        writer.close();
    }

    if(firstWrite) firstWrite = false;
}
void Output::write_snapshot_per_body_compact_nbody(double t, double tcpu, vector<Body> &bodies) {
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        string outputFile = dir + bodies[i].name + ".out";

        create_output_file(outputFile, firstWrite);

        ofstream writer;
        writer.open(outputFile, ios::app);

        if (!writer) {
            cerr << "Cannot open output file: " << outputFile << endl;
            exit(1);
        }

        writer << setprecision(16);
        writer << scientific;

        if(firstWrite) {
            string s = "#";
            writer << s;
            
            s =  "t";
            writer << setw(n_width-1) << s + time_unit_output;

            s = "x";
            writer << setw(n_width) << s + length_unit_output;
            s = "y";
            writer << setw(n_width) << s + length_unit_output;
            s = "z";
            writer << setw(n_width) << s + length_unit_output;

            s = "vx";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vy";
            writer << setw(n_width) << s + speed_unit_output;
            s = "vz";
            writer << setw(n_width) << s + speed_unit_output;
            
            writer << endl;
        }

        writer << setw(n_width) << units.convert_time_from_code_to_output(t, time_unit_output);

        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_length_from_code_to_output(bodies[i].r[j], length_unit_output); 
        for(int j=0; j<3; j++) writer << setw(n_width) << units.convert_speed_from_code_to_output(bodies[i].v[j], speed_unit_output); 
        writer << endl;

        writer.close();
    }

    if(firstWrite) firstWrite = false;
}

void Output::save_to_binary(double &t, int &N, double &tcpu, double &dt_prev, int &num_integration_steps, int &collision_flag, int &roche_flag, int &breakup_flag, vector<Body> &bodies, double &dt_snapshot, double &dt0_log, double &fmul_log, int &num_snapshot) {
    string bin_name = file_out + ".bin" + to_string(bin_number);
    bin_number = (bin_number + 1) % 2;

    fstream fb(bin_name, ios::in | ios::out | ios::binary | ios::trunc);

    // header
    fb.write((char*)&t, sizeof (double));
    fb.write((char*)&N, sizeof (int));
    fb.write((char*)&tcpu, sizeof (double));
    fb.write((char*)&dt_prev, sizeof (double));
    fb.write((char*)&num_integration_steps, sizeof (int));
     
    fb.write((char*)&collision_flag, sizeof (int));
    fb.write((char*)&roche_flag, sizeof (int));
    fb.write((char*)&breakup_flag, sizeof (int));     
 
    fb.write((char*)&dt_snapshot, sizeof (double));
    fb.write((char*)&dt0_log, sizeof (double));
    fb.write((char*)&fmul_log, sizeof (double));

    fb.write((char*)&num_snapshot, sizeof (int));
 
    // id
    vector<int> id(N);
    for(int i=0; i<N; i++) id[i] = bodies[i].id;
    
    fb.write((char*)&id[0], sizeof (int)*N);
    
    // name
    for(int i=0; i<N; i++) {
        string name = bodies[i].name;
        int name_len = name.length();                
        fb.write((char*)&name_len, sizeof (int));
        fb.write((char*)&name[0], sizeof (char)*name_len);        
    }
    
    // internal props
    vector<double> m(N), R(N), xi(N);
    for(int i=0; i<N; i++) {
        m[i]  = bodies[i].m;
        R[i]  = bodies[i].R;
        xi[i] = bodies[i].xi;
    }

    fb.write((char*)&m[0], sizeof (double)*N);
    fb.write((char*)&R[0], sizeof (double)*N);
    fb.write((char*)&xi[0], sizeof (double)*N);   
   
    // Tidal parameters
    vector<double> kf(N), tau(N);
    for(int i=0; i<N; i++) {
        kf[i]  = bodies[i].kf;
        tau[i] = bodies[i].tau;
    }

    fb.write((char*)&kf[0], sizeof (double)*N);
    fb.write((char*)&tau[0], sizeof (double)*N);
    
    // Other physical ingredients
    vector<double> a_mb(N);
    for(int i=0; i<N; i++) {
        a_mb[i]  = bodies[i].a_mb;
    }

    fb.write((char*)&a_mb[0], sizeof (double)*N);

    // Inertia tensors
    vector<double> Ixx(N), Ixy(N), Ixz(N), Iyy(N), Iyz(N), Izz(N);
    for(int i=0; i<N; i++) {
        Ixx[i] = bodies[i].I[0];
        Ixy[i] = bodies[i].I[1];
        Ixz[i] = bodies[i].I[2];
        Iyy[i] = bodies[i].I[3];
        Iyz[i] = bodies[i].I[4];
        Izz[i] = bodies[i].I[5];
    }

    fb.write((char*)&Ixx[0], sizeof (double)*N);
    fb.write((char*)&Ixy[0], sizeof (double)*N);
    fb.write((char*)&Ixz[0], sizeof (double)*N);
    fb.write((char*)&Iyy[0], sizeof (double)*N);
    fb.write((char*)&Iyz[0], sizeof (double)*N);
    fb.write((char*)&Izz[0], sizeof (double)*N);

    vector<double> Ixx_p(N), Ixy_p(N), Ixz_p(N), Iyy_p(N), Iyz_p(N), Izz_p(N);
    for(int i=0; i<N; i++) {
        Ixx_p[i] = bodies[i].I_p[0];
        Ixy_p[i] = bodies[i].I_p[1];
        Ixz_p[i] = bodies[i].I_p[2];
        Iyy_p[i] = bodies[i].I_p[3];
        Iyz_p[i] = bodies[i].I_p[4];
        Izz_p[i] = bodies[i].I_p[5];
    }

    fb.write((char*)&Ixx_p[0], sizeof (double)*N);
    fb.write((char*)&Ixy_p[0], sizeof (double)*N);
    fb.write((char*)&Ixz_p[0], sizeof (double)*N);
    fb.write((char*)&Iyy_p[0], sizeof (double)*N);
    fb.write((char*)&Iyz_p[0], sizeof (double)*N);
    fb.write((char*)&Izz_p[0], sizeof (double)*N);

    vector<double> Ixx_n(N), Ixy_n(N), Ixz_n(N), Iyy_n(N), Iyz_n(N), Izz_n(N);
    for(int i=0; i<N; i++) {
        Ixx_n[i] = bodies[i].I_n[0];
        Ixy_n[i] = bodies[i].I_n[1];
        Ixz_n[i] = bodies[i].I_n[2];
        Iyy_n[i] = bodies[i].I_n[3];
        Iyz_n[i] = bodies[i].I_n[4];
        Izz_n[i] = bodies[i].I_n[5];
    }

    fb.write((char*)&Ixx_n[0], sizeof (double)*N);
    fb.write((char*)&Ixy_n[0], sizeof (double)*N);
    fb.write((char*)&Ixz_n[0], sizeof (double)*N);
    fb.write((char*)&Iyy_n[0], sizeof (double)*N);
    fb.write((char*)&Iyz_n[0], sizeof (double)*N);
    fb.write((char*)&Izz_n[0], sizeof (double)*N);

    vector<double> Ixx_inv(N), Ixy_inv(N), Ixz_inv(N), Iyy_inv(N), Iyz_inv(N), Izz_inv(N);
    for(int i=0; i<N; i++) {
        Ixx_inv[i] = bodies[i].I_inv[0];
        Ixy_inv[i] = bodies[i].I_inv[1];
        Ixz_inv[i] = bodies[i].I_inv[2];
        Iyy_inv[i] = bodies[i].I_inv[3];
        Iyz_inv[i] = bodies[i].I_inv[4];
        Izz_inv[i] = bodies[i].I_inv[5];
    }

    fb.write((char*)&Ixx_inv[0], sizeof (double)*N);
    fb.write((char*)&Ixy_inv[0], sizeof (double)*N);
    fb.write((char*)&Ixz_inv[0], sizeof (double)*N);
    fb.write((char*)&Iyy_inv[0], sizeof (double)*N);
    fb.write((char*)&Iyz_inv[0], sizeof (double)*N);
    fb.write((char*)&Izz_inv[0], sizeof (double)*N);

    vector<double> Ixx_e_r(N), Ixy_e_r(N), Ixz_e_r(N), Iyy_e_r(N), Iyz_e_r(N), Izz_e_r(N);
    for(int i=0; i<N; i++) {
        Ixx_e_r[i] = bodies[i].I_e_r[0];
        Ixy_e_r[i] = bodies[i].I_e_r[1];
        Ixz_e_r[i] = bodies[i].I_e_r[2];
        Iyy_e_r[i] = bodies[i].I_e_r[3];
        Iyz_e_r[i] = bodies[i].I_e_r[4];
        Izz_e_r[i] = bodies[i].I_e_r[5];
    }

    fb.write((char*)&Ixx_e_r[0], sizeof (double)*N);
    fb.write((char*)&Ixy_e_r[0], sizeof (double)*N);
    fb.write((char*)&Ixz_e_r[0], sizeof (double)*N);
    fb.write((char*)&Iyy_e_r[0], sizeof (double)*N);
    fb.write((char*)&Iyz_e_r[0], sizeof (double)*N);
    fb.write((char*)&Izz_e_r[0], sizeof (double)*N);

    vector<double> dIxx_e_r(N), dIxy_e_r(N), dIxz_e_r(N), dIyy_e_r(N), dIyz_e_r(N), dIzz_e_r(N);
    for(int i=0; i<N; i++) {
        dIxx_e_r[i] = bodies[i].dI_e_r[0];
        dIxy_e_r[i] = bodies[i].dI_e_r[1];
        dIxz_e_r[i] = bodies[i].dI_e_r[2];
        dIyy_e_r[i] = bodies[i].dI_e_r[3];
        dIyz_e_r[i] = bodies[i].dI_e_r[4];
        dIzz_e_r[i] = bodies[i].dI_e_r[5];
    }

    fb.write((char*)&dIxx_e_r[0], sizeof (double)*N);
    fb.write((char*)&dIxy_e_r[0], sizeof (double)*N);
    fb.write((char*)&dIxz_e_r[0], sizeof (double)*N);
    fb.write((char*)&dIyy_e_r[0], sizeof (double)*N);
    fb.write((char*)&dIyz_e_r[0], sizeof (double)*N);
    fb.write((char*)&dIzz_e_r[0], sizeof (double)*N);

    vector<double> Ixx_e_w(N), Ixy_e_w(N), Ixz_e_w(N), Iyy_e_w(N), Iyz_e_w(N), Izz_e_w(N);
    for(int i=0; i<N; i++) {
        Ixx_e_w[i] = bodies[i].I_e_w[0];
        Ixy_e_w[i] = bodies[i].I_e_w[1];
        Ixz_e_w[i] = bodies[i].I_e_w[2];
        Iyy_e_w[i] = bodies[i].I_e_w[3];
        Iyz_e_w[i] = bodies[i].I_e_w[4];
        Izz_e_w[i] = bodies[i].I_e_w[5];
    }

    fb.write((char*)&Ixx_e_w[0], sizeof (double)*N);
    fb.write((char*)&Ixy_e_w[0], sizeof (double)*N);
    fb.write((char*)&Ixz_e_w[0], sizeof (double)*N);
    fb.write((char*)&Iyy_e_w[0], sizeof (double)*N);
    fb.write((char*)&Iyz_e_w[0], sizeof (double)*N);
    fb.write((char*)&Izz_e_w[0], sizeof (double)*N);

    vector<double> Ixx_e(N), Ixy_e(N), Ixz_e(N), Iyy_e(N), Iyz_e(N), Izz_e(N);
    for(int i=0; i<N; i++) {
        Ixx_e[i] = bodies[i].I_e[0];
        Ixy_e[i] = bodies[i].I_e[1];
        Ixz_e[i] = bodies[i].I_e[2];
        Iyy_e[i] = bodies[i].I_e[3];
        Iyz_e[i] = bodies[i].I_e[4];
        Izz_e[i] = bodies[i].I_e[5];
    }

    fb.write((char*)&Ixx_e[0], sizeof (double)*N);
    fb.write((char*)&Ixy_e[0], sizeof (double)*N);
    fb.write((char*)&Ixz_e[0], sizeof (double)*N);
    fb.write((char*)&Iyy_e[0], sizeof (double)*N);
    fb.write((char*)&Iyz_e[0], sizeof (double)*N);
    fb.write((char*)&Izz_e[0], sizeof (double)*N);

    vector<double> dIxx_e(N), dIxy_e(N), dIxz_e(N), dIyy_e(N), dIyz_e(N), dIzz_e(N);
    for(int i=0; i<N; i++) {
        dIxx_e[i] = bodies[i].dI_e[0];
        dIxy_e[i] = bodies[i].dI_e[1];
        dIxz_e[i] = bodies[i].dI_e[2];
        dIyy_e[i] = bodies[i].dI_e[3];
        dIyz_e[i] = bodies[i].dI_e[4];
        dIzz_e[i] = bodies[i].dI_e[5];
    }

    fb.write((char*)&dIxx_e[0], sizeof (double)*N);
    fb.write((char*)&dIxy_e[0], sizeof (double)*N);
    fb.write((char*)&dIxz_e[0], sizeof (double)*N);
    fb.write((char*)&dIyy_e[0], sizeof (double)*N);
    fb.write((char*)&dIyz_e[0], sizeof (double)*N);
    fb.write((char*)&dIzz_e[0], sizeof (double)*N);

    vector<double> dIxx_n(N), dIxy_n(N), dIxz_n(N), dIyy_n(N), dIyz_n(N), dIzz_n(N);
    for(int i=0; i<N; i++) {
        dIxx_n[i] = bodies[i].dI_n[0];
        dIxy_n[i] = bodies[i].dI_n[1];
        dIxz_n[i] = bodies[i].dI_n[2];
        dIyy_n[i] = bodies[i].dI_n[3];
        dIyz_n[i] = bodies[i].dI_n[4];
        dIzz_n[i] = bodies[i].dI_n[5];
    }

    fb.write((char*)&dIxx_n[0], sizeof (double)*N);
    fb.write((char*)&dIxy_n[0], sizeof (double)*N);
    fb.write((char*)&dIxz_n[0], sizeof (double)*N);
    fb.write((char*)&dIyy_n[0], sizeof (double)*N);
    fb.write((char*)&dIyz_n[0], sizeof (double)*N);
    fb.write((char*)&dIzz_n[0], sizeof (double)*N);
    
    vector<double> Ixx_e_prev(N), Ixy_e_prev(N), Ixz_e_prev(N), Iyy_e_prev(N), Iyz_e_prev(N), Izz_e_prev(N);
    for(int i=0; i<N; i++) {
        Ixx_e_prev[i] = bodies[i].I_e_prev[0];
        Ixy_e_prev[i] = bodies[i].I_e_prev[1];
        Ixz_e_prev[i] = bodies[i].I_e_prev[2];
        Iyy_e_prev[i] = bodies[i].I_e_prev[3];
        Iyz_e_prev[i] = bodies[i].I_e_prev[4];
        Izz_e_prev[i] = bodies[i].I_e_prev[5];
    }

    fb.write((char*)&Ixx_e_prev[0], sizeof (double)*N);
    fb.write((char*)&Ixy_e_prev[0], sizeof (double)*N);
    fb.write((char*)&Ixz_e_prev[0], sizeof (double)*N);
    fb.write((char*)&Iyy_e_prev[0], sizeof (double)*N);
    fb.write((char*)&Iyz_e_prev[0], sizeof (double)*N);
    fb.write((char*)&Izz_e_prev[0], sizeof (double)*N);
    
    vector<double> Ixx_n_prev(N), Ixy_n_prev(N), Ixz_n_prev(N), Iyy_n_prev(N), Iyz_n_prev(N), Izz_n_prev(N);
    for(int i=0; i<N; i++) {
        Ixx_n_prev[i] = bodies[i].I_n_prev[0];
        Ixy_n_prev[i] = bodies[i].I_n_prev[1];
        Ixz_n_prev[i] = bodies[i].I_n_prev[2];
        Iyy_n_prev[i] = bodies[i].I_n_prev[3];
        Iyz_n_prev[i] = bodies[i].I_n_prev[4];
        Izz_n_prev[i] = bodies[i].I_n_prev[5];
    }

    fb.write((char*)&Ixx_n_prev[0], sizeof (double)*N);
    fb.write((char*)&Ixy_n_prev[0], sizeof (double)*N);
    fb.write((char*)&Ixz_n_prev[0], sizeof (double)*N);
    fb.write((char*)&Iyy_n_prev[0], sizeof (double)*N);
    fb.write((char*)&Iyz_n_prev[0], sizeof (double)*N);
    fb.write((char*)&Izz_n_prev[0], sizeof (double)*N);

    vector<double> Ixx_e_prev_bu(N), Ixy_e_prev_bu(N), Ixz_e_prev_bu(N), Iyy_e_prev_bu(N), Iyz_e_prev_bu(N), Izz_e_prev_bu(N);
    for(int i=0; i<N; i++) {
        Ixx_e_prev_bu[i] = bodies[i].I_e_prev_bu[0];
        Ixy_e_prev_bu[i] = bodies[i].I_e_prev_bu[1];
        Ixz_e_prev_bu[i] = bodies[i].I_e_prev_bu[2];
        Iyy_e_prev_bu[i] = bodies[i].I_e_prev_bu[3];
        Iyz_e_prev_bu[i] = bodies[i].I_e_prev_bu[4];
        Izz_e_prev_bu[i] = bodies[i].I_e_prev_bu[5];
    }

    fb.write((char*)&Ixx_e_prev_bu[0], sizeof (double)*N);
    fb.write((char*)&Ixy_e_prev_bu[0], sizeof (double)*N);
    fb.write((char*)&Ixz_e_prev_bu[0], sizeof (double)*N);
    fb.write((char*)&Iyy_e_prev_bu[0], sizeof (double)*N);
    fb.write((char*)&Iyz_e_prev_bu[0], sizeof (double)*N);
    fb.write((char*)&Izz_e_prev_bu[0], sizeof (double)*N);

    vector<double> Ixx_e_rh(N), Ixy_e_rh(N), Ixz_e_rh(N), Iyy_e_rh(N), Iyz_e_rh(N), Izz_e_rh(N);
    for(int i=0; i<N; i++) {
        Ixx_e_rh[i] = bodies[i].I_e_rh[0];
        Ixy_e_rh[i] = bodies[i].I_e_rh[1];
        Ixz_e_rh[i] = bodies[i].I_e_rh[2];
        Iyy_e_rh[i] = bodies[i].I_e_rh[3];
        Iyz_e_rh[i] = bodies[i].I_e_rh[4];
        Izz_e_rh[i] = bodies[i].I_e_rh[5];
    }

    fb.write((char*)&Ixx_e_rh[0], sizeof (double)*N);
    fb.write((char*)&Ixy_e_rh[0], sizeof (double)*N);
    fb.write((char*)&Ixz_e_rh[0], sizeof (double)*N);
    fb.write((char*)&Iyy_e_rh[0], sizeof (double)*N);
    fb.write((char*)&Iyz_e_rh[0], sizeof (double)*N);
    fb.write((char*)&Izz_e_rh[0], sizeof (double)*N);

    // Spin and angular momenta
    vector<double> wx(N), wy(N), wz(N);
    vector<double> Lx(N), Ly(N), Lz(N);

    for(int i=0; i<N; i++) {
        wx[i] = bodies[i].w[0];
        wy[i] = bodies[i].w[1];
        wz[i] = bodies[i].w[2];
        Lx[i] = bodies[i].L[0];
        Ly[i] = bodies[i].L[1];
        Lz[i] = bodies[i].L[2];
    }

    fb.write((char*)&wx[0], sizeof (double)*N);
    fb.write((char*)&wy[0], sizeof (double)*N);
    fb.write((char*)&wz[0], sizeof (double)*N);
    fb.write((char*)&Lx[0], sizeof (double)*N);
    fb.write((char*)&Ly[0], sizeof (double)*N);
    fb.write((char*)&Lz[0], sizeof (double)*N);

    vector<double> Tx(N), Ty(N), Tz(N);

    for(int i=0; i<N; i++) {
        Tx[i] = bodies[i].T[0];
        Ty[i] = bodies[i].T[1];
        Tz[i] = bodies[i].T[2];
    }

    fb.write((char*)&Tx[0], sizeof (double)*N);
    fb.write((char*)&Ty[0], sizeof (double)*N);
    fb.write((char*)&Tz[0], sizeof (double)*N);
    
    // pos and vel
    vector<double> x(N), y(N), z(N);
    vector<double> vx(N), vy(N), vz(N);

    for(int i=0; i<N; i++) {
        x[i]  = bodies[i].r[0];
        y[i]  = bodies[i].r[1];
        z[i]  = bodies[i].r[2];
        vx[i] = bodies[i].v[0];
        vy[i] = bodies[i].v[1];
        vz[i] = bodies[i].v[2];
    }

    fb.write((char*)&x[0], sizeof (double)*N);
    fb.write((char*)&y[0], sizeof (double)*N);
    fb.write((char*)&z[0], sizeof (double)*N);
    fb.write((char*)&vx[0], sizeof (double)*N);
    fb.write((char*)&vy[0], sizeof (double)*N);
    fb.write((char*)&vz[0], sizeof (double)*N);

    vector<double> ax(N), ay(N), az(N);

    for(int i=0; i<N; i++) {
        ax[i]  = bodies[i].a[0];
        ay[i]  = bodies[i].a[1];
        az[i]  = bodies[i].a[2];
    }

    fb.write((char*)&ax[0], sizeof (double)*N);
    fb.write((char*)&ay[0], sizeof (double)*N);
    fb.write((char*)&az[0], sizeof (double)*N);
    
    // Aux quantities
    vector<double> R5(N), R5_3(N), kf_R5(N), kf_R5_3(N), tau_inv(N);
    for(int i=0; i<N; i++) {
        R5[i]      = bodies[i].R5;
        R5_3[i]    = bodies[i].R5_3;
        kf_R5[i]   = bodies[i].kf_R5;
        kf_R5_3[i] = bodies[i].kf_R5_3;
        tau_inv[i] = bodies[i].tau_inv;
    }

    fb.write((char*)&R5[0], sizeof (double)*N);
    fb.write((char*)&R5_3[0], sizeof (double)*N);
    fb.write((char*)&kf_R5[0], sizeof (double)*N);   
    fb.write((char*)&kf_R5_3[0], sizeof (double)*N);   
    fb.write((char*)&tau_inv[0], sizeof (double)*N);   
    
    vector<double> vvx(N), vvy(N), vvz(N);

    for(int i=0; i<N; i++) {
        vvx[i]  = bodies[i].vv[0];
        vvy[i]  = bodies[i].vv[1];
        vvz[i]  = bodies[i].vv[2];
    }

    fb.write((char*)&vvx[0], sizeof (double)*N);
    fb.write((char*)&vvy[0], sizeof (double)*N);
    fb.write((char*)&vvz[0], sizeof (double)*N);

    // Auxiliary variables
    vector<double> Jxx_n(N), Jxy_n(N), Jxz_n(N), Jyy_n(N), Jyz_n(N), Jzz_n(N);
    for(int i=0; i<N; i++) {
        Jxx_n[i] = bodies[i].J_n[0];
        Jxy_n[i] = bodies[i].J_n[1];
        Jxz_n[i] = bodies[i].J_n[2];
        Jyy_n[i] = bodies[i].J_n[3];
        Jyz_n[i] = bodies[i].J_n[4];
        Jzz_n[i] = bodies[i].J_n[5];
    }

    fb.write((char*)&Jxx_n[0], sizeof (double)*N);
    fb.write((char*)&Jxy_n[0], sizeof (double)*N);
    fb.write((char*)&Jxz_n[0], sizeof (double)*N);
    fb.write((char*)&Jyy_n[0], sizeof (double)*N);
    fb.write((char*)&Jyz_n[0], sizeof (double)*N);
    fb.write((char*)&Jzz_n[0], sizeof (double)*N);

    vector<double> Jxx(N), Jxy(N), Jxz(N), Jyy(N), Jyz(N), Jzz(N);
    for(int i=0; i<N; i++) {
        Jxx[i] = bodies[i].J[0];
        Jxy[i] = bodies[i].J[1];
        Jxz[i] = bodies[i].J[2];
        Jyy[i] = bodies[i].J[3];
        Jyz[i] = bodies[i].J[4];
        Jzz[i] = bodies[i].J[5];
    }

    fb.write((char*)&Jxx[0], sizeof (double)*N);
    fb.write((char*)&Jxy[0], sizeof (double)*N);
    fb.write((char*)&Jxz[0], sizeof (double)*N);
    fb.write((char*)&Jyy[0], sizeof (double)*N);
    fb.write((char*)&Jyz[0], sizeof (double)*N);
    fb.write((char*)&Jzz[0], sizeof (double)*N);

    vector<double> Jxx_inv(N), Jxy_inv(N), Jxz_inv(N), Jyy_inv(N), Jyz_inv(N), Jzz_inv(N);
    for(int i=0; i<N; i++) {
        Jxx_inv[i] = bodies[i].J_inv[0];
        Jxy_inv[i] = bodies[i].J_inv[1];
        Jxz_inv[i] = bodies[i].J_inv[2];
        Jyy_inv[i] = bodies[i].J_inv[3];
        Jyz_inv[i] = bodies[i].J_inv[4];
        Jzz_inv[i] = bodies[i].J_inv[5];
    }

    fb.write((char*)&Jxx_inv[0], sizeof (double)*N);
    fb.write((char*)&Jxy_inv[0], sizeof (double)*N);
    fb.write((char*)&Jxz_inv[0], sizeof (double)*N);
    fb.write((char*)&Jyy_inv[0], sizeof (double)*N);
    fb.write((char*)&Jyz_inv[0], sizeof (double)*N);
    fb.write((char*)&Jzz_inv[0], sizeof (double)*N);

    vector<double> Kx(N), Ky(N), Kz(N);

    for(int i=0; i<N; i++) {
        Kx[i] = bodies[i].K[0];
        Ky[i] = bodies[i].K[1];
        Kz[i] = bodies[i].K[2];
    }

    fb.write((char*)&Kx[0], sizeof (double)*N);
    fb.write((char*)&Ky[0], sizeof (double)*N);
    fb.write((char*)&Kz[0], sizeof (double)*N);

    // Output parameters
    fb.write((char*)&firstWrite, sizeof (bool));
    fb.write((char*)&firstWriteDiag, sizeof (bool));
    fb.write((char*)&file_counter, sizeof (int));
    fb.write((char*)&bin_number, sizeof (int));
        
    fb.close();
}
void Output::read_binary_backup(double &t, int &N, double &tcpu, double &dt_prev, int &num_integration_steps, vector<Body> &bodies, int &collision_flag, int &roche_flag, int &breakup_flag, double &dt_snapshot, double &dt0_log, double &fmul_log, int &num_snapshot) {
    string bin0_name = file_out + ".bin0";

    bool isComplete0;
    double t0, tcpu0, dt_prev0;
    int N0;
    vector<Body> bodies0;
    bool firstWrite0;
    bool firstWriteDiag0;
    int file_counter0;
    int bin_number0;
    int num_integration_steps0, num_snapshot0;
    int collision_flag0, roche_flag0, breakup_flag0;
    double dt_snapshot0, dt0_log0, fmul_log0;
    
    try {
        read_from_binary(t0, N0, tcpu0, dt_prev0, num_integration_steps0, bodies0, bin0_name, firstWrite0, firstWriteDiag0, file_counter0, bin_number0, collision_flag0, roche_flag0, breakup_flag0, dt_snapshot0, dt0_log0, fmul_log0, num_snapshot0); 
        isComplete0 = true;
    }
    catch(const std::exception&) {
        isComplete0 = false;    
    }

    string bin1_name = file_out + ".bin1";

    bool isComplete1;
    double t1, tcpu1, dt_prev1;
    int N1;
    vector<Body> bodies1;
    bool firstWrite1;
    bool firstWriteDiag1;
    int file_counter1;
    int bin_number1;
    int num_integration_steps1, num_snapshot1;
    int collision_flag1, roche_flag1, breakup_flag1;
    double dt_snapshot1, dt0_log1, fmul_log1;
    
    try {
        read_from_binary(t1, N1, tcpu1, dt_prev1, num_integration_steps1, bodies1, bin1_name, firstWrite1, firstWriteDiag1, file_counter1, bin_number1, collision_flag1, roche_flag1, breakup_flag1, dt_snapshot1, dt0_log1, fmul_log1, num_snapshot1); 
        isComplete1 = true;
    }
    catch(const std::exception&) {
        isComplete1 = false;    
    }

    if(isComplete0) {
        if(isComplete1) {
            if(t1 > t0) {
                t = t1;
                N = N1;
                tcpu = tcpu1;
                dt_prev = dt_prev1;
                bodies = bodies1;
                firstWrite = firstWrite1;
                firstWriteDiag = firstWriteDiag1;
                file_counter = file_counter1;
                bin_number = bin_number1;
                num_integration_steps = num_integration_steps1;
                collision_flag = collision_flag1;
                roche_flag = roche_flag1;
                breakup_flag = breakup_flag1;
                dt_snapshot = dt_snapshot1;
                dt0_log = dt0_log1;
                fmul_log = fmul_log1;
                num_snapshot = num_snapshot1;
            }
            else {
                t = t0;
                N = N0;
                tcpu = tcpu0;
                dt_prev = dt_prev0;
                bodies = bodies0;            
                firstWrite = firstWrite0;
                firstWriteDiag = firstWriteDiag0;
                file_counter = file_counter0;
                bin_number = bin_number0;
                num_integration_steps = num_integration_steps0;
                collision_flag = collision_flag0;
                roche_flag = roche_flag0;
                breakup_flag = breakup_flag0;
                dt_snapshot = dt_snapshot0;
                dt0_log = dt0_log0;
                fmul_log = fmul_log0;
                num_snapshot = num_snapshot0;
            }
        }
        else {
            t = t0;
            N = N0;
            tcpu = tcpu0;
            dt_prev = dt_prev0;
            bodies = bodies0;   
            firstWrite = firstWrite0;
            firstWriteDiag = firstWriteDiag0;
            file_counter = file_counter0;
            bin_number = bin_number0; 
            num_integration_steps = num_integration_steps0;                    
            collision_flag = collision_flag0;
            roche_flag = roche_flag0;
            breakup_flag = breakup_flag0;
            dt_snapshot = dt_snapshot0;
            dt0_log = dt0_log0;
            fmul_log = fmul_log0;           
            num_snapshot = num_snapshot0; 
        }
    }
    else {
        if(isComplete1) {
            t = t1;
            N = N1;
            tcpu = tcpu1;
            dt_prev = dt_prev1;
            bodies = bodies1;
            firstWrite = firstWrite1;
            firstWriteDiag = firstWriteDiag1;
            file_counter = file_counter1;
            bin_number = bin_number1;       
            num_integration_steps = num_integration_steps1;    
            collision_flag = collision_flag1;
            roche_flag = roche_flag1;
            breakup_flag = breakup_flag1;             
            dt_snapshot = dt_snapshot1;
            dt0_log = dt0_log1;
            fmul_log = fmul_log1;
            num_snapshot = num_snapshot1;
        }
        else {
            cerr << " " << endl;
            cerr << "Incomplete/invalid binary files: cannot continue simulation." << endl;
            cerr << " " << endl;            
            exit(1);
        }    
    }
}
void Output::read_from_binary(double &t, int &N, double &tcpu, double &dt_prev, int &num_integration_steps, vector<Body> &bodies, string &bin_name, bool &myfirstWrite, bool &myfirstWriteDiag, int &myfile_counter, int &mybin_number, int &collision_flag, int &roche_flag, int &breakup_flag, double &dt_snapshot, double &dt0_log, double &fmul_log, int &num_snapshot) { 
    fstream fb(bin_name, ios::in | ios::out | ios::binary);

    // header
    fb.read((char*)&t, sizeof (double));
    fb.read((char*)&N, sizeof (int));
    fb.read((char*)&tcpu, sizeof (double));
    fb.read((char*)&dt_prev, sizeof (double));
    fb.read((char*)&num_integration_steps, sizeof (int));

    fb.read((char*)&collision_flag, sizeof (int));
    fb.read((char*)&roche_flag, sizeof (int));
    fb.read((char*)&breakup_flag, sizeof (int));
    
    fb.read((char*)&dt_snapshot, sizeof (double));
    fb.read((char*)&dt0_log, sizeof (double));
    fb.read((char*)&fmul_log, sizeof (double));

    fb.read((char*)&num_snapshot, sizeof (int));
    
    bodies.resize(N);

    // id
    vector<int> id(N);
    fb.read((char*)&id[0], sizeof (int)*N);    

    for(int i=0; i<N; i++) bodies[i].id = id[i];

    // name
    for(int i=0; i<N; i++) {
        int name_len;        
        string name;
        fb.read((char*)&name_len, sizeof (int));
        name.resize(name_len);
        fb.read((char*)&name[0], sizeof (char)*name_len);
        bodies[i].name = name;
    }

    // internal props
    vector<double> m(N), R(N), xi(N);
    fb.read((char*)&m[0], sizeof (double)*N);    
    fb.read((char*)&R[0], sizeof (double)*N);    
    fb.read((char*)&xi[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].m = m[i];
        bodies[i].R    = R[i];
        bodies[i].xi   = xi[i];
    }    
    
    // Tidal parameters
    vector<double> kf(N), tau(N);
    fb.read((char*)&kf[0], sizeof (double)*N);    
    fb.read((char*)&tau[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].kf  = kf[i];
        bodies[i].tau = tau[i];
    }    
    
    // Other physical ingredients
    vector<double> a_mb(N);
    fb.read((char*)&a_mb[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].a_mb = a_mb[i];
    }    

    // Inertia tensors
    vector<double> Ixx(N), Ixy(N), Ixz(N), Iyy(N), Iyz(N), Izz(N);

    fb.read((char*)&Ixx[0], sizeof (double)*N);    
    fb.read((char*)&Ixy[0], sizeof (double)*N);    
    fb.read((char*)&Ixz[0], sizeof (double)*N);    
    fb.read((char*)&Iyy[0], sizeof (double)*N);    
    fb.read((char*)&Iyz[0], sizeof (double)*N);    
    fb.read((char*)&Izz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I[0] = Ixx[i];
        bodies[i].I[1] = Ixy[i];
        bodies[i].I[2] = Ixz[i];
        bodies[i].I[3] = Iyy[i];
        bodies[i].I[4] = Iyz[i];
        bodies[i].I[5] = Izz[i];
    }    

    vector<double> Ixx_p(N), Ixy_p(N), Ixz_p(N), Iyy_p(N), Iyz_p(N), Izz_p(N);

    fb.read((char*)&Ixx_p[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_p[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_p[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_p[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_p[0], sizeof (double)*N);    
    fb.read((char*)&Izz_p[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_p[0] = Ixx_p[i];
        bodies[i].I_p[1] = Ixy_p[i];
        bodies[i].I_p[2] = Ixz_p[i];
        bodies[i].I_p[3] = Iyy_p[i];
        bodies[i].I_p[4] = Iyz_p[i];
        bodies[i].I_p[5] = Izz_p[i];
    }    

    vector<double> Ixx_n(N), Ixy_n(N), Ixz_n(N), Iyy_n(N), Iyz_n(N), Izz_n(N);

    fb.read((char*)&Ixx_n[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_n[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_n[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_n[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_n[0], sizeof (double)*N);    
    fb.read((char*)&Izz_n[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_n[0] = Ixx_n[i];
        bodies[i].I_n[1] = Ixy_n[i];
        bodies[i].I_n[2] = Ixz_n[i];
        bodies[i].I_n[3] = Iyy_n[i];
        bodies[i].I_n[4] = Iyz_n[i];
        bodies[i].I_n[5] = Izz_n[i];
    }    

    vector<double> Ixx_inv(N), Ixy_inv(N), Ixz_inv(N), Iyy_inv(N), Iyz_inv(N), Izz_inv(N);

    fb.read((char*)&Ixx_inv[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_inv[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_inv[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_inv[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_inv[0], sizeof (double)*N);    
    fb.read((char*)&Izz_inv[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_inv[0] = Ixx_inv[i];
        bodies[i].I_inv[1] = Ixy_inv[i];
        bodies[i].I_inv[2] = Ixz_inv[i];
        bodies[i].I_inv[3] = Iyy_inv[i];
        bodies[i].I_inv[4] = Iyz_inv[i];
        bodies[i].I_inv[5] = Izz_inv[i];
    }    

    vector<double> Ixx_e_r(N), Ixy_e_r(N), Ixz_e_r(N), Iyy_e_r(N), Iyz_e_r(N), Izz_e_r(N);

    fb.read((char*)&Ixx_e_r[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_e_r[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_e_r[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_e_r[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_e_r[0], sizeof (double)*N);    
    fb.read((char*)&Izz_e_r[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_e_r[0] = Ixx_e_r[i];
        bodies[i].I_e_r[1] = Ixy_e_r[i];
        bodies[i].I_e_r[2] = Ixz_e_r[i];
        bodies[i].I_e_r[3] = Iyy_e_r[i];
        bodies[i].I_e_r[4] = Iyz_e_r[i];
        bodies[i].I_e_r[5] = Izz_e_r[i];
    }    

    vector<double> dIxx_e_r(N), dIxy_e_r(N), dIxz_e_r(N), dIyy_e_r(N), dIyz_e_r(N), dIzz_e_r(N);

    fb.read((char*)&dIxx_e_r[0], sizeof (double)*N);    
    fb.read((char*)&dIxy_e_r[0], sizeof (double)*N);    
    fb.read((char*)&dIxz_e_r[0], sizeof (double)*N);    
    fb.read((char*)&dIyy_e_r[0], sizeof (double)*N);    
    fb.read((char*)&dIyz_e_r[0], sizeof (double)*N);    
    fb.read((char*)&dIzz_e_r[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].dI_e_r[0] = dIxx_e_r[i];
        bodies[i].dI_e_r[1] = dIxy_e_r[i];
        bodies[i].dI_e_r[2] = dIxz_e_r[i];
        bodies[i].dI_e_r[3] = dIyy_e_r[i];
        bodies[i].dI_e_r[4] = dIyz_e_r[i];
        bodies[i].dI_e_r[5] = dIzz_e_r[i];
    }    

    vector<double> Ixx_e_w(N), Ixy_e_w(N), Ixz_e_w(N), Iyy_e_w(N), Iyz_e_w(N), Izz_e_w(N);

    fb.read((char*)&Ixx_e_w[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_e_w[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_e_w[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_e_w[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_e_w[0], sizeof (double)*N);    
    fb.read((char*)&Izz_e_w[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_e_w[0] = Ixx_e_w[i];
        bodies[i].I_e_w[1] = Ixy_e_w[i];
        bodies[i].I_e_w[2] = Ixz_e_w[i];
        bodies[i].I_e_w[3] = Iyy_e_w[i];
        bodies[i].I_e_w[4] = Iyz_e_w[i];
        bodies[i].I_e_w[5] = Izz_e_w[i];
    }    

    vector<double> Ixx_e(N), Ixy_e(N), Ixz_e(N), Iyy_e(N), Iyz_e(N), Izz_e(N);

    fb.read((char*)&Ixx_e[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_e[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_e[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_e[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_e[0], sizeof (double)*N);    
    fb.read((char*)&Izz_e[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_e[0] = Ixx_e[i];
        bodies[i].I_e[1] = Ixy_e[i];
        bodies[i].I_e[2] = Ixz_e[i];
        bodies[i].I_e[3] = Iyy_e[i];
        bodies[i].I_e[4] = Iyz_e[i];
        bodies[i].I_e[5] = Izz_e[i];
    }    

    vector<double> dIxx_e(N), dIxy_e(N), dIxz_e(N), dIyy_e(N), dIyz_e(N), dIzz_e(N);

    fb.read((char*)&dIxx_e[0], sizeof (double)*N);    
    fb.read((char*)&dIxy_e[0], sizeof (double)*N);    
    fb.read((char*)&dIxz_e[0], sizeof (double)*N);    
    fb.read((char*)&dIyy_e[0], sizeof (double)*N);    
    fb.read((char*)&dIyz_e[0], sizeof (double)*N);    
    fb.read((char*)&dIzz_e[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].dI_e[0] = dIxx_e[i];
        bodies[i].dI_e[1] = dIxy_e[i];
        bodies[i].dI_e[2] = dIxz_e[i];
        bodies[i].dI_e[3] = dIyy_e[i];
        bodies[i].dI_e[4] = dIyz_e[i];
        bodies[i].dI_e[5] = dIzz_e[i];
    }    

    vector<double> dIxx_n(N), dIxy_n(N), dIxz_n(N), dIyy_n(N), dIyz_n(N), dIzz_n(N);

    fb.read((char*)&dIxx_n[0], sizeof (double)*N);    
    fb.read((char*)&dIxy_n[0], sizeof (double)*N);    
    fb.read((char*)&dIxz_n[0], sizeof (double)*N);    
    fb.read((char*)&dIyy_n[0], sizeof (double)*N);    
    fb.read((char*)&dIyz_n[0], sizeof (double)*N);    
    fb.read((char*)&dIzz_n[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].dI_n[0] = dIxx_n[i];
        bodies[i].dI_n[1] = dIxy_n[i];
        bodies[i].dI_n[2] = dIxz_n[i];
        bodies[i].dI_n[3] = dIyy_n[i];
        bodies[i].dI_n[4] = dIyz_n[i];
        bodies[i].dI_n[5] = dIzz_n[i];
    }    

    vector<double> Ixx_e_prev(N), Ixy_e_prev(N), Ixz_e_prev(N), Iyy_e_prev(N), Iyz_e_prev(N), Izz_e_prev(N);

    fb.read((char*)&Ixx_e_prev[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_e_prev[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_e_prev[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_e_prev[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_e_prev[0], sizeof (double)*N);    
    fb.read((char*)&Izz_e_prev[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_e_prev[0] = Ixx_e_prev[i];
        bodies[i].I_e_prev[1] = Ixy_e_prev[i];
        bodies[i].I_e_prev[2] = Ixz_e_prev[i];
        bodies[i].I_e_prev[3] = Iyy_e_prev[i];
        bodies[i].I_e_prev[4] = Iyz_e_prev[i];
        bodies[i].I_e_prev[5] = Izz_e_prev[i];
    }    

    vector<double> Ixx_n_prev(N), Ixy_n_prev(N), Ixz_n_prev(N), Iyy_n_prev(N), Iyz_n_prev(N), Izz_n_prev(N);

    fb.read((char*)&Ixx_n_prev[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_n_prev[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_n_prev[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_n_prev[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_n_prev[0], sizeof (double)*N);    
    fb.read((char*)&Izz_n_prev[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_n_prev[0] = Ixx_n_prev[i];
        bodies[i].I_n_prev[1] = Ixy_n_prev[i];
        bodies[i].I_n_prev[2] = Ixz_n_prev[i];
        bodies[i].I_n_prev[3] = Iyy_n_prev[i];
        bodies[i].I_n_prev[4] = Iyz_n_prev[i];
        bodies[i].I_n_prev[5] = Izz_n_prev[i];
    }    

    vector<double> Ixx_e_prev_bu(N), Ixy_e_prev_bu(N), Ixz_e_prev_bu(N), Iyy_e_prev_bu(N), Iyz_e_prev_bu(N), Izz_e_prev_bu(N);

    fb.read((char*)&Ixx_e_prev_bu[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_e_prev_bu[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_e_prev_bu[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_e_prev_bu[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_e_prev_bu[0], sizeof (double)*N);    
    fb.read((char*)&Izz_e_prev_bu[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_e_prev_bu[0] = Ixx_e_prev_bu[i];
        bodies[i].I_e_prev_bu[1] = Ixy_e_prev_bu[i];
        bodies[i].I_e_prev_bu[2] = Ixz_e_prev_bu[i];
        bodies[i].I_e_prev_bu[3] = Iyy_e_prev_bu[i];
        bodies[i].I_e_prev_bu[4] = Iyz_e_prev_bu[i];
        bodies[i].I_e_prev_bu[5] = Izz_e_prev_bu[i];
    }    

    vector<double> Ixx_e_rh(N), Ixy_e_rh(N), Ixz_e_rh(N), Iyy_e_rh(N), Iyz_e_rh(N), Izz_e_rh(N);

    fb.read((char*)&Ixx_e_rh[0], sizeof (double)*N);    
    fb.read((char*)&Ixy_e_rh[0], sizeof (double)*N);    
    fb.read((char*)&Ixz_e_rh[0], sizeof (double)*N);    
    fb.read((char*)&Iyy_e_rh[0], sizeof (double)*N);    
    fb.read((char*)&Iyz_e_rh[0], sizeof (double)*N);    
    fb.read((char*)&Izz_e_rh[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].I_e_rh[0] = Ixx_e_rh[i];
        bodies[i].I_e_rh[1] = Ixy_e_rh[i];
        bodies[i].I_e_rh[2] = Ixz_e_rh[i];
        bodies[i].I_e_rh[3] = Iyy_e_rh[i];
        bodies[i].I_e_rh[4] = Iyz_e_rh[i];
        bodies[i].I_e_rh[5] = Izz_e_rh[i];
    }    

    // Spin and angular momenta
    vector<double> wx(N), wy(N), wz(N);
    vector<double> Lx(N), Ly(N), Lz(N);

    fb.read((char*)&wx[0], sizeof (double)*N);    
    fb.read((char*)&wy[0], sizeof (double)*N);    
    fb.read((char*)&wz[0], sizeof (double)*N);    
    fb.read((char*)&Lx[0], sizeof (double)*N);    
    fb.read((char*)&Ly[0], sizeof (double)*N);    
    fb.read((char*)&Lz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].w[0] = wx[i];
        bodies[i].w[1] = wy[i];
        bodies[i].w[2] = wz[i];
        bodies[i].L[0] = Lx[i];
        bodies[i].L[1] = Ly[i];
        bodies[i].L[2] = Lz[i];
    }    

    vector<double> Tx(N), Ty(N), Tz(N);

    fb.read((char*)&Tx[0], sizeof (double)*N);    
    fb.read((char*)&Ty[0], sizeof (double)*N);    
    fb.read((char*)&Tz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].T[0] = Tx[i];
        bodies[i].T[1] = Ty[i];
        bodies[i].T[2] = Tz[i];
    }    
    
    // pos and vel
    vector<double> x(N), y(N), z(N);
    vector<double> vx(N), vy(N), vz(N);

    fb.read((char*)&x[0], sizeof (double)*N);    
    fb.read((char*)&y[0], sizeof (double)*N);    
    fb.read((char*)&z[0], sizeof (double)*N);    
    fb.read((char*)&vx[0], sizeof (double)*N);    
    fb.read((char*)&vy[0], sizeof (double)*N);    
    fb.read((char*)&vz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].r[0] = x[i];
        bodies[i].r[1] = y[i];
        bodies[i].r[2] = z[i];
        bodies[i].v[0] = vx[i];
        bodies[i].v[1] = vy[i];
        bodies[i].v[2] = vz[i];
    }    

    vector<double> ax(N), ay(N), az(N);

    fb.read((char*)&ax[0], sizeof (double)*N);    
    fb.read((char*)&ay[0], sizeof (double)*N);    
    fb.read((char*)&az[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].a[0] = ax[i];
        bodies[i].a[1] = ay[i];
        bodies[i].a[2] = az[i];
    }    
    
    // Aux quantities
    vector<double> R5(N), R5_3(N), kf_R5(N), kf_R5_3(N), tau_inv(N);
    fb.read((char*)&R5[0], sizeof (double)*N);    
    fb.read((char*)&R5_3[0], sizeof (double)*N);    
    fb.read((char*)&kf_R5[0], sizeof (double)*N);    
    fb.read((char*)&kf_R5_3[0], sizeof (double)*N);    
    fb.read((char*)&tau_inv[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].R5  = R5[i];
        bodies[i].R5_3  = R5_3[i];
        bodies[i].kf_R5 = kf_R5[i];
        bodies[i].kf_R5_3 = kf_R5_3[i];
        bodies[i].tau_inv = tau_inv[i];
    }    

    vector<double> vvx(N), vvy(N), vvz(N);

    fb.read((char*)&vvx[0], sizeof (double)*N);    
    fb.read((char*)&vvy[0], sizeof (double)*N);    
    fb.read((char*)&vvz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].vv[0] = vvx[i];
        bodies[i].vv[1] = vvy[i];
        bodies[i].vv[2] = vvz[i];
    }    
    
    // Auxiliary variables
    vector<double> Jxx_n(N), Jxy_n(N), Jxz_n(N), Jyy_n(N), Jyz_n(N), Jzz_n(N);

    fb.read((char*)&Jxx_n[0], sizeof (double)*N);    
    fb.read((char*)&Jxy_n[0], sizeof (double)*N);    
    fb.read((char*)&Jxz_n[0], sizeof (double)*N);    
    fb.read((char*)&Jyy_n[0], sizeof (double)*N);    
    fb.read((char*)&Jyz_n[0], sizeof (double)*N);    
    fb.read((char*)&Jzz_n[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].J_n[0] = Jxx_n[i];
        bodies[i].J_n[1] = Jxy_n[i];
        bodies[i].J_n[2] = Jxz_n[i];
        bodies[i].J_n[3] = Jyy_n[i];
        bodies[i].J_n[4] = Jyz_n[i];
        bodies[i].J_n[5] = Jzz_n[i];
    }    

    vector<double> Jxx(N), Jxy(N), Jxz(N), Jyy(N), Jyz(N), Jzz(N);

    fb.read((char*)&Jxx[0], sizeof (double)*N);    
    fb.read((char*)&Jxy[0], sizeof (double)*N);    
    fb.read((char*)&Jxz[0], sizeof (double)*N);    
    fb.read((char*)&Jyy[0], sizeof (double)*N);    
    fb.read((char*)&Jyz[0], sizeof (double)*N);    
    fb.read((char*)&Jzz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].J[0] = Jxx[i];
        bodies[i].J[1] = Jxy[i];
        bodies[i].J[2] = Jxz[i];
        bodies[i].J[3] = Jyy[i];
        bodies[i].J[4] = Jyz[i];
        bodies[i].J[5] = Jzz[i];
    }    

    vector<double> Jxx_inv(N), Jxy_inv(N), Jxz_inv(N), Jyy_inv(N), Jyz_inv(N), Jzz_inv(N);

    fb.read((char*)&Jxx_inv[0], sizeof (double)*N);    
    fb.read((char*)&Jxy_inv[0], sizeof (double)*N);    
    fb.read((char*)&Jxz_inv[0], sizeof (double)*N);    
    fb.read((char*)&Jyy_inv[0], sizeof (double)*N);    
    fb.read((char*)&Jyz_inv[0], sizeof (double)*N);    
    fb.read((char*)&Jzz_inv[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].J_inv[0] = Jxx_inv[i];
        bodies[i].J_inv[1] = Jxy_inv[i];
        bodies[i].J_inv[2] = Jxz_inv[i];
        bodies[i].J_inv[3] = Jyy_inv[i];
        bodies[i].J_inv[4] = Jyz_inv[i];
        bodies[i].J_inv[5] = Jzz_inv[i];
    }    

    vector<double> Kx(N), Ky(N), Kz(N);

    fb.read((char*)&Kx[0], sizeof (double)*N);    
    fb.read((char*)&Ky[0], sizeof (double)*N);    
    fb.read((char*)&Kz[0], sizeof (double)*N);    

    for(int i=0; i<N; i++) {
        bodies[i].K[0] = Kx[i];
        bodies[i].K[1] = Ky[i];
        bodies[i].K[2] = Kz[i];
    }    
    
    // Output parameters
    fb.read((char*)&myfirstWrite, sizeof (bool));
    fb.read((char*)&myfirstWriteDiag, sizeof (bool));
    fb.read((char*)&myfile_counter, sizeof (int));
    fb.read((char*)&mybin_number, sizeof (int));
        
    fb.close();
}

// Diagnostics log
void Output::write_log(int argc, char* argv[], double t, double tcpu, int &num_integration_steps, int N0, int N1, double dr, double dv, double dLx, double dLy, double dLz, double dLx_rel, double dLy_rel, double dLz_rel, double dL_abs, double dL_rel, double dE_abs, double dE_rel, double x0, int outcome_type, int collision_flag, int roche_flag, int breakup_flag) {
    string outputFile = file_out + ".log";

    create_output_file(outputFile, true);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << endl;
        cerr << "Cannot open output file: " << outputFile << endl;
        cerr << endl;
        exit(1);
    }

    writer << endl;
    writer << "Diagnostic summary: " << endl;
    
    writer << "Command        = ";
    for(int i=0; i<argc; i++) writer << argv[i] << " "; 
    writer << endl;
    
    writer << "T              = " << units.convert_time_from_code_to_output(t, time_unit_output) << " " << time_unit_output << endl;
    writer << "t_cpu          = " << tcpu << " [s]" << endl;
    writer << "Time steps     = " << num_integration_steps << endl;

    writer << "dr_cm_abs      = " << units.convert_length_from_code_to_output(dr, length_unit_output) << " " << length_unit_output << endl;
    writer << "dv_cm_abs      = " << units.convert_speed_from_code_to_output(dv, speed_unit_output) << " " << speed_unit_output << endl;

    writer << "dLx_abs        = " << units.convert_angular_momentum_from_code_to_output(dLx, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    writer << "dLy_abs        = " << units.convert_angular_momentum_from_code_to_output(dLy, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    writer << "dLz_abs        = " << units.convert_angular_momentum_from_code_to_output(dLz, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    
    writer << "dLx_rel        = " << dLx_rel << endl;
    writer << "dLy_rel        = " << dLy_rel << endl;
    writer << "dLz_rel        = " << dLz_rel << endl;

    writer << "dL_abs         = " << units.convert_angular_momentum_from_code_to_output(dL_abs, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    writer << "dL_rel         = " << dL_rel << endl;

    writer << "dE_abs         = " << units.convert_energy_from_code_to_output(dE_abs, energy_unit_output) << " " << energy_unit_output << endl;    
    writer << "dE_rel         = " << dE_rel << endl;

    writer << "N init         = " << N0 << endl;
    writer << "N final        = " << N1 << endl;

    writer << "Outcome        = " << outcome_type << endl;
    writer << "collision flag = " << collision_flag << endl;
    writer << "roche flag     = " << roche_flag << endl;
    writer << "breakup flag   = " << breakup_flag << endl;

    writer << "x0             = " << x0 << endl;
    writer << endl;

    writer.close();
}
void Output::print_log(int argc, char* argv[], double t, double tcpu, int &num_integration_steps, int N0, int N1, double dr, double dv, double dLx, double dLy, double dLz, double dLx_rel, double dLy_rel, double dLz_rel, double dL_abs, double dL_rel, double dE_abs, double dE_rel, double x0, int outcome_type, int collision_flag, int roche_flag, int breakup_flag) {
    cout << endl;
    cout << "Diagnostic summary: " << endl;
    
    cout << "Command        = ";
    for(int i=0; i<argc; i++) cout << argv[i] << " "; 
    cout << endl;
    
    cout << "T              = " << units.convert_time_from_code_to_output(t, time_unit_output) << " " << time_unit_output << endl;
    cout << "t_cpu          = " << tcpu << " [s]" << endl;
    cout << "Time steps     = " << num_integration_steps << endl;

    cout << "dr_cm_abs      = " << units.convert_length_from_code_to_output(dr, length_unit_output) << " " << length_unit_output << endl;
    cout << "dv_cm_abs      = " << units.convert_speed_from_code_to_output(dv, speed_unit_output) << " " << speed_unit_output << endl;

    cout << "dLx_abs        = " << units.convert_angular_momentum_from_code_to_output(dLx, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    cout << "dLy_abs        = " << units.convert_angular_momentum_from_code_to_output(dLy, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    cout << "dLz_abs        = " << units.convert_angular_momentum_from_code_to_output(dLz, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    
    cout << "dLx_rel        = " << dLx_rel << endl;
    cout << "dLy_rel        = " << dLy_rel << endl;
    cout << "dLz_rel        = " << dLz_rel << endl;

    cout << "dL_abs         = " << units.convert_angular_momentum_from_code_to_output(dL_abs, angular_momentum_unit_output) << " " << angular_momentum_unit_output << endl;
    cout << "dL_rel         = " << dL_rel << endl;

    cout << "dE_abs         = " << units.convert_energy_from_code_to_output(dE_abs, energy_unit_output) << " " << energy_unit_output << endl;    
    cout << "dE_rel         = " << dE_rel << endl;

    cout << "N init         = " << N0 << endl;
    cout << "N final        = " << N1 << endl;
 
    cout << "Outcome        = " << outcome_type << endl;
    cout << "collision flag = " << collision_flag << endl;
    cout << "roche flag     = " << roche_flag << endl;
    cout << "breakup flag   = " << breakup_flag << endl;

    cout << "x0             = " << x0 << endl;
}

void Output::write_diag(double &t, double &t_cpu, int &num_integration_steps, int &N, array<double, 3> &r0, array<double, 3> &v0, array<double, 3> &L0_orb, array<double, 3> &L0_spin, double &Ekin0_orb, double Epot0, double &Ekin0_spin) {
    string outputFile = file_out + ".diag";

    create_output_file(outputFile, firstWriteDiag);

    ofstream writer;
    writer.open(outputFile, ios::app);

    if (!writer) {
        cerr << endl;
        cerr << "Cannot open output file: " << outputFile << endl;
        cerr << endl;
        exit(1);
    }

    if(firstWriteDiag) {
        string str = "# t" + time_unit_output;
        writer << setw(n_width) << str;
        str = "t_cpu[s] ";
        writer << setw(n_width) << str;
        str = "num_step ";
        writer << setw(n_width) << str;
        str = "N ";
        writer << setw(n_width) << str;
        str = "xcm" + length_unit_output + " ";
        writer << setw(n_width) << str;
        str = "ycm" + length_unit_output + " ";
        writer << setw(n_width) << str;
        str = "zcm" + length_unit_output + " ";
        writer << setw(n_width) << str;
        str = "vxcm" + speed_unit_output + " ";
        writer << setw(n_width) << str;
        str = "vycm" + speed_unit_output + " ";
        writer << setw(n_width) << str;
        str = "vzcm" + speed_unit_output + " ";
        writer << setw(n_width) << str;
        str = "Lx_orb" + angular_momentum_unit_output + " ";
        writer << setw(n_width) << str;
        str = "Ly_orb" + angular_momentum_unit_output + " ";
        writer << setw(n_width) << str;
        str = "Lz_orb" + angular_momentum_unit_output + " ";
        writer << setw(n_width) << str;
        str = "Lx_spin" + angular_momentum_unit_output + " ";
        writer << setw(n_width) << str;
        str = "Ly_spin" + angular_momentum_unit_output + " ";
        writer << setw(n_width) << str;
        str = "Lz_spin" + angular_momentum_unit_output + " ";
        writer << setw(n_width) << str;
        str = "E_kin_orb" + energy_unit_output + " ";
        writer << setw(n_width) << str;
        str = "E_pot" + energy_unit_output + " ";
        writer << setw(n_width) << str;
        str = "E_spin" + energy_unit_output;
        writer << setw(n_width) << str << endl;
    }

    writer << setprecision(16);
    writer << scientific;

    writer << setw(n_width) << units.convert_time_from_code_to_output(t, time_unit_output);
    writer << setw(n_width) << t_cpu;
    writer << setw(n_width) << num_integration_steps;
    writer << setw(n_width) << N;
    for(int i=0; i<3; i++) writer << setw(n_width) << units.convert_length_from_code_to_output(r0[i], length_unit_output);
    for(int i=0; i<3; i++) writer << setw(n_width) << units.convert_speed_from_code_to_output(v0[i], speed_unit_output);
    for(int i=0; i<3; i++) writer << setw(n_width) << units.convert_angular_momentum_from_code_to_output(L0_orb[i], angular_momentum_unit_output);
    for(int i=0; i<3; i++) writer << setw(n_width) << units.convert_angular_momentum_from_code_to_output(L0_spin[i], angular_momentum_unit_output);
    writer << setw(n_width) << units.convert_energy_from_code_to_output(Ekin0_orb, energy_unit_output);
    writer << setw(n_width) << units.convert_energy_from_code_to_output(Epot0, energy_unit_output);
    writer << setw(n_width) << units.convert_energy_from_code_to_output(Ekin0_spin, energy_unit_output);
    writer << endl;     

    writer.close();
    
    if(firstWriteDiag) firstWriteDiag = false;
} 
    
void Output::backup_file(string f1) {
    string f2 = "";

    int N_str = f1.length();
    for(int k=0; k<N_str; k++) {
        f2 = f1.substr(N_str-1-k, 1+k);            

        if(f1.at(N_str-1-k) == '/') {
            f2 = f1.substr(N_str-1-k+1, 1+k-1);            
            break;
        }
    }
    
    f2 = dir + f2;

    copy_file(f1, f2);
}

void Output::copy_file(string f1, string f2) {
    create_output_file(f2, true);

    ifstream str1;
    str1.open(f1.c_str());
    if(!str1) {
        cerr << "Cannot open " << f1 << "!" << endl;
        exit(1);
    }
        
    ofstream str2;
    str2.open(f2.c_str(), ios::app);    
    if(!str2) {
        cerr << "Cannot open " << f2 << "!" << endl;
        exit(1);
    }

    string line;
    getline(str1, line);

    while(!str1.eof()) {
        str2 << line << endl;
        getline(str1, line);
    }

    str1.close();
    str2.close();
}



    
    
