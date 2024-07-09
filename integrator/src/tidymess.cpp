#include "Timer.h"
#include "Banner.h"

#include "Initializer.h"
#include "Output.h"

#include "Tidy.h"

#include "Collision.h"
#include "Breakup.h"

int main(int argc, char* argv[]) {

    // Parse user arguments
    Initializer init;    

    init.parse_arguments(argc, argv);

    // Set up tidymess
    Tidy tidymess;

    tidymess.set_model_time(init.t_begin);
    tidymess.set_particles(init.name, init.id, init.data_ic);

    tidymess.set_dt_mode(init.dt_mode);
    tidymess.set_dt_const(init.dt_const);
    tidymess.set_eta(init.eta);

    tidymess.set_tidal_model(init.tidal_model);

    tidymess.set_pn_order(init.pn_order);
    tidymess.set_speed_of_light(init.speed_of_light);

    tidymess.set_collision_mode(init.collisions);
    tidymess.set_roche_mode(init.roche_limit);

    tidymess.set_magnetic_braking(init.magnetic_braking);

    tidymess.set_n_iter(init.n_iter);

    // Set initial shapes and angular momenta
    if(tidymess.get_tidal_model() > 0) {
        switch(init.initial_shape) {
            case 0:
                tidymess.set_to_spherical_shape();
                break;
            case 1:
                tidymess.set_to_equilibrium_shape();
                break;
        }

        tidymess.update_angular_momentum();    
    }
    
    tidymess.commit_parameters();
                            
    // Set up collision detection
    Collision collision;
    
    collision.set_collision_mode(init.collisions); 
    collision.set_roche_mode(init.roche_limit);     
    collision.setup();
    
    // Set up rotational breakup detection
    Breakup breakup;
    
    breakup.set_breakup_mode(init.breakup_speed);    
    breakup.setup();
    
    // Set up output
    Output output;

    output.set_dir(init.output_dir);
    output.set_overwrite(init.overwrite);
    output.set_file_out(init.file_out);        
    
    output.set_format(init.output_format);
    output.set_info(init.output_info);
    output.set_coor(init.output_coor);

    output.set_units(init.physical_units, init.mass_unit, init.length_unit, init.time_unit, init.speed_unit);
    output.set_conversion_factors(init.get_Cm(), init.get_Cr(), init.get_Cv(), init.get_Ct());

    output.set_tidal_model(init.tidal_model);
    output.determine_output_mode();        
            
    // Set up simulation and diagnostics
    double t           = tidymess.get_model_time();
    double t_end       = init.t_end;

    int snapshot_mode  = init.snapshot_mode;
    int N_snapshot     = init.n_snapshot;
    double t1_log      = init.t1_code;
    
    vector<Body> bodies = tidymess.get_particles();
    int N_init = bodies.size();

    //---------------------------------------------------------------------

    bool dt_pos;
    int dt_sgn;
    if(t_end > t) {
        dt_pos = true;
        dt_sgn = 1;
    }        
    else {
        dt_pos = false;
        dt_sgn = -1;
    }

    tidymess.set_dt_sgn(dt_sgn);

    int num_snapshot = 0;
    int num_integration_step = 0;

    double t_begin = t;
    double dt_snapshot = 0;

    double dt0_log = 0;
    double fmul_log = 0;

    if(snapshot_mode == 0) { // Constant time interval
        double t_sim = t_end - t_begin;
        dt_snapshot = t_sim / N_snapshot;        
    }
    else if(snapshot_mode == 1) { // Logarithmic time interval
        // Assume t_end > t_begin >= 0
        if(t_begin == 0) {
            double tmin = t1_log;
            double tmax = t_end;
            double logtmin = log(tmin);
            double logtmax = log(tmax);
               
            double dlogt = (logtmax-logtmin)/N_snapshot;
                
            //logtmax = log(t_end);
            //logtmin = logtmax - (N_snapshot+1)*dlogt;
            logtmin = logtmin - dlogt;
                
            dt0_log = logtmin;
            fmul_log = dlogt;

        }
        else {
            double tmin = t_begin;
            double tmax = t_end;
            double logtmin = log(tmin);
            double logtmax = log(tmax);
                
            double dlogt = (logtmax-logtmin)/N_snapshot;
                                
            dt0_log = logtmin;
            fmul_log = dlogt;

        }
    }
    
    //---------------------------------------------------------------------

    array<double, 3> r0 = tidymess.get_center_of_mass();
    array<double, 3> v0 = tidymess.get_center_of_mass_velocity();
    array<double, 3> L0_orb = tidymess.get_orbital_angular_momentum();
    array<double, 3> L0_spin = tidymess.get_spin_angular_momentum();
    double Ekin0_orb = tidymess.get_orbital_kinetic_energy();
    double Ekin0_spin = tidymess.get_spin_kinetic_energy();
    double Epot0 = tidymess.get_potential_energy();

    //---------------------------------------------------------------------

    int N;
    double tcpu_offset, t_cpu;

    double dt_prev = tidymess.get_dt_prev();

    int collision_flag = 0;
    int roche_flag = 0;
    int breakup_flag = 0;

    if(init.to_continue) {
        double t_bin;
        int N_bin;
        double tcpu_bin;
        double dt_prev_bin, t_end_bin;
        vector<Body> bodies_bin;
        int num_snapshot_bin;
        int num_integration_step_bin;
        int collision_flag_bin, roche_flag_bin, breakup_flag_bin;
        double dt_snapshot_bin, dt0_log_bin, fmul_log_bin;

        output.read_binary_backup(t_bin, N_bin, tcpu_bin, dt_prev_bin, t_end_bin, num_integration_step_bin, bodies_bin, collision_flag_bin, roche_flag_bin, breakup_flag_bin, dt_snapshot_bin, dt0_log_bin, fmul_log_bin, num_snapshot_bin);

        t = t_bin;
        N = N_bin;
        tcpu_offset = tcpu_bin;
        dt_prev = dt_prev_bin;
        num_integration_step = num_integration_step_bin;
        collision_flag = collision_flag_bin;
        roche_flag = roche_flag_bin;
        breakup_flag = breakup_flag_bin;
        
        bodies = bodies_bin;
        //for(int i=0; i<N; i++) {
        //    bodies[i].update_aux_properties();
        //}        
        
        tidymess.set_model_time(t);
        tidymess.set_particles(bodies);
        
        tidymess.set_dt_prev(dt_prev);
        tidymess.set_num_integration_step(num_integration_step);

        num_snapshot = 0;        
        t_begin = t; //init.t_begin; // t;
        dt0_log = log(t_begin);

        if(t_end == t_end_bin) { // Finish a simulation
            dt_snapshot = dt_snapshot_bin;
            fmul_log = fmul_log_bin;
        }
        else { // Extend a simulation
            //double t_sim = t_end - t_begin;
            //dt_snapshot = t_sim / N_snapshot;        

            dt_snapshot = dt_snapshot_bin;

            //double tmin = t_begin; //t;
            //double tmax = t_end;
            //double logtmin = log(tmin);
            //double logtmax = log(tmax);
                
            //double dlogt = (logtmax-logtmin)/N_snapshot;                  
            //fmul_log = dlogt;

            fmul_log = fmul_log_bin;
        }
    }
    else {
        N = bodies.size();
        
        tcpu_offset = 0;
        t_cpu = tcpu_offset;

        output.write_snapshot(t, t_cpu, bodies);

        output.backup_file(init.file_par);
        output.backup_file(init.file_ic);

        output.save_to_binary(t, N, t_cpu, dt_prev, t_end, num_integration_step, collision_flag, roche_flag, breakup_flag, bodies, dt_snapshot, dt0_log, fmul_log, num_snapshot);

        if(init.output_diag) {
            output.write_diag(t, t_cpu, num_integration_step, N, r0, v0, L0_orb, L0_spin, Ekin0_orb, Epot0, Ekin0_spin); 
        }
    }
                
    if(t == 0 && init.snapshot_mode == 1) N_snapshot++;            
                
    //---------------------------------------------------------------------
    
    Banner banner;

    if(init.output_terminal) {
        banner.print_stars();
        banner.print_banner();
        banner.print_stars();
    
        banner.print_intro(init.to_continue);
    }

    //---------------------------------------------------------------------

    double max_cpu_time = init.max_cpu_time;
    bool check_max_cpu = (max_cpu_time == 0) ? false : true;

    // 0 = end time reached
    // 1 = maximum cpu time reached
    // 2 = halted due to collision
    // 3 = halted due to roche limit
    // 4 = single body left in the system
    // 5 = body reached breakup speed
    int stopping_condition_type = -1;

    Timer timer;
    timer.start();

    // Time integration
    while(num_snapshot < N_snapshot) {

        num_snapshot++;

        // Integration step
        if(snapshot_mode == 0) {
            t = t_begin + num_snapshot * dt_snapshot;
            /*
            if(dt_pos) {
                if(t > t_end) {
                    t = t_end;
                    num_snapshot--;
                }
            }
            else {
                if(t < t_end) {
                    t = t_end;
                    num_snapshot--;
                }
            }   
            */            
            tidymess.evolve_model(t);
        }
        else if(snapshot_mode == 1) {
            //double logt = dt0_log + num_snapshot * fmul_log;    
            dt0_log += fmul_log;
            
            double logt = dt0_log;
            t = exp(logt);    
                        
            tidymess.evolve_model(t);
        }
        else {
            tidymess.evolve_model(N_snapshot);

            t = tidymess.get_model_time();
        }

        t = tidymess.get_model_time();
        bodies = tidymess.get_particles();

        // Store a snapshot
        t_cpu = tcpu_offset + timer.read();        
        output.write_snapshot(t, t_cpu, bodies);

        // Save snapshot to binary file
        N = bodies.size();
        dt_prev = tidymess.get_dt_prev();
        num_integration_step = tidymess.get_num_integration_step();

        output.save_to_binary(t, N, t_cpu, dt_prev, t_end, num_integration_step, collision_flag, roche_flag, breakup_flag, bodies, dt_snapshot, dt0_log, fmul_log, num_snapshot);

        // Store diagnostics
        if(init.output_diag) {
            array<double, 3> r1 = tidymess.get_center_of_mass();
            array<double, 3> v1 = tidymess.get_center_of_mass_velocity();
            array<double, 3> L1_orb = tidymess.get_orbital_angular_momentum();
            array<double, 3> L1_spin = tidymess.get_spin_angular_momentum();
            double Ekin1_orb = tidymess.get_orbital_kinetic_energy();
            double Ekin1_spin = tidymess.get_spin_kinetic_energy();
            double Epot1 = tidymess.get_potential_energy();
        
            output.write_diag(t, t_cpu, num_integration_step, N, r1, v1, L1_orb, L1_spin, Ekin1_orb, Epot1, Ekin1_spin); 
        }
                
        if(init.output_terminal) {
            cout << "Simulation progress:\t";
            cout << scientific<< setw(12) << setprecision(6) << 100*t/t_end << " %";
            cout << "\tCPU time:\t";
            cout << scientific<< setw(12) << setprecision(6) << t_cpu << " [s]";
            cout << "\tETA:\t";
            cout << scientific << setw(12) << setprecision(6) << (t_end-t) * t_cpu/(t-init.t_begin) << " [s]" << endl;
        }

        // Stopping condition end time
        if(dt_pos) {
            if(t >= t_end) {
                stopping_condition_type = 0;
                break;
            }
        }
        else {
            if(t <= t_end) {
                stopping_condition_type = 0;
                break;
            }        
        }
        
        // Collision handling
        if(collision.to_detect_collision) {
            if(tidymess.is_collision_detected() == true) {
                collision_flag = tidymess.get_collision_flag();        
                
                if(collision.collision_mode == 1) {
                    if(init.output_terminal) cout << "Collision detected. Ignoring it. " << endl;
                }
                else if(collision.collision_mode == 2) { // exception = End program
                    stopping_condition_type = 2;
            
                    if(init.output_terminal) cout << "Collision detected. Ending simulation." << endl;
                    break;
                }
                else if(collision.collision_mode == 3) { // sticky sphere approximation
                    if(init.output_terminal) cout << "Collision detected. Replace by center of mass particle." << endl;

                    vector< array<int, 2> > cindex = tidymess.get_collision_indices(); 
                    collision.replace(bodies, cindex);

                    tidymess.set_particles(bodies);
                    tidymess.commit_particles();

                    t_cpu = tcpu_offset + timer.read();   
                    output.write_snapshot(t, t_cpu, bodies);

                    if(bodies.size() == 1) {
                        stopping_condition_type = 4;
                
                        if(init.output_terminal) cout << "Only one body left. Ending simulation." << endl;
                        break;
                    }
                }
            }
        }
        
        // Roche limit handling
        if(collision.to_detect_roche) {                             
            if(tidymess.is_roche_detected() == true) {
                roche_flag = tidymess.get_roche_flag();
            
                if(collision.roche_mode == 1) {
                    if(init.output_terminal) cout << "Roche limit breached. Ignoring it. " << endl;
                }
                else if(collision.roche_mode == 2) {
                    stopping_condition_type = 3;
            
                    if(init.output_terminal) cout << "Roche limit breached. Ending simulation." << endl;
                    break;
                }
            }
        }

        // Breakup speed handling
        if(breakup.to_detect) {
            bool rotational_breakup = breakup.detect_breakup(bodies);           
            if(rotational_breakup) {
                breakup_flag = 1;
            
                if(breakup.mode == 1) {
                    if(init.output_terminal) cout << "Rotational breakup detected. Ignoring it. " << endl;
                }
                else if(breakup.mode == 2) {
                    stopping_condition_type = 5;
            
                    if(init.output_terminal) cout << "Rotational breakup detected. Ending simulation." << endl;
                    break;
                }
            } 
        }
        
        // Maximum CPU running time
        if(check_max_cpu) {
            if(t_cpu > max_cpu_time) {
                stopping_condition_type = 1;

                if(init.output_terminal) {
                    cout << "t_cpu        = " << t_cpu << " [s]" << endl;
                    cout << "max_cpu_time = " << max_cpu_time << " [s]" << endl; 
                    cout << "Halting simulation!" << endl;
                }
                break;
            }
        }
    }

    //---------------------------------------------------------------------
        
    // Save snapshot to binary file
    N = bodies.size();
    t_cpu = tcpu_offset + timer.read();
    dt_prev = tidymess.get_dt_prev();
    num_integration_step = tidymess.get_num_integration_step();
 
    output.save_to_binary(t, N, t_cpu, dt_prev, t_end, num_integration_step, collision_flag, roche_flag, breakup_flag, bodies, dt_snapshot, dt0_log, fmul_log, num_snapshot);

    // Clean up
    timer.stop();
    
    // Final log data
    int N_final = bodies.size();
    
    array<double, 3> r1 = tidymess.get_center_of_mass();
    array<double, 3> v1 = tidymess.get_center_of_mass_velocity();
    array<double, 3> L1_orb = tidymess.get_orbital_angular_momentum();
    array<double, 3> L1_spin = tidymess.get_spin_angular_momentum();
    double Ekin1_orb = tidymess.get_orbital_kinetic_energy();
    double Ekin1_spin = tidymess.get_spin_kinetic_energy();
    double Epot1 = tidymess.get_potential_energy();

    if(init.output_diag) {
        output.write_diag(t, t_cpu, num_integration_step, N, r1, v1, L1_orb, L1_spin, Ekin1_orb, Epot1, Ekin1_spin); 
    }

    double dr = sqrt(pow(r1[0]-r0[0], 2) + pow(r1[1]-r0[1], 2) + pow(r1[2]-r0[2], 2)); 
    double dv = sqrt(pow(v1[0]-v0[0], 2) + pow(v1[1]-v0[1], 2) + pow(v1[2]-v0[2], 2));

    vector<double> L0(3);
    L0[0] = L0_orb[0] + L0_spin[0];
    L0[1] = L0_orb[1] + L0_spin[1];
    L0[2] = L0_orb[2] + L0_spin[2];

    vector<double> L1(3);
    L1[0] = L1_orb[0] + L1_spin[0];
    L1[1] = L1_orb[1] + L1_spin[1];
    L1[2] = L1_orb[2] + L1_spin[2];

    double dLx = L1[0]-L0[0];
    double dLy = L1[1]-L0[1];
    double dLz = L1[2]-L0[2];

    double L0mag = sqrt(L0[0]*L0[0] + L0[1]*L0[1] + L0[2]*L0[2]);
    double L1mag = sqrt(L1[0]*L1[0] + L1[1]*L1[1] + L1[2]*L1[2]);
    double dLmag = sqrt(dLx*dLx + dLy*dLy + dLz*dLz);

    double dLx_rel = dLx;
    double dLy_rel = dLy;
    double dLz_rel = dLz;
    if(L0[0] != 0) dLx_rel /= L0[0];
    if(L0[1] != 0) dLy_rel /= L0[1];
    if(L0[2] != 0) dLz_rel /= L0[2];

    double dL_rel = dLmag;
    if(L0mag != 0) dL_rel /= (L0mag);

    double Etot0 = Ekin0_orb + Epot0 + Ekin0_spin;
    double Etot1 = Ekin1_orb + Epot1 + Ekin1_spin;

    double dE_abs = Etot1-Etot0;
    double dE_rel = dE_abs;
    if(Etot0 != 0) dE_rel /= Etot0;

    output.write_log(argc, argv, t, t_cpu, num_integration_step, N_init, N_final, dr, dv, dLx, dLy, dLz, dLx_rel, dLy_rel, dLz_rel, dLmag, dL_rel, dE_abs, dE_rel, bodies[0].r[0], stopping_condition_type, collision_flag, roche_flag, breakup_flag); 

    if(init.output_terminal) {
        output.print_log(argc, argv, t, t_cpu, num_integration_step, N_init, N_final, dr, dv, dLx, dLy, dLz, dLx_rel, dLy_rel, dLz_rel, dLmag, dL_rel, dE_abs, dE_rel, bodies[0].r[0], stopping_condition_type, collision_flag, roche_flag, breakup_flag);

        banner.print_outro(init.output_dir);

        banner.print_stars();
        banner.print_banner();
        banner.print_reference();
        banner.print_stars();
        cout << endl;
    }

    return 1;
}


