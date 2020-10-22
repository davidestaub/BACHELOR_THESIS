#include <iostream>
#include <array>
#include <vector>
#include"particles.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <chrono>
#include <random>
#include <math.h>
#include <utility>
#include <functional>
#include <set>
#include <fstream>
#include <iomanip>

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::Matrix2d;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;

class Simulation{
public:
    double mscale = std::pow(10,-6);
    Electrode lower_electrode = Electrode(1.0,{100.0* mscale,0.0* mscale, -500.0 * mscale}, 10.0,1000.0, 1000.0,0.0,40,240);
    Electrode upper_electrode = Electrode(-1.0,{100.0* mscale,0.0* mscale, 500.0 * mscale},10.0,1000.0,1000.0,0.0,40,240);
    int num_steps = 1;
    double time_step = Particle::time_step;
    int time_index = 0;
    std::vector<double> total_energy;
    std::vector<double> kinetic_energy;
    std::vector<double> potential_energy;
    VectorXd spring_force_vector;
    VectorXd friction_vector;
    VectorXd electricalfield_force_vector;
    MatrixXd spring_force_matrix_dx;
    MatrixXd spring_force_matrix_dv;
    MatrixXd friction_matrix_dv;
    MatrixXd mass_matrix;
    VectorXd velocity_vector;
    VectorXd position_vector;
    VectorXd position_vector_minus_1;
    VectorXd position_vector_minus_2;
    VectorXd initial_position_vector;
    bool first_itteration = true;
    bool dimers_only = false;
    bool tracer_on = false;
    bool dimer_tracers = false;
    bool erase_dimer_tracers = false;
    VectorXd inverse_number_of_neighbours;
    VectorXd drag_occlusion_factor_vector;



    double beta_ehd = 0.03;
    double K_plus_charge = 1.0;
    double CL_minus_charge = -1.0;
    double elementary_charge = 1.60217662 * std::pow(10,-19);
    double molarity = std::pow(10,-5) * std::pow(10,-3); // convert from mol/L to mol/m^3
    double NA = 6.02214076  * std::pow(10,23);
    double D_K_plus = 1.960 * std::pow(10,-9);
    double D_Cl_minus = 2.030 * std::pow(10,-9);

    bool safety = false;


    // all units in micro scale, (10^-6)

    double eta = 1.0 * std::pow(10,-3); //viscosity of water in PA * seconds
    double dynamic_viscosity =  8.9 * std::pow(10,-4); //dynamic vsicosity of water in PA * s
    double viscosity_multiplier = 1.0;

    double gas_constant = 8.31446261815324;
    double faraday_constant = 9.64853321233100184 * std::pow(10.0,4);

    double vacuum_permittivity = 8.8541878176 * std::pow(10.0,-12);
    double relative_permittivity = 80.2; //water at room temperatur

    double conductivity = 0.0015 * std::pow(10.0,-1);
    double density = 997.0;

    double T = 293; //room temperatur
    double kB = 1.38 * std::pow(10.0,-23); //Boltzmann constant
    Vector2d force_damper;

    double noise_damper = 4;
    double drag_damper = 0;
    double brownian_motion_multiplier = 0.0;

    Vector3d Boxcenter = {0,0,0};
    Vector3d Boxsize = {1500,1500,1500};
    double xbuffer, ybuffer, zbuffer = 5.0;

    //SPRING STUFF//
    double stiffnes_constant = 4000.0;
    //double rest_length = 2 * Particle::radius_for_spring;
    double damping_constant = 50.0;


    std::vector<double> velocities_over_time1_in_x;
    std::vector<double> velocities_over_time1_in_y;
    std::vector<double> e_vec;
    std::vector<double> spring_force_vec_over_time_x;
    std::vector<double> spring_force_vec_over_time_y;
    std::vector<double> spring_force_derivative_x_in_x;
    std::vector<double> A_B_DISTANCE;
    std::vector<double> position_vec_over_time_in_x;
    std::vector<double> position_vec_over_time_in_y;
    std::vector<double> position_vec_over_time_in_x_2;
    std::vector<double> position_vec_over_time_in_y_2;
    std::vector<double> position_vec_over_time_in_x_3;
    std::vector<double> position_vec_over_time_in_y_3;
    std::vector<double> position_vec_over_time_in_x_4;
    std::vector<double> position_vec_over_time_in_y_4;


    std::vector<double> rotator_position_vec_over_time_in_x;
    std::vector<double> rotator_position_vec_over_time_in_y;
    std::vector<double> rotator_position_vec_over_time_in_x_2;
    std::vector<double> rotator_position_vec_over_time_in_y_2;

    std::vector<double> triangle_position_vec_over_time_in_x;
    std::vector<double> triangle_position_vec_over_time_in_y;
    std::vector<double> triangle_position_vec_over_time_in_x_2;
    std::vector<double> triangle_position_vec_over_time_in_y_2;
    std::vector<double> triangle_position_vec_over_time_in_x_3;
    std::vector<double> triangle_position_vec_over_time_in_y_3;

    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> z_values1;
    std::vector<double> z_values2;
    std::vector<double> z_values;
    std::vector<double> best_found;
    std::vector<double> forces_FD;
    std::vector<double> energy_FD;
    std::vector<double> force_jacobian_FD;
    std::vector<double> dFx_dx_FD;
    std::vector<double> F_ma;
    std::vector<double> b_force_x;
    std::vector<double> b_force_y;
    std::vector<double> b_force_z;
    std::vector<double> ehd_force_x;
    std::vector<double> ehd_force_y;
    std::vector<double> spring_force_x;
    std::vector<double> spring_force_y;

    std::vector<double> e_x;
    std::vector<double> e_y;
    std::vector<double> e_z;

    std::vector<Vector3d> tracer;
    std::vector<std::vector<Vector3d>> dimer_tracers_vector;


    std::vector<double> velocities;
    std::vector<double> velocities_2;
    std::vector<double> frequencies;
    std::vector<double> Efields;

    std::vector<double> amount_x;
    std::vector<double> amount_y;
    std::vector<double> amount_z;
    std::vector<double> direction_norm;
    double last_frequency = 0.0;
    double last_peak_voltage = 0.0;

    std::vector<double> drag_x, drag_y, drag_px, drag_py;

    std::vector<double> friction_force_over_time_x;

    Vector3d ex = {1.0,0.0,0.0};

    bool test_bool = true;

    int max_iterations = 1000;

    int size = -1;

    VectorXd current_thermal_force;
    VectorXd current_rotational_thermal_force;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0.0,1.0};



    Simulation() {
    }






    //INTERFACES AND NETWONS METHOD

    void run_simulation_for_connected_Particles(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        std::cout<<"begin ========================"<<std::endl;
        for(int i = 0; i<num_steps;i++) {
            total_energy.push_back(0.0);
            kinetic_energy.push_back(0.0);
            potential_energy.push_back(0.0);

           UPDATE_SYSTEM_NEW(connected_particles);


        }
        std::cout<<"end ========================"<<std::endl;
    }

    void UPDATE_SYSTEM_NEW(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

      //  std::cout<<" test 1"<<std::endl;

        current_thermal_force = safety  * get_current_thermal_force(connected_particles); //make sure to only call it once a time step due to the random nature of the force.
        //std::cout<<"current thermal force x = "<<current_thermal_force[0]<<std::endl<<"current thermal force x2 = "<<current_thermal_force[3]<<std::endl;
        current_rotational_thermal_force = safety * get_current_rotational_thermal_force(connected_particles);
        b_force_x.push_back(current_thermal_force[0]);
        b_force_y.push_back(current_thermal_force[3]);
        //std::cout<<" test 2"<<std::endl;

        fill_spring_force_and_derivatives_NEW(connected_particles);
        //std::cout<<" test 3"<<std::endl;

        first_itteration = false;

        VectorXd x_t_plus_1_init = position_vector + (position_vector - position_vector_minus_1);
        VectorXd x_t_plus_1_old = x_t_plus_1_init;

        //std::cout<<" test 4"<<std::endl;


        //int maxiter, double beta, double tol, VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles
        THE_FINAL_NEWTON_with_drag(max_iterations,0.5,std::pow(10,(-10)), x_t_plus_1_old, connected_particles);

       // std::cout<<" test 5"<<std::endl;

        //Newtons_Method_NEW_with_force_update(100000,std::pow(10,(-8)), 0.25, 0.5,  x_t_plus_1_old,connected_particles);


        ehd_force_x.push_back(get_current_EHD_force(connected_particles,x_t_plus_1_old)[0]);
       // std::cout<<" test 6"<<std::endl;
        ehd_force_y.push_back(get_current_EHD_force(connected_particles,x_t_plus_1_old)[1]);
       // std::cout<<" test 7"<<std::endl;

        spring_force_x.push_back(get_current_spring_force(connected_particles,x_t_plus_1_old)[0]);
       // std::cout<<" test 8"<<std::endl;
        spring_force_y.push_back(get_current_spring_force(connected_particles,x_t_plus_1_old)[1]);
       // std::cout<<" test 9"<<std::endl;

        VectorXd tmp = position_vector_minus_1;
        //std::cout<<" test 10"<<std::endl;
        if(tracer_on) {
            //std::cout<<" test 11"<<std::endl;
            Vector3d com = Eigen::VectorXd::Zero(3);
            for(int i = 0; i<x_t_plus_1_old.size(); i+=3){
                Vector3d current = {x_t_plus_1_old[i], x_t_plus_1_old[i+1],x_t_plus_1_old[i+2]};
                com = com + current;
            }
            com = com/(x_t_plus_1_old.size()/3);
            tracer.push_back({com[0], com[1],com[2]});
        }else{
           // std::cout<<" test 12"<<std::endl;
            tracer.erase(tracer.begin(), tracer.end());
        }
        //std::cout<<" test 13"<<std::endl;

        if(dimer_tracers){
           // std::cout<<" test 14"<<std::endl;
            for(int i = 0; i<connected_particles.size(); i++){
                //dimer_tracers_vector.push_back(sdt::vector<Vector3d>)
                std::cout<<"test 1";
                std::tuple<Particle&, Particle&,double> particle_pair = connected_particles[i];
                std::cout<<"test 2";
                Vector3d com = Eigen::VectorXd::Zero(3);
                std::cout<<"test 3";
                com = (std::get<0>(particle_pair).position + std::get<1>(particle_pair).position)/2.0;
                std::cout<<"test 4";
                dimer_tracers_vector[i].push_back(com);
                std::cout<<"test 5";


            }
           // std::cout<<" test 15"<<std::endl;
        }if(erase_dimer_tracers){
            for(int i = 0; i<connected_particles.size(); i++){
                //dimer_tracers_vector.erase(dimer_tracers_vector.begin(), dimer_tracers_vector.end());
                dimer_tracers_vector[i].erase(dimer_tracers_vector[i].begin(), dimer_tracers_vector[i].end());
            }
        }
       // std::cout<<" test 16"<<std::endl;


        Vector3d current_amount_EHD = get_amount(std::get<0>(connected_particles[0]), std::get<1>(connected_particles[0]));
        amount_x.push_back(current_amount_EHD[0] * safety);
        amount_y.push_back(current_amount_EHD[1] * safety);
        amount_z.push_back(current_amount_EHD[2] * safety);
        direction_norm.push_back(current_amount_EHD.norm() * safety);
        VectorXd current_drag = get_current_friction_force(connected_particles,x_t_plus_1_old);
       // std::cout<<" test 17"<<std::endl;

        drag_x.push_back(current_drag[0] * safety);
        drag_y.push_back(current_drag[1] * safety);
        drag_px.push_back(get_stokes_friction_force(std::get<0>(connected_particles[0]))[0] * safety);
        drag_py.push_back(get_stokes_friction_force(std::get<0>(connected_particles[0]))[1] * safety);
       // std::cout<<" test 18"<<std::endl;


       // std::cout<<" Position at t-1= \n"<<tmp<<"\n Positoin at t = \n"<<position_vector<<"\n Position at t+1 = \n"<<x_t_plus_1_old<<"\n \n \n";
        //energy_FD.push_back(get_dE_dx(connected_particles,position_vector,x_t_plus_1_old));
        //forces_FD.push_back(get_F_x_FD(connected_particles));
        //force_jacobian_FD.push_back(get_F_x_d_x(connected_particles));
        //std::cout<<" test 20"<<std::endl;

        //double current_epot = get_Potential_Energy(connected_particles);
        //double current_ekin = get_Kinetic_Energy(connected_particles);
        //potential_energy.push_back(current_epot);
        //kinetic_energy.push_back(current_ekin);
        //total_energy.push_back(current_epot + current_ekin);
        //std::cout<<" test 21"<<std::endl;

        VectorXd F_x = get_current_spring_force(connected_particles,x_t_plus_1_old) - safety *get_current_friction_force(connected_particles,x_t_plus_1_old) + safety *electricalfield_force_vector + safety *get_current_EHD_force(connected_particles,x_t_plus_1_old) + safety *current_thermal_force;
        VectorXd M_a = mass_matrix * (x_t_plus_1_old  - 2.0 * position_vector + position_vector_minus_1);
        VectorXd F_minus_M_time_a = F_x - M_a;
        F_ma.push_back(F_minus_M_time_a.norm());
        std::cout<<" test 22"<<std::endl;



        VectorXd tmp_tmp = position_vector;
        std::vector<std::tuple <Particle&,Particle&,double> > tmp_connected_particles = connected_particles;
       // std::cout<<" test 23"<<std::endl;
       // std::cout<<"pos old"<<std::get<0>(connected_particles[0]).position<<std::endl;
        update_positions_NEW(connected_particles,x_t_plus_1_old);
      //  std::cout<<" test 24"<<std::endl;
       // std::cout<<"pos new"<<std::get<0>(connected_particles[0]).position<<std::endl;
       // std::cout<<"time step"<<time_step<<std::endl;
        dFx_dx_FD.push_back( get_dFx_dx_FD(tmp_connected_particles, connected_particles, tmp, tmp_tmp, x_t_plus_1_old));
       // std::cout<<" test 25"<<std::endl;

        double e = get_e(x_t_plus_1_old,connected_particles);

       // std::cout<<" e = "<<std::endl<<e<<std::endl;
        e_vec.push_back(e);

       // add_brownian_motion(connected_particles);
      //  std::cout<<" test 26"<<std::endl;

        update_rest_length(connected_particles);
       // std::cout<<" test 27"<<std::endl;

        generate_data_for_v_w_plot(connected_particles);
        //std::cout<<" test 28"<<std::endl;

        apply_wind_shadow(connected_particles);

       // std::cout<<" test 29"<<std::endl;

       // dissolve(connected_particles,0.1);

    }

    void update_positions_NEW(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles,VectorXd x_t_plus_1_new){
        position_vector_minus_2 = position_vector_minus_1;
        position_vector_minus_1 = position_vector;
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited){
                velocities_over_time1_in_x.push_back(A.velocity[0]);
                velocities_over_time1_in_y.push_back(A.velocity[1]);

                A.velocity[0] = (x_t_plus_1_new(A.index) -A.position[0])/time_step;
               //A.velocity = A.get_velocity(x_t_plus_1_new);

                A.velocity[1] = (x_t_plus_1_new(A.index+1) -A.position[1])/time_step;

                A.velocity[2] = (x_t_plus_1_new(A.index+2) -A.position[2])/time_step;

                A.position[0] = x_t_plus_1_new(A.index);
                A.position[1] = x_t_plus_1_new(A.index+1);
                A.position[2] = x_t_plus_1_new(A.index+2);
                A.visited=true;
            }
            if(!B.visited){
               // B.velocity = B.get_velocity(x_t_plus_1_new);

                B.velocity[0] = (x_t_plus_1_new(B.index) - B.position[0])/time_step;

                B.velocity[1] = (x_t_plus_1_new(B.index+1) - B.position[1])/time_step;

                B.velocity[2] = (x_t_plus_1_new(B.index+2) - B.position[2])/time_step;

                B.position[0] = x_t_plus_1_new(B.index);
                B.position[1] = x_t_plus_1_new(B.index+1);
                B.position[2] = x_t_plus_1_new(B.index+2);
                B.visited=true;
            }
        }
        position_vector = x_t_plus_1_new;

        reset_flags(connected_particles);
    }

    VectorXd initialize_position_vector_minus_1(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited){
                position_vector_minus_1(A.index) = A.position[0];
                position_vector_minus_1(A.index+1) = A.position[1];
                position_vector_minus_1(A.index+2) = A.position[2];
                A.visited = true;
            }
            if(!B.visited){
                position_vector_minus_1(B.index) = B.position[0];
                position_vector_minus_1(B.index+1) = B.position[1];
                position_vector_minus_1(B.index+2) = B.position[2];
                B.visited = true;
            }

        }

        reset_flags(connected_particles);
        position_vector_minus_2 = position_vector_minus_1;
        return position_vector_minus_1;
    }

    void THE_FINAL_NEWTON_with_drag(int maxiter, double beta, double tol, VectorXd& x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double t;
        double r;
        VectorXd x_i = x;
        int i;

        for(i=0;i<maxiter;i++){
            t=1.0;
            r=0.0;
            //     std::cout<<"iteration = "<<i<<std::endl;

            VectorXd g = get_g_with_drag(x_i,connected_particles);
            //     std::cout<<" g = "<<g<<std::endl;
            MatrixXd g_dx = get_g_dx(x_i,connected_particles);


            // MatrixXd g_dv = get_g_dv(x_i,connected_particles);
            //   std::cout<<"g_dx = "<<g_dx<<std::endl;
            //std::cout<<"g_dv = "<<g_dv<<std::endl;

            MatrixXd H = g_dx;
            // std::cout<<"H = "<<H<<std::endl;



            if(g.transpose() * (-H.inverse() * g) <= 0){ //< or <=
                //    std::cout<<"good search direction"<<std::endl;
                //   std::cout<<" gT delta x = "<<g.transpose() * (-g_dx * g)<<std::endl;
            }else{
                //    std::cout<<"wrong search direaction"<<std::endl;
                //   std::cout<<" gT delta x = "<<g.transpose() * (-g_dx * g)<<std::endl;
                r += 0.01;
                //     std::cout<<" r = "<<r<<std::endl;
                // g_dx = g_dx + Eigen::MatrixXd::Identity(g_dx.rows(),g_dx.cols()) *r;
                H = H + Eigen::MatrixXd::Identity(H.rows(),H.cols()) * r;
                //   std::cout<<" gT delta x = "<<g.transpose() * (-g_dx * g)<<std::endl;

                while(g.transpose() * (-H.inverse() * g) > 0){ //> or >=
                    //      std::cout<<"inside while loop"<<std::endl;
                    r *=10;
                    //         std::cout<<" r now is "<<r<<std::endl;
                    //g_dx = g_dx + Eigen::MatrixXd::Identity(g_dx.rows(),g_dx.cols()) *r;
                    H = H + Eigen::MatrixXd::Identity(H.rows(),H.cols()) * r;
                    //     std::cout<<" gT delta x = "<<g.transpose() * (-g_dx * g)<<std::endl;
                }
            }

            VectorXd search_direction = -1.0 * H.inverse() * g;


            // VectorXd x_i_plus_1_tmp = x_i - t * g_dx.inverse() * g;
            VectorXd x_i_plus_1_tmp = x_i + t * search_direction;
            //   std::cout<<"the new value ="<<std::endl<<x_i_plus_1_tmp<<std::endl;
            //    std::cout<<" H inverse = "<<std::endl<<g_dx.inverse()<<std::endl;
            //    std::cout<<"search direction = "<<std::endl<<g_dx.inverse() * g<<std::endl;

            //  std::cout<<" x i +1="<<std::endl<<x_i_plus_1_tmp<<std::endl<<" e at x i+1 = "<<std::endl<<get_e(x_i_plus_1_tmp,connected_particles);
            //   std::cout<<" x i ="<<std::endl<<x_i<<std::endl<<" e at x i = "<<std::endl<<get_e(x_i,connected_particles);

            while(get_g_with_drag(x_i_plus_1_tmp,connected_particles).norm() > get_g_with_drag(x_i,connected_particles).norm()/* get_e(x_i_plus_1_tmp,connected_particles) > get_e(x_i,connected_particles)*/){
                //       std::cout<<" inside second while loop, e(x_i+1) = "<<get_e(x_i_plus_1_tmp,connected_particles)<<" and e(xi) = "<<get_e(x_i,connected_particles)<<std::endl;
                t *= beta;
                //      std::cout<<"t = "<<t<<std::endl;
                x_i_plus_1_tmp = x_i + t * search_direction;
            }
            VectorXd x_i_plus_1 = x_i_plus_1_tmp;
            //  std::cout<<" e after while "<<std::endl<<get_e(x_i_plus_1,connected_particles)<<std::endl;

            //  std::cout<<" x i "<<std::endl<<x_i<<std::endl;
            //  std::cout<<" x i +1 "<<std::endl<<x_i_plus_1<<std::endl;



            x_i = x_i_plus_1;

            if(get_g_with_drag(x_i,connected_particles).norm()<tol){
                //   std::cout<<"Newton's Method took "<<i<<" iterations to converge and the value of e at the end is " <<get_e(x_i,connected_particles)<<std::endl;
                x = x_i;
                break;

            }else{
                //      std::cout<<" residual = "<<get_g(x_i,connected_particles).norm()<<std::endl;
            }


        }
        x = x_i;
        std::cout<<" \n \n \n";
        std::cout<<" NEwtons method terminated after "<<i<<" iterations "<<" the residual is "<<get_g_with_drag(x_i,connected_particles).norm()<<std::endl;
        std::cout<<" \n \n \n";
    }

    void reset_flags(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            A.visited = false;
            B.visited = false;
        }
    }

    void assign_index(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        int count = 0;
        for(auto& particle_pair : connected_particles) {
            //lexicographic order
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited){
                A.visited = true;
                A.index = count;
                count+=3;
            }
            if(!B.visited){
                B.visited = true;
                B.index = count;
                count+=3;
            }
        }
        size = count;
        spring_force_vector = VectorXd::Zero(size);
        friction_vector = VectorXd::Zero(size);
        spring_force_matrix_dx = MatrixXd::Zero(size,size);
        spring_force_matrix_dv = MatrixXd::Zero(size,size);
        friction_matrix_dv = MatrixXd::Zero(size,size);
        mass_matrix = MatrixXd::Zero(size,size);
        velocity_vector = VectorXd::Zero(size);
        position_vector = VectorXd::Zero(size);
        position_vector_minus_1 = VectorXd::Zero(size);
        electricalfield_force_vector = VectorXd::Zero(size);
        inverse_number_of_neighbours = VectorXd::Zero(size);
        drag_occlusion_factor_vector = VectorXd::Ones(size);
        dimer_tracers_vector.resize(size/6);
        reset_flags(connected_particles);
    }

    void count_neighbours(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            A.number_of_neighbours += 1.0;
            Particle &B = std::get<1>(particle_pair);
            B.number_of_neighbours += 1.0;
        }
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            inverse_number_of_neighbours[A.index] = 1.0/A.number_of_neighbours;
            inverse_number_of_neighbours[A.index+1] = 1.0/A.number_of_neighbours;
            inverse_number_of_neighbours[A.index+2] = 1.0/A.number_of_neighbours;

            inverse_number_of_neighbours[B.index] = 1.0/B.number_of_neighbours;
            inverse_number_of_neighbours[B.index+1] = 1.0/B.number_of_neighbours;
            inverse_number_of_neighbours[B.index+2] = 1.0/B.number_of_neighbours;

        }
    }

    void connect(Particle &A, Particle &B, double restlength, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        std::tuple<Particle&,Particle&,double> particle_pair(A,B,restlength);
        connected_particles.push_back(particle_pair);
    }
    void connect_new(Particle &A, Particle &B, double restlength, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, int i){
        std::tuple<Particle&,Particle&,double> particle_pair(A,B,restlength);
        connected_particles[i] = particle_pair;
    }

    void update_rest_length(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        if(dimers_only) {
            for (auto &particle_pair : connected_particles) {
                std::get<2>(particle_pair) = set_rest_length_rad_plus_rad(std::get<0>(particle_pair),
                                                                          std::get<1>(particle_pair));
            }
        }
    }

    void reset_simulation(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

        spring_force_vector.setZero();
        spring_force_matrix_dx.setZero();
        spring_force_matrix_dv.setZero();
        friction_vector.setZero();
        friction_matrix_dv.setZero();
        velocity_vector.setZero();
        position_vector.setZero();
        tracer.erase(tracer.begin(),tracer.end());
        for(auto particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited) {
                A.position[0] = initial_position_vector(A.index);
                A.position[1] = initial_position_vector(A.index+1);
                A.position[2] = initial_position_vector(A.index+2);
                A.velocity = A.initial_velocity;
                position_vector(A.index) = A.position[0];
                position_vector(A.index+1) = A.position[1];
                position_vector(A.index+2) = A.position[2];
                A.visited = true;
            }
            if(!B.visited){
                B.position[0] = initial_position_vector(B.index);
                B.position[1] = initial_position_vector(B.index+1);
                B.position[2] = initial_position_vector(B.index+2);
                B.velocity = B.initial_velocity;
                position_vector(B.index) = B.position[0];
                position_vector(B.index+1) = B.position[1];
                position_vector(B.index+2) = B.position[2];
                B.visited = true;
            }
        }
        reset_flags(connected_particles);
        position_vector_minus_1 = position_vector;

    }

    double set_rest_length_rad_plus_rad(Particle& A, Particle& B){
        return A.radius + B.radius;
    }

    void generate_data_for_v_w_plot(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        auto pair = connected_particles[0];
        Particle& A = std::get<0>(pair);
        double velocity = A.velocity[0];

        if(!(last_frequency == lower_electrode.frequency) && lower_electrode.peak_voltage>= 0.0001 && safety){
            velocities.push_back(velocity);
            frequencies.push_back(lower_electrode.frequency);
        }
        last_frequency = lower_electrode.frequency;

        if(!(last_peak_voltage == lower_electrode.peak_voltage) && lower_electrode.peak_voltage >= 0.0001){
            velocities_2.push_back(velocity);
            Efields.push_back(get_Eo());
        }
        last_peak_voltage = lower_electrode.peak_voltage;
    }








    //CALCULATION FOR GLOBAL STUFF DEPENDANT ON CURRENT POSITION TO USE INSIDE NEWTONS METHOD

    MatrixXd get_current_spring_force_jacobian(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_matrix_dx.rows();
        int  m = spring_force_matrix_dx.cols();

        MatrixXd spring_force_jacobian(size,size);
        spring_force_jacobian.setZero();

        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = current_positions[A.index];
            A_tmp.position[1] = current_positions[A.index+1];
            A_tmp.position[2] = current_positions[A.index+2];
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];
            B_tmp.position[2] = current_positions[B.index+2];

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

            Matrix3d current_spring_forceA_dx = NEW_get_damped_spring_force_dXA_without_damping(A_tmp,B_tmp,std::get<2>(particle_pair));
           /* Matrix3d current_spring_forceB_dx = NEW_get_damped_spring_force_dXA_without_damping(B_tmp,A_tmp,std::get<2>(particle_pair));

            spring_force_jacobian.block<3,3>(A.index,A.index) += current_spring_forceA_dx;
            spring_force_jacobian.block<3,3>(B.index,B.index) += current_spring_forceB_dx;
            spring_force_jacobian.block<3,3>(A.index,B.index) += -1.0 *current_spring_forceA_dx;
            spring_force_jacobian.block<3,3>(B.index,A.index) += -1.0 *current_spring_forceB_dx;*/

            spring_force_jacobian.block<3,3>(A.index,A.index) += current_spring_forceA_dx;
            spring_force_jacobian.block<3,3>(B.index,B.index) += current_spring_forceA_dx;
            spring_force_jacobian.block<3,3>(A.index,B.index) += (-1.0 * current_spring_forceA_dx);
            spring_force_jacobian.block<3,3>(B.index,A.index) += (-1.0 * current_spring_forceA_dx);


        }

        return spring_force_jacobian;
    }

    VectorXd get_current_spring_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_vector.size();
       // std::cout<<" current positions = "<<std::endl<<current_positions<<std::endl;


        VectorXd current_spring_force(n);
        current_spring_force.setZero();

        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = current_positions[A.index];
            A_tmp.position[1] = current_positions[A.index+1];
            A_tmp.position[2] = current_positions[A.index+2];
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];
            B_tmp.position[2] = current_positions[B.index+2];

           // std::cout<< "in here the positions are:"<<std::endl<<A_tmp.position[0]<<std::endl<<A_tmp.position[1]<<std::endl<<B_tmp.position[0]<<std::endl<<B_tmp.position[1]<<std::endl;

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

           Vector3d current_spring_force_tmp = NEW_get_damped_spring_force_without_damping(particle_pair_tmp);
           //std::cout<<"current spring force  tmp = "<<std::endl<<current_spring_force_tmp<<std::endl;

            current_spring_force[A.index] += current_spring_force_tmp[0];
            current_spring_force[A.index + 1] += current_spring_force_tmp[1];
            current_spring_force[A.index + 2] += current_spring_force_tmp[2];
            current_spring_force[B.index] += (-1.0 * current_spring_force_tmp[0]);
            current_spring_force[B.index + 1] += (-1.0 *current_spring_force_tmp[1]);
            current_spring_force[B.index + 2] += (-1.0 *current_spring_force_tmp[2]);



        }
       // std::cout<<"current spring force  = "<<std::endl<<current_spring_force<<std::endl;

        return current_spring_force;
    }

    VectorXd get_current_friction_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_vector.size();
        VectorXd current_friction_force(n);
        current_friction_force.setZero();
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited){
                Particle A_tmp = A;
                A_tmp.position[0] = current_positions[A.index];
                A_tmp.position[1] = current_positions[A.index+1];
                A_tmp.position[2] = current_positions[A.index+2];
                A_tmp.radius = A.radius;
                Vector3d pos_at_A;
                pos_at_A[0] = position_vector[A.index];
                pos_at_A[1] = position_vector[A.index+1];
                pos_at_A[2] = position_vector[A.index+2];
                A_tmp.velocity = ( A_tmp.position - pos_at_A )/time_step;
                Vector3d current_friction_at_A = get_stokes_friction_force(A_tmp);
                current_friction_force[A.index] = current_friction_at_A[0];
                current_friction_force[A.index + 1] = current_friction_at_A[1];
                current_friction_force[A.index + 2] = current_friction_at_A[2];
                A.visited=true;
            }

            if(!B.visited){
                Particle B_tmp = B;
                B_tmp.position[0] = current_positions[B.index];
                B_tmp.position[1] = current_positions[B.index+1];
                B_tmp.position[2] = current_positions[B.index+2];
                Vector3d pos_at_B;
                B_tmp.radius = B.radius;
                pos_at_B[0] = position_vector[B.index];
                pos_at_B[1] = position_vector[B.index+1];
                pos_at_B[2] = position_vector[B.index+2];
                B_tmp.velocity = ( B_tmp.position -pos_at_B )/time_step;
                Vector3d current_friction_at_B = get_stokes_friction_force(B_tmp);
                current_friction_force[B.index] = current_friction_at_B[0];
                current_friction_force[B.index + 1] =current_friction_at_B[1];
                current_friction_force[B.index + 2] =current_friction_at_B[2];
                B.visited=true;
            }
        }
        reset_flags(connected_particles);
        return current_friction_force/*.cwiseProduct(drag_occlusion_factor_vector)*/;
    }

    MatrixXd get_current_friction_force_jacobian(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        MatrixXd result(friction_matrix_dv.rows(),friction_matrix_dv.cols());
        result.setZero();
        double constant = (6.0 * M_PI * dynamic_viscosity * viscosity_multiplier)/time_step;
        for(auto& particle_pair : connected_particles) {

            Particle &A = std::get<0>(particle_pair);
            Particle &B = std::get<1>(particle_pair);

            if(!A.visited){
                result(A.index,A.index) = constant * A.radius;
                result(A.index+1,A.index+1) = constant * A.radius;
                result(A.index+2,A.index+2) = constant * A.radius;
                A.visited=true;
            }
            if(!B.visited){
                result(B.index,B.index) = constant * B.radius;
                result(B.index+1,B.index+1) = constant * B.radius;
                result(B.index+2,B.index+2) = constant * B.radius;
                B.visited = true;
            }

        }
        reset_flags(connected_particles);
        return result /** drag_occlusion_factor_vector.asDiagonal()*/;
    }

    MatrixXd get_current_electrical_force_jacobian(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        MatrixXd result(friction_matrix_dv.rows(),friction_matrix_dv.cols());
        result.setZero();
        return result;
    }

    void fill_spring_force_and_derivatives_NEW(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

        spring_force_vector.setZero();
        spring_force_matrix_dx.setZero();
        spring_force_matrix_dv.setZero();
        friction_vector.setZero();
        friction_matrix_dv.setZero();
        //mass_matrix.setZero();
        velocity_vector.setZero();



        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
          //  std::cout<<" Velocity  of particle"<<A.index<<"= "<<A.velocity[0]  << " micrometer/s"<<" and frequency is = "<<lower_electrode.frequency<<std::endl<<" time step = "<<time_step<<std::endl;
          //  std::cout<<"A.mass = "<<A.mass<<std::endl;
          //  std::cout<<"A.radius = "<<A.radius<<std::endl;
           // std::cout<<std::endl<<std::endl<<"POSITION"<<std::endl<<A.position<<std::endl;
            Particle& B = std::get<1>(particle_pair);
          //  std::cout<<" Velocity  of particle"<<B.index<<"= "<<B.velocity[0]  << " micrometer/s"<<" and frequency is = "<<lower_electrode.frequency<<std::endl<<" time step = "<<time_step<<std::endl;

         // std::cout<<"A.index = "<<A.index<<"A.radius = "<<A.radius<<std::endl;
         // std::cout<<" B index = "<<B.index<<"B.radius = "<<B.radius<<std::endl;

            Vector3d current_spring_force = NEW_get_damped_spring_force_without_damping(particle_pair);
            Matrix3d current_spring_forceA_dx = NEW_get_damped_spring_force_dXA_without_damping(A,B,std::get<2>(particle_pair));
            Matrix3d current_spring_forceB_dx = NEW_get_damped_spring_force_dXA_without_damping(B,A,std::get<2>(particle_pair));
           // std::cout<<"Spring Force Vector Before Filling \n"<<spring_force_vector<<"\n";
            spring_force_vector[A.index] += current_spring_force[0];
            spring_force_vector[A.index + 1] += current_spring_force[1];
            spring_force_vector[A.index + 2] += current_spring_force[2];
            spring_force_vector[B.index] += -1.0 * current_spring_force[0];
            spring_force_vector[B.index + 1] += -1.0 *current_spring_force[1];
            spring_force_vector[B.index + 2] += -1.0 *current_spring_force[2];
          //  std::cout<<"Spring Force Vector After Filling at positions "<<A.index<<" and "<<B.index<<"\n"<<spring_force_vector<<"\n";
            spring_force_vec_over_time_x.push_back(spring_force_vector[A.index]);
           // std::cout<<spring_force_vec_over_time_x[0]<<std::endl;
            spring_force_vec_over_time_y.push_back(spring_force_vector[A.index+1]);
           // std::cout<<"Spring Force DX Before Filling \n"<<spring_force_matrix_dx<<"\n";
            spring_force_matrix_dx.block<3,3>(A.index,A.index) += current_spring_forceA_dx;
            spring_force_matrix_dx.block<3,3>(B.index,B.index) += current_spring_forceB_dx;
            spring_force_matrix_dx.block<3,3>(A.index,B.index) += -1.0 * current_spring_forceA_dx;
            spring_force_matrix_dx.block<3,3>(B.index,A.index) += -1.0 * current_spring_forceB_dx;

            //friction_vector[A.index] = get_stokes_friction_force(A)[0];
            //friction_vector[A.index+1] = get_stokes_friction_force(A)[1];
            //friction_vector[A.index+2] = get_stokes_friction_force(A)[2];
            //friction_vector[B.index] = get_stokes_friction_force(B)[0];
            //friction_vector[B.index+1] = get_stokes_friction_force(B)[1];
            //friction_vector[B.index+2] = get_stokes_friction_force(B)[2];

            //friction_force_over_time_x.push_back(friction_vector[A.index]);

            friction_matrix_dv.block<3,3>(A.index,A.index) = get_friction_force_dv(A);
            friction_matrix_dv.block<3,3>(B.index,B.index) = get_friction_force_dv(B);

            electricalfield_force_vector[A.index] = get_electrical_force(A)[0];
            electricalfield_force_vector[A.index +1] = get_electrical_force(A)[1];
            electricalfield_force_vector[A.index +2] = get_electrical_force(A)[2];
            electricalfield_force_vector[B.index] = get_electrical_force(B)[0];
            electricalfield_force_vector[B.index+1] = get_electrical_force(B)[1];
            electricalfield_force_vector[B.index+2] = get_electrical_force(B)[2];
           // std::cout<<" e force = "<<electricalfield_force_vector[A.index +2]<<std::endl;



           // std::cout<<"Spring Force DX After Filling at positions "<<A.index<<" and "<<B.index<<"\n"<<spring_force_matrix_dx<<"\n";
            spring_force_derivative_x_in_x.push_back(spring_force_matrix_dx(A.index+1,A.index+1));

            //spring_force_matrix_dv.block<2,2>(A.index,A.index) += current_spring_force_dv;
            //spring_force_matrix_dv.block<2,2>(B.index,B.index) += -1.0 * current_spring_force_dv;

          //  std::cout<<"Mass MAtrix BEfore filling \n"<<mass_matrix<<"\n";
            mass_matrix(A.index,A.index) = A.mass;
            mass_matrix(A.index+1,A.index+1) = A.mass;
            mass_matrix(A.index+2,A.index+2) = A.mass;

            mass_matrix(B.index,B.index) = B.mass;
            mass_matrix(B.index+1,B.index+1) = B.mass;
            mass_matrix(B.index+2,B.index+2) = B.mass;
           // std::cout<<"Mass MAtrix After filling \n"<<mass_matrix<<"\n";

           // std::cout<<"Position Vector BEfore filling \n"<<position_vector<<"\n";
            position_vector(A.index) = A.position[0];
            position_vector(A.index+1) = A.position[1];
            position_vector(A.index+2) = A.position[2];
            position_vector(B.index) = B.position[0];
            position_vector(B.index+1) = B.position[1];
            position_vector(B.index+2) = B.position[2];
            if(first_itteration){
                initial_position_vector = position_vector;
            }
           // std::cout<<"Position Vector After filling \n"<<position_vector<<"\n";
           if(A.index == 0 && !A.visited) {
           //    std::cout<<" true 0";

               position_vec_over_time_in_x.push_back(A.position[0]);
               position_vec_over_time_in_y.push_back(A.position[1]);

               rotator_position_vec_over_time_in_x.push_back(A.position[0]);
               rotator_position_vec_over_time_in_y.push_back(A.position[1]);

               triangle_position_vec_over_time_in_x.push_back(A.position[0]);
               triangle_position_vec_over_time_in_y.push_back(A.position[1]);
               A.visited = true;
           }
           if(B.index == 3 && !B.visited){
               triangle_position_vec_over_time_in_x_2.push_back(B.position[0]);
               triangle_position_vec_over_time_in_y_2.push_back(B.position[1]);
               B.visited = true;
           }

            if(B.index == 6 && !B.visited){
                triangle_position_vec_over_time_in_x_3.push_back(B.position[0]);
                triangle_position_vec_over_time_in_y_3.push_back(B.position[1]);
                B.visited = true;
            }


            if(A.index == 6 && !A.visited) {
           //     std::cout<<" true 6";

                position_vec_over_time_in_x_2.push_back(A.position[0]);
                position_vec_over_time_in_y_2.push_back(A.position[1]);
                A.visited = true;
            }


            if(A.index == 12 && !A.visited) {
            //    std::cout<<"true 12";

                position_vec_over_time_in_x_3.push_back(A.position[0]);
                position_vec_over_time_in_y_3.push_back(A.position[1]);
                A.visited = true;
            }

            if(B.index == 9 && !B.visited) {
             //   std::cout<<"true 9";

                rotator_position_vec_over_time_in_x_2.push_back(B.position[0]);
                rotator_position_vec_over_time_in_y_2.push_back(B.position[1]);
                B.visited = true;
            }


            if(A.index == 18 && !A.visited) {
            //    std::cout<<"true 18"<<std::endl;
                position_vec_over_time_in_x_4.push_back(A.position[0]);
                position_vec_over_time_in_y_4.push_back(A.position[1]);


                A.visited = true;
            }


            velocity_vector(A.index) = A.velocity[0];
            velocity_vector(A.index+1) = A.velocity[1];
            velocity_vector(A.index+2) = A.velocity[2];
            velocity_vector(B.index) = B.velocity[0];
            velocity_vector(B.index+1) = B.velocity[1];
            velocity_vector(B.index+2) = B.velocity[2];

        }
        reset_flags(connected_particles);


    }

    VectorXd get_current_EHD_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_vector.size();
        VectorXd current_EHD_force(n);
        current_EHD_force.setZero();

        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = current_positions[A.index];
            A_tmp.position[1] = current_positions[A.index+1];
            A_tmp.position[2] = current_positions[A.index+2];
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];
            B_tmp.position[2] = current_positions[B.index+2];

            // std::cout<< "in here the positions are:"<<std::endl<<A_tmp.position[0]<<std::endl<<A_tmp.position[1]<<std::endl<<B_tmp.position[0]<<std::endl<<B_tmp.position[1]<<std::endl;

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

            Vector3d local_EHD_force = get_my_EHD_force(particle_pair_tmp);
            current_EHD_force[A.index] += local_EHD_force[0];
            current_EHD_force[A.index+1] += local_EHD_force[1];
            current_EHD_force[A.index+2] += local_EHD_force[2];
            current_EHD_force[B.index] += local_EHD_force[0];
            current_EHD_force[B.index+1] += local_EHD_force[1];
            current_EHD_force[B.index+2] += local_EHD_force[2];

        }
      //  std::cout<<" current full ehd force = "<< current_EHD_force<<std::endl;
        return current_EHD_force/*.cwiseProduct(inverse_number_of_neighbours)*/;
    }

  /*  VectorXd get_current_brownian_motion(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        int n = spring_force_vector.size();
        VectorXd current_Bm(n);
        current_Bm.setZero();
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            Vector3d local_Bm_A = brownian_motion(std::get<0>(particle_pair));
            Vector3d local_Bm_B = brownian_motion(std::get<1>(particle_pair));
            current_Bm[A.index] = local_Bm_A[0];
            current_Bm[A.index+1] = local_Bm_A[1];
            current_Bm[A.index+2] = local_Bm_A[2];
            current_Bm[B.index] = local_Bm_B[0];
            current_Bm[B.index+1] = local_Bm_B[1];
            current_Bm[B.index+2] = local_Bm_B[2];

        }
        return current_Bm;
    }
*/
    //from thermal noise paper
    VectorXd get_current_thermal_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        int n = spring_force_vector.size();
        VectorXd current_thermal_force(n);
        current_thermal_force.setZero();
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            if(!A.visited){

                double x_rand = d(gen);
                double y_rand = d(gen);
                if(x_rand == y_rand){

                    std::cout<<"error occured, random number generator failed!"<<std::endl;
                }
                Vector3d local_tf_A = thermal_force(A,x_rand,y_rand);
                current_thermal_force[A.index] = local_tf_A[0];
                current_thermal_force[A.index+1] = local_tf_A[1];
                current_thermal_force[A.index+2] = local_tf_A[2];
                A.visited=true;
            }
            if(!B.visited){

                double x_rand_2 = d(gen);
                double y_rand_2 = d(gen);
                if(x_rand_2 == y_rand_2){
                    std::cout<<"error occured, random number generator failed!"<<std::endl;
                }
                Vector3d local_tf_B = thermal_force(B,x_rand_2,y_rand_2);
                current_thermal_force[B.index] = local_tf_B[0];
                current_thermal_force[B.index+1] = local_tf_B[1];
                current_thermal_force[B.index+2] = local_tf_B[2];
                B.visited=true;
            }
        }
        reset_flags(connected_particles);
        return current_thermal_force;
    }

    VectorXd get_current_rotational_thermal_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        int n = spring_force_vector.size();
        VectorXd current_rotational_thermal_force(n);
        current_rotational_thermal_force.setZero();
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            if(!A.visited){
                std::random_device rd{};
                std::mt19937 gen{rd()};
                std::normal_distribution<> d{0.0,1.0};
                double x_rand = d(gen);
                apply_rotational_noise(A,x_rand);
                A.visited=true;
            }
            if(!B.visited){
                std::random_device rd_2{};
                std::mt19937 gen_2{rd_2()};
                std::normal_distribution<> d_2{0.0,1.0};
                double x_rand_2 = d_2(gen_2);
                apply_rotational_noise(B,x_rand_2);
                B.visited=true;
            }

            std::vector<Vector3d> local_rotational_thermal_force = rotational_thermal_force(particle_pair);
            current_rotational_thermal_force[A.index] += local_rotational_thermal_force[0][0];
            current_rotational_thermal_force[A.index+1] += local_rotational_thermal_force[0][1];
            current_rotational_thermal_force[A.index+2] += local_rotational_thermal_force[0][2];
            current_rotational_thermal_force[B.index] += local_rotational_thermal_force[1][0];
            current_rotational_thermal_force[B.index+1] += local_rotational_thermal_force[1][1];
            current_rotational_thermal_force[B.index+2] += local_rotational_thermal_force[1][2];

        }
        reset_flags(connected_particles);
        return current_rotational_thermal_force;
    }

    MatrixXd get_current_EHD_jacobian(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_matrix_dx.rows();
        int  m = spring_force_matrix_dx.cols();

        MatrixXd EHD_force_jacobian(n,m);
        EHD_force_jacobian.setZero();

        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = current_positions[A.index];
            A_tmp.position[1] = current_positions[A.index+1];
            A_tmp.position[2] = current_positions[A.index+2];
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];
            B_tmp.position[2] = current_positions[B.index+2];

            //std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

            //Matrix3d current_spring_force_dx = NEW_get_damped_spring_force_dXA_without_damping(particle_pair_tmp);
            Matrix3d current_EHD_force_A_dx = get_EHD_force_dx(A_tmp,B_tmp);
            //Matrix3d current_EHD_force_B_dx = get_EHD_force_dx(B_tmp,A_tmp);
          /*  Matrix3d current_EHD_force_B_dx = -1.0 * current_EHD_force_A_dx;

            EHD_force_jacobian.block<3,3>(A.index,A.index) += current_EHD_force_A_dx;
            EHD_force_jacobian.block<3,3>(B.index,B.index) += current_EHD_force_B_dx;
            EHD_force_jacobian.block<3,3>(A.index,B.index) += -1.0 *current_EHD_force_A_dx;
            EHD_force_jacobian.block<3,3>(B.index,A.index) += -1.0 *current_EHD_force_B_dx;*/

            EHD_force_jacobian.block<3,3>(A.index,A.index) += current_EHD_force_A_dx;
            EHD_force_jacobian.block<3,3>(B.index,B.index) += (-1.0 *current_EHD_force_A_dx);
            EHD_force_jacobian.block<3,3>(A.index,B.index) += (-1.0 *current_EHD_force_A_dx);
            EHD_force_jacobian.block<3,3>(B.index,A.index) += current_EHD_force_A_dx;

        }
       // std::cout<<" EHD JAcobian = "<<EHD_force_jacobian<<std::endl;
        return EHD_force_jacobian;
    }

    void apply_wind_shadow(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        drag_occlusion_factor_vector = Eigen::VectorXd::Ones(size);

        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!((A.velocity[0] == 0.0) && (A.velocity[1] ==0.0))) {

                //A ocluded by B

                Vector2d vel_A_from_xa = {A.velocity[0] + A.position[0], A.velocity[1] + A.position[1]};
               // std::cout << "vel_A_from_xa " << vel_A_from_xa << std::endl;

                Vector2d line_vector;
                if(A.velocity[0] == 0.0){
                    line_vector = {1.0,0.0};
                }
                else if(A.velocity[1] == 0.0){
                    line_vector = {0.0,1.0};
                }
                else {
                    double orthogonal_y = 1.0;
                    double orthogonal_x = (-1.0 * A.velocity[1] * orthogonal_y) / A.velocity[0];
                    line_vector = {orthogonal_x, orthogonal_y};
                }


              //  std::cout << "line_vector" << line_vector << std::endl;
                Vector2d point_on_line_1 = {line_vector[0] + A.position[0], line_vector[1] + A.position[1]};
              //  std::cout << "point_on_line_1" << point_on_line_1 << std::endl;
                Vector2d point_on_line_2 = {A.position[0], A.position[1]};
              //  std::cout << "point_on_line_2" << point_on_line_2 << std::endl;
                double value_of_known_point =
                        (vel_A_from_xa[0] - point_on_line_1[0]) * (point_on_line_2[1] - point_on_line_1[1]) -
                        (vel_A_from_xa[1] - point_on_line_1[1]) * (point_on_line_2[0] - point_on_line_1[0]);
              //  std::cout << "value_of_known_point" << value_of_known_point << std::endl;
                int sign_of_known_point = std::copysign(1.0, value_of_known_point);
               // std::cout << "sign_of_known_point" << sign_of_known_point << std::endl;

                double value_of_other_center =
                        (B.position[0] - point_on_line_1[0]) * (point_on_line_2[1] - point_on_line_1[1]) -
                        (B.position[1] - point_on_line_1[1]) * (point_on_line_2[0] - point_on_line_1[0]);
               // std::cout << "value_of_other_center" << value_of_other_center << std::endl;
                int sign_of_other_center = std::copysign(1.0, value_of_other_center);
               // std::cout << "sign_of_other_center" << sign_of_other_center << std::endl;

                if (sign_of_known_point == sign_of_other_center) {
                  //  std::cout << "true" << std::endl;
                    Vector2d A_to_B = {B.position[0] - A.position[0], B.position[1] - A.position[1]};
                    //A_to_B.normalize();
                 //   std::cout << "A_to_B" << A_to_B << std::endl;
                    //line_vector.normalize();
                 //   std::cout<<"line vector = "<<line_vector<<std::endl;
                    double norm_prod = A_to_B.norm() * line_vector.norm();
                  //  std::cout << "norm_prod" << norm_prod << std::endl;
                  //  std::cout<<"thing inside acos = "<<(line_vector.dot(A_to_B)) / norm_prod<<std::endl;
                    double angle = std::acos((line_vector.dot(A_to_B)) / norm_prod); //should be bewteen 0 and pi
                  //  std::cout << "angle" << angle << std::endl;
                    double factor = std::sin(angle); //between 0 and 1
                  //  std::cout << "factor" << factor << std::endl;
                    double reduce_by = (1.0 - ((2.0 * B.radius * factor) / (2.0 * A.radius)));
                  //  std::cout << "reduce_by" << reduce_by << std::endl;
                    drag_occlusion_factor_vector[A.index] *= reduce_by;
                    drag_occlusion_factor_vector[A.index + 1] *= reduce_by;
                    drag_occlusion_factor_vector[A.index + 2] *= reduce_by;

                    /*Vector2d reflected_vector = A_to_B - 2.0 * (A_to_B.dot(line_vector.normalized())) * line_vector.normalized();
                    double alpha = std::pow(std::sin(2.0 * angle),2);
                    reflected_vector.normalize();
                    reflected_vector *= 6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius * A.velocity.norm();
                    drag_occlusion_factor_vector[A.index] *= (1.0 - alpha);
                    drag_occlusion_factor_vector[A.index+1] *= (1.0 - alpha);
                    drag_occlusion_factor_vector[A.index+2] *= (1.0 - alpha);*/

                    //drag_occlusion_direction_change_vector[A.index]
                }
            }
            if(!((B.velocity[0] == 0.0) && (B.velocity[1] == 0.0))) {

                Vector2d vel_B_from_xb = {B.velocity[0] + B.position[0], B.velocity[1] + B.position[1]};


                Vector2d line_vector_B;
                if(B.velocity[0] == 0.0){
                    line_vector_B = {1.0,0.0};
                }
                else if(B.velocity[1] == 0.0){
                    line_vector_B = {0.0,1.0};
                }
                else{
                    double orthogonal_y_B = 1.0;
                    double orthogonal_x_B = (-1.0 * B.velocity[1] * orthogonal_y_B) / B.velocity[0];
                    line_vector_B = {orthogonal_x_B, orthogonal_y_B};
                }

                Vector2d point_on_line_1_B = {line_vector_B[0] + B.position[0], line_vector_B[1] + B.position[1]};
                Vector2d point_on_line_2_B = {B.position[0], B.position[1]};
                double value_of_known_point_B =
                        (vel_B_from_xb[0] - point_on_line_1_B[0]) * (point_on_line_2_B[1] - point_on_line_1_B[1]) -
                        (vel_B_from_xb[1] - point_on_line_1_B[1]) * (point_on_line_2_B[0] - point_on_line_1_B[0]);
                int sign_of_known_point_B = std::copysign(1.0, value_of_known_point_B);

                double value_of_other_center_B =
                        (A.position[0] - point_on_line_1_B[0]) * (point_on_line_2_B[1] - point_on_line_1_B[1]) -
                        (A.position[1] - point_on_line_1_B[1]) * (point_on_line_2_B[0] - point_on_line_1_B[0]);
                int sign_of_other_center_B = std::copysign(1.0, value_of_other_center_B);

                if (sign_of_known_point_B == sign_of_other_center_B) {
                    Vector2d B_to_A = {A.position[0] - B.position[0], A.position[1] - B.position[1]};
                    B_to_A.normalize();
                    line_vector_B.normalize();
                    double norm_prod_B = B_to_A.norm() * line_vector_B.norm();
                    double angle_B = std::acos((line_vector_B.dot(B_to_A)) / norm_prod_B); //should be bewteen 0 and pi
                    double factor_B = std::sin(angle_B); //between 0 and 1
                    double reduce_by_B = (1.0 - ((2.0 * A.radius * factor_B) / (2.0 * B.radius)));
                    drag_occlusion_factor_vector[B.index] *= reduce_by_B;
                    drag_occlusion_factor_vector[B.index + 1] *= reduce_by_B;
                    drag_occlusion_factor_vector[B.index + 2] *= reduce_by_B;
                }

            }
        }
       // std::cout<<"drag_occlusion_factor_vector"<<drag_occlusion_factor_vector<<std::endl<<std::endl;

    }




    //CALCULATIONS FOR e, g AND g/dx

    double get_e(VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        VectorXd a = (x - 2.0 * position_vector + position_vector_minus_1);
        double x_t_M_x = 0.5 * a.transpose() * mass_matrix * a;

        double E_x =  get_current_potential_energy(x,connected_particles)+  safety *get_E_drag_indefinit_integral(x,connected_particles)+ safety * get_E_electrical(x) + safety *get_E_EHD(x,connected_particles) + safety * current_thermal_force.transpose() * x;

        return   x_t_M_x +  std::pow(time_step,2) * E_x; //check again
    }

    VectorXd get_g_with_drag(VectorXd x, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

        VectorXd F_x = get_current_spring_force(connected_particles,x) - safety *get_current_friction_force(connected_particles,x)/* + safety *electricalfield_force_vector */+ safety *get_current_EHD_force(connected_particles,x) + safety *current_thermal_force /*+ safety * current_rotational_thermal_force*/;



      //  std::cout<<"FS = \n"<<get_current_spring_force(connected_particles,x)<<" \n Fd ="<<get_current_friction_force(connected_particles,x)<<"\n Fel = "<<electricalfield_force_vector;

      //  std::cout<<"EHD FORCE ="<<get_current_EHD_force(connected_particles,x)<<std::endl;
        //  std::cout<<" F = "<<std::endl<<F_x<<std::endl;
      //   std::cout<<"mass matrix = "<<std::endl<<mass_matrix<<std::endl;
      //  std::cout<<" M * (x-y) = \n"<<mass_matrix * (x - 2.0 * position_vector + position_vector_minus_1)<<std::endl;

      //std::cout<<" Spring force = "<<get_current_spring_force(connected_particles,x)<<" drag force = "<<-get_current_friction_force(connected_particles,x)<<" electrical_force = "<<electricalfield_force_vector<<" EHD force = "<<get_current_EHD_force(connected_particles,x)<<" thermal force ="<<current_thermal_force<<std::endl;

        return mass_matrix * (x  - 2.0 * position_vector + position_vector_minus_1) - std::pow(time_step, 2) * F_x;
    }

    MatrixXd get_g_dx(VectorXd x, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        MatrixXd dF_dx = get_current_spring_force_jacobian(connected_particles,x) + safety  * -1.0 *get_current_friction_force_jacobian(connected_particles,x) + safety *get_current_electrical_force_jacobian(connected_particles,x) + safety * get_current_EHD_jacobian(connected_particles,x);
      //  std::cout<<" cuurent spring force jacobian = "<<std::endl<<dF_dx<<std::endl;
      //std::cout<<"dfdx = \n"<<dF_dx;

        MatrixXd result = mass_matrix - (std::pow(time_step,2) * dF_dx);
      //  std::cout<<"get g dx returns \n"<<result<<std::endl;
        return result;
    }

    double get_current_potential_energy(VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double potential_energy = 0.0;
        for(auto particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = x[A.index];
            A_tmp.position[1] = x[A.index+1];
            A_tmp.position[2] = x[A.index + 2];
            B_tmp.position[0] = x[B.index];
            B_tmp.position[1] = x[B.index+1];
            B_tmp.position[2] = x[B.index+2];

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

            Vector3d displacement =  std::get<0>(particle_pair_tmp).position - std::get<1>(particle_pair_tmp).position;
            double displcacement_norm = displacement.norm();
            //new
            double current_spring_lenght = displcacement_norm;
            potential_energy += 0.5 * stiffnes_constant * std::pow((current_spring_lenght-std::get<2>(particle_pair_tmp)),2);
        }
        return potential_energy;
    }

    //the one to use acording to Stelian
    double get_E_drag_indefinit_integral(VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double E = 0.0;
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            if(!A.visited){
                Vector3d new_pos;
                new_pos[0] = x[A.index];
                new_pos[1] = x[A.index+1];
                new_pos[2] = x[A.index+2];

                Vector3d old_pos;
                old_pos[0] = position_vector[A.index];
                old_pos[1] = position_vector[A.index+1];
                old_pos[2] = position_vector[A.index+2];

                double t1 = new_pos.dot(new_pos);

                double t2 = old_pos.dot(new_pos);


                E += (-6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius)/time_step * (t1/2.0 - t2);

                A.visited = true;
            }
            if(!B.visited){
                Vector3d new_pos;
                new_pos[0] = x[B.index];
                new_pos[1] = x[B.index+1];
                new_pos[2] = x[B.index+2];

                Vector3d old_pos;
                old_pos[0] = position_vector[B.index];
                old_pos[1] = position_vector[B.index+1];
                old_pos[2] = position_vector[B.index+2];

                double t1 = new_pos.dot(new_pos);

                double t2 = old_pos.dot(new_pos);


                E += (6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * B.radius)/time_step * (t1/2.0 - t2);

                B.visited = true;
            }
        }
        reset_flags(connected_particles);

        return E;
    }

    double get_E_electrical(VectorXd x){
        return electricalfield_force_vector.transpose() * x;
    }

    //not sure about this
    double get_E_EHD(VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        return get_current_EHD_force(connected_particles,x).transpose() * x;
    }







    //CALCULATION FOR SINGLE PARTICLE OR PAIR

    Vector3d NEW_get_damped_spring_force_without_damping(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        Vector3d xa = A.position;
        Vector3d xb = B.position;
        Vector3d xat = {position_vector[A.index], position_vector[A.index +1], position_vector[A.index +2]};
        Vector3d xbt = {position_vector[B.index], position_vector[B.index +1], position_vector[B.index +2]};
        Vector3d va = (xa - xat)/time_step;
        Vector3d vb = (xb - xbt)/time_step;

        Vector3d direction = (A.position - B.position)/(A.position - B.position).norm();

        double damping_term = damping_constant * (va -vb).dot(direction);
        double spring_term = stiffnes_constant * ( (A.position - B.position).norm() - rest_length);
        Vector3d result_ = -1.0 *(spring_term /*+ damping_term*/) * (direction);

       return result_;

    }


    Matrix3d NEW_get_damped_spring_force_dXA_without_damping(Particle myself, Particle other, double rest_length){

        Matrix3d Identity = MatrixXd::Identity(3,3);
        Particle A = myself;
        Particle B = other;
        Vector3d xa = A.position;
        Vector3d xb = B.position;
        Vector3d xat = {position_vector[A.index], position_vector[A.index +1], position_vector[A.index +2]};
        Vector3d xbt = {position_vector[B.index], position_vector[B.index +1], position_vector[B.index +2]};
        Vector3d va = (xa - xat)/time_step;
        Vector3d vb = (xb - xbt)/time_step;

        Vector3d t0 = xa - xb;

        double t1 = t0.norm();

        Matrix3d T2 = t0 * t0.transpose();

       double t3 = stiffnes_constant * (t1 - rest_length);

       double t2 = 1.0/time_step;
       Matrix3d T3 = T2;
       double t4 = std::pow(t1,2);
       Vector3d t5 = t2 * (xa - xat) - t2 * (xb - xbt);
       double t6 = damping_constant/t4;
       double t7 = t5.transpose() * t0;

         Matrix3d result = -1.0 * ((stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2 + (t3/t1) * Identity);

         Matrix3d damping_gradient = -1.0 * (damping_constant/(time_step * t4) * T3 - (2.0 * damping_constant * t7)/std::pow(t1,4) * T3 + t6 * t0 * t5.transpose() + t6 * t7 * Identity );


        return (result /*+ damping_gradient*/);
    }

    Vector3d get_stokes_friction_force(Particle particle){
        double Area = 4.0/3.0 * M_PI * std::pow(particle.radius,3);
        //prev
      //  std::cout<<" particle.velocity from particle "<<particle.index<<" = "<<particle.velocity<<std::endl;
         return 6.0 * M_PI * viscosity_multiplier*dynamic_viscosity * particle.radius * particle.velocity;


        //new
        //Vector3d squared_velocity = particle.velocity.array().square;
        //return 0.5 * density * squared_velocity * Area * drag_coefficient

    }

    Matrix3d get_friction_force_dv(Particle particle){ //should it be minus
        return Eigen::MatrixXd::Identity(3,3) * 6.0 * M_PI * viscosity_multiplier * dynamic_viscosity * particle.radius;
    }

    //electric stuff from particle and electric fiel paper
    Vector3d get_electrical_force(Particle particle){
        return get_coulomb_force(particle) /*+ get_image_force(particle) + get_gradient_force(particle)*/;
    }

    Vector3d get_coulomb_force(Particle particle){
        return get_electrical_field() * particle.charge;
    }

    Vector3d get_electrical_field(){
        Vector3d E_field;
        E_field[1] = 0.0;
        E_field[0] = 0.0;
        E_field[2] = (lower_electrode.voltage - upper_electrode.voltage)/(std::abs(lower_electrode.position[2] - upper_electrode.position[2]));
        return E_field;
    }

    Vector3d brownian_motion(Particle &particle){

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::normal_distribution<double> distribution (0.0,1.0);

        unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator2(seed2);

        double x_rand = distribution(generator);
        double y_rand = distribution(generator2);

        //from active brownian particles paper

        long double D = (kB * T* brownian_motion_multiplier)/(6 * M_PI * dynamic_viscosity * (particle.radius) ); //not sure

        long double k = std::sqrt(2 * D * time_step);


        double dx = x_rand * k;
        double dy = y_rand * k;

        Vector3d result;
        result[0] = dx;
        result[1] = dy;
        result[2] = 0.0;
        return result;
    }

    Vector3d thermal_force(Particle &particle,double x_rand, double y_rand){
       /* unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::normal_distribution<double> distribution (0.0,1.0);

        std::default_random_engine generator2;
        std::normal_distribution<double> distribution2(0.0,1.0);*/



        //from active brownian particles paper

        double constant = std::sqrt(2 * kB * T * 6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * particle.radius) * brownian_motion_multiplier;

        double alpha = 3.0 * M_PI * dynamic_viscosity * viscosity_multiplier * particle.radius* 2.0;
        //std::cout<<"alpha = "<<alpha<<std::endl;
        //std::cout<<"particle.mass = "<<particle.mass;
        double new_beta = alpha/particle.mass;
        //std::cout<<"new_beta = "<<new_beta<<std::endl;
        //double R_s = (std::pow(new_beta,2) * std::pow(time_step,2) + 2.0 * std::exp(-new_beta * time_step) - std::exp(-2.0*new_beta * time_step) -1)/(2.0 * new_beta * time_step + 4.0 * std::exp(-new_beta * time_step) - std::exp(-2.0 * new_beta * time_step) - 3.0);
        //std::cout<<"Rs = "<<R_s<<std::endl;
        double R_s = 1.5;
        double effective_mass = particle.mass + 0.5 * density * (M_PI * std::pow(particle.radius,3));
        //std::cout<<"effective_mass = "<<effective_mass<<std::endl;
        double intensity = R_s * (2.0 * kB * T * new_beta)/(effective_mass);
        //std::cout<<"intensity = "<<intensity<<std::endl;
       // double beta_newnew = 3.0 * M_PI * dynamic_viscosity * viscosity_multiplier * particle.radius * 2.0 /(particle.mass);
        //double intensity_2 = 2.0 * kB * T * beta_newnew / (M_PI * particle.mass);

       // double newest_constant = std::sqrt(M_PI * intensity_2/time_step);

        double new_constant = std::sqrt(M_PI * intensity/time_step) * brownian_motion_multiplier;
        //std::cout<<"new_constant = "<<new_constant<<std::endl;

        double new_force = new_constant * x_rand * particle.mass;
       // std::cout<<"new_force = "<<new_force<<std::endl;
       // std::cout<<" new force = "<<new_force<<std::endl;
       // std::cout<<" new constant = "<<new_constant<<std::endl;
       // std::cout<<" newest constant = "<<newest_constant<<std::endl;
       // std::cout<<" constant = "<<constant<<std::endl;
        Vector3d thermal_force = {new_constant * x_rand * particle.mass, new_constant * y_rand * particle.mass , 0.0};
        //std::cout<<" thermal force = "<<thermal_force<<std::endl;
        //std::cout<<"x rand inside  = "<<x_rand<<std::endl<<"y rand inside = "<<y_rand<<std::endl;
        return thermal_force;
    }

    std::vector<Vector3d> rotational_thermal_force(std::tuple <Particle&,Particle&,double> particle_pair){
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        double radi_sum = std::get<2>(particle_pair);
        Vector3d trig_angle_A = {std::cos(B.theta), std::sin(B.theta), 0.0};
        Vector3d F_rot_A = radi_sum * trig_angle_A * A.mass/std::pow(time_step,2);
        Vector3d trig_angle_B = {std::cos(A.theta), std::sin(A.theta), 0.0};
        Vector3d F_rot_B = radi_sum * trig_angle_B * B.mass/std::pow(time_step,2);
        std::vector<Vector3d> rot_result(2);
        rot_result[0] = F_rot_A;
        rot_result[1] = F_rot_B;
        return rot_result;
    }

    void apply_rotational_noise(Particle& A ,double theta_rand_1){
        double rotational_diffusion_constant_A = (kB * T)/(8.0 * M_PI * dynamic_viscosity * viscosity_multiplier * std::pow(A.radius,3));
       // std::cout<<" Dr = "<<rotational_diffusion_constant_A<<std::endl;
       // std::cout<<"A theta before = "<<A.theta<<std::endl;
        A.theta += std::sqrt(2.0 * rotational_diffusion_constant_A * time_step) * theta_rand_1;
       // std::cout<<" A theta middle = "<<A.theta<<std::endl;
        A.theta = fmod(A.theta,(2.0 * M_PI));
       // std::cout<<" A theta test = "<<A.theta<<std::endl;
        if(A.theta <0){
            A.theta += 2.0 * M_PI;
        }
       // std::cout<<" A theta after = "<<A.theta<<std::endl;





    }

    void add_brownian_motion(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(auto particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited) {
                A.position = A.position + (brownian_motion_multiplier * brownian_motion(A));
                A.visited = true;
            }
            if(!B.visited){
                B.position = B.position + (brownian_motion_multiplier * brownian_motion(B));
                B.visited = true;
            }
        }
        reset_flags(connected_particles);
    }

    Vector3d get_my_EHD_force_old(std::tuple <Particle&,Particle&,double> particle_pair){
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);

        double r_A;
        double r_B;
        double r_A_y;
        double r_B_y;

        Vector3d direction = get_amount(A,B);

        double x_amount = direction[0];
        double y_amount = direction[1];

        Vector3d displacement = A.position - B.position;
        r_B   = displacement[0];
        r_B_y = displacement[1];
        r_A   = -1.0 * displacement[0];
        r_A_y = -1.0 * displacement[1];


        double U_EHD_A = get_U_EHD(A,r_A) * x_amount;
        double U_EHD_A_y = get_U_EHD(A,r_A_y) * y_amount;
        double U_EHD_B = get_U_EHD(B,r_B) * x_amount;
        double U_EHD_B_y = get_U_EHD(B,r_B_y) * y_amount;
        Vector3d U_EHD_A_vec = {U_EHD_A, U_EHD_A_y,0.0};
        Vector3d U_EHD_B_vec = {U_EHD_B,U_EHD_B_y,0.0};

        //double U_test = (U_EHD_A * x_amount * A.radius + U_EHD_B * x_amount* B.radius) / (A.radius + B.radius);
        //double U_test_Y = (U_EHD_A_y * y_amount * A.radius + U_EHD_B_y * y_amount* B.radius) / (A.radius + B.radius);

        Vector3d F_EHD_A = {6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius * U_EHD_A,6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius * U_EHD_A_y, 0.0};
        Vector3d F_EHD_B = {6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * B.radius * U_EHD_B,6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * B.radius * U_EHD_B_y, 0.0};

        double radius_of_imaginary_sphere = A.radius + B.radius;
        double viscosity_constant = 6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * radius_of_imaginary_sphere;
        //return (U_EHD_A_vec + U_EHD_B_vec) * viscosity_constant;

        return (F_EHD_A + F_EHD_B);
    }

    Vector3d get_amount(Particle myself, Particle other){
        /* Vector3d direction = myself.position - other.position;
         for(int i = 0; i<direction.size();i++){
             direction[i] = std::abs(direction[i]);
         }
         double norm = direction.norm();
         for(int i = 0; i<direction.size();i++){
             direction[i] = direction[i]/norm;
         }
         double dprod = direction.dot(ex);
         double theta_angle = std::acos(dprod);
         double x_amount = (-2.0/M_PI) * theta_angle + 1.0;
         double y_amount = 1- x_amount;
         Vector3d amount = {x_amount, y_amount, 0.0};
         return amount;*/

        Vector3d direction = myself.position - other.position;
        double norm = direction.norm();
        Vector3d abs_direction = direction.array().abs();
        Vector3d norm_dir = abs_direction/norm;
        return norm_dir;
    }

    Vector3d get_my_EHD_force(std::tuple <Particle&,Particle&,double> particle_pair){
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);

        Vector3d direction_ = A.position - B.position;
        double norm = direction_.norm();
        Vector3d abs_direction = direction_.array().abs();
        Vector3d norm_dir = abs_direction/norm;

        Vector3d direction = norm_dir;

        // a value from 0 to 1 to determine the x -y force split
        double x_amount = direction[0];
        double y_amount = direction[1];

        double magnitude_A = get_U_EHD(A, A.radius + B.radius);
        double magnitude_B = get_U_EHD(B, A.radius + B.radius);

        Vector3d displacement = A.position - B.position;
        double sign_A_x = -1.0 * std::copysign(1.0,displacement[0]);
        double sign_A_y = -1.0 * std::copysign(1.0, displacement[1]);
        double sign_B_x = -1.0 * sign_A_x;
        double sign_B_y = -1.0 * sign_A_y;

        double U_EHD_A_x = magnitude_A * sign_A_x * x_amount;
        double U_EHD_A_y = magnitude_A * sign_A_y * y_amount;
        double U_EHD_B_x = magnitude_B * sign_B_x * x_amount;
        double U_EHD_B_y = magnitude_B * sign_B_y * y_amount;


        Vector3d F_EHD_A = {6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius * U_EHD_A_x,6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius * U_EHD_A_y, 0.0};
        Vector3d F_EHD_B = {6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * B.radius * U_EHD_B_x,6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * B.radius * U_EHD_B_y, 0.0};

        return F_EHD_A + F_EHD_B;
    }


    Matrix3d get_EHD_force_dx(Particle& myself, Particle& other){
        Vector3d xa = myself.position;
        Vector3d xb = other.position;
        double magnitude_A = get_U_EHD(myself, myself.radius + other.radius);
        double magnitude_B = get_U_EHD(other, myself.radius + other.radius);
        Vector3d magB = {magnitude_B,magnitude_B,0.0};
        Vector3d magA = {magnitude_A,magnitude_A,0.0};
        double rB = other.radius;
        double rA = myself.radius;
       /* double c1 = 6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * myself.radius;
        double c2 = ((beta_ehd * get_c(myself))/(dynamic_viscosity * viscosity_multiplier));
        Matrix3d amount_matrix;
        Vector3d amount = get_amount(myself,other);
        amount_matrix.col(0) = amount;
        amount_matrix.col(1) = amount;
        amount_matrix.col(2) = Eigen::VectorXd::Zero(3);
      //  std::cout<<" amount matrix = "<<amount_matrix<<std::endl;

        Matrix3d f_r_matrix;
        Vector3d f_r = get_f_r(myself,other);
        f_r_matrix.col(0) = f_r;
        f_r_matrix.col(1) = f_r;
        f_r_matrix.col(2) = Eigen::VectorXd::Zero(3);
       // std::cout<<" fr matrix = "<<f_r_matrix<<std::endl;
       // std::cout<<" c1 = "<<c1<<std::endl;
       // std::cout<<" c2 = "<<c2<<std::endl;
       // std::cout<<" fr dx = "<<get_f_rA_dxA(myself,other)<<std::endl;
       // std::cout<<" d amount dx = "<<get_amount_dx(myself,other)<<std::endl;


            return c1 * (c2 * amount_matrix.cwiseProduct(get_f_rA_dxA(myself,other))  + c2  * f_r_matrix.cwiseProduct(get_amount_dx(myself,other)));*/


       Vector3d t0 = xa -xb;
       Vector3d t1 = t0.array().sign();
       double t2 = t0.norm();
       Vector3d t3 = magB.cwiseProduct(t0);
       double t4 = 6.0 * M_PI * dynamic_viscosity * viscosity_multiplier;
       double t5 = std::pow(t2,3);
       Vector3d t6 = magA.cwiseProduct(t1);
       Vector3d t7 = t0.array().abs();
       Matrix3d mat1 = (t3.cwiseProduct(t1)).asDiagonal();
       Matrix3d mat2 = (t6.cwiseProduct(t1)).asDiagonal();

       Matrix3d dF_A_B_dA;
        dF_A_B_dA =   ((rB * t4)/t2) * mat1
                    - ((rB * t4)/t5) * (t3.cwiseProduct(t7)) * t0.transpose()
                    - ((rA * t4)/t2) * mat2
                    + ((rA * t4)/t5) * (t6.cwiseProduct(t7)) * t0.transpose();

        return dF_A_B_dA;


    }

    Matrix3d get_f_rA_dxA(Particle myself, Particle other){
       /* Matrix3d d_fr_rA_d_xA;
        Vector3d xa = myself.position;
        std::cout<<" xa = "<<xa<<std::endl;
        Vector3d xb = other.position;
        std::cout<<" xb = "<<xb<<std::endl;
        Vector3d t0 = xb - xa;
        std::cout<<" t0 = "<<t0<<std::endl;
        Vector3d o_o_z = {1.0,1.0,0.0};
        std::cout<<" ooz = "<<o_o_z<<std::endl;
        Vector3d t1 = o_o_z + 1.0/myself.radius  * t0;
        std::cout<<" t1 = "<<t1<<std::endl;
        Vector3d t2 = t1.array().square();
        std::cout<<" t2 = "<<t2<<std::endl;
        std::cout<<" middle term = "<<t0.cwiseProduct(t2)<<std::endl;
        std::cout<<" after. array()"<<(t0.cwiseProduct(t2)).array()<<std::endl;
        std::cout<<" after pow = "<<(t0.cwiseProduct(t2)).array().pow(-(1.0 + 5.0/2.0))<<std::endl;
        Vector3d t3 = (t0.cwiseProduct(t2)).array().pow(-(1.0 + 5.0/2.0));
        Vector3d alternative_t3 = {std::pow((t0.cwiseProduct(t2)).array()[0], -(1.0 + 5.0/2.0)), std::pow((t0.cwiseProduct(t2)).array()[1], -(1.0 + 5.0/2.0)), 0.0};
        std::cout<<" alt t3 = "<<alternative_t3<<std::endl;
        std::cout<<" t3 = "<<t3<<std::endl;
        Vector3d term1 = t3.cwiseProduct(t2);
        std::cout<<" term1 = "<<term1<<std::endl;
        Vector3d term2 = (t3.cwiseProduct(t0)).cwiseProduct(t1);
        std::cout<<" term2 = "<<term2<<std::endl;
        Matrix3d mat1 = term1.asDiagonal();
        std::cout<<" mat1 = "<<mat1<<std::endl;
        Matrix3d mat2 = term2.asDiagonal();
        std::cout<<" mat2 = "<<mat2<<std::endl;
        Matrix3d result = ((15.0/(4.0 * myself.radius)) * mat1) + ((30.0/(std::pow(myself.radius,2) * 4.0)) * mat2);
        std::cout<<" result = "<<result<<std::endl;
        return result;*/

        Matrix3d d_fr_rA_d_xA;
        Vector3d xa = myself.position;
      //  std::cout<<" xa = "<<xa<<std::endl;
        Vector3d xb = other.position;
      //  std::cout<<" xb = "<<xb<<std::endl;
        Vector3d t0 = xb - xa;
     //   std::cout<<" t0 = "<<t0<<std::endl;
        Vector3d o_o_z = {1.0,1.0,0.0};
     //   std::cout<<" ooz = "<<o_o_z<<std::endl;
        Vector3d rhs = (1.0/myself.radius  * t0).array().square();
        Vector3d t1 = o_o_z + rhs;
      //  std::cout<<" t1 = "<<t1<<std::endl;

        Vector3d  middle_term = {std::pow(t1[0], -(1.0 + 5.0/2.0)),std::pow(t1[1], -(1.0 + 5.0/2.0)),0.0};;
     //   std::cout<<"middl term = "<<middle_term<<std::endl;

        Matrix3d mat1 = ( (t0.cwiseProduct(middle_term)).cwiseProduct(t0) ).asDiagonal();
     //   std::cout<<"mat1 = "<<mat1<<std::endl;
        Vector3d term =  {std::pow(t1[0], -(5.0/2.0)),std::pow(t1[1], -(5.0/2.0)),0.0};
     //   std::cout<<"term = "<<term<<std::endl;
        Matrix3d mat2 = term.asDiagonal();
      //  std::cout<<"mat2= "<<mat2<<std::endl;
        Matrix3d result = 30.0/(std::pow(myself.radius,3) * 4) * mat1 - (3.0/(2.0 * myself.radius)) * mat2;

      //  std::cout<<" result = "<<result<<std::endl;
        return result;
    }

    //dfr_A_dxb = - dfr_A_dxA = - get_f_r_A_dxA(A,B)
    //dfr_B_dxb = get_f_r_A_dxA(B,A)
    //dfr_B_dxA = - dfr_B_dxb



    Matrix3d get_amount_dx(Particle myself, Particle other){
        Vector3d e_x = {1.0,0.0,0.0};
        Vector3d xa = myself.position;
        Vector3d xb = other.position;
        Vector3d t0 = xa - xb;
        double t1 = t0.norm();
     //   std::cout<<" t1= "<<t1<<std::endl;
       /* Vector3d t2 = {std::abs(t0[0]), std::abs(t0[1]),0.0};
        std::cout<<"t2 = "<<t2<<std::endl;
        double tmp = t2[0];
        double t3 = ex.transpose() * t2;
        std::cout<<" t3= "<<t3<<std::endl;
        double t4 = -0.5;
        double test2 = (t2.dot(ex))/t1;
        std::cout<<"test2="<<test2<<std::endl;
        Vector3d sign_of_t0 = {std::copysign(1.0,t0[0]), std::copysign(1.0,t0[1]), 0.0};
        std::cout<<" sign of t0="<<sign_of_t0<<std::endl;
        Vector3d d_amount_x_d_xA = (-2.0 * std::pow((1.0 - std::pow(t3/t1,2)),t4) * t3)/(std::pow(t1,3) * M_PI) * t0
                + (( std::pow( 1.0 - std::pow( (( t2.dot(e_x))/t1) ,2.0) ,t4 ) ) * 2.0 )/(M_PI * t1) * e_x.cwiseProduct(sign_of_t0);
        std::cout<<"damount = "<<d_amount_x_d_xA<<std::endl;
        std::cout<<" t3/t1 = "<<t3/t1<<std::endl;
        std::cout<<" first term="<<(-2.0 * std::pow((1.0 - std::pow(t3/t1,2)),t4) * t3)/(std::pow(t1,3)* M_PI)<<std::endl;
        std::cout<<" second term = "<<(( std::pow( 1.0 - std::pow( (( t2.dot(e_x))/t1) ,2.0) ,t4 ) ) * 2.0 )/(M_PI * t1) * e_x.cwiseProduct(sign_of_t0)<<std::endl;
        std::cout<<" first power = "<<std::pow((1.0 - std::pow(t3/t1,2)),t4)<<std::endl;
        std::cout<<" second power ="<<std::pow(t1,3)<<std::endl;
        std::cout<<" third power = "<<std::pow( 1.0 - std::pow( (( t2.dot(e_x))/t1) ,2.0) ,t4 )<<std::endl;
        std::cout<<" cwise = "<<e_x.cwiseProduct(sign_of_t0)<<std::endl;
        */


       /* Matrix3d amount_dx;
        amount_dx.row(0) = d_amount_x_d_xA;
        amount_dx.row(1) = -1.0 *d_amount_x_d_xA;
        amount_dx.row(2) = Eigen::VectorXd::Zero(3);
        return amount_dx;*/
        Vector3d sign_t0 = {std::copysign(1.0,t0[0]), std::copysign(1.0,t0[1]), 0.0};
        Matrix3d mat1 = sign_t0.asDiagonal();
        Vector3d abs_t0 = {std::abs(t0[0]),std::abs(t0[1]),0.0};
       return (1.0/t1) * mat1 - (1.0/std::pow(t1,3)) * abs_t0 * t0.transpose();
    }
    //amountdxA = get_amount_dx(A,B)
    //amount dXb = - amountdXA = get_amount_dx(B,A)





    double get_c(Particle particle){
        double w_bar = (lower_electrode.frequency * get_H()) /( (1.0 /get_debye_length()) * get_ion_diffusivity());
        // std::cout<<"W bar = "<<w_bar<<std::endl;

        double w_bar_2 = std::pow(w_bar,2);
        return  relative_permittivity * vacuum_permittivity *(1.0/get_debye_length()) * get_H() * std::pow((lower_electrode.peak_voltage/(2.0 * get_H())),2)
                    * ((get_K1(particle) + w_bar * get_K2(particle))/(1.0 + w_bar_2) );
    }

    Vector3d get_f_r(Particle myself, Particle other){
        Vector3d ones = Eigen::VectorXd::Ones(3);
        Vector3d term = ((other.position - myself.position)/(myself.radius)).array().square();
        Vector3d term2 = (ones + term).array().pow(-5.0/2.0);
        return (3.0/(2.0 * myself.radius)) * (other.position - myself.position).cwiseProduct(term2);
    }

    double get_U_EHD(Particle particle,double r){ //deleted vacuum permittivity

        //physrevlet paper
         double w_bar = (lower_electrode.frequency * get_H()) /( (1.0 /get_debye_length()) * get_ion_diffusivity());
        // std::cout<<"W bar = "<<w_bar<<std::endl;

         double w_bar_2 = std::pow(w_bar,2);
         double C =  relative_permittivity * vacuum_permittivity *(1.0/get_debye_length()) * get_H() * std::pow((lower_electrode.peak_voltage/(2.0 * get_H())),2)
                 * ((get_K1(particle) + w_bar * get_K2(particle))/(1.0 + w_bar_2) );
        // std::cout<<" C = "<<C<<std::endl;
         double f_r = (3.0 * (r/particle.radius))
                 / (  2.0 * std::pow((1.0 + std::pow((r/particle.radius),2) ),(5.0/2.0))   );
        // std::cout<<" U EHD for particle "<<particle.index<<" with radius "<<particle.radius<<" = "<<beta * (C / (dynamic_viscosity* viscosity_multiplier)) * f_r<<std::endl;
       // std::cout<<"fr ="<<f_r<<std::endl;
       // std::cout<<" beta ehd  ="<<beta_ehd<<std::endl;
      //  std::cout<<" U ehd = "<<beta_ehd * (C / (dynamic_viscosity* viscosity_multiplier)) * f_r<<std::endl;
         return  beta_ehd * (C / (dynamic_viscosity* viscosity_multiplier)) * f_r;

        /* double alpha = (lower_electrode.frequency * get_H()) /( (1.0 /get_debye_length()) * get_ion_diffusivity());
         double alpha_2 = std::pow(alpha,2);
         double Vp = lower_electrode.peak_voltage;
         double H = get_H();
         double C = relative_permittivity * vacuum_permittivity * std::pow( (Vp/(2.0 * H)) ,2) *
                 (alpha_2/(1.0 + alpha_2)) *
                 (1.0/get_debye_length()) * get_ion_diffusivity() /lower_electrode.frequency;
         double up = 3.0 * (r/particle.radius);
         double down = 2.0 * std::pow(1.0 + std::pow((r/particle.radius),2),(5.0/2.0));
         double f_r = up/down;
         return (C * get_K2(particle) / (dynamic_viscosity * viscosity_multiplier) ) * f_r; */
    }

    double get_H(){
       // std::cout<<"H = "<<0.5 * std::abs(lower_electrode.position[2] - upper_electrode.position[2])<<std::endl;
        return 0.5 * std::abs(lower_electrode.position[2] - upper_electrode.position[2]);
   }

    double get_omega_bar(Particle particle){
        return (lower_electrode.frequency * get_H()) /(get_debye_length() * get_ion_diffusivity());
    }


    double get_ion_diffusivity(){
      //  std::cout<<" D = "<<(kB * T)/( 6.0 * M_PI * (dynamic_viscosity * viscosity_multiplier) * particle.radius)<<std::endl;
        //return (kB * T)/( 6.0 * M_PI * (dynamic_viscosity * viscosity_multiplier) * particle.radius);
       // return 1.47 * std::pow(10,-9);
       return 0.5 * (D_K_plus + D_Cl_minus);
    }
    double get_ionic_strength(){
        return 0.5 * (molarity * std::pow(K_plus_charge,2) + molarity * std::pow(CL_minus_charge,2));
    }


    double get_debye_length(){
        // std::cout<<" debye length = "<<2.0 * particle.charge * std::sqrt((M_PI * particle.density)/(kB * T))<<std::endl;
        //  return 2.0 * particle.charge * std::sqrt((M_PI * particle.density)/(kB * T));
      //  std::cout<<" debye length = "<< std::sqrt((vacuum_permittivity * relative_permittivity * kB * T)/(2.0 * NA * std::pow(elementary_charge,2) * get_ionic_strength()))<<std::endl;
      //  std::cout<<" new debye length with higher concentration = "<< std::sqrt((vacuum_permittivity * relative_permittivity * kB * T)/(2.0 * NA * std::pow(elementary_charge,2) * get_ionic_strength())) * std::pow(10,5)<<std::endl;
        return std::sqrt((vacuum_permittivity * relative_permittivity * kB * T)/(2.0 * NA * std::pow(elementary_charge,2) * get_ionic_strength()));
       // return 0.6 * std::pow(10,-9);
        // return 150 * std::pow(10,-9);
    }

    //from conductivity paper
    double get_diffuse_layer_conductance(Particle particle){
        double m = std::pow((gas_constant * T) / faraday_constant, 2) * (2.0 * relative_permittivity * vacuum_permittivity)/(3.0 * dynamic_viscosity * get_ion_diffusivity());
      //  std::cout<<"m= "<<m<<std::endl;
        double t1 = (4.0 * std::pow(faraday_constant ,2)* molarity * get_ion_diffusivity() * (1.0 + 3.0 * m))/(gas_constant * T * (1.0/get_debye_length()));
      //  std::cout<<"t1= "<<t1<<std::endl;
      //  std::cout<<" cosh = "<<cosh((elementary_charge * particle.zeta_potential * std::pow(10,-3))/(2.0 * kB * T))<<std::endl;
        double t2 = cosh((elementary_charge * particle.zeta_potential* std::pow(10,-3))/(2.0 * kB * T)) - 1.0;
       // std::cout<<" zeta = "<<particle.zeta_potential<<std::endl;
      //  std::cout<<"t2= "<<t2<<std::endl;

        return t1 * t2;

    }

    double get_conductivity(Particle particle){
     //   std::cout<<" get conductivity = "<<(particle.stern_layer_conductance / particle.radius) + 2.0 * (get_diffuse_layer_conductance(particle)/particle.radius)<<std::endl;
     double diffusive_layer_conductance = get_diffuse_layer_conductance(particle);
     double stern_layer_conductance = particle.stern_layer_conductance;
     double particle_radius = particle.radius;
    // std::cout<<" diffuse layer conductance = "<<diffusive_layer_conductance<<std::endl;
   //  std::cout<<" stern layer conductance = "<<stern_layer_conductance<<std::endl;
   //  std::cout<<" mass = "<<particle.mass<<std::endl;
   //  std::cout<<" conductivity w d = "<<2.0 *  (particle.stern_layer_conductance / particle.radius) + 2.0 * (get_diffuse_layer_conductance(particle)/particle.radius)<<std::endl;
   //  std::cout<<" conductivity w.o d = "<<2.0 *  (particle.stern_layer_conductance / particle.radius)<<std::endl;
        return 2.0 *  (particle.stern_layer_conductance / particle.radius) + 2.0 * (get_diffuse_layer_conductance(particle)/particle.radius);
    }

    double get_Eo(){
        return (lower_electrode.peak_voltage/std::sqrt(2.0))/(2.0 * get_H());
    }




    // bith from clausius mossotti paper
    double get_K1(Particle particle){
        //double w_bar = (lower_electrode.frequency * get_H()) /( (1.0 /get_debye_length()) * get_ion_diffusivity());
      //  std::cout<<"W bar = "<<w_bar<<std::endl;
        double  w = lower_electrode.frequency;
        double w_2 = std::pow(w,2);
       // double w_2 = w_bar_2;
        double sigma_m = conductivity;
        double sigma_p = get_conductivity(particle);
        double epsilon_p = particle.permittivity * vacuum_permittivity;
        double epsilon_m = relative_permittivity* vacuum_permittivity;
        double up = w_2 * (epsilon_p - epsilon_m) * (epsilon_p + 2.0 * epsilon_m) + (sigma_p - sigma_m) * (sigma_p + 2.0 * sigma_m);
        double down = w_2 * std::pow((epsilon_p + 2.0 * epsilon_m),2) + std::pow((sigma_p+ 2.0 * sigma_m),2);
    //    std::cout<<" K1 = "<<up/down<<std::endl;
        return up/down;


    }

    double get_K2(Particle particle){
       // double w_bar = (lower_electrode.frequency * get_H()) /( (1.0 /get_debye_length()) * get_ion_diffusivity());
      //  std::cout<<"W bar = "<<w_bar<<std::endl;

       // double w_bar_2 = std::pow(w_bar,2);
       // double w_2 = w_bar_2;
       // double w = w_bar; // ?????;
       double w = lower_electrode.frequency;
       double w_2 = std::pow(w,2);

        double sigma_m = conductivity;
        double sigma_p = get_conductivity(particle);
        double epsilon_p = particle.permittivity* vacuum_permittivity;
        double epsilon_m = relative_permittivity* vacuum_permittivity;
        double up = w * (sigma_m - sigma_p) * (epsilon_p + 2.0*epsilon_m) - (epsilon_p - epsilon_m) * (sigma_p + 2.0 * sigma_m);
        double down = w_2 * std::pow((epsilon_p + 2.0 * epsilon_m),2) + std::pow((sigma_p + 2.0 * sigma_m),2);
        //return ( w * (sigma_m - sigma_p)*(epsilon_p + 2.0 * epsilon_m) - (epsilon_p - epsilon_m)*(sigma_p + 2.0 * sigma_m)   )/(w_2 * std::pow((epsilon_p + 2.0 * epsilon_m),2) + std::pow((sigma_p + 2.0 * sigma_m),2));
     //   std::cout<<" K2 = "<<up/down<<std::endl;
        return up/down;
    }










    //EXPERIMENTAL STUFF FOR LATER

  /*  void check_Boundaries(Particle &particle){

            if(particle.position[0] + particle.radius >= Boxcenter[0] + Boxsize[0]/2 - xbuffer){
                particle.velocity[0] = -1 *particle.velocity[0];
            }

            if(particle.position[1]+ particle.radius >= Boxcenter[1] + Boxsize[1]/2 - ybuffer){
                particle.velocity[1] = -1 *particle.velocity[1];
            }

            if(particle.position[0]- particle.radius <= Boxcenter[0] - Boxsize[0]/2 + xbuffer){
                particle.velocity[0] = -1 *particle.velocity[0];
            }

            if(particle.position[1]- particle.radius <= Boxcenter[1] - Boxsize[1]/2 + ybuffer){
                particle.velocity[1] = -1 *particle.
                [1];
            }

    }

    void remove_connected_element(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, Particle& particle_to_be_removed){
        for(auto it = std::begin(connected_particles); it != std::end(connected_particles); ++it){
            if(std::get<0>(*it).index == particle_to_be_removed.index ||  std::get<1>(*it).index == particle_to_be_removed.index){
                connected_particles.erase(it);
            }
        }
    }

    void dissolve(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles,double dissolve_rate){
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited){
                if(A.dissolve){
                    A.dissolve_particle(dissolve_rate);
                    if(A.radius < std::pow(10,-1)){
                        remove_connected_element(connected_particles,A);

                    }
                }
                A.visited = true;
            }
            if(!B.visited){
                if(B.dissolve){
                    B.dissolve_particle(dissolve_rate);
                    if(B.radius < std::pow(10,-2)){
                        remove_connected_element(connected_particles,B);

                    }
                }
                B.visited = true;
            }
        }
        reset_flags(connected_particles);
    }

*/



    //CHECKING GRADIENT WITH FINITE DIFFERENCE STUFF (DOES NOT WORK PROPERLY BUT NOT A PRIORITY RIGHT NOW)

    //only works for a system of two particles
    double get_dE_dx (std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd x_old, VectorXd x_new){
        std::vector<std::tuple <Particle&,Particle&,double> > connected_particles_new;
        for(auto particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            double rest_length = std::get<2>(particle_pair);
            Particle A_new = A;
            A_new.position[0] = x_new(A.index);
            Particle B_new = B;
            B_new.position[0] = x_new(B.index);
            connected_particles_new.push_back(std::tuple<Particle&,Particle&,double>(A_new,B_new,rest_length));
        }
        double E_at_x = get_Potential_Energy(connected_particles);
        double E_at_x_plus_delta_x = get_Potential_Energy(connected_particles_new);
        double delta_x = (std::abs(x_new[0] - x_new[2]) - std::get<2>(connected_particles_new[0])) - (std::abs(x_old[0] - x_old[2]) - std::get<2>(connected_particles[0]));

        return (E_at_x_plus_delta_x - E_at_x)/delta_x;
    }

    double get_dFx_dx_FD(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles_new, VectorXd x_old, VectorXd x_new, VectorXd x_new_new){
        double second =  get_dE_dx(connected_particles,x_old,x_new);
        double first = get_dE_dx(connected_particles_new,x_new,x_new_new);
        double delta_x = (std::abs(x_new[0] - x_new[2]) - std::get<2>(connected_particles_new[0])) - (std::abs(x_old[0] - x_old[2]) - std::get<2>(connected_particles[0]));

        return (first - second)/delta_x;
    }

    //not generic only for testing
    double get_F_x_FD(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        return spring_force_vector[0];
    }

    double get_F_x_d_x(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double result = spring_force_matrix_dx(0,0);
        // std::cout<<std::endl<<"Fdx  HERE IS"<<std::endl<<result<<std::endl;
        return result;
    }

    double get_Potential_Energy(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double potential_energy = 0;
        for(auto particle_pair : connected_particles) {
            Vector3d displacement =  std::get<0>(particle_pair).position - std::get<1>(particle_pair).position;
            double displcacement_norm = displacement.norm();
            potential_energy += 0.5 * stiffnes_constant * std::pow((displcacement_norm-std::get<2>(particle_pair)),2);
        }
        return potential_energy;
    }

    double get_Kinetic_Energy(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double kinetic_energy = 0;
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited_ekin) {
                kinetic_energy += 0.5 * A.mass * (A.velocity).transpose() * (A.velocity);
                A.visited_ekin = true;
            }
            if(!B.visited_ekin) {
                kinetic_energy += 0.5 * B.mass * (B.velocity).transpose() * (B.velocity);
                B.visited_ekin = true;
            }
        }
        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            A.visited_ekin = false;
            B.visited_ekin = false;
        }
        return kinetic_energy;
    }




};