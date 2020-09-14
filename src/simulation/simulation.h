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

class Simulation{
public:
    double mscale = std::pow(10,-6);
    Electrode lower_electrode = Electrode(1.0,{-100.0* mscale,-500.0* mscale}, 10.0,1000.0,240.0,40,240);
    Electrode upper_electrode = Electrode(-1.0,{100.0* mscale,500.0* mscale},10.0,1000.0,0.0,40,240);
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
    double micro_time_step = time_step * std::pow(10,6);


    // all units in micro scale, (10^-6)

    double eta = 1.0 * std::pow(10,-3); //viscosity of water in PA * seconds
    double dynamic_viscosity =  8.9 * std::pow(10,-4); //dynamic vsicosity of water in PA * s
    double viscosity_multiplier = 1.0;

    double vacuum_permittivity = 8.8541878176 * std::pow(10,-12);
    double relative_permittivity = 80.2; //water at room temperatur

    double conductivity = 0.05 * std::pow(10,-6);


    double T = 293; //room temperatur
    double kB = 1.38 * std::pow(10,-23); //Boltzmann constant
    Vector2d force_damper;

    double noise_damper = 4;
    double drag_damper = 0;
    double brownian_motion_multiplier = 0.0;

    Vector2d Boxcenter = {0,0};
    Vector2d Boxsize = {1500,1500};
    double xbuffer, ybuffer = 5.0;

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

    std::vector<double> e_x;
    std::vector<double> e_y;
    std::vector<double> e_z;

    std::vector<double> friction_force_over_time_x;

    bool test_bool = true;

    int max_iterations = 1000;

    int size = -1;



    Simulation() {}






    //INTERFACES AND NETWONS METHOD

    void run_simulation_for_connected_Particles(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(int i = 0; i<num_steps;i++) {
            total_energy.push_back(0.0);
            kinetic_energy.push_back(0.0);
            potential_energy.push_back(0.0);

           UPDATE_SYSTEM_NEW(connected_particles);


        }
    }

    void UPDATE_SYSTEM_NEW(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

        fill_spring_force_and_derivatives_NEW(connected_particles);

        first_itteration = false;

        VectorXd x_t_plus_1_init = position_vector + (position_vector - position_vector_minus_1);
        VectorXd x_t_plus_1_old = x_t_plus_1_init;

        //int maxiter, double beta, double tol, VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles
        THE_FINAL_NEWTON_with_drag(max_iterations,0.5,std::pow(10,(-9)), x_t_plus_1_old, connected_particles);

        //Newtons_Method_NEW_with_force_update(100000,std::pow(10,(-8)), 0.25, 0.5,  x_t_plus_1_old,connected_particles);


        VectorXd tmp = position_vector_minus_1;

       // std::cout<<" Position at t-1= \n"<<tmp<<"\n Positoin at t = \n"<<position_vector<<"\n Position at t+1 = \n"<<x_t_plus_1_old<<"\n \n \n";
        energy_FD.push_back(get_dE_dx(connected_particles,position_vector,x_t_plus_1_old));
        forces_FD.push_back(get_F_x_FD(connected_particles));
        force_jacobian_FD.push_back(get_F_x_d_x(connected_particles));

        double current_epot = get_Potential_Energy(connected_particles);
        double current_ekin = get_Kinetic_Energy(connected_particles);
        potential_energy.push_back(current_epot);
        kinetic_energy.push_back(current_ekin);
        total_energy.push_back(current_epot + current_ekin);

        VectorXd tmp_tmp = position_vector;
        std::vector<std::tuple <Particle&,Particle&,double> > tmp_connected_particles = connected_particles;
        update_positions_NEW(connected_particles,x_t_plus_1_old);
        dFx_dx_FD.push_back( get_dFx_dx_FD(tmp_connected_particles, connected_particles, tmp, tmp_tmp, x_t_plus_1_old));

        double e = get_e(x_t_plus_1_old,connected_particles);

       // std::cout<<" e = "<<std::endl<<e<<std::endl;
        e_vec.push_back(e);

        add_brownian_motion(connected_particles);

        update_rest_length(connected_particles);

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
                A.velocity[1] = (x_t_plus_1_new(A.index+1) -A.position[1])/time_step;
                A.position[0] = x_t_plus_1_new(A.index);
                A.position[1] = x_t_plus_1_new(A.index+1);
                A.visited=true;
            }
            if(!B.visited){
                B.velocity[0] = (x_t_plus_1_new(B.index) - B.position[0])/time_step;
                B.velocity[1] = (x_t_plus_1_new(B.index+1) - B.position[1])/time_step;
                B.position[0] = x_t_plus_1_new(B.index);
                B.position[1] = x_t_plus_1_new(B.index+1);
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
                A.visited = true;
            }
            if(!B.visited){
                position_vector_minus_1(B.index) = B.position[0];
                position_vector_minus_1(B.index+1) = B.position[1];
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

            while(get_g_with_drag(x_i_plus_1_tmp,connected_particles).norm() > get_g_with_drag(x_i,connected_particles).norm()){
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
                count+=2;
            }
            if(!B.visited){
                B.visited = true;
                B.index = count;
                count+=2;
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
        reset_flags(connected_particles);
    }

    void connect(Particle &A, Particle &B, double restlength, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        std::tuple<Particle&,Particle&,double> particle_pair(A,B,restlength);
        connected_particles.push_back(particle_pair);
    }

    void update_rest_length(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(auto& particle_pair : connected_particles){
            std::get<2>(particle_pair) = set_rest_length_rad_plus_rad(std::get<0>(particle_pair), std::get<1>(particle_pair));
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
        for(auto particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited) {
                A.position[0] = initial_position_vector(A.index);
                A.position[1] = initial_position_vector(A.index+1);
                B.velocity = B.initial_velocity;
                position_vector(A.index) = A.position[0];
                position_vector(A.index+1) = A.position[1];
                A.visited = true;
            }
            if(!B.visited){
                B.position[0] = initial_position_vector(B.index);
                B.position[1] = initial_position_vector(B.index+1);
                B.velocity = B.initial_velocity;
                position_vector(B.index) = B.position[0];
                position_vector(B.index+1) = B.position[1];
                B.visited = true;
            }
        }
        reset_flags(connected_particles);
        position_vector_minus_1 = position_vector;

    }

    double set_rest_length_rad_plus_rad(Particle& A, Particle& B){
        return A.radius + B.radius;
    }








    //CALCULATION FOR GLOBAL STUFF DEPENDANT ON CURRENT POSITION TO USE INSIDE NEWTONS METHOD

    MatrixXd get_current_spring_force_jacobian(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_matrix_dx.rows();
        int  m = spring_force_matrix_dx.cols();

        MatrixXd spring_force_jacobian(n,m);
        spring_force_jacobian.setZero();

        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = current_positions[A.index];
            A_tmp.position[1] = current_positions[A.index+1];
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

            Matrix2d current_spring_force_dx = NEW_get_damped_spring_force_dXA_without_damping(particle_pair_tmp);

            spring_force_jacobian.block<2,2>(A.index,A.index) += current_spring_force_dx;
            spring_force_jacobian.block<2,2>(B.index,B.index) += current_spring_force_dx;
            spring_force_jacobian.block<2,2>(A.index,B.index) += -1.0 *current_spring_force_dx;
            spring_force_jacobian.block<2,2>(B.index,A.index) += -1.0 *current_spring_force_dx;

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
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];

           // std::cout<< "in here the positions are:"<<std::endl<<A_tmp.position[0]<<std::endl<<A_tmp.position[1]<<std::endl<<B_tmp.position[0]<<std::endl<<B_tmp.position[1]<<std::endl;

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

           Vector2d current_spring_force_tmp = NEW_get_damped_spring_force_without_damping(particle_pair_tmp);
           //std::cout<<"current spring force  tmp = "<<std::endl<<current_spring_force_tmp<<std::endl;

            current_spring_force[A_tmp.index] += current_spring_force_tmp[0];
            current_spring_force[A_tmp.index + 1] += current_spring_force_tmp[1];
            current_spring_force[B_tmp.index] += -1.0 * current_spring_force_tmp[0];
            current_spring_force[B_tmp.index + 1] += -1.0 *current_spring_force_tmp[1];



        }
       // std::cout<<"current spring force  = "<<std::endl<<current_spring_force<<std::endl;

        return current_spring_force;
    }

    VectorXd get_current_friction_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = friction_vector.size();
        // std::cout<<" current positions = "<<std::endl<<current_positions<<std::endl;


        VectorXd current_friction_force(n);
        current_friction_force.setZero();

        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Particle A_tmp = A;
            Particle B_tmp = B;
            A_tmp.position[0] = current_positions[A.index];
            A_tmp.position[1] = current_positions[A.index+1];
            B_tmp.position[0] = current_positions[B.index];
            B_tmp.position[1] = current_positions[B.index+1];

            Vector2d pos_at_A;
            pos_at_A[0] = position_vector[A.index];
            pos_at_A[1] = position_vector[A.index+1];

            Vector2d pos_at_B;
            pos_at_B[0] = position_vector[B.index];
            pos_at_B[1] = position_vector[B.index+1];

            A_tmp.velocity = (A_tmp.position - pos_at_A)/time_step;
            B_tmp.velocity = (B_tmp.position - pos_at_B)/time_step;

            Vector2d current_friction_at_A = get_stokes_friction_force(A_tmp);
            Vector2d current_friction_at_B = get_stokes_friction_force(B_tmp);

            current_friction_force[A_tmp.index] = current_friction_at_A[0];
            current_friction_force[A_tmp.index + 1] = current_friction_at_A[1];
            current_friction_force[B_tmp.index] = current_friction_at_B[0];
            current_friction_force[B_tmp.index + 1] =current_friction_at_B[1];

        }

        return current_friction_force;
    }

    MatrixXd get_current_friction_force_jacobian(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        MatrixXd result(friction_matrix_dv.rows(),friction_matrix_dv.cols());
        result.setZero();
        double constant = 6.0 * M_PI * dynamic_viscosity * viscosity_multiplier;
        for(auto& particle_pair : connected_particles) {

            Particle &A = std::get<0>(particle_pair);
            Particle &B = std::get<1>(particle_pair);

            result(A.index,A.index) = constant * A.radius;
            result(A.index+1,A.index+1) = constant * A.radius;
            result(B.index,B.index) = constant * B.radius;
            result(B.index+1,B.index+1) = constant * B.radius;

        }
        return result;
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
           // std::cout<<std::endl<<std::endl<<"POSITION"<<std::endl<<A.position<<std::endl;
            Particle& B = std::get<1>(particle_pair);

            Vector2d current_spring_force = NEW_get_damped_spring_force_without_damping(particle_pair);
            Matrix2d current_spring_force_dx = NEW_get_damped_spring_force_dXA_without_damping(particle_pair);

           // std::cout<<"Spring Force Vector Before Filling \n"<<spring_force_vector<<"\n";
            spring_force_vector[A.index] += current_spring_force[0];
            spring_force_vector[A.index + 1] += current_spring_force[1];
            spring_force_vector[B.index] += -1.0 * current_spring_force[0];
            spring_force_vector[B.index + 1] += -1.0 *current_spring_force[1];
          //  std::cout<<"Spring Force Vector After Filling at positions "<<A.index<<" and "<<B.index<<"\n"<<spring_force_vector<<"\n";
            spring_force_vec_over_time_x.push_back(spring_force_vector[A.index]);
           // std::cout<<spring_force_vec_over_time_x[0]<<std::endl;
            spring_force_vec_over_time_y.push_back(spring_force_vector[A.index+1]);
           // std::cout<<"Spring Force DX Before Filling \n"<<spring_force_matrix_dx<<"\n";
            spring_force_matrix_dx.block<2,2>(A.index,A.index) += current_spring_force_dx;
            spring_force_matrix_dx.block<2,2>(B.index,B.index) += current_spring_force_dx;
            spring_force_matrix_dx.block<2,2>(A.index,B.index) += -1.0 * current_spring_force_dx;
            spring_force_matrix_dx.block<2,2>(B.index,A.index) += -1.0 * current_spring_force_dx;

            friction_vector[A.index] = get_stokes_friction_force(A)[0];
            friction_vector[A.index+1] = get_stokes_friction_force(A)[1];
            friction_vector[B.index] = get_stokes_friction_force(B)[0];
            friction_vector[B.index+1] = get_stokes_friction_force(B)[1];

            friction_force_over_time_x.push_back(friction_vector[A.index]);

            friction_matrix_dv.block<2,2>(A.index,A.index) = get_friction_force_dv(A);
            friction_matrix_dv.block<2,2>(B.index,B.index) = get_friction_force_dv(B);

            electricalfield_force_vector[A.index] = get_electrical_force(A)[0];
            electricalfield_force_vector[A.index +1] = get_electrical_force(A)[1];
            electricalfield_force_vector[B.index] = get_electrical_force(B)[0];
            electricalfield_force_vector[B.index+1] = get_electrical_force(B)[1];
            std::cout<<" e force = "<<electricalfield_force_vector[A.index +1]<<std::endl;



           // std::cout<<"Spring Force DX After Filling at positions "<<A.index<<" and "<<B.index<<"\n"<<spring_force_matrix_dx<<"\n";
            spring_force_derivative_x_in_x.push_back(spring_force_matrix_dx(A.index+1,A.index+1));

            //spring_force_matrix_dv.block<2,2>(A.index,A.index) += current_spring_force_dv;
            //spring_force_matrix_dv.block<2,2>(B.index,B.index) += -1.0 * current_spring_force_dv;

          //  std::cout<<"Mass MAtrix BEfore filling \n"<<mass_matrix<<"\n";
            mass_matrix(A.index,A.index) = A.mass;
            mass_matrix(A.index+1,A.index+1) = A.mass;

            mass_matrix(B.index,B.index) = B.mass;
            mass_matrix(B.index+1,B.index+1) = B.mass;
           // std::cout<<"Mass MAtrix After filling \n"<<mass_matrix<<"\n";

           // std::cout<<"Position Vector BEfore filling \n"<<position_vector<<"\n";
            position_vector(A.index) = A.position[0];
            position_vector(A.index+1) = A.position[1];
            position_vector(B.index) = B.position[0];
            position_vector(B.index+1) = B.position[1];
            if(first_itteration){
                initial_position_vector = position_vector;
            }
           // std::cout<<"Position Vector After filling \n"<<position_vector<<"\n";
            position_vec_over_time_in_x.push_back(A.position[0]);
            position_vec_over_time_in_y.push_back(A.position[1]);


            velocity_vector(A.index) = A.velocity[0];
            velocity_vector(A.index+1) = A.velocity[1];
            velocity_vector(B.index) = B.velocity[0];
            velocity_vector(B.index+1) = B.velocity[1];

        }


    }

    VectorXd get_current_EHD_force(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles, VectorXd current_positions){
        int n = spring_force_vector.size();
        VectorXd current_EHD_force(n);
        current_EHD_force.setZero();

        for(auto& particle_pair : connected_particles){
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            Vector2d local_EHD_force = get_my_EHD_force(particle_pair);
            current_EHD_force[A.index] = local_EHD_force[0];
            current_EHD_force[A.index+1] = local_EHD_force[1];
            current_EHD_force[B.index] = local_EHD_force[0];
            current_EHD_force[B.index+1] = local_EHD_force[1];

        }
        return current_EHD_force;
    }



    //CALCULATIONS FOR e, g AND g/dx

    double get_e(VectorXd x,std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        VectorXd a = (x - 2.0 * position_vector + position_vector_minus_1);
        double x_t_M_x = 0.5 * a.transpose() * mass_matrix * a;

        double E_x =  get_current_potential_energy(x,connected_particles)+ get_E_drag_indefinit_integral(x,connected_particles)+ get_E_electrical(x) + get_E_EHD(x,connected_particles);

        return   x_t_M_x +  std::pow(time_step,2) * E_x; //check again
    }

    VectorXd get_g_with_drag(VectorXd x, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        VectorXd F_x = get_current_spring_force(connected_particles,x) - get_current_friction_force(connected_particles,x) + electricalfield_force_vector + get_current_EHD_force(connected_particles,x);
       // std::cout<<"FS = \n"<<get_current_spring_force(connected_particles,x)<<"Fd = \n"<<get_current_friction_force(connected_particles,x)<<" Fel = \n "<<electricalfield_force_vector;
       //   std::cout<<" F = "<<std::endl<<F_x<<std::endl;
        // std::cout<<"mass matrix = "<<std::endl<<mass_matrix<<std::endl;
      //  std::cout<<" M * (x-y) = \n"<<mass_matrix * (x - 2.0 * position_vector + position_vector_minus_1)<<std::endl;
        return mass_matrix * (x - 2.0 * position_vector + position_vector_minus_1) - std::pow(time_step, 2) * F_x;
    }

    MatrixXd get_g_dx(VectorXd x, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        MatrixXd dF_dx = get_current_spring_force_jacobian(connected_particles,x) - get_current_friction_force_jacobian(connected_particles,x) + get_current_electrical_force_jacobian(connected_particles,x);
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
            B_tmp.position[0] = x[B.index];
            B_tmp.position[1] = x[B.index+1];

            std::tuple<Particle&,Particle&,double> particle_pair_tmp(A_tmp,B_tmp,std::get<2>(particle_pair));

            Vector2d displacement =  std::get<0>(particle_pair_tmp).position - std::get<1>(particle_pair_tmp).position;
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
                Vector2d new_pos;
                new_pos[0] = x[A.index];
                new_pos[1] = x[A.index+1];

                Vector2d old_pos;
                old_pos[0] = position_vector[A.index];
                old_pos[1] = position_vector[A.index+1];

                double t1 = new_pos.dot(new_pos);

                double t2 = old_pos.dot(new_pos);


                E += (6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius)/time_step * (t1/2.0 - t2);
                A.visited = true;
            }
            if(!B.visited){
                Vector2d new_pos;
                new_pos[0] = x[B.index];
                new_pos[1] = x[B.index+1];

                Vector2d old_pos;
                old_pos[0] = position_vector[B.index];
                old_pos[1] = position_vector[B.index+1];

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

    Vector2d NEW_get_damped_spring_force_without_damping(std::tuple<Particle,Particle,double> particle_pair){
       // std::cout<<" ======= inside =========="<<std::endl;
        double rest_length = std::get<2>(particle_pair);
      //  std::cout<<"rest length = "<<rest_length<<std::endl;
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
       // std::cout<<" positions are = "<<std::endl<<A.position[0]<<std::endl<<A.position[1]<<std::endl<<B.position[0]<<std::endl<<B.position[1]<<std::endl;
        //changed plus form damping to minus
        A_B_DISTANCE.push_back((A.position - B.position).norm());
        Vector2d direction = (A.position - B.position)/(A.position - B.position).norm();
      //  std::cout<<"direction = "<<direction<<std::endl;
        long double norm_ = (A.position - B.position).norm();
        long double radi_ = A.radius + B.radius;
        long double rest_length_ = rest_length;
        long double test = (norm_ - radi_) - rest_length_;
        //std::cout << std::fixed;
        //std::cout << std::setprecision(15);
      //  std::cout<<" norm = "<<norm_<<"radi = "<<radi_<<" rest_length = "<<rest_length_<<"sum = "<<test<<std::endl;
      //  std::cout<<" magnitude = "<<((A.position - B.position).norm()- (A.radius + B.radius)) - rest_length<<" = "<<(A.position - B.position).norm()<<" - "<<(A.radius + B.radius) <<" - "<<rest_length<<std::endl;
        Vector2d result_ = -1.0 *(stiffnes_constant * ( ((A.position - B.position).norm()) - rest_length) ) * (direction);
       // std::cout<<std::endl<<"FORCE = "<<std::endl<<result<<std::endl;
      //  std::cout<<" results = "<<std::endl<<result_<<std::endl;
      //  std::cout<<" ======= exit =========="<<std::endl;
       return result_;

    }

    Matrix2d NEW_get_damped_spring_force_dXA_without_damping(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        Vector2d xa = A.position;
        Vector2d xb = B.position;
        Vector2d va = A.velocity;
        Vector2d vb = B.velocity;

        Vector2d t0 = xa - xb;
      //  std::cout<<" t0 = \n"<<t0<<"\n;";
        double t1 = t0.norm();
       // std::cout<<" t1 = \n"<<t1<<"\n";
        Matrix2d T2 = t0 * t0.transpose();
        //std::cout<<std::endl<<"T2="<<std::endl<<T2<<std::endl;

       // std::cout<<" T2 = \n"<<T2<<"\n";
       //double t3 = stiffnes_constant * (t1 - rest_length);
       double t3 = stiffnes_constant * (t1 - rest_length);
       // std::cout<<" t3 = \n"<<t3<<"\n";
       // std::cout<<"SPRING FORCE DX = \n"<<-1.0<<" * "<< stiffnes_constant<<" / "<<std::pow(t1,2)<<" * "<<T2<<" - "<<t3<<" / "<<std::pow(t1,3)<<" * "<<T2<<" + "<<t3<<" / "<<t1<<" * "<<Identity<<"\n";
       // std::cout<<"Spring force Dx = \n"<<-1.0<<" * "<<(stiffnes_constant/std::pow(t1,2)) * T2<<" - "<<(t3/std::pow(t1,3)) * T2<<" + "<<(t3/t1) * Identity<<"\n";
      //  std::cout<<std::endl<<"term 1 = "<<(stiffnes_constant/std::pow(t1,2)) * T2<<std::endl<<"Term 2"<<(t3/std::pow(t1,3)) * T2<<std::endl<<" 1 - 2"<<std::endl<<(stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2<<std::endl<<"Term 3"<<std::endl<<(t3/t1) * Identity<<std::endl;
        Matrix2d result = -1.0 * ((stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2 + (t3/t1) * Identity);


        //OG

        //Vector2d dir = (xa - xb)/(xa -xb).norm();
        //double length = (xa-xb).norm();
        //return ((Identity - dir * dir.transpose() ) * std::min(1.0,length/rest_length) - Identity) * (-1.0 * stiffnes_constant);
        //TRYING TO AVOID SINGULARITY

        //NEWONE
       /* Matrix2d result = -stiffnes_constant * (
                (1.0 - rest_length/(xa-xb).norm()) * (Identity - (t0/t1) * (t0/t1).transpose())
                + (t0/t1) * (t0/t1).transpose()
                );
                */
       // std::cout<<std::endl<<std::endl<<result(0,0)<<std::endl<<std::endl;
       // std::cout<<result<<std::endl;
        //double x = result(0,0);
        //double y = result(1,1);
        //result(0,0) = y;
        //result(1,1)=x;
        return result;
    }

    Vector2d get_stokes_friction_force(Particle particle){
        return 6.0 * M_PI * viscosity_multiplier*dynamic_viscosity * particle.radius * particle.velocity;
    }

    Matrix2d get_friction_force_dv(Particle particle){ //should it be minus
        return Eigen::MatrixXd::Identity(2,2) * 6.0 * M_PI * viscosity_multiplier * dynamic_viscosity * particle.radius;
    }

    //electric stuff from particle and electric fiel paper
    Vector2d get_electrical_force(Particle particle){
        return get_coulomb_force(particle) /*+ get_image_force(particle) + get_gradient_force(particle)*/;
    }

    Vector2d get_coulomb_force(Particle particle){
        return get_electrical_field() * particle.charge;
    }

    Vector2d get_electrical_field(){
        Vector2d E_field;
        E_field[1] =  (lower_electrode.voltage - upper_electrode.voltage)/(std::abs(lower_electrode.position[1] - upper_electrode.position[1]));
        E_field[0] = 0.0;
        return E_field;
    }

    Vector2d brownian_motion(Particle &particle){

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution (0.0,1.0);

        double x_rand = distribution(generator);
        double y_rand = distribution(generator);
        long double D = (kB * T)/(3 * M_PI * eta * (particle.radius*2)*std::pow(10,-6) ); //not sure
        long double k = std::sqrt(2 * D * time_step);

        double dx = x_rand * k;
        double dy = y_rand * k;
        Vector2d result;
        result[0] = dx;
        result[1] = dy;
        return result;
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

    double get_debye_length(Particle particle){
        return 2.0 * particle.charge * std::sqrt((M_PI * particle.density)/(kB * T));
    }

    Vector2d get_my_EHD_force(std::tuple <Particle&,Particle&,double> particle_pair){
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);

        double r_A;
        double r_B;

        if(A.position[0] > B.position[0]){
            r_A = 2.0* A.radius;
            r_B = -2.0 * B.radius;
        }else{
            r_A = -2.0* A.radius;
            r_B = 2.0 * B.radius;
        }
        double U_EHD_A = get_U_EHD(A,1.0,r_A);
        double U_EHD_B = get_U_EHD(B,1.0,r_B);

        Vector2d F_EHD_A = {6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * A.radius * U_EHD_A,0.0};
        Vector2d F_EHD_B = {6.0 * M_PI * dynamic_viscosity * viscosity_multiplier * B.radius * U_EHD_B,0.0};

        std::cout<<"FA = "<<F_EHD_A<<std::endl<<"FB = "<<F_EHD_B<<std::endl<<"FTOT = "<<F_EHD_A + F_EHD_B<<std::endl;

        return F_EHD_A + F_EHD_B;
    }

    double get_U_EHD(Particle particle, double beta,double r){ //deleted vacuum permittivity
        /*double t1 = beta * relative_permittivity * get_debye_length(particle) * get_H() / dynamic_viscosity;
        double t2 = std::pow(lower_electrode.peak_voltage/(2.0 * get_H()),2);
        double t3 = (get_K1() + get_omega_bar(particle) * get_K2(particle))/(1.0 + std::pow(get_omega_bar(particle),2));
        double t4 = (3.0 * (r/particle.radius))/(2.0 * std::pow(1.0 + std::pow(r/particle.radius,2),(5.0/2.0)));
        //std::cout<<" U EHD for Particle "<<particle.index<<" = "<<t1 * t2 * t3 * t4<<std::endl;
        return  t1 * t2 * t3 * t4;*/

        double alpha = (lower_electrode.frequency * get_H()) /( get_debye_length(particle) * get_ion_diffusivity(particle));
        double alpha_2 = std::pow(alpha,2);
        double C = relative_permittivity * std::pow(lower_electrode.peak_voltage/(2.0 * get_H()),2) * (alpha_2/(1 + alpha_2)) * (get_debye_length(particle) * get_ion_diffusivity(particle))/lower_electrode.frequency;
        return (C * get_K2(particle) / dynamic_viscosity) * ( (3.0 * (r/particle.radius))/( 2.0 * std::pow(1 + std::pow(r/particle.radius, 2),(5.0/2.0))) );
    }

    double get_H(){

        return 0.5 * std::abs(lower_electrode.position[1] - upper_electrode.position[1]);
   }

    double get_omega_bar(Particle particle){
        return (lower_electrode.frequency * get_H()) /(get_debye_length(particle) * get_ion_diffusivity(particle));
    }

    double get_ion_diffusivity(Particle particle){
        return (kB * T)/( 6.0 * M_PI * dynamic_viscosity * particle.radius);
    }



    double get_K1(){
        return 1.0;
    }

    double get_K2(Particle particle){
        double w = lower_electrode.frequency;
        double w_2 = std::pow(w,2);
        double sigma_m = conductivity;
        double sigma_p = particle.conductivity;
        double epsilon_p = particle.permittivity;
        double epsilon_m = relative_permittivity;
        return ( w * (sigma_m - sigma_p)*(epsilon_p + 2.0 * epsilon_m) - (epsilon_p - epsilon_m)*(sigma_p + 2.0 * sigma_m)   )/(w_2 * std::pow((epsilon_p + 2.0 * epsilon_m),2) + std::pow((sigma_p + 2.0 * sigma_m),2));
    }










    //EXPERIMENTAL STUFF FOR LATER

    void check_Boundaries(Particle &particle){

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
                particle.velocity[1] = -1 *particle.velocity[1];
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
            Vector2d displacement =  std::get<0>(particle_pair).position - std::get<1>(particle_pair).position;
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