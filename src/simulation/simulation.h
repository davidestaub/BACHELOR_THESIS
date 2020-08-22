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

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::Matrix2d;
using Eigen::Matrix;
using Eigen::MatrixXd;

class Simulation{
public:
    int num_steps = 1;
    double time_step = Particle::time_step;
    int time_index = 0;
    std::vector<double> total_energy;
    std::vector<double> kinetic_energy;
    std::vector<double> potential_energy;
    VectorXd spring_force_vector;
    MatrixXd spring_force_matrix_dx;
    MatrixXd spring_force_matrix_dv;
    MatrixXd mass_matrix;
    VectorXd velocity_vector;
    double micro_time_step = time_step * std::pow(10,6);

    // all units in micro scale, (10^-6)

    double eta = 1.0 * std::pow(10,-3); //viscosity of water in PA * seconds

    double T = 293; //room temperatur
    double kB = 1.38 * std::pow(10,-23); //Boltzmann constant
    Vector2d force_damper;

    double noise_damper = 4;
    double drag_damper = 0;

    Vector2d Boxcenter = {0,0};
    Vector2d Boxsize = {1500,1500};
    double xbuffer, ybuffer = 5.0;

    //SPRING STUFF//
    double stiffnes_constant = 1000.0;
    //double rest_length = 2 * Particle::radius_for_spring;
    double damping_constant = 200.0;
    int newton_counter = 0;

    std::vector<double> velocities_over_time1_in_x;
    std::vector<double> velocities_over_time1_in_y;

    int size = -1;

    //capture data for ploting









    Simulation() {}


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

    void connect(Particle &A, Particle &B, double restlength, std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        std::tuple<Particle&,Particle&,double> particle_pair(A,B,restlength);
        connected_particles.push_back(particle_pair);
    }

    Vector2d eval_F(Vector2d delta_v_init_A, double h, double six_Pi_mu_r_A, double h2, double d_f_spring_x_dx_A, double d_f_spring_y_dx_A, double d_f_spring_x_dy_A, double d_f_spring_y_dy_A, double d_f_circle_y_dx, double d_f_circle_x_dy, Vector2d F_A, Vector2d v_k_minus_1_A ){
        Vector2d F_Newton_A;
        F_Newton_A[0] = delta_v_init_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_x_dx_A) - delta_v_init_A[1] * h2 * (d_f_circle_x_dy + d_f_spring_x_dy_A) - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * d_f_spring_x_dx_A + v_k_minus_1_A[1] * (d_f_circle_x_dy + d_f_spring_x_dy_A)) )  );
        F_Newton_A[1] = delta_v_init_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_y_dy_A) - delta_v_init_A[0] * h2 * (d_f_circle_y_dx + d_f_spring_y_dx_A) - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * d_f_spring_y_dy_A + v_k_minus_1_A[0] * (d_f_circle_y_dx + d_f_spring_y_dx_A)) )  );
        return F_Newton_A;
    }

    Vector2d NEW_get_damped_spring_force(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        //changed plus form damping to minus
        return -1.0 *(stiffnes_constant * ((A.position - B.position).norm() - rest_length) + (damping_constant * (A.velocity - B.velocity).transpose()) * ((A.position - B.position)/(A.position - B.position).norm()) ) * ((A.position - B.position)/(A.position - B.position).norm());
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
        spring_force_matrix_dx = MatrixXd::Zero(size,size);
        spring_force_matrix_dv = MatrixXd::Zero(size,size);
        mass_matrix = MatrixXd::Zero(size,size);
        velocity_vector = VectorXd::Zero(size);
        reset_flags(connected_particles);
    }

    void update_velocities_and_positions(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles,VectorXd delta_v){
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);
            if(!A.visited){
                A.velocity[0] += delta_v(A.index);
                A.velocity[1] += delta_v(A.index+1);
                A.position += time_step * A.velocity;
                A.visited = true;
            }
            if(!B.visited){
                B.velocity[0] += delta_v(B.index);
                B.velocity[1] += delta_v(B.index+1);
                B.position += time_step * B.velocity;
                B.visited = true;
            }
            std::cout<<A.visited<<std::endl<<B.visited<<std::endl;
        }

        reset_flags(connected_particles);
    }

    void fill_spring_force_and_derivatives(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

        spring_force_vector.setZero();
        spring_force_matrix_dx.setZero();
        spring_force_matrix_dv.setZero();
        mass_matrix.setZero();
        velocity_vector.setZero();


        for(auto& particle_pair : connected_particles) {

            Particle& A = std::get<0>(particle_pair);
            Particle& B = std::get<1>(particle_pair);

            Vector2d current_spring_force = NEW_get_damped_spring_force(particle_pair);
            Matrix2d current_spring_force_dx = NEW_get_damped_spring_force_dXA(particle_pair);
            Matrix2d current_spring_force_dv = NEW_get_damped_spring_force_dVA(particle_pair);

            spring_force_vector[A.index] += current_spring_force[0];
            spring_force_vector[A.index + 1] += current_spring_force[1];
            spring_force_vector[B.index] += -1.0 * current_spring_force[0];
            spring_force_vector[B.index + 1] += -1.0 *current_spring_force[1];

            spring_force_matrix_dx.block<2,2>(A.index,A.index) += current_spring_force_dx;
            spring_force_matrix_dx.block<2,2>(B.index,B.index) += -1.0 * current_spring_force_dx;

            spring_force_matrix_dv.block<2,2>(A.index,A.index) += current_spring_force_dv;
            spring_force_matrix_dv.block<2,2>(B.index,B.index) += -1.0 * current_spring_force_dv;

            mass_matrix(A.index,A.index) = A.mass;
            mass_matrix(A.index+1,A.index+1) = A.mass;

            mass_matrix(B.index,B.index) = B.mass;
            mass_matrix(B.index+1,B.index+1) = B.mass;

            velocity_vector(A.index) = A.velocity[0];
            velocity_vector(A.index+1) = A.velocity[1];
            velocity_vector(B.index) = B.velocity[0];
            velocity_vector(B.index+1) = B.velocity[1];
        }
    }

    Vector2d get_damped_spring_force(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        //changed plus form damping to minus
        return (stiffnes_constant * ((A.position - B.position).norm() - rest_length) + (damping_constant * (A.velocity - B.velocity).transpose()) * ((A.position - B.position)/(A.position - B.position).norm()) ) * ((A.position - B.position)/(A.position - B.position).norm());
    }

    Matrix2d get_damped_spring_force_dXA(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        Vector2d xa = A.position;
        Vector2d xb = B.position;
        Vector2d va = A.position;
        Vector2d vb = B.position;

        Vector2d t0 = xa -xb;
        double t1 = t0.norm();
        double t2 = t0.squaredNorm();
        Vector2d t3 = va -vb;
        Matrix2d T4 = t0 * t0.transpose();
        double t5 = t3.transpose() * t0;
        double t6 = stiffnes_constant * (t1 - rest_length) + (damping_constant*t5)/t1;

        return -1.0* (  stiffnes_constant/t2 * T4 + damping_constant/t2 * t0 * t3.transpose() - (damping_constant * t5)/std::pow(t1,4)*T4 - t6/std::pow(t1,3) * T4 + t6/t1 * Identity );


    }

    Matrix2d NEW_get_damped_spring_force_dXA(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        Vector2d xa = A.position;
        Vector2d xb = B.position;
        Vector2d va = A.position;
        Vector2d vb = B.position;

        Vector2d t0 = xa -xb;
        double t1 = t0.norm();
        double t2 = std::pow(t1,2);
        Vector2d t3 = va -vb;
        Matrix2d T4 = t0 * t0.transpose();
        double t5 = t3.transpose() * t0;
        double t6 = stiffnes_constant * (t1 - rest_length) + (damping_constant*t5)/t1;

        return -1.0 * ((stiffnes_constant/t2) * T4 + (damping_constant/t2) * t0 * t3.transpose() - ((damping_constant * t5)/std::pow(t1,4)) * T4 - (t6/std::pow(t1,3)) * T4 + (t6/t1) * Identity);


    }



    Matrix2d get_damped_spring_force_dVA(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        Vector2d t0 = A.position - B.position;
        return (damping_constant)/t0.squaredNorm() * t0 * t0.transpose();
    }

    Matrix2d NEW_get_damped_spring_force_dVA(std::tuple<Particle,Particle,double> particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle A = std::get<0>(particle_pair);
        Particle B = std::get<1>(particle_pair);
        Vector2d t0 = A.position - B.position;
        return ((-1.0 *damping_constant)/t0.squaredNorm()) * t0 * t0.transpose();
    }

    Matrix2d evaluate_Jacobian(Matrix2d d_f_dv_A, Particle& A, Matrix2d d_f_dx_A){
        return Eigen::MatrixXd::Identity(2,2) - d_f_dv_A * (time_step/A.mass) - d_f_dx_A * (std::pow(time_step,2)/A.mass);
    }

    MatrixXd evaluate_Jacobian_ALL(MatrixXd d_f_dv, MatrixXd d_f_dx){
        return Eigen::MatrixXd::Identity(size,size) - d_f_dv * (time_step* mass_matrix.inverse()) - d_f_dx * (std::pow(time_step,2) * mass_matrix.inverse());
    }

    Matrix2d get_df_dx_B(std::tuple<Particle,Particle,double> particle_pair){
        Particle B =  std::get<1>(particle_pair);
        Matrix2d d_f_spring_dxB;
        d_f_spring_dxB =-1.0 * get_damped_spring_force_dXA(particle_pair);
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
        return d_fcircle_dxy_B + d_f_spring_dxB;
    }

    Matrix2d get_df_dv_B(std::tuple<Particle,Particle,double> particle_pair){
        Particle B = std::get<1>(particle_pair);
        Matrix2d d_f_dv_B;
        d_f_dv_B << -B.six_Pi_mu_r,0,0,-B.six_Pi_mu_r;
        Matrix2d d_f_spring_dvA = get_damped_spring_force_dVA(particle_pair);
        Matrix2d d_f_spring_dvB = -1.0 * d_f_spring_dvA;
        d_f_dv_B += d_f_spring_dvB;
        return d_f_dv_B;
    }

    Vector2d get_F_B(std::tuple<Particle,Particle,double> particle_pair){
        Particle B = std::get<1>(particle_pair);
        return B.circle_force2(B) +get_drag(B) + (-1.0 *get_damped_spring_force(particle_pair));
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

    double get_current_velocities_1_in_x(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
       double result = 0;
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
           result = A.velocity[0];
        }
        return result;
    }

    double get_current_velocities_1_in_y(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        double result = 0;
        for(auto& particle_pair : connected_particles) {
            Particle& A = std::get<0>(particle_pair);
            result = A.velocity[1];
        }
        return result;
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


    void Newtons_Method(int maxiter,double tol, double alpha, double beta, VectorXd delta_v_init, VectorXd& delta_v_old, VectorXd v_k_minus_1_A ){

        double h2 = std::pow(time_step,2);
        double h = time_step;
        int count =0;

        VectorXd F_Newton;
        MatrixXd Jacobian;
        VectorXd F_Newton_new;

        for(int i = 0; i<maxiter;i++){
            double t = 1.0;
            F_Newton = evaluate_F_ALL(spring_force_matrix_dv,spring_force_matrix_dx, delta_v_old,spring_force_vector);
            Jacobian = evaluate_Jacobian_ALL( spring_force_matrix_dv, spring_force_matrix_dx);
            double old_eval = F_Newton.norm();

            VectorXd delta_v_new = delta_v_old-  t *Jacobian.inverse() * F_Newton;
            F_Newton_new = evaluate_F_ALL(spring_force_matrix_dv,spring_force_matrix_dx,delta_v_new,spring_force_vector);
            double new_eval = F_Newton_new.norm();
            std::cout<<"new eval = "<<new_eval<<std::endl;

            while(new_eval > old_eval){ //decrement t until step gives better result
                t = beta * t;
                delta_v_new = delta_v_old - t * Jacobian.inverse() * F_Newton;
                F_Newton_new = evaluate_F_ALL(spring_force_matrix_dv,spring_force_matrix_dx,delta_v_new,spring_force_vector);
                new_eval = F_Newton_new.norm();
            }

            delta_v_new = delta_v_old -  t *Jacobian.inverse() * F_Newton; // Newton update
            delta_v_old = delta_v_new;

            F_Newton = evaluate_F_ALL(spring_force_matrix_dv,spring_force_matrix_dx, delta_v_old,spring_force_vector);

            count++;

            if((Jacobian.inverse() * F_Newton).norm() < tol ){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init<<" and delta_v_end = "<<delta_v_old<<std::endl;
                i = maxiter;
                break;

            }else{
                std::cout<<std::endl<<" value is = "<<(Jacobian.inverse() * F_Newton).norm()<<std::endl;
            }
        }
    }



    //The latest version
    void run_simulation_for_connected_Particles(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(int i = 0; i<num_steps;i++) {
            total_energy.push_back(0.0);
            kinetic_energy.push_back(0.0);
            potential_energy.push_back(0.0);
            velocities_over_time1_in_x.push_back(0.0);
            velocities_over_time1_in_y.push_back(0.0);
            /*for (auto& particle_pair : connected_particles) {
                UPDATE_damped_spring_together(particle_pair);
            }*/
            UPDATE_SYSTEM(connected_particles);


        }
    }

    void UPDATE_SYSTEM(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){

        fill_spring_force_and_derivatives(connected_particles);

        MatrixXd Full_Identity = MatrixXd::Identity(size,size);
        MatrixXd A_init = Full_Identity - (time_step * mass_matrix.inverse()) * spring_force_matrix_dv - std::pow(time_step,2) * mass_matrix.inverse() * spring_force_matrix_dx;
        VectorXd b_init = time_step * mass_matrix.inverse() * (spring_force_vector + time_step * spring_force_matrix_dx * velocity_vector);
        VectorXd delta_v_init = A_init.colPivHouseholderQr().solve(b_init);

        VectorXd delta_v_old = delta_v_init;

        Newtons_Method(100000,std::pow(10,(-12)), 0.25, 0.5, delta_v_init,  delta_v_old, velocity_vector);
        newton_counter++;

        update_velocities_and_positions(connected_particles,delta_v_old);
    }

    void UPDATE_damped_spring_together(std::tuple<Particle&,Particle&,double> &particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        /*add_brownian_motion(A);
        add_brownian_motion(B);*/
        std::cout<<"added brownian motion"<<std::endl;
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;



        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);

        Vector2d f_drag_A =  get_drag(A);
        Vector2d f_drag_B =  get_drag(B);


        Vector2d f_spring_A = NEW_get_damped_spring_force(particle_pair);
        Vector2d f_spring_B = -1.0 * f_spring_A;

        Vector2d F_A = /*f_c_A   - f_drag_A+*/ f_spring_A;
        std::cout<<"F_A = "<<F_A<<" with f_circle = "<<f_c_A<<" and f_drag= "<<f_drag_A<<"and f_spring = "<<f_spring_A<<std::endl;


        Vector2d F_B = /*f_c_B - f_drag_B+*/ f_spring_B;
        std::cout<<"F_B = "<<F_B<<" with f_circle = "<<f_c_B<<" and f_drag= "<<f_drag_B<<"and f_spring = "<<f_spring_B<<std::endl;


        Matrix2d d_f_spring_dxA;
        d_f_spring_dxA =NEW_get_damped_spring_force_dXA(particle_pair);

        Matrix2d d_f_spring_dxB;        //Check again
        d_f_spring_dxB= -1.0 * d_f_spring_dxA;


        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;




        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        //d_f_dv_A << -A.six_Pi_mu_r,0,0,-A.six_Pi_mu_r;
        Matrix2d d_f_spring_dvA = NEW_get_damped_spring_force_dVA(particle_pair);
        d_f_dv_A = d_f_spring_dvA; //Check if this is correct
        //std::cout<<" d F/dvA = "<<d_f_dv_A<<" with d_fdrag/dv= "<<(d_f_dv_A-d_f_spring_dvA)<<" and d_fspring/dvA = "<<d_f_spring_dvA<<std::endl;
        Matrix2d d_f_dv_B;
        Matrix2d d_f_spring_dvB = -1.0 * d_f_spring_dvA;
        //d_f_dv_B << -B.six_Pi_mu_r,0,0,-B.six_Pi_mu_r;
        d_f_dv_B = d_f_spring_dvB; //Check if this is correct
       // std::cout<<" d F/dvB = "<<d_f_dv_B<<" with d_fdrag/dv= "<<(d_f_dv_B-d_f_spring_dvB)<<" and d_fspring/dvB = "<<d_f_spring_dvB<<std::endl;

        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
        //d_fcircle_dxy << 0,-1.0,1.0,0;
        d_f_dx_A = /*d_fcircle_dxy_A +*/ d_f_spring_dxA;
        //std::cout<<" d_f/dxA = "<<d_f_dx_A<<" with dfcircle/dx = "<<d_fcircle_dxy_A<<" and dfspring/dx = "<<d_f_spring_dxA<<std::endl;
        d_f_dx_B = /*d_fcircle_dxy_B +*/d_f_spring_dxB;
       // std::cout<<" d_f/dxB = "<<d_f_dx_B<<" with dfcircle/dx = "<<d_fcircle_dxy_B<<" and dfspring/dx = "<<d_f_spring_dxB<<std::endl;
        std::cout<<std::endl<<"=======================CHECKS BEFORE INIT==============================="<<std::endl<<"iteration= "<<std::endl<<time_index<<std::endl<<std::endl<<std::endl<<std::endl<<"Identity="<<std::endl<<Identity<<std::endl<<"time_step="<<std::endl<<time_step<<std::endl<<"A.mass="<<std::endl<<A.mass<<std::endl<<"B.mass="<<std::endl<<B.mass<<std::endl<<"A.velocity="<<std::endl<<A.velocity<<std::endl<<"B.velocity="<<std::endl<<B.velocity<<std::endl<<"d_f_dv_A="<<std::endl<<d_f_dv_A<<std::endl<<"d_f_dv_B="<<std::endl<<d_f_dv_B<<std::endl<<"d_f_dx_A="<<std::endl<<d_f_dx_A<<std::endl<<"d_f_dx_B="<<std::endl<<d_f_dx_B<<std::endl;



        std::cout<<std::endl<<"F_A = "<<std::endl<<F_A<<std::endl;
        std::cout<<std::endl<<"F_B = "<<std::endl<<F_B<<std::endl;

        Matrix2d A_init_A = Identity - (time_step/A.mass) * d_f_dv_A - (std::pow(time_step,2)/A.mass) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step/B.mass) * d_f_dv_B - (std::pow(time_step,2)/B.mass) * d_f_dx_B;
        std::cout<<std::endl<<"A_init_A = "<<std::endl<<A_init_A<<std::endl;
        std::cout<<std::endl<<"A_init_B = "<<std::endl<<A_init_B<<std::endl;

        Vector2d b_init_A = (time_step/A.mass) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step/B.mass) * ( F_B  + time_step * d_f_dx_B * B.velocity);

        std::cout<<std::endl<<"b_init_A = "<<std::endl<<b_init_A<<std::endl;
        std::cout<<std::endl<<"b_init_B = "<<std::endl<<b_init_B<<std::endl;

        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);

        std::cout<<"delta v init A = "<<std::endl<<delta_v_init_A<<std::endl;
        std::cout<<"delta v init B = "<<std::endl<<delta_v_init_B<<std::endl;
        std::cout<<std::endl<<"========================================================================"<<std::endl;
        // end of get initial delta v

        //NEWTONS METHOD //

        int maxiter = 100000;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        Vector2d v_k_minus_1_A = A.velocity; //changed
        Vector2d v_k_minus_1_B =  B.velocity;
        int count =0;
        double tol = std::pow(10,(-12));
        double alpha = 0.25;
        double beta = 0.5;
        Matrix2d Jacobian_A;
        Matrix2d Jacobian_B;
        Vector2d F_Newton_A_new;
        Vector2d F_Newton_B_new;



        std::cout<<"AAAAAAAaaaaaaaaaaaaaaa"<<std::endl;
        for(int i = 0; i<maxiter;i++){
            double t_A = 1.0;
            double t_B = 1.0;
            F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A,v_k_minus_1_A);
            F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B,v_k_minus_1_B);
            Jacobian_A = evaluate_Jacobian( d_f_dv_A, A,  d_f_dx_A);
            Jacobian_B = evaluate_Jacobian( d_f_dv_B, B,  d_f_dx_B);
            double old_eval_A = F_Newton_A.norm();
            double old_eval_B = F_Newton_B.norm();
            std::cout<<"old eval_A before while = "<<old_eval_A<<std::endl;
            std::cout<<"old eval_B before while = "<<old_eval_B<<std::endl;

            Vector2d delta_v_new_A = delta_v_old_A -  t_A *Jacobian_A.inverse() * F_Newton_A;
            Vector2d delta_v_new_B = delta_v_old_B -  t_B *Jacobian_B.inverse() * F_Newton_B;

            F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);
            F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);

            double new_eval_A = F_Newton_A_new.norm();
            double new_eval_B = F_Newton_B_new.norm();
            std::cout<<"new_eval_A before while = "<<new_eval_A<<std::endl;

            while(new_eval_A > old_eval_A || new_eval_B > old_eval_B /*+ alpha * t * F_Newton_A.transpose() *(-Jacobian_A.inverse() * F_Newton_A)*/){ //decrement t until step gives better result
                t_A = beta * t_A;
                delta_v_new_A = delta_v_old_A - t_A * Jacobian_A.inverse() * F_Newton_A;
                F_Newton_A_new = evaluate_F(A, d_f_dv_A, d_f_dx_A, delta_v_new_A, F_A, v_k_minus_1_A);
                new_eval_A = F_Newton_A_new.norm();
                std::cout << " t_B = " << t_A << "new_eval_A= " << new_eval_A << " || ";
                t_B = beta * t_B;
                delta_v_new_B = delta_v_old_B - t_B * Jacobian_B.inverse() * F_Newton_B;
                F_Newton_B_new = evaluate_F(B, d_f_dv_B, d_f_dx_B, delta_v_new_B, F_B, v_k_minus_1_B);
                new_eval_B = F_Newton_B_new.norm();
                std::cout << " t_B = " << t_B << "new_eval_B= " << new_eval_B << " || ";
            }

            delta_v_new_A = delta_v_old_A -  t_A *Jacobian_A.inverse() * F_Newton_A; // Newton update
            delta_v_new_B = delta_v_old_B -  t_B *Jacobian_B.inverse() * F_Newton_B; // Newton update

            delta_v_old_A = delta_v_new_A;
            delta_v_old_B = delta_v_new_B;

            F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A,v_k_minus_1_A);
            F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B,v_k_minus_1_B);
            count++;

            if((Jacobian_A.inverse() * F_Newton_A).norm() < tol && (Jacobian_B.inverse() * F_Newton_B).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_A<<" and delta_v_end = "<<delta_v_old_A<<std::endl<<"with delta v init_B ="<<delta_v_init_B<< "and delta_v_end ="<<delta_v_old_B<<std::endl;

                i = maxiter;
                break;

            }else{
                std::cout<<std::endl;
                std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<" B value is = "<<(Jacobian_B.inverse() * F_Newton_B).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }
        newton_counter++;
        std::cout<<"vel 1"<<A.velocity<<std::endl<<"pos 1"<<A.position<<std::endl;
        Vector2d delta_v_A = delta_v_old_A;
        A.velocity += delta_v_A;
        std::cout<<"delta v A after Newton ="<<std::endl<<delta_v_A<<std::endl;
        A.position += time_step * A.velocity;
        Vector2d delta_v_B = delta_v_old_B;
        B.velocity += delta_v_B;
        std::cout<<"delta v B after Newton ="<<std::endl<<delta_v_B<<std::endl;
        B.position += time_step * B.velocity;
        std::cout<<"vel 2"<<A.velocity<<std::endl<<"pos 2"<<A.position<<std::endl;
        std::cout<<"vel 3"<<std::get<0>(particle_pair).velocity<<std::endl<<"pos 3"<<std::get<0>(particle_pair).position<<std::endl;
    }

    void UPDATE_Separate_but_new(std::tuple<Particle&,Particle&,double> &particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        add_brownian_motion(A);
        add_brownian_motion(B);

        //A.velocity *=0.1;
       // B.velocity *= 0.1;

        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);

        Vector2d f_drag_A =  get_drag(A);
        Vector2d f_drag_B =  get_drag(B);

        Vector2d f_spring_A = get_damped_spring_force(particle_pair);
        Vector2d f_spring_B = -1.0 * f_spring_A;


        Vector2d F_A = f_c_A   - f_drag_A+ f_spring_A;

        Vector2d F_B = f_c_B - f_drag_B+ f_spring_B;


        Matrix2d d_f_spring_dA;
        d_f_spring_dA = get_damped_spring_force_dXA(particle_pair); // new
        Matrix2d d_f_spring_dB;
        d_f_spring_dB= -1.0 * d_f_spring_dA;

        double d_f_spring_x_dx_A = d_f_spring_dA(0,0);
        double d_f_spring_x_dx_B = d_f_spring_dB(0,0);
        double d_f_spring_x_dy_A = d_f_spring_dA(0,1);
        double d_f_spring_x_dy_B = d_f_spring_dB(0,1);
        double d_f_spring_y_dx_A = d_f_spring_dA(1,0);
        double d_f_spring_y_dx_B = d_f_spring_dB(1,0);
        double d_f_spring_y_dy_A = d_f_spring_dA(1,1);
        double d_f_spring_y_dy_B = d_f_spring_dB(1,1);


        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;


        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        d_f_dv_A << -A.six_Pi_mu_r,0,0,-A.six_Pi_mu_r;
        Matrix2d d_f_spring_dvA = get_damped_spring_force_dVA(particle_pair); //new
        d_f_dv_A += d_f_spring_dvA;
        Matrix2d d_f_dv_B;
        Matrix2d d_f_spring_dvB = -1.0 * d_f_spring_dvA;
        d_f_dv_B << -B.six_Pi_mu_r,0,0,-B.six_Pi_mu_r;
        d_f_dv_B += d_f_spring_dvB;



        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
        //d_fcircle_dxy << 0,-1.0,1.0,0;
        d_f_dx_A = d_fcircle_dxy_A + d_f_spring_dA;
        d_f_dx_B = d_fcircle_dxy_B +d_f_spring_dB;


        Matrix2d A_init_A = Identity - (time_step/A.mass) * d_f_dv_A - (std::pow(time_step,2)/A.mass) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step/B.mass) * d_f_dv_B - (std::pow(time_step,2)/B.mass) * d_f_dx_B;

        Vector2d b_init_A = (time_step/A.mass) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step/B.mass) * ( F_B   + time_step * d_f_dx_B * B.velocity);



        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);





        // end of get initial delta v

        //NEWTONS METHOD //

        int maxiter = 100;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        Vector2d v_k_minus_1_A = delta_v_init_A - A.velocity;
        Vector2d v_k_minus_1_B = delta_v_init_B - B.velocity;
        int count =0;
        double tol = std::pow(10,(-10));
        double alpha = 0.25;
        double beta = 0.5;


        std::cout<<"AAAAAAAaaaaaaaaaaaaaaa"<<std::endl; /// STILL WRONG; maybe delta_v init should change in line 268, abruch kriterium doesnt change and is wrong aswell
        for(int i = 0; i<maxiter;i++){
            double t = 1.0;
            F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A,v_k_minus_1_A);

            Matrix2d Jacobian_A;
            //Jacobian_A = Identity - d_f_dv_A * time_step/A.mass - d_f_dx_A * std::pow(time_step,2)/A.mass;
            Jacobian_A = evaluate_Jacobian( d_f_dv_A, A,  d_f_dx_A);

            double old_eval = F_Newton_A.norm();

            Vector2d F_Newton_A_new;

            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);


            double new_eval = F_Newton_A_new.norm();

            while(new_eval > old_eval + alpha * t * F_Newton_A.transpose() *(-Jacobian_A.inverse() * F_Newton_A)){
                t = beta * t;
                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

                F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);

                new_eval = F_Newton_A_new.norm();
            }

            delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            delta_v_old_A = delta_v_new_A;
            count++;

            if((Jacobian_A.inverse() * F_Newton_A).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_A<<" and delta_v_end = "<<delta_v_old_A<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
                //  std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }
        Vector2d delta_v_A = delta_v_old_A;
        A.velocity += delta_v_A;
        A.position += time_step * A.velocity;

        std::cout<<"Bbbbbbbbbbbbbbbbbbb"<<std::endl;

        for(int i = 0; i<maxiter;i++){
            double t = 1.0;
            F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B,v_k_minus_1_B);
            std::cout<<" F Newton = "<<F_Newton_B<<std::endl;
            Matrix2d Jacobian_B;
            //Jacobian_B = Identity - d_f_dv_B * time_step/B.mass - d_f_dx_B * std::pow(time_step,2)/B.mass;
            Jacobian_B = evaluate_Jacobian( d_f_dv_B, B,  d_f_dx_B);
            double old_eval = F_Newton_B.norm();
            Vector2d delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
            Vector2d F_Newton_B_new;
            F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);
            double new_eval = F_Newton_B_new.norm();

            int while_counter = 0;
            while(new_eval > old_eval + alpha * t * F_Newton_B_new.transpose() *(-Jacobian_B.inverse() * F_Newton_B)){
                t = beta * t;
                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
                F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);
                new_eval = F_Newton_B_new.norm();
                std::cout<<"t = "<<t<<std::endl;
                while_counter++;
                std::cout<<"while counter= "<<while_counter<<std::endl;
            }
            delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
            delta_v_old_B = delta_v_new_B;
            count++;
            if((Jacobian_B.inverse() * F_Newton_B).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_B<<" and delta_v_end = "<<delta_v_old_B<<std::endl;
                break;
            }
        }
        Vector2d delta_v_B = delta_v_old_B;
        B.velocity += delta_v_B;
        B.position += time_step * B.velocity;
    }

    void UPDATE(std::tuple<Particle&,Particle&,double> &particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        add_brownian_motion(A);
        add_brownian_motion(B);

        //A.velocity *=0.1;
       // B.velocity *= 0.1;

        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);
        std::cout<<"f c A =  " <<f_c_A<<std::endl;
        std::cout<<"f c B = " <<f_c_B<<std::endl;
        double six_Pi_mu_r_A = ( 6.0 * M_PI * eta * (A.radius*std::pow(10,-6))  * std::pow(10,6) );
        double six_Pi_mu_r_B = ( 6.0 * M_PI * eta * (B.radius*std::pow(10,-6))  * std::pow(10,6) );
        std::cout<<"six pi" <<six_Pi_mu_r_A<<std::endl;
        Vector2d f_drag_A =  A.velocity *  six_Pi_mu_r_A;
        Vector2d f_drag_B =  B.velocity *  six_Pi_mu_r_B;
        std::cout<<"f drag " <<f_drag_A<<std::endl;
        Vector2d f_spring_A = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) * ((A.position - B.position)/(A.position - B.position).norm());
       // f_spring_A;
        Vector2d f_spring_B = -1.0 * f_spring_A;
        std::cout<<"f spring " <<f_spring_A<<std::endl;
        Vector2d F_A = f_c_A   - f_drag_A+ f_spring_A;
        std::cout<<"F A = "<<F_A<<std::endl;
        Vector2d F_B = f_c_B - f_drag_B+ f_spring_B;
        std::cout<<"F B = "<<F_B<<std::endl;
        double root_for_spring = (A.position - B.position).norm();

        Vector2d t0 = A.position - B.position;
        double t1 = t0.norm();
        Matrix2d T2 = t0 * t0.transpose();
        double t3 = stiffnes_constant * (t1 - rest_length);
        Matrix2d d_f_spring_dA;
        d_f_spring_dA = -1.0 * ((stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2 + (t3/t1)* Identity );
        Matrix2d d_f_spring_dB;
        d_f_spring_dB= -1.0 * d_f_spring_dA;

        double d_f_spring_x_dx_A = d_f_spring_dA(0,0);
        double d_f_spring_x_dx_B = d_f_spring_dB(0,0);
        double d_f_spring_x_dy_A = d_f_spring_dA(0,1);
        double d_f_spring_x_dy_B = d_f_spring_dB(0,1);
        double d_f_spring_y_dx_A = d_f_spring_dA(1,0);
        double d_f_spring_y_dx_B = d_f_spring_dB(1,0);
        double d_f_spring_y_dy_A = d_f_spring_dA(1,1);
        double d_f_spring_y_dy_B = d_f_spring_dB(1,1);



       // double d_f_spring_x_dx_A = stiffnes_constant * std::pow(A.position[0] - B.position[0],2) * ( root_for_spring - rest_length) / std::pow(root_for_spring,3)
                         //        - stiffnes_constant * (root_for_spring - rest_length) / root_for_spring - stiffnes_constant * std::pow(A.position[0] - B.position[0],2) / std::pow(root_for_spring,2);
        //double d_f_spring_x_dx_B = -1.0 * d_f_spring_x_dx_A;
        //double d_f_spring_x_dy_A  = stiffnes_constant * rest_length * (A.position[0] - B.position[0]) * (A.position[1] - B.position[1]) / std::pow(root_for_spring,3);
        //double d_f_spring_x_dy_B =  -1.0 * d_f_spring_x_dy_A;
        //double d_f_spring_y_dx_A = stiffnes_constant * rest_length * (A.position[0] - B.position[0]) * (A.position[1] - B.position[1]) / std::pow(root_for_spring,3);
        //double d_f_spring_y_dx_B = -1.0 * d_f_spring_y_dx_A;
        //double d_f_spring_y_dy_A = stiffnes_constant * std::pow(A.position[1] - B.position[1],2) * ( root_for_spring - rest_length) / std::pow(root_for_spring,3)
                            //     - stiffnes_constant * (root_for_spring - rest_length) / root_for_spring - stiffnes_constant * std::pow(A.position[1] - B.position[1],2) / std::pow(root_for_spring,2);
        //double d_f_spring_y_dy_B = -1.0 * d_f_spring_y_dy_A;



        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;



        //tmp
/*
        Vector2d F_on_a;
        Vector2d F_on_b;

        F_on_a = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) *
                 ((A.position - B.position)/(A.position - B.position).norm());
        // std::cout<<"F on A = "<<F_on_a<<std::endl;
        F_on_b = -1.0 *F_on_a;
        // std::cout<<"F on B = "<<F_on_b<<std::endl;
        Vector2d t0_ = A.position - B.position;
        double t1_ = t0.norm();
        double t2_ = stiffnes_constant * (rest_length - t1);
        Matrix2d T3_ = t0 * t0.transpose();
        Matrix2d dF_on_a_dx;

        dF_on_a_dx = t2_/t1_ * Identity -
                     ( (stiffnes_constant/std::pow(t1_,2)) * T3_ + (t2_/std::pow(t1_,3))*T3_ );

        // std::cout<<"DF on A DX = "<<dF_on_a_dx<<std::endl;

        Matrix2d dF_on_b_dx;
        dF_on_b_dx = (stiffnes_constant/std::pow(t1_,2)) * T3_
                     + (t2_/std::pow(t1_,3))*T3_
                     - t2_/t1_*Identity;

        F_A.setZero();
        F_B.setZero();
        F_A =f_c_A   - f_drag_A + F_on_a;
        F_B =f_c_B - f_drag_B+ F_on_b;
        */

        //endtmp


        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        d_f_dv_A << -six_Pi_mu_r_A,0,0,-six_Pi_mu_r_A;
        Matrix2d d_f_dv_B;
        d_f_dv_B << -six_Pi_mu_r_B,0,0,-six_Pi_mu_r_B;
       // d_f_dv_A;
       // d_f_dv_B;

        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
                //d_fcircle_dxy << 0,-1.0,1.0,0;
        d_f_dx_A = d_fcircle_dxy_A + d_f_spring_dA;
        d_f_dx_B = d_fcircle_dxy_B +d_f_spring_dB;
      //  d_f_dx_A;
      //  d_f_dx_B;
       // d_f_dx_A << d_f_spring_x_dx_A, d_f_circle_x_dy + d_f_spring_x_dy_A, d_f_circle_y_dx + d_f_spring_y_dx_A, d_f_spring_y_dy_A;
       // d_f_dx_B << d_f_spring_x_dx_B, d_f_circle_x_dy + d_f_spring_x_dy_B, d_f_circle_y_dx + d_f_spring_y_dx_B, d_f_spring_y_dy_B;




        /*Matrix2d A_init_A = Identity - time_step * d_f_dv_A - std::pow(time_step,2) * d_f_dx_A;
        Matrix2d A_init_B = Identity - time_step * d_f_dv_B - std::pow(time_step,2) * d_f_dx_B;
        Vector2d b_init_A = time_step * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = time_step * ( F_B   + time_step * d_f_dx_B * B.velocity);*/

        Matrix2d A_init_A = Identity - (time_step/A.mass) * d_f_dv_A - (std::pow(time_step,2)/A.mass) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step/B.mass) * d_f_dv_B - (std::pow(time_step,2)/B.mass) * d_f_dx_B;

        Vector2d b_init_A = (time_step/A.mass) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step/B.mass) * ( F_B   + time_step * d_f_dx_B * B.velocity);



        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);





        // end of get initial delta v

        //NEWTONS METHOD //

        int maxiter = 100;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        Vector2d v_k_minus_1_A = delta_v_init_A - A.velocity;
        Vector2d v_k_minus_1_B = delta_v_init_B - B.velocity;
        int count =0;
        double tol = std::pow(10,(-10));
        double alpha = 0.25;
        double beta = 0.5;


        std::cout<<"AAAAAAAaaaaaaaaaaaaaaa"<<std::endl; /// STILL WRONG; maybe delta_v init should change in line 268, abruch kriterium doesnt change and is wrong aswell
        for(int i = 0; i<maxiter;i++){
            double t = 1.0;

           /* F_Newton_A[0] = delta_v_old_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_old_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                        - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
            F_Newton_A[1] = delta_v_old_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_old_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                          - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );
            */

            F_Newton_A[0] = delta_v_old_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_old_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
            F_Newton_A[1] = delta_v_old_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_old_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );



            // F_Newton_A = (Identity - d_f_dv_A * time_step- d_f_dx_A * std::pow(time_step,2)) * delta_v_old_A - time_step * (F_A + d_f_dx_A * time_step * A.velocity);

            Matrix2d Jacobian_A;

            //Jacobian_A << 1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_x_dx_A + h* six_Pi_mu_r_A, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_A),
            //        -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_A)  ,  1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_y_dy_A + h* six_Pi_mu_r_A;

            Jacobian_A << 1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(0,0) + (h/A.mass)* six_Pi_mu_r_A, -(h2/A.mass) * (d_f_circle_x_dy + d_f_spring_dA(0,1)), -(h2/A.mass) * (d_f_circle_y_dx + d_f_spring_dA(1,0))  ,  1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(1,1) + (h/A.mass)* six_Pi_mu_r_A;


          //  Vector2d old_eval_F = eval_F(delta_v_old_A,  h, six_Pi_mu_r_A, h2, d_f_spring_x_dx_A, d_f_spring_y_dx_A,  d_f_spring_x_dy_A,  d_f_spring_y_dy_A, d_f_circle_y_dx, d_f_circle_x_dy, F_A, v_k_minus_1_A);

            Matrix2d A_MAT = Identity - time_step * d_f_dv_A - std::pow(time_step,2) * d_f_dx_A;
            Vector2d b_vec = time_step * ( F_A  + time_step * d_f_dx_A * A.velocity);
            Vector2d delta_v_init_A = A_MAT.colPivHouseholderQr().solve(b_vec);
            //double old_eval = (A_MAT * delta_v_old_A - b_vec).norm();

            double old_eval = F_Newton_A.norm();
            //here

            Vector2d F_Newton_A_new;




            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            /*F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                            - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
            F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                            - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );
*/
            F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
            F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );

            //double new_eval = (A_MAT * delta_v_new_A - b_vec).norm();
            double new_eval = F_Newton_A_new.norm();

            while(new_eval > old_eval + alpha * t * /*old_eval_F*/F_Newton_A.transpose() *(-Jacobian_A.inverse() * F_Newton_A)){
                t = beta * t;
                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
              /*  F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                                    - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
                F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                                    - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );
*/
                F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
                F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );


                new_eval = F_Newton_A_new.norm();
                //new_eval = (A_MAT * delta_v_new_A - b_vec).norm();

            }

            delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            delta_v_old_A = delta_v_new_A;
            count++;

            if((Jacobian_A.inverse() * F_Newton_A).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_A<<" and delta_v_end = "<<delta_v_old_A<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
              //  std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }
        Vector2d delta_v_A = delta_v_old_A;
        A.velocity += delta_v_A;
        A.position += time_step * A.velocity;

        std::cout<<"Bbbbbbbbbbbbbbbbbbb"<<std::endl;

        for(int i = 0; i<maxiter;i++){
            double t = 1.0;

            /*F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                            - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
            F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                            - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );
*/
            F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
            F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );

            std::cout<<" F Newton = "<<F_Newton_B<<std::endl;

            Matrix2d Jacobian_B;
            //Jacobian_B << 1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_x_dx_B + h* six_Pi_mu_r_B, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_B), -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_B)  ,  1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_y_dy_B + h* six_Pi_mu_r_B;
            Jacobian_B << 1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(0,0) + h/B.mass * six_Pi_mu_r_B, -h2/B.mass * (d_f_circle_x_dy + d_f_spring_dB(0,1)), -h2/B.mass * (d_f_circle_y_dx + d_f_spring_dB(1,0))  ,  1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(1,1) + h/B.mass * six_Pi_mu_r_B;



            //Vector2d old_eval_F = eval_F(delta_v_old_B,  h, six_Pi_mu_r_B, h2, d_f_spring_x_dx_B, d_f_spring_y_dx_B,  d_f_spring_x_dy_B,  d_f_spring_y_dy_B, d_f_circle_y_dx, d_f_circle_x_dy, F_B, v_k_minus_1_B);
            Vector2d old_eval_F_B = F_Newton_B;
            Matrix2d A_MAT = Identity - time_step * d_f_dv_B - std::pow(time_step,2) * d_f_dx_B;
            Vector2d b_vec = time_step * ( F_B  + time_step * d_f_dx_B * B.velocity);
            Vector2d delta_v_init_B = A_MAT.colPivHouseholderQr().solve(b_vec);
            //double old_eval = (A_MAT * delta_v_old_B - b_vec).norm();


            double old_eval = F_Newton_B.norm();
            Vector2d delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

            Vector2d F_Newton_B_new;

           /* F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                            - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
            F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                            - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );
*/
            F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
            F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );

            // double new_eval = (A_MAT * delta_v_new_B - b_vec).norm();

           double new_eval = F_Newton_B_new.norm();

                int while_counter = 0;
            while(new_eval > (old_eval + alpha * t * F_Newton_B.transpose() *(-Jacobian_B.inverse() * F_Newton_B))){
                t = beta * t;
                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

               /* F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                                    - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
                F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                                    - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );
*/
                F_Newton_B[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
                F_Newton_B[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
                //!!!!!! should be F_Newton_B_new but that does not work and i dont get it anymore



                //new_eval = (A_MAT * delta_v_new_B - b_vec).norm();
                new_eval = F_Newton_B_new.norm();
                std::cout<<"t = "<<t<<std::endl;
                while_counter++;
                std::cout<<"while counter= "<<while_counter<<std::endl;
            }



            delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

            delta_v_old_B = delta_v_new_B;
            count++;




            if((Jacobian_B.inverse() * F_Newton_B).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_B<<" and delta_v_end = "<<delta_v_old_B<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
               // std::cout<<" value is = "<<(Jacobian_B.inverse() * F_Newton_B).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }



       // Vector2d delta_v_A = delta_v_old_A;
        Vector2d delta_v_B = delta_v_old_B;
      /*  std::cout<<" delta v A 2= "<<delta_v_A<<std::endl;
        std::cout<<" delta v B 2= "<<delta_v_B<<std::endl;

        //END OF NEWTONSMETHOD

        std::cout<<" A vel pre = "<<A.velocity<<std::endl;
        std::cout<<" B vel pre = "<<B.velocity<<std::endl;
        A.velocity += delta_v_A;
        std::cout<<" A vel post = "<<A.velocity<<std::endl;
*/
        B.velocity += delta_v_B;
       /* std::cout<<" B vel post = "<<B.velocity<<std::endl;
        std::cout<<" A pos pre = "<<A.position<<std::endl;
        std::cout<<" B pos pre = "<<B.position<<std::endl;
        //A.velocity *= 0.1;
       // B.velocity *= 0.1;
        A.position += time_step * A.velocity;
        */
        B.position += time_step * B.velocity;
      /*  std::cout<<" A pos post = "<<A.position<<std::endl;
        std::cout<<" B pos post = "<<B.position<<std::endl;
*/


        //TMP

        /*
         Vector2d F_on_a;
        Vector2d F_on_b;

        F_on_a = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) *
                 ((A.position - B.position)/(A.position - B.position).norm());
        // std::cout<<"F on A = "<<F_on_a<<std::endl;
        F_on_b = -1.0 *F_on_a;
        // std::cout<<"F on B = "<<F_on_b<<std::endl;
        Vector2d t0 = A.position - B.position;
        double t1 = t0.norm();
        double t2 = stiffnes_constant * (rest_length - t1);
        Matrix2d T3 = t0 * t0.transpose();
        Matrix2d dF_on_a_dx;

        dF_on_a_dx = t2/t1 * Identity -
                     ( (stiffnes_constant/std::pow(t1,2)) * T3 + (t2/std::pow(t1,3))*T3 );

        // std::cout<<"DF on A DX = "<<dF_on_a_dx<<std::endl;

        Matrix2d dF_on_b_dx;
        dF_on_b_dx = (stiffnes_constant/std::pow(t1,2)) * T3
                     + (t2/std::pow(t1,3))*T3
                     - t2/t1*Identity;

        // std::cout<<"DF on B DX = "<<dF_on_b_dx<<std::endl;

        Matrix2d dF_on_a_dv;
        dF_on_a_dv <<0,0,0,0;
        Matrix2d dF_on_b_dv = dF_on_a_dv;


        Matrix2d A_MAT;

        A_MAT= Identity - time_step * dF_on_a_dv - std::pow(time_step, 2) * dF_on_a_dx;
        Vector2d A_b;
        A_b = time_step * (F_on_a +time_step * dF_on_a_dx * A.velocity);
        Vector2d delta_v_a;
        delta_v_a =  A_MAT.colPivHouseholderQr().solve(A_b);
        // std::cout<<"delta v A = "<<delta_v_a<<std::endl;




        Matrix2d B_MAT;
        B_MAT= Identity - time_step * dF_on_b_dv - std::pow(time_step, 2) * dF_on_b_dx;
        Vector2d B_b;
        B_b = time_step * (F_on_b +time_step * dF_on_b_dx * B.velocity);
        Vector2d delta_v_b;
        delta_v_b =  B_MAT.colPivHouseholderQr().solve(B_b);

        //  std::cout<<"delta v B = "<<delta_v_b<<std::endl;



        A.velocity += delta_v_a;
        B.velocity += delta_v_b;

        // std::cout<<"A vel = " <<A.velocity<<std::endl;
        // std::cout<<"B vel = " <<B.velocity<<std::endl;

        //std::cout<<"A pos prev = " <<A.position<<std::endl;
        // std::cout<<"B pos  prev= " <<B.position<<std::endl;
        // std::cout<<"time step = "<<time_step<<std::endl;

        A.position += time_step * A.velocity;
        B.position += time_step * B.velocity;

        //END TMP
         */








    }

    void UPDATE_OG(std::tuple<Particle&,Particle&,double> &particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        add_brownian_motion(A);
        add_brownian_motion(B);

         A.velocity *=0.1;
        B.velocity *= 0.1;

        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);
        std::cout<<"f c A =  " <<f_c_A<<std::endl;
        std::cout<<"f c B = " <<f_c_B<<std::endl;
        double six_Pi_mu_r_A = ( 6.0 * M_PI * eta * (A.radius*std::pow(10,-6))  * std::pow(10,6) );
        double six_Pi_mu_r_B = ( 6.0 * M_PI * eta * (B.radius*std::pow(10,-6))  * std::pow(10,6) );
        std::cout<<"six pi" <<six_Pi_mu_r_A<<std::endl;
        Vector2d f_drag_A =  A.velocity *  six_Pi_mu_r_A;
        Vector2d f_drag_B =  B.velocity *  six_Pi_mu_r_B;
        std::cout<<"f drag " <<f_drag_A<<std::endl;
        Vector2d f_spring_A = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) * ((A.position - B.position)/(A.position - B.position).norm());
        // f_spring_A;
        Vector2d f_spring_B = -1.0 * f_spring_A;
        std::cout<<"f spring " <<f_spring_A<<std::endl;
        Vector2d F_A = f_c_A   - f_drag_A+ f_spring_A;
        std::cout<<"F A = "<<F_A<<std::endl;
        Vector2d F_B = f_c_B - f_drag_B+ f_spring_B;
        std::cout<<"F B = "<<F_B<<std::endl;
        double root_for_spring = (A.position - B.position).norm();

        Vector2d t0 = A.position - B.position;
        double t1 = t0.norm();
        Matrix2d T2 = t0 * t0.transpose();
        double t3 = stiffnes_constant * (t1 - rest_length);
        Matrix2d d_f_spring_dA;
        d_f_spring_dA = -1.0 * ((stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2 + (t3/t1)* Identity );
        Matrix2d d_f_spring_dB;
        d_f_spring_dB= -1.0 * d_f_spring_dA;

        double d_f_spring_x_dx_A = d_f_spring_dA(0,0);
        double d_f_spring_x_dx_B = d_f_spring_dB(0,0);
        double d_f_spring_x_dy_A = d_f_spring_dA(0,1);
        double d_f_spring_x_dy_B = d_f_spring_dB(0,1);
        double d_f_spring_y_dx_A = d_f_spring_dA(1,0);
        double d_f_spring_y_dx_B = d_f_spring_dB(1,0);
        double d_f_spring_y_dy_A = d_f_spring_dA(1,1);
        double d_f_spring_y_dy_B = d_f_spring_dB(1,1);



        // double d_f_spring_x_dx_A = stiffnes_constant * std::pow(A.position[0] - B.position[0],2) * ( root_for_spring - rest_length) / std::pow(root_for_spring,3)
        //        - stiffnes_constant * (root_for_spring - rest_length) / root_for_spring - stiffnes_constant * std::pow(A.position[0] - B.position[0],2) / std::pow(root_for_spring,2);
        //double d_f_spring_x_dx_B = -1.0 * d_f_spring_x_dx_A;
        //double d_f_spring_x_dy_A  = stiffnes_constant * rest_length * (A.position[0] - B.position[0]) * (A.position[1] - B.position[1]) / std::pow(root_for_spring,3);
        //double d_f_spring_x_dy_B =  -1.0 * d_f_spring_x_dy_A;
        //double d_f_spring_y_dx_A = stiffnes_constant * rest_length * (A.position[0] - B.position[0]) * (A.position[1] - B.position[1]) / std::pow(root_for_spring,3);
        //double d_f_spring_y_dx_B = -1.0 * d_f_spring_y_dx_A;
        //double d_f_spring_y_dy_A = stiffnes_constant * std::pow(A.position[1] - B.position[1],2) * ( root_for_spring - rest_length) / std::pow(root_for_spring,3)
        //     - stiffnes_constant * (root_for_spring - rest_length) / root_for_spring - stiffnes_constant * std::pow(A.position[1] - B.position[1],2) / std::pow(root_for_spring,2);
        //double d_f_spring_y_dy_B = -1.0 * d_f_spring_y_dy_A;



        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;



        //tmp
/*
        Vector2d F_on_a;
        Vector2d F_on_b;

        F_on_a = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) *
                 ((A.position - B.position)/(A.position - B.position).norm());
        // std::cout<<"F on A = "<<F_on_a<<std::endl;
        F_on_b = -1.0 *F_on_a;
        // std::cout<<"F on B = "<<F_on_b<<std::endl;
        Vector2d t0_ = A.position - B.position;
        double t1_ = t0.norm();
        double t2_ = stiffnes_constant * (rest_length - t1);
        Matrix2d T3_ = t0 * t0.transpose();
        Matrix2d dF_on_a_dx;

        dF_on_a_dx = t2_/t1_ * Identity -
                     ( (stiffnes_constant/std::pow(t1_,2)) * T3_ + (t2_/std::pow(t1_,3))*T3_ );

        // std::cout<<"DF on A DX = "<<dF_on_a_dx<<std::endl;

        Matrix2d dF_on_b_dx;
        dF_on_b_dx = (stiffnes_constant/std::pow(t1_,2)) * T3_
                     + (t2_/std::pow(t1_,3))*T3_
                     - t2_/t1_*Identity;

        F_A.setZero();
        F_B.setZero();
        F_A =f_c_A   - f_drag_A + F_on_a;
        F_B =f_c_B - f_drag_B+ F_on_b;
        */

        //endtmp


        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        d_f_dv_A << -six_Pi_mu_r_A,0,0,-six_Pi_mu_r_A;
        Matrix2d d_f_dv_B;
        d_f_dv_B << -six_Pi_mu_r_B,0,0,-six_Pi_mu_r_B;
        // d_f_dv_A;
        // d_f_dv_B;

        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
        //d_fcircle_dxy << 0,-1.0,1.0,0;
        d_f_dx_A = d_fcircle_dxy_A + d_f_spring_dA;
        d_f_dx_B = d_fcircle_dxy_B +d_f_spring_dB;
        //  d_f_dx_A;
        //  d_f_dx_B;
        // d_f_dx_A << d_f_spring_x_dx_A, d_f_circle_x_dy + d_f_spring_x_dy_A, d_f_circle_y_dx + d_f_spring_y_dx_A, d_f_spring_y_dy_A;
        // d_f_dx_B << d_f_spring_x_dx_B, d_f_circle_x_dy + d_f_spring_x_dy_B, d_f_circle_y_dx + d_f_spring_y_dx_B, d_f_spring_y_dy_B;




        /*Matrix2d A_init_A = Identity - time_step * d_f_dv_A - std::pow(time_step,2) * d_f_dx_A;
        Matrix2d A_init_B = Identity - time_step * d_f_dv_B - std::pow(time_step,2) * d_f_dx_B;
        Vector2d b_init_A = time_step * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = time_step * ( F_B   + time_step * d_f_dx_B * B.velocity);*/

        Matrix2d A_init_A = Identity - (time_step) * d_f_dv_A - (std::pow(time_step,2)) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step) * d_f_dv_B - (std::pow(time_step,2)) * d_f_dx_B;

        Vector2d b_init_A = (time_step) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step) * ( F_B   + time_step * d_f_dx_B * B.velocity);



        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);





        // end of get initial delta v

        //NEWTONS METHOD //

        int maxiter = 100;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        Vector2d v_k_minus_1_A = delta_v_init_A - A.velocity;
        Vector2d v_k_minus_1_B = delta_v_init_B - B.velocity;
        int count =0;
        double tol = std::pow(10,(-10));
        double alpha = 0.25;
        double beta = 0.5;


        std::cout<<"AAAAAAAaaaaaaaaaaaaaaa"<<std::endl; /// STILL WRONG; maybe delta_v init should change in line 268, abruch kriterium doesnt change and is wrong aswell
        for(int i = 0; i<maxiter;i++){
            double t = 1.0;

             F_Newton_A[0] = delta_v_old_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_old_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                         - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
             F_Newton_A[1] = delta_v_old_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_old_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                           - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );


           // F_Newton_A[0] = delta_v_old_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_old_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
           // F_Newton_A[1] = delta_v_old_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_old_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );
//this is not OG but new


            // F_Newton_A = (Identity - d_f_dv_A * time_step- d_f_dx_A * std::pow(time_step,2)) * delta_v_old_A - time_step * (F_A + d_f_dx_A * time_step * A.velocity);
            //this is not OG but new new

            Matrix2d Jacobian_A;

            Jacobian_A << 1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_x_dx_A + h* six_Pi_mu_r_A, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_A),
                   -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_A)  ,  1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_y_dy_A + h* six_Pi_mu_r_A;

            //Jacobian_A << 1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(0,0) + (h/A.mass)* six_Pi_mu_r_A, -(h2/A.mass) * (d_f_circle_x_dy + d_f_spring_dA(0,1)), -(h2/A.mass) * (d_f_circle_y_dx + d_f_spring_dA(1,0))  ,  1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(1,1) + (h/A.mass)* six_Pi_mu_r_A;


              Vector2d old_eval_F = eval_F(delta_v_old_A,  h, six_Pi_mu_r_A, h2, d_f_spring_x_dx_A, d_f_spring_y_dx_A,  d_f_spring_x_dy_A,  d_f_spring_y_dy_A, d_f_circle_y_dx, d_f_circle_x_dy, F_A, v_k_minus_1_A);

            Matrix2d A_MAT = Identity - time_step * d_f_dv_A - std::pow(time_step,2) * d_f_dx_A;
            Vector2d b_vec = time_step * ( F_A  + time_step * d_f_dx_A * A.velocity);
            Vector2d delta_v_init_A = A_MAT.colPivHouseholderQr().solve(b_vec);
            //double old_eval = (A_MAT * delta_v_old_A - b_vec).norm();

            double old_eval = F_Newton_A.norm();
            //here

            Vector2d F_Newton_A_new;




            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                            - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
            F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                            - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );

            //F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
            //F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );

            //double new_eval = (A_MAT * delta_v_new_A - b_vec).norm();
            double new_eval = F_Newton_A_new.norm();

            while(new_eval > old_eval + alpha * t * old_eval_F.transpose() *(-Jacobian_A.inverse() * F_Newton_A)){
                t = beta * t;
                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
                  F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                                      - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
                  F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                                      - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );

              //  F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
              //  F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );


                new_eval = F_Newton_A_new.norm();
                //new_eval = (A_MAT * delta_v_new_A - b_vec).norm();

            }

            delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            delta_v_old_A = delta_v_new_A;
            count++;

            if((Jacobian_A.inverse() * F_Newton_A).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_A<<" and delta_v_end = "<<delta_v_old_A<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
                //  std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }


        std::cout<<"Bbbbbbbbbbbbbbbbbbb"<<std::endl;

        for(int i = 0; i<maxiter;i++){
            double t = 1.0;

            F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                            - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
            F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                            - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );

         //   F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
          //  F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );

            std::cout<<" F Newton = "<<F_Newton_B<<std::endl;

            Matrix2d Jacobian_B;
            Jacobian_B << 1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_x_dx_B + h* six_Pi_mu_r_B, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_B), -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_B)  ,  1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_y_dy_B + h* six_Pi_mu_r_B;
            //Jacobian_B << 1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(0,0) + h/B.mass * six_Pi_mu_r_B, -h2/B.mass * (d_f_circle_x_dy + d_f_spring_dB(0,1)), -h2/B.mass * (d_f_circle_y_dx + d_f_spring_dB(1,0))  ,  1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(1,1) + h/B.mass * six_Pi_mu_r_B;



            Vector2d old_eval_F = eval_F(delta_v_old_B,  h, six_Pi_mu_r_B, h2, d_f_spring_x_dx_B, d_f_spring_y_dx_B,  d_f_spring_x_dy_B,  d_f_spring_y_dy_B, d_f_circle_y_dx, d_f_circle_x_dy, F_B, v_k_minus_1_B);
            Vector2d old_eval_F_B = F_Newton_B;
            Matrix2d A_MAT = Identity - time_step * d_f_dv_B - std::pow(time_step,2) * d_f_dx_B;
            Vector2d b_vec = time_step * ( F_B  + time_step * d_f_dx_B * B.velocity);
            Vector2d delta_v_init_B = A_MAT.colPivHouseholderQr().solve(b_vec);
            //double old_eval = (A_MAT * delta_v_old_B - b_vec).norm();


            double old_eval = F_Newton_B.norm();
            Vector2d delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

            Vector2d F_Newton_B_new;

             F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                             - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
             F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                             - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );

          //  F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
           // F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );

            // double new_eval = (A_MAT * delta_v_new_B - b_vec).norm();

            double new_eval = F_Newton_B_new.norm();

            int while_counter = 0;
            while(new_eval > old_eval + alpha * t * old_eval_F_B.transpose() *(-Jacobian_B.inverse() * F_Newton_B)){
                t = beta * t;
                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;


                 F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                                     - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
                 F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                                     - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );

               // F_Newton_B[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
                //F_Newton_B[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );

                //new_eval = (A_MAT * delta_v_new_B - b_vec).norm();
                new_eval = F_Newton_B_new.norm();
                std::cout<<"t = "<<t<<std::endl;
                while_counter++;
                std::cout<<"while counter= "<<while_counter<<std::endl;
            }



            delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

            delta_v_old_B = delta_v_new_B;
            count++;




            if((Jacobian_B.inverse() * F_Newton_B).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_B<<" and delta_v_end = "<<delta_v_old_B<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
                // std::cout<<" value is = "<<(Jacobian_B.inverse() * F_Newton_B).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }

        Vector2d delta_v_A = delta_v_old_A;
        A.velocity += delta_v_A;
        A.position += time_step * A.velocity;

        // Vector2d delta_v_A = delta_v_old_A;
        Vector2d delta_v_B = delta_v_old_B;
        /*  std::cout<<" delta v A 2= "<<delta_v_A<<std::endl;
          std::cout<<" delta v B 2= "<<delta_v_B<<std::endl;

          //END OF NEWTONSMETHOD

          std::cout<<" A vel pre = "<<A.velocity<<std::endl;
          std::cout<<" B vel pre = "<<B.velocity<<std::endl;
          A.velocity += delta_v_A;
          std::cout<<" A vel post = "<<A.velocity<<std::endl;
  */
        B.velocity += delta_v_B;
        /* std::cout<<" B vel post = "<<B.velocity<<std::endl;
         std::cout<<" A pos pre = "<<A.position<<std::endl;
         std::cout<<" B pos pre = "<<B.position<<std::endl;
         //A.velocity *= 0.1;
        // B.velocity *= 0.1;
         A.position += time_step * A.velocity;
         */
        B.position += time_step * B.velocity;
        /*  std::cout<<" A pos post = "<<A.position<<std::endl;
          std::cout<<" B pos post = "<<B.position<<std::endl;
  */


        //TMP

        /*
         Vector2d F_on_a;
        Vector2d F_on_b;

        F_on_a = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) *
                 ((A.position - B.position)/(A.position - B.position).norm());
        // std::cout<<"F on A = "<<F_on_a<<std::endl;
        F_on_b = -1.0 *F_on_a;
        // std::cout<<"F on B = "<<F_on_b<<std::endl;
        Vector2d t0 = A.position - B.position;
        double t1 = t0.norm();
        double t2 = stiffnes_constant * (rest_length - t1);
        Matrix2d T3 = t0 * t0.transpose();
        Matrix2d dF_on_a_dx;

        dF_on_a_dx = t2/t1 * Identity -
                     ( (stiffnes_constant/std::pow(t1,2)) * T3 + (t2/std::pow(t1,3))*T3 );

        // std::cout<<"DF on A DX = "<<dF_on_a_dx<<std::endl;

        Matrix2d dF_on_b_dx;
        dF_on_b_dx = (stiffnes_constant/std::pow(t1,2)) * T3
                     + (t2/std::pow(t1,3))*T3
                     - t2/t1*Identity;

        // std::cout<<"DF on B DX = "<<dF_on_b_dx<<std::endl;

        Matrix2d dF_on_a_dv;
        dF_on_a_dv <<0,0,0,0;
        Matrix2d dF_on_b_dv = dF_on_a_dv;


        Matrix2d A_MAT;

        A_MAT= Identity - time_step * dF_on_a_dv - std::pow(time_step, 2) * dF_on_a_dx;
        Vector2d A_b;
        A_b = time_step * (F_on_a +time_step * dF_on_a_dx * A.velocity);
        Vector2d delta_v_a;
        delta_v_a =  A_MAT.colPivHouseholderQr().solve(A_b);
        // std::cout<<"delta v A = "<<delta_v_a<<std::endl;




        Matrix2d B_MAT;
        B_MAT= Identity - time_step * dF_on_b_dv - std::pow(time_step, 2) * dF_on_b_dx;
        Vector2d B_b;
        B_b = time_step * (F_on_b +time_step * dF_on_b_dx * B.velocity);
        Vector2d delta_v_b;
        delta_v_b =  B_MAT.colPivHouseholderQr().solve(B_b);

        //  std::cout<<"delta v B = "<<delta_v_b<<std::endl;



        A.velocity += delta_v_a;
        B.velocity += delta_v_b;

        // std::cout<<"A vel = " <<A.velocity<<std::endl;
        // std::cout<<"B vel = " <<B.velocity<<std::endl;

        //std::cout<<"A pos prev = " <<A.position<<std::endl;
        // std::cout<<"B pos  prev= " <<B.position<<std::endl;
        // std::cout<<"time step = "<<time_step<<std::endl;

        A.position += time_step * A.velocity;
        B.position += time_step * B.velocity;

        //END TMP
         */








    }

    void UPDATE_damped_spring(std::tuple<Particle&,Particle&,double> &particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        add_brownian_motion(A);
        add_brownian_motion(B);

       //  A.velocity *=0.1;
        //B.velocity *= 0.1;

        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);

        Vector2d f_drag_A =  get_drag(A);
        Vector2d f_drag_B =  get_drag(B);


        Vector2d f_spring_A = get_damped_spring_force(particle_pair);
        Vector2d f_spring_B = -1.0 * f_spring_A;

        Vector2d F_A = f_c_A   - f_drag_A+ f_spring_A;

        Vector2d F_B = f_c_B - f_drag_B+ f_spring_B;

        Matrix2d d_f_spring_dxA;
        d_f_spring_dxA =get_damped_spring_force_dXA(particle_pair);
        Matrix2d d_f_spring_dxB;
        d_f_spring_dxB= -1.0 * d_f_spring_dxA;

        double d_f_spring_x_dx_A = d_f_spring_dxA(0,0);
        double d_f_spring_x_dx_B = d_f_spring_dxB(0,0);
        double d_f_spring_x_dy_A = d_f_spring_dxA(0,1);
        double d_f_spring_x_dy_B = d_f_spring_dxB(0,1);
        double d_f_spring_y_dx_A = d_f_spring_dxA(1,0);
        double d_f_spring_y_dx_B = d_f_spring_dxB(1,0);
        double d_f_spring_y_dy_A = d_f_spring_dxA(1,1);
        double d_f_spring_y_dy_B = d_f_spring_dxB(1,1);

        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;


        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        d_f_dv_A << -A.six_Pi_mu_r,0,0,-A.six_Pi_mu_r;
        Matrix2d d_f_spring_dvA = get_damped_spring_force_dVA(particle_pair);
        d_f_dv_A += d_f_spring_dvA; //Check if this is correct
        Matrix2d d_f_dv_B;
        Matrix2d d_f_spring_dvB = -1.0 * d_f_spring_dvA;
        d_f_dv_B << -B.six_Pi_mu_r,0,0,-B.six_Pi_mu_r;
        d_f_dv_B += d_f_spring_dvB; //Check if this is correct


        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
        //d_fcircle_dxy << 0,-1.0,1.0,0;
        d_f_dx_A = d_fcircle_dxy_A + d_f_spring_dxA;
        d_f_dx_B = d_fcircle_dxy_B +d_f_spring_dxB;

        Matrix2d A_init_A = Identity - (time_step/A.mass) * d_f_dv_A - (std::pow(time_step,2)/A.mass) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step/B.mass) * d_f_dv_B - (std::pow(time_step,2)/B.mass) * d_f_dx_B;

        Vector2d b_init_A = (time_step/A.mass) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step/B.mass) * ( F_B   + time_step * d_f_dx_B * B.velocity);

        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);
        // end of get initial delta v

        //NEWTONS METHOD //

        int maxiter = 100;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        Vector2d v_k_minus_1_A = delta_v_init_A - A.velocity;
        Vector2d v_k_minus_1_B = delta_v_init_B - B.velocity;
        int count =0;
        double tol = std::pow(10,(-10));
        double alpha = 0.25;
        double beta = 0.5;


        std::cout<<"AAAAAAAaaaaaaaaaaaaaaa"<<std::endl;
        for(int i = 0; i<maxiter;i++){
            double t = 1.0;

           /* F_Newton_A[0] = delta_v_old_A[0] * (1.0 + h/A.mass * A.six_Pi_mu_r - h2/A.mass * (d_f_spring_dxA(0,0) + d_fcircle_dxy_A(0,0)))
                    - delta_v_old_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dxA(0,1))
                    -( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dxA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dxA(0,1))) )  );
            F_Newton_A[1] = delta_v_old_A[1] * (1.0 + h/A.mass * A.six_Pi_mu_r - h2/A.mass * (d_f_spring_dxA(1,1) + d_fcircle_dxy_A(1,1)))
                    - delta_v_old_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dxA(1,0))-
                    ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dxA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dxA(1,0))) )  );

*/ //integrate spring dv terms into this later if other version is to slow
            F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A,v_k_minus_1_A);


            // F_Newton_A = (Identity - d_f_dv_A * time_step- d_f_dx_A * std::pow(time_step,2)) * delta_v_old_A - time_step * (F_A + d_f_dx_A * time_step * A.velocity);

            Matrix2d Jacobian_A;

            //Jacobian_A << 1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_x_dx_A + h* six_Pi_mu_r_A, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_A),
            //        -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_A)  ,  1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_y_dy_A + h* six_Pi_mu_r_A;

           // Jacobian_A << 1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(0,0) + (h/A.mass)* six_Pi_mu_r_A, -(h2/A.mass) * (d_f_circle_x_dy + d_f_spring_dA(0,1)), -(h2/A.mass) * (d_f_circle_y_dx + d_f_spring_dA(1,0))  ,  1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(1,1) + (h/A.mass)* six_Pi_mu_r_A;

           Jacobian_A = evaluate_Jacobian( d_f_dv_A, A,  d_f_dx_A);

            double old_eval = F_Newton_A.norm();
            //here

            Vector2d F_Newton_A_new;




            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            /*F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                            - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
            F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                            - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );
*/
            //F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
            //F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );

            F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);
            //double new_eval = (A_MAT * delta_v_new_A - b_vec).norm();
            double new_eval = F_Newton_A_new.norm();

            while(new_eval > old_eval + alpha * t * F_Newton_A.transpose() *(-Jacobian_A.inverse() * F_Newton_A)){
                t = beta * t;
                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
                /*  F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                                      - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
                  F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                                      - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );
  */
                //F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
                //F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );
                F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);

                new_eval = F_Newton_A_new.norm();
                //new_eval = (A_MAT * delta_v_new_A - b_vec).norm();

            }

            delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            delta_v_old_A = delta_v_new_A;
            count++;

            if((Jacobian_A.inverse() * F_Newton_A).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_A<<" and delta_v_end = "<<delta_v_old_A<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
                //  std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }
        Vector2d delta_v_A = delta_v_old_A;
        A.velocity += delta_v_A;
        A.position += time_step * A.velocity;

        std::cout<<"Bbbbbbbbbbbbbbbbbbb"<<std::endl;

        for(int i = 0; i<maxiter;i++){
            double t = 1.0;

            /*F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                            - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
            F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                            - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );
*/
            //F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
            //F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
            F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B,v_k_minus_1_B);
            std::cout<<" F Newton = "<<F_Newton_B<<std::endl;

            Matrix2d Jacobian_B;
            //Jacobian_B << 1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_x_dx_B + h* six_Pi_mu_r_B, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_B), -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_B)  ,  1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_y_dy_B + h* six_Pi_mu_r_B;
            //Jacobian_B << 1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(0,0) + h/B.mass * six_Pi_mu_r_B, -h2/B.mass * (d_f_circle_x_dy + d_f_spring_dB(0,1)), -h2/B.mass * (d_f_circle_y_dx + d_f_spring_dB(1,0))  ,  1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(1,1) + h/B.mass * six_Pi_mu_r_B;
            Jacobian_B = evaluate_Jacobian( d_f_dv_B, B,  d_f_dx_B);


            //Vector2d old_eval_F = eval_F(delta_v_old_B,  h, six_Pi_mu_r_B, h2, d_f_spring_x_dx_B, d_f_spring_y_dx_B,  d_f_spring_x_dy_B,  d_f_spring_y_dy_B, d_f_circle_y_dx, d_f_circle_x_dy, F_B, v_k_minus_1_B);

            double old_eval = F_Newton_B.norm();
            Vector2d delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

            Vector2d F_Newton_B_new;

            /* F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                             - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
             F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                             - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );
 */
           // F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
           // F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
            F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);
            // double new_eval = (A_MAT * delta_v_new_B - b_vec).norm();

            double new_eval = F_Newton_B_new.norm();

            int while_counter = 0;
            while(new_eval > old_eval + alpha * t * F_Newton_B.transpose() *(-Jacobian_B.inverse() * F_Newton_B)){
                t = beta * t;
                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

                /* F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                                     - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
                 F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                                     - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );
 */
                //F_Newton_B[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
               // F_Newton_B[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
                F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);
                //new_eval = (A_MAT * delta_v_new_B - b_vec).norm();
                new_eval = F_Newton_B_new.norm();
                std::cout<<"t = "<<t<<std::endl;
                while_counter++;
                std::cout<<"while counter= "<<while_counter<<std::endl;
            }



            delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

            delta_v_old_B = delta_v_new_B;
            count++;




            if((Jacobian_B.inverse() * F_Newton_B).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_B<<" and delta_v_end = "<<delta_v_old_B<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
                // std::cout<<" value is = "<<(Jacobian_B.inverse() * F_Newton_B).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }



        // Vector2d delta_v_A = delta_v_old_A;
        Vector2d delta_v_B = delta_v_old_B;
        /*  std::cout<<" delta v A 2= "<<delta_v_A<<std::endl;
          std::cout<<" delta v B 2= "<<delta_v_B<<std::endl;

          //END OF NEWTONSMETHOD

          std::cout<<" A vel pre = "<<A.velocity<<std::endl;
          std::cout<<" B vel pre = "<<B.velocity<<std::endl;
          A.velocity += delta_v_A;
          std::cout<<" A vel post = "<<A.velocity<<std::endl;
  */
        B.velocity += delta_v_B;
        /* std::cout<<" B vel post = "<<B.velocity<<std::endl;
         std::cout<<" A pos pre = "<<A.position<<std::endl;
         std::cout<<" B pos pre = "<<B.position<<std::endl;
         //A.velocity *= 0.1;
        // B.velocity *= 0.1;
         A.position += time_step * A.velocity;
         */
        B.position += time_step * B.velocity;
        /*  std::cout<<" A pos post = "<<A.position<<std::endl;
          std::cout<<" B pos post = "<<B.position<<std::endl;
  */


        //TMP

        /*
         Vector2d F_on_a;
        Vector2d F_on_b;

        F_on_a = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) *
                 ((A.position - B.position)/(A.position - B.position).norm());
        // std::cout<<"F on A = "<<F_on_a<<std::endl;
        F_on_b = -1.0 *F_on_a;
        // std::cout<<"F on B = "<<F_on_b<<std::endl;
        Vector2d t0 = A.position - B.position;
        double t1 = t0.norm();
        double t2 = stiffnes_constant * (rest_length - t1);
        Matrix2d T3 = t0 * t0.transpose();
        Matrix2d dF_on_a_dx;

        dF_on_a_dx = t2/t1 * Identity -
                     ( (stiffnes_constant/std::pow(t1,2)) * T3 + (t2/std::pow(t1,3))*T3 );

        // std::cout<<"DF on A DX = "<<dF_on_a_dx<<std::endl;

        Matrix2d dF_on_b_dx;
        dF_on_b_dx = (stiffnes_constant/std::pow(t1,2)) * T3
                     + (t2/std::pow(t1,3))*T3
                     - t2/t1*Identity;

        // std::cout<<"DF on B DX = "<<dF_on_b_dx<<std::endl;

        Matrix2d dF_on_a_dv;
        dF_on_a_dv <<0,0,0,0;
        Matrix2d dF_on_b_dv = dF_on_a_dv;


        Matrix2d A_MAT;

        A_MAT= Identity - time_step * dF_on_a_dv - std::pow(time_step, 2) * dF_on_a_dx;
        Vector2d A_b;
        A_b = time_step * (F_on_a +time_step * dF_on_a_dx * A.velocity);
        Vector2d delta_v_a;
        delta_v_a =  A_MAT.colPivHouseholderQr().solve(A_b);
        // std::cout<<"delta v A = "<<delta_v_a<<std::endl;




        Matrix2d B_MAT;
        B_MAT= Identity - time_step * dF_on_b_dv - std::pow(time_step, 2) * dF_on_b_dx;
        Vector2d B_b;
        B_b = time_step * (F_on_b +time_step * dF_on_b_dx * B.velocity);
        Vector2d delta_v_b;
        delta_v_b =  B_MAT.colPivHouseholderQr().solve(B_b);

        //  std::cout<<"delta v B = "<<delta_v_b<<std::endl;



        A.velocity += delta_v_a;
        B.velocity += delta_v_b;

        // std::cout<<"A vel = " <<A.velocity<<std::endl;
        // std::cout<<"B vel = " <<B.velocity<<std::endl;

        //std::cout<<"A pos prev = " <<A.position<<std::endl;
        // std::cout<<"B pos  prev= " <<B.position<<std::endl;
        // std::cout<<"time step = "<<time_step<<std::endl;

        A.position += time_step * A.velocity;
        B.position += time_step * B.velocity;

        //END TMP
         */








    }

    void UPDATE_damped_spring_with_interupdate(std::tuple<Particle&,Particle&,double> &particle_pair){
        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        add_brownian_motion(A);
        add_brownian_motion(B);
        std::cout<<"added brownian motion"<<std::endl;
        std::cout<< " A position = "<<A.position<<"  B position = "<<B.position<<std::endl<<"  A velocity = "<<A.velocity<<"  B.velocity = "<<B.velocity<<std::endl;

        // A.velocity *= 0.1;
       // B.velocity *= 0.1;

        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);

        Vector2d f_drag_A =  get_drag(A);
        Vector2d f_drag_B =  get_drag(B);


        Vector2d f_spring_A = get_damped_spring_force(particle_pair);
        Vector2d f_spring_B = -1.0 * f_spring_A;

        Vector2d F_A = f_c_A   - f_drag_A+ f_spring_A;
        std::cout<<"F_A = "<<F_A<<" with f_circle = "<<f_c_A<<" and f_drag= "<<f_drag_A<<"and f_spring = "<<f_spring_A<<std::endl;


        Vector2d F_B = f_c_B - f_drag_B+ f_spring_B;
        std::cout<<"F_B = "<<F_B<<" with f_circle = "<<f_c_B<<" and f_drag= "<<f_drag_B<<"and f_spring = "<<f_spring_B<<std::endl;


        Matrix2d d_f_spring_dxA;
        d_f_spring_dxA =get_damped_spring_force_dXA(particle_pair);

        Matrix2d d_f_spring_dxB;
        d_f_spring_dxB= -1.0 * d_f_spring_dxA;


        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;




        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        d_f_dv_A << -A.six_Pi_mu_r,0,0,-A.six_Pi_mu_r;
        Matrix2d d_f_spring_dvA = get_damped_spring_force_dVA(particle_pair);
        d_f_dv_A += d_f_spring_dvA; //Check if this is correct
        std::cout<<" d F/dvA = "<<d_f_dv_A<<" with d_fdrag/dv= "<<(d_f_dv_A-d_f_spring_dvA)<<" and d_fspring/dvA = "<<d_f_spring_dvA<<std::endl;
        Matrix2d d_f_dv_B;
        Matrix2d d_f_spring_dvB = -1.0 * d_f_spring_dvA;
        d_f_dv_B << -B.six_Pi_mu_r,0,0,-B.six_Pi_mu_r;
        d_f_dv_B += d_f_spring_dvB; //Check if this is correct
        std::cout<<" d F/dvB = "<<d_f_dv_B<<" with d_fdrag/dv= "<<(d_f_dv_B-d_f_spring_dvB)<<" and d_fspring/dvB = "<<d_f_spring_dvB<<std::endl;

        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;
        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);
        //d_fcircle_dxy << 0,-1.0,1.0,0;
        d_f_dx_A = d_fcircle_dxy_A + d_f_spring_dxA;
        std::cout<<" d_f/dxA = "<<d_f_dx_A<<" with dfcircle/dx = "<<d_fcircle_dxy_A<<" and dfspring/dx = "<<d_f_spring_dxA<<std::endl;
        d_f_dx_B = d_fcircle_dxy_B +d_f_spring_dxB;
        std::cout<<" d_f/dxB = "<<d_f_dx_B<<" with dfcircle/dx = "<<d_fcircle_dxy_B<<" and dfspring/dx = "<<d_f_spring_dxB<<std::endl;
        Matrix2d A_init_A = Identity - (time_step/A.mass) * d_f_dv_A - (std::pow(time_step,2)/A.mass) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step/B.mass) * d_f_dv_B - (std::pow(time_step,2)/B.mass) * d_f_dx_B;

        Vector2d b_init_A = (time_step/A.mass) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step/B.mass) * ( F_B   + time_step * d_f_dx_B * B.velocity);

        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        std::cout<<"delta v init A = "<<delta_v_init_A<<std::endl;
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);
        std::cout<<"delta v init B = "<<delta_v_init_B<<std::endl;
        // end of get initial delta v

        //NEWTONS METHOD //

        int maxiter = 100;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        Vector2d v_k_minus_1_A = delta_v_init_A - A.velocity;
        Vector2d v_k_minus_1_B = delta_v_init_B - B.velocity;
        int count =0;
        double tol = std::pow(10,(-10));
        double alpha = 0.25;
        double beta = 0.5;



        std::cout<<"AAAAAAAaaaaaaaaaaaaaaa"<<std::endl;
        for(int i = 0; i<maxiter;i++){
            double t = 1.0;
            F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A,v_k_minus_1_A);
            Matrix2d Jacobian_A;
            Jacobian_A = evaluate_Jacobian( d_f_dv_A, A,  d_f_dx_A);
            double old_eval = F_Newton_A.norm();
            std::cout<<"old eval before while = "<<old_eval<<std::endl;
            Vector2d F_Newton_A_new;
            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
            F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);
            double new_eval = F_Newton_A_new.norm();
            std::cout<<"new_eval before while = "<<new_eval<<std::endl;

            while(new_eval > old_eval /*+ alpha * t * F_Newton_A.transpose() *(-Jacobian_A.inverse() * F_Newton_A)*/){ //decrement t until step gives better result
                t = beta * t;
                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
                 F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A,v_k_minus_1_A);
                new_eval = F_Newton_A_new.norm();
                std::cout<<" t = "<<t<<"new_eval= "<<new_eval<<" || ";
            }

            delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A; // Newton update
            delta_v_old_A = delta_v_new_A;
            F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A,v_k_minus_1_A);
            count++;

            if((Jacobian_A.inverse() * F_Newton_A).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_A<<" and delta_v_end = "<<delta_v_old_A<<std::endl;
               i = maxiter;
                break;

            }else{
                std::cout<<std::endl;
                  std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<std::endl;
            }

        }
        std::cout<<"vel 1"<<A.velocity<<std::endl<<"pos 1"<<A.position<<std::endl;
        Vector2d delta_v_A = delta_v_old_A;
        A.velocity += delta_v_A;
        A.position += time_step * A.velocity;
        std::cout<<"vel 2"<<A.velocity<<std::endl<<"pos 2"<<A.position<<std::endl;
        std::cout<<"vel 3"<<std::get<0>(particle_pair).velocity<<std::endl<<"pos 3"<<std::get<0>(particle_pair).position<<std::endl;

        F_B = get_F_B(particle_pair);
        d_f_dv_B = get_df_dv_B(particle_pair);
        d_f_dx_B = get_df_dx_B(particle_pair);
        newton_counter++;
        std::cout<<"newton_counter = "<<newton_counter<<std::endl;

        for(int i = 0; i<maxiter;i++){
            double t = 1.0;
            F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B,v_k_minus_1_B);
            Matrix2d Jacobian_B;
             Jacobian_B = evaluate_Jacobian( d_f_dv_B, B,  d_f_dx_B);
             double old_eval = F_Newton_B.norm();
            Vector2d delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
            Vector2d F_Newton_B_new;
            F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);
            double new_eval = F_Newton_B_new.norm();

            int while_counter = 0;
            while(new_eval > old_eval /*+ alpha * t * F_Newton_B.transpose() *(-Jacobian_B.inverse() * F_Newton_B)*/){
                t = beta * t;
                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
                F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B,v_k_minus_1_B);
                new_eval = F_Newton_B_new.norm();
                std::cout<<"t = "<<t<<std::endl;
                while_counter++;
                std::cout<<"while counter= "<<while_counter<<std::endl;
            }

            delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
            delta_v_old_B = delta_v_new_B;
            F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B,v_k_minus_1_B);
            count++;

            if((Jacobian_B.inverse() * F_Newton_B).norm() < tol){
                std::cout<<" Newtons method converged after "<<count<<" itterations with delta v init = "<<delta_v_init_B<<" and delta_v_end = "<<delta_v_old_B<<std::endl;
                break;

            }else{
                std::cout<<std::endl;
                 std::cout<<" value is = "<<(Jacobian_B.inverse() * F_Newton_B).norm()<<std::endl;
                std::cout<<std::endl;
            }
        }

        Vector2d delta_v_B = delta_v_old_B;
        B.velocity += delta_v_B;
        B.position += time_step * B.velocity;
    }

    Vector2d get_drag(Particle& particle){
        return particle.velocity * ( 6.0 * M_PI * eta * (particle.radius*std::pow(10,-6))  * std::pow(10,6) );
    }

    Vector2d get_spring_A(Vector2d A_position, Vector2d B_position,double rest_length){
        Vector2d f_spring_A = (stiffnes_constant * (rest_length - (A_position - B_position).norm()) ) * ((A_position - B_position)/(A_position - B_position).norm());
        return f_spring_A;
    }
    Vector2d get_spring_B(Vector2d A_position, Vector2d B_position,double rest_length){
        Vector2d f_spring_B = (stiffnes_constant * (rest_length - (A_position - B_position).norm()) ) * ((A_position - B_position)/(A_position - B_position).norm());
        return -1.0 *f_spring_B;
    }

    Matrix2d get_d_f_spring_dA(Vector2d A_position, Vector2d B_position,double rest_length){
        Matrix2d Identity = Eigen::MatrixXd::Identity(2,2);
        double root_for_spring = (A_position - B_position).norm();

        Vector2d t0 = A_position - B_position;
        double t1 = t0.norm();
        Matrix2d T2 = t0 * t0.transpose();
        double t3 = stiffnes_constant * (t1 - rest_length);
        Matrix2d d_f_spring_dA;
        d_f_spring_dA = -1.0 * ((stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2 + (t3/t1)* Identity );
        return d_f_spring_dA;
    }

    Matrix2d get_d_f_spring_dB(Vector2d A_position, Vector2d B_position,double rest_length){
        Matrix2d Identity = Eigen::MatrixXd::Identity(2,2);
        double root_for_spring = (A_position - B_position).norm();

        Vector2d t0 = A_position - B_position;
        double t1 = t0.norm();
        Matrix2d T2 = t0 * t0.transpose();
        double t3 = stiffnes_constant * (t1 - rest_length);
        Matrix2d d_f_spring_dA;
        d_f_spring_dA = -1.0 * ((stiffnes_constant/std::pow(t1,2)) * T2 - (t3/std::pow(t1,3)) * T2 + (t3/t1)* Identity );
        return -1.0 * d_f_spring_dA;
    }

    VectorXd evaluate_F_ALL(MatrixXd d_f_dv, MatrixXd d_f_dx, VectorXd delta_v, VectorXd forces){
        VectorXd v_k_minus_1 =  velocity_vector;
        MatrixXd Identity = MatrixXd::Identity(size,size);
        MatrixXd A_init = Identity - ((time_step *mass_matrix.inverse()) * d_f_dv) - ((std::pow(time_step,2) * mass_matrix.inverse()) * d_f_dx);
        VectorXd b_init = (time_step * mass_matrix.inverse()) * ( forces  + time_step * d_f_dx * v_k_minus_1);
        VectorXd F;
        F = A_init * delta_v  - b_init;
        return F;
    }


    Vector2d evaluate_F(Particle A, Matrix2d d_f_dv_A, Matrix2d d_f_dx_A, Vector2d delta_v, Vector2d forces, Vector2d delta_v_init){
        // Changed from v_k_minus_1 = delta_v_init - A.velocity; to v_k_minus_1 =  A.velocity;
        Vector2d v_k_minus_1 =  A.velocity;
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Matrix2d A_init_A = Identity - ((time_step/A.mass) * d_f_dv_A) - ((std::pow(time_step,2)/A.mass) * d_f_dx_A);
        Vector2d b_init_A = (time_step/A.mass) * ( forces  + time_step * d_f_dx_A * v_k_minus_1);
        Vector2d F;
        F = A_init_A * delta_v  - b_init_A;
        return F;
     /*  double h = time_step;
       double h2 = std::pow(time_step,2);
        F_Newton_B[0] = delta_v[0] * (1.0 + h/A.mass * A.six_Pi_mu_r - h2/A.mass * (d_f_dx_A(0,0))) - delta_v[1] * h2/A.mass * (d_f_dx_A(0,1))- ( h/A.mass * (  forces[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
        F_Newton_B[1] = delta_v[1] * (1.0 + h/A.mass * A.six_Pi_mu_r - h2/A.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
*/


    }

    //stupid name change later
    void UPDATE_AB_TOGETHER(std::tuple<Particle&,Particle&,double> &particle_pair){

        double rest_length = std::get<2>(particle_pair);
        Matrix2d Identity = MatrixXd::Identity(2,2);
        Particle& A = std::get<0>(particle_pair);
        Particle& B = std::get<1>(particle_pair);


        add_brownian_motion(A);
        add_brownian_motion(B);

       // A.velocity *=0.1;
       // B.velocity *= 0.1;
        ///WHYYYY???????




        // COMPUTATIONS FOR FIRST PARTICLE A //

        Vector2d f_c_A = A.circle_force2(A);
        Vector2d f_c_B = B.circle_force2(B);

        double six_Pi_mu_r_A = ( 6.0 * M_PI * eta * (A.radius*std::pow(10,-6))  * std::pow(10,6) );
        double six_Pi_mu_r_B = ( 6.0 * M_PI * eta * (B.radius*std::pow(10,-6))  * std::pow(10,6) );

        Vector2d f_drag_A =  A.velocity *  six_Pi_mu_r_A;
        Vector2d f_drag_B =  B.velocity *  six_Pi_mu_r_B;

        Vector2d f_spring_A = get_spring_A(A.position,B.position,rest_length);
        Vector2d f_spring_B = -1.0 * f_spring_A;




        Vector2d F_A = f_c_A   - f_drag_A+ f_spring_A;

        Vector2d F_B = f_c_B - f_drag_B+ f_spring_B;



        Matrix2d d_f_spring_dA;
        d_f_spring_dA = get_d_f_spring_dA(A.position,B.position,rest_length);
        Matrix2d d_f_spring_dB;
        d_f_spring_dB= -1.0 * d_f_spring_dA;
        double d_f_spring_x_dx_A = d_f_spring_dA(0,0);
        double d_f_spring_x_dx_B = d_f_spring_dB(0,0);
        double d_f_spring_x_dy_A = d_f_spring_dA(0,1);
        double d_f_spring_x_dy_B = d_f_spring_dB(0,1);
        double d_f_spring_y_dx_A = d_f_spring_dA(1,0);
        double d_f_spring_y_dx_B = d_f_spring_dB(1,0);
        double d_f_spring_y_dy_A = d_f_spring_dA(1,1);
        double d_f_spring_y_dy_B = d_f_spring_dB(1,1);

        double d_f_circle_x_dy = -1.0;
        double d_f_circle_y_dx = 1.0;

        // Get initial delta v to use for newtons method
        Matrix2d d_f_dv_A;
        d_f_dv_A << -six_Pi_mu_r_A,0,0,-six_Pi_mu_r_A;
        Matrix2d d_f_dv_B;
        d_f_dv_B << -six_Pi_mu_r_B,0,0,-six_Pi_mu_r_B;

        Matrix2d d_f_dx_A;
        Matrix2d d_f_dx_B;

        Matrix2d d_fcircle_dxy_A = A.circle_dforce_dx2(A);
        Matrix2d d_fcircle_dxy_B = B.circle_dforce_dx2(B);

        d_f_dx_A = d_fcircle_dxy_A + d_f_spring_dA;
        d_f_dx_B = d_fcircle_dxy_B +d_f_spring_dB;

        Matrix2d A_init_A = Identity - (time_step/A.mass) * d_f_dv_A - (std::pow(time_step,2)/A.mass) * d_f_dx_A;
        Matrix2d A_init_B = Identity - (time_step/B.mass) * d_f_dv_B - (std::pow(time_step,2)/B.mass) * d_f_dx_B;

        Vector2d b_init_A = (time_step/A.mass) * ( F_A  + time_step * d_f_dx_A * A.velocity);
        Vector2d b_init_B = (time_step/B.mass) * ( F_B   + time_step * d_f_dx_B * B.velocity);

        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);
       // delta_v_init_A *=0.1;
        //delta_v_init_B *=0.1;

        //NEWTONS METHOD //

        std::cout<<"mass = "<<A.mass<<" and "<<B.mass<<std::endl;
        int maxiter = 100;
        double h2 = std::pow(time_step,2);
        double h = time_step;
        Vector2d F_Newton_A;
        Vector2d F_Newton_B;
        Matrix2d Jacobian_A;
        Matrix2d Jacobian_B;
        Vector2d delta_v_old_A = delta_v_init_A;
        Vector2d delta_v_old_B = delta_v_init_B;
        double testnorm = delta_v_old_A.norm();
        double testnorm2 = delta_v_old_B.norm();
        Vector2d tmp_position_A = A.position;
        Vector2d tmp_position_B = B.position;


        Vector2d v_k_minus_1_A = delta_v_init_A - A.velocity;
        Vector2d v_k_minus_1_B = delta_v_init_B - B.velocity;
        int count =0;
        double tol = std::pow(10,(-14));
        double alpha = 0.25;
        double beta = 0.5;

       for(int i = 0; i<maxiter;i++){
            double t = 1.0;
            //Vector2d evaluate_F(Particle& A, Vector2d d_f_dv_A, Vector2d d_f_dx_A, Vector2d delta_v, Vector2d forces){
           // F_Newton_A = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_old_A,F_A);
           // F_Newton_B = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_old_B,F_B);


          //  d_f_spring_dA = get_d_f_spring_dA(tmp_position_A,tmp_position_B,rest_length);


            F_Newton_A[0] = delta_v_old_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_old_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
            F_Newton_A[1] = delta_v_old_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_old_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );
           // Jacobian_A = Identity - d_f_dv_A * time_step/A.mass - d_f_dx_A * std::pow(time_step,2)/A.mass;
           // Jacobian_B = Identity - d_f_dv_B * time_step/B.mass - d_f_dx_B * std::pow(time_step,2)/B.mass;
            Jacobian_A << 1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(0,0) + (h/A.mass)* six_Pi_mu_r_A, -(h2/A.mass) * (d_f_circle_x_dy + d_f_spring_dA(0,1)), -(h2/A.mass) * (d_f_circle_y_dx + d_f_spring_dA(1,0))  ,  1.0 + (h/A.mass) * six_Pi_mu_r_A - (h2/A.mass) * d_f_spring_dA(1,1) + (h/A.mass)* six_Pi_mu_r_A;
            Vector2d old_eval_F_A = F_Newton_A;
            double old_eval_A = F_Newton_A.norm();
            Vector2d F_Newton_A_new;
            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
          //  tmp_position_A = A.position + time_step * (delta_v_new_A + A.velocity);


          //  d_f_spring_dB= get_d_f_spring_dB(tmp_position_A,tmp_position_B,rest_length);

            F_Newton_B[0] = delta_v_old_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_old_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
            F_Newton_B[1] = delta_v_old_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_old_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
            Jacobian_B << 1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(0,0) + h/B.mass * six_Pi_mu_r_B, -h2/B.mass * (d_f_circle_x_dy + d_f_spring_dB(0,1)), -h2/B.mass * (d_f_circle_y_dx + d_f_spring_dB(1,0))  ,  1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * d_f_spring_dB(1,1) + h/B.mass * six_Pi_mu_r_B;
            Vector2d old_eval_F_B = F_Newton_B;
            double old_eval_B = F_Newton_B.norm();
            Vector2d F_Newton_B_new;
            Vector2d delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
            //tmp_position_B = B.position + time_step * (delta_v_new_B + B.velocity);

           // F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A);
           // F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B);

           // d_f_spring_dB= get_d_f_spring_dB(tmp_position_A,tmp_position_B,rest_length);
           // d_f_spring_dA = get_d_f_spring_dA(tmp_position_A,tmp_position_B,rest_length);


            F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
            F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );

            F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
            F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );
            std::cout<<" F Newton = "<<F_Newton_B<<std::endl;


            double new_eval_A = F_Newton_A_new.norm();
            double new_eval_B = F_Newton_B_new.norm();

            while((new_eval_A > old_eval_A + alpha * t * old_eval_F_A.transpose() *(-Jacobian_A.inverse() * F_Newton_A)) || (new_eval_B > old_eval_B + alpha * t * old_eval_F_B.transpose() *(-Jacobian_B.inverse() * F_Newton_B))){
                t = beta * t;

                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
                //tmp_position_A = A.position + time_step * (delta_v_new_A + A.velocity);
               // d_f_spring_dA = get_d_f_spring_dA(tmp_position_A,tmp_position_B,rest_length);

              //  F_Newton_A_new = evaluate_F(A,d_f_dv_A,d_f_dx_A,delta_v_new_A,F_A);

                F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2/A.mass * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))- ( h/A.mass * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_dA(0,0) + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_dA(0,1))) )  );
                F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h/A.mass * six_Pi_mu_r_A - h2/A.mass * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2/A.mass * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))- ( h/A.mass * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_dA(1,1) + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_dA(1,0))) )  );

                new_eval_A = F_Newton_A_new.norm();

                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
               // tmp_position_B = B.position + time_step * (delta_v_new_B + B.velocity);
               // d_f_spring_dB = get_d_f_spring_dB(tmp_position_A,tmp_position_B,rest_length);
               // F_Newton_B_new = evaluate_F(B,d_f_dv_B,d_f_dx_B,delta_v_new_B,F_B);

                F_Newton_B[0] = delta_v_new_B[0] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2/B.mass * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))- ( h/B.mass * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_dB(0,0) + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_dB(0,1))) )  );
                F_Newton_B[1] = delta_v_new_B[1] * (1.0 + h/B.mass * six_Pi_mu_r_B - h2/B.mass * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2/B.mass * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))- ( h/B.mass * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_dB(1,1) + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_dB(1,0))) )  );

                new_eval_B = F_Newton_B_new.norm();

                std::cout<<"still in here"<<std::endl;

            }

            delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
           // tmp_position_A = A.position + time_step * (delta_v_new_A + A.velocity);
            delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;
           // tmp_position_B = B.position + time_step * (delta_v_new_B + B.velocity);

            delta_v_old_A = delta_v_new_A;
            delta_v_old_B = delta_v_new_B;
            count++;
            std::cout<<"count "<<count<<std::endl;

            if( ((Jacobian_A.inverse() * F_Newton_A).norm() < tol) && ((Jacobian_B.inverse() * F_Newton_B).norm() < tol) ){

                break;

            }else{

                std::cout<<std::endl;
                  std::cout<<" value is = "<<(Jacobian_A.inverse() * F_Newton_A).norm()<<std::endl;
                std::cout<<std::endl;

            }
        }

        Vector2d delta_v_A = delta_v_old_A;
        Vector2d delta_v_B = delta_v_old_B;
        //END OF NEWTONSMETHOD

        A.velocity += delta_v_A;
        B.velocity += delta_v_B;

        A.position += time_step * A.velocity;
        B.position += time_step * B.velocity;
    }




    void compute_spring_force_and_update_Particles_simple(std::tuple<Particle&,Particle&,double> &particle_pair){



        Particle& A = std::get<0>(particle_pair);

        Particle& B = std::get<1>(particle_pair);

/*
        std::cout<<"A vel = " <<A.velocity<<std::endl;
        std::cout<<"B vel = " <<B.velocity<<std::endl;


        std::cout<<"A pos = " <<A.position<<std::endl;
        std::cout<<"B pos = " <<B.position<<std::endl;*/

//only for noww



    A.velocity *= 0.1;
    B.velocity *=0.1;


        double rest_length = std::get<2>(particle_pair);
        Vector2d F_on_a;
        Vector2d F_on_b;
        Matrix2d Identity = MatrixXd::Identity(2,2);
        F_on_a = (stiffnes_constant * (rest_length - (A.position - B.position).norm()) ) *
                ((A.position - B.position)/(A.position - B.position).norm());
       // std::cout<<"F on A = "<<F_on_a<<std::endl;
        F_on_b = -1.0 *F_on_a;
       // std::cout<<"F on B = "<<F_on_b<<std::endl;
        Vector2d t0 = A.position - B.position;
        double t1 = t0.norm();
        double t2 = stiffnes_constant * (rest_length - t1);
        Matrix2d T3 = t0 * t0.transpose();
        Matrix2d dF_on_a_dx;

        dF_on_a_dx = t2/t1 * Identity -
                ( (stiffnes_constant/std::pow(t1,2)) * T3 + (t2/std::pow(t1,3))*T3 );

       // std::cout<<"DF on A DX = "<<dF_on_a_dx<<std::endl;

        Matrix2d dF_on_b_dx;
        dF_on_b_dx = (stiffnes_constant/std::pow(t1,2)) * T3
                + (t2/std::pow(t1,3))*T3
                - t2/t1*Identity;

       // std::cout<<"DF on B DX = "<<dF_on_b_dx<<std::endl;

        Matrix2d dF_on_a_dv;
        dF_on_a_dv <<0,0,0,0;
        Matrix2d dF_on_b_dv = dF_on_a_dv;


        Matrix2d A_MAT;

        A_MAT= Identity - time_step * dF_on_a_dv - std::pow(time_step, 2) * dF_on_a_dx;
        Vector2d A_b;
        A_b = time_step * (F_on_a +time_step * dF_on_a_dx *A.velocity);
        Vector2d delta_v_a;
        delta_v_a =  A_MAT.colPivHouseholderQr().solve(A_b);
       // std::cout<<"delta v A = "<<delta_v_a<<std::endl;




        Matrix2d B_MAT;
        B_MAT= Identity - time_step * dF_on_b_dv - std::pow(time_step, 2) * dF_on_b_dx;
        Vector2d B_b;
        B_b = time_step * (F_on_b +time_step * dF_on_b_dx  *B.velocity);
        Vector2d delta_v_b;
        delta_v_b =  B_MAT.colPivHouseholderQr().solve(B_b);

      //  std::cout<<"delta v B = "<<delta_v_b<<std::endl;



        A.velocity += delta_v_a;
        B.velocity += delta_v_b;

       // std::cout<<"A vel = " <<A.velocity<<std::endl;
       // std::cout<<"B vel = " <<B.velocity<<std::endl;

        //std::cout<<"A pos prev = " <<A.position<<std::endl;
       // std::cout<<"B pos  prev= " <<B.position<<std::endl;
       // std::cout<<"time step = "<<time_step<<std::endl;

        A.position += time_step * A.velocity;
        B.position += time_step * B.velocity;

      //  std::cout<<"A pos post= " <<A.position<<std::endl;
      //  std::cout<<"B pos post= " <<B.position<<std::endl;

        check_Boundaries(A);
         update_particle(A);
        add_viscous_drag(A);
        add_brownian_motion(A);

        check_Boundaries(B);
         update_particle(B);
        add_viscous_drag(B);
        add_brownian_motion(B);




    }





    //An OLD version that isnt simplified and produces instable/wrong results
    void connect_particles(Particle &A, Particle &B){
        double rest_length = A.radius + B.radius;
        Vector2d xA = A.position;
        Vector2d xB = B.position;
        Vector2d vA = A.velocity;
        Vector2d vB = B.velocity;
        Vector2d force_on_A_by_B;
        Vector2d force_on_B_by_A;
        Matrix2d Identity = MatrixXd::Identity(2,2);
        double first_term = stiffnes_constant * ((xA - xB).norm() - rest_length);
        Vector2d second_term = damping_constant * (vA - vB);
        Vector2d third_term = (xA - xB);
        Vector2d last_term = (xA - xB)/(xA - xB).norm();
        force_on_A_by_B = - 1.0 * (first_term + second_term.dot(third_term) ) * last_term;
        //force_on_A_by_B = -1.0 * ( stiffnes_constant * ((xA - xB).norm() - rest_length) + (damping_constant * ((vA - vB) * (xA - xB).transpose())) ) * (xA - xB)/(xA - xB).norm();
        force_on_B_by_A = -force_on_A_by_B;

        Vector2d t0 = xA -xB;
        double t1 = t0.norm();
        double t2 = std::pow(t1,2);
        Vector2d t3 = vA - vB;
        Matrix2d T4 = t0 * t0.transpose();
        double t5 = t3.transpose() * t0;
        double t6 = stiffnes_constant * (t1 - rest_length) + (damping_constant * t5)/t1;
        Matrix2d dfA_dxA = -1.0 *(stiffnes_constant/t2 * T4
                + damping_constant/t2 * t0 * t3.transpose()
                - ((damping_constant * t5)/std::pow(t1,4) ) * T4
                - (t6/(std::pow(t1,3))) * T4 + t6/t1 * Identity);

        Matrix2d dfA_dvA = -1.0 *((damping_constant/std::pow(t0.norm(),2)) * t0 * t0.transpose());
        Matrix2d dfB_dxB = -1.0 * dfA_dxA;
        Matrix2d dfB_dvB = -1.0 * dfA_dvA;


        Eigen::Matrix2d MAT_A;

        MAT_A = Identity - time_step * dfA_dvA - std::pow(time_step, 2) * dfA_dxA;
        Vector2d b;
        b = time_step * force_on_A_by_B +time_step * dfA_dxA * A.velocity;
        Vector2d delta_v_A;

        delta_v_A = MAT_A.colPivHouseholderQr().solve(b);
       // std::cout<<"delta v A = "<<delta_v_A<<std::endl;



        Eigen::Matrix2d MAT_B;

        MAT_B = Identity - time_step * dfB_dvB - std::pow(time_step, 2) * dfB_dxB;
        Vector2d b_B;
        b_B = time_step * force_on_B_by_A +time_step * dfB_dxB * B.velocity;
        Vector2d delta_v_B;

        delta_v_B = MAT_B.colPivHouseholderQr().solve(b_B);
       // std::cout<<"delta v b = "<<delta_v_B<<std::endl;



        A.velocity += delta_v_A;
        B.velocity += delta_v_B;
        A.position += (time_step * A.velocity);
        B.position += (time_step * B.velocity);



    }

    void run_simulation_with_spring_connected_Particles(Particle &A, Particle &B, std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
           // connect_particles(A,B);
            for (int j = 0; j<particles.size(); j++) {
                check_Boundaries(particles[j]);
                update_particle(particles[j]);
                add_viscous_drag(particles[j]);
                add_brownian_motion(particles[j]);
               /* std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;*/
            }
            connect_particles(A,B);
        }
    }
    //doesnt work
    void run_simulation_with_spring_connected_Particles_averaged_movement(Particle &A, Particle &B, std::vector<Particle> &particles) {

        for (int i = 0; i < num_steps; i++) {
            connect_particles(A,B);
            for (int j = 0; j<particles.size(); j++) {
                check_Boundaries(particles[j]);

                compute_update_particle(particles[j]);
                compute_viscous_drag(particles[j]);
                compute_brownian_motion(particles[j]);
                std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
            for(int j = 0; j<particles.size(); j++){
                std::cout<<"A.tmp_drag = "<<A.tmp_drag<<std::endl<<"B.tmp_drag = "<<B.tmp_drag<<std::endl<<" A.tmp_noise = "<<A.tmp_noise<<std::endl<<"B.tmp_noise = " <<B.tmp_noise<<std::endl<<" A.tmp_update = "<<A.tmp_update<<std::endl<<" B.tmp_update = "<<B.tmp_update<<std::endl;

                particles[j].velocity += ((A.tmp_drag + B.tmp_drag)/2.0 + (A.tmp_noise + B.tmp_noise)/2.0 + (A.tmp_update + B.tmp_update)/2.0);
                particles[j].position = particles[j].velocity * time_step;
            }
            connect_particles(A,B);
        }
    }



     void update_particle(Particle &particle) {
        double tmpinit = particle.position[0];
        double tmpinit2 = particle.position[1];
        double velinit = particle.velocity[0];
        double vel2init = particle.velocity[1];

        particle.velocity += particle.delta_v_circle2(particle);
        force_damper[0] = 1.0; //1.0 / (10.0 * log(std::abs(particle.velocity[0])));//0.001; //(1.0/0.5 *std::abs(particle.velocity[0]));
        force_damper[1] =1.0;// 1.0 / (10.0 * log(std::abs(particle.velocity[1]))); //0.001; //(1.0/0.5 *std::abs(particle.velocity[1]));

        particle.velocity[0] = particle.velocity[0] * force_damper[0];
        particle.velocity[1] = particle.velocity[1] * force_damper[1];

        particle.position += (time_step * particle.velocity);
        double vel = particle.velocity[0];
        double vel2 = particle.velocity[1];
        double tmp = particle.position[0];
        double tmp2 = particle.position[1];

    }



    void add_viscous_drag(Particle &particle){
        double tmp = 6.0 * M_PI * eta * (particle.radius*std::pow(10,-6)) * (time_step* std::pow(10,-6));
        tmp = (tmp * std::pow(10,6));
        Vector2d delta_v_drag = particle.velocity * (tmp/(1.0 - tmp));
        particle.velocity += delta_v_drag;
       // std::cout<<"change in velocity = "<<delta_v_drag<<std::endl;
        particle.position += particle.velocity * time_step;
    }



    void add_brownian_motion(Particle &particle){

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution (0.0,1.0);

        double x_rand = distribution(generator);
        double y_rand = distribution(generator);
        long double D = (kB * T)/(3 * M_PI * eta * (particle.radius*2)*std::pow(10,-6) ); //not sure
        long double k = std::sqrt(2 * D * time_step);
        k = noise_damper *k * std::pow(10,6); // because we are in micron scale /
        double dx = x_rand * k;
        double dy = y_rand * k;
        particle.position[0] += dx;
        particle.position[1] += dy;
    }

    void run_simulation(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                check_Boundaries(particles[j]);
                update_particle(particles[j]);
                check_Boundaries(particles[j]);
                add_viscous_drag(particles[j]);
                check_Boundaries(particles[j]);
                add_brownian_motion(particles[j]);
                check_Boundaries(particles[j]);
               // std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
               // std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                //std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                //std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
        }
    }

    void run_simulation_with_forcefield_only(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                update_particle(particles[j]);
                std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
        }
    }

    void run_simulation_with_drag(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {

                update_particle(particles[j]);
                add_viscous_drag(particles[j]); // assumption add the new velocities together
                std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
        }
    }

    void run_simulation_with_brownian_motion(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                update_particle(particles[j]);
                add_brownian_motion(particles[j]);
                std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
        }
    }

    void run_simulation_with_drag_only(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                add_viscous_drag(particles[j]);
                std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
        }
    }

    void run_simulation_with_brownian_motion_only(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                add_brownian_motion(particles[j]);
                std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
            }
        }
    }

    void run_simulation_with_brownian_motion_and_drag_only(std::vector<Particle> &particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                add_viscous_drag(particles[j]);
                add_brownian_motion(particles[j]);
                std::cout<<"start particle "<<j<<std::endl;
                 std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                 std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                 std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                 std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;
                 std::cout<<"end particle "<<j<<std::endl;
            }
        }
    }

    void reset_simulation(std::vector<Particle> &particles){
        for (int j = 0; j<particles.size(); j++) {
            particles[j].position = particles[j].initial_position;
            particles[j].velocity = particles[j].initial_velocity;
        }
    }

    //COMPUTE AND NOT UPDATE FUNCTIONS:

    void compute_viscous_drag(Particle &particle){
        double tmp = 6.0 * M_PI * eta * (particle.radius*std::pow(10,-6)) * (time_step* std::pow(10,-6));
        tmp = (tmp * std::pow(10,6));
        particle.tmp_drag = particle.velocity * (tmp/(1.0 - tmp)); //tmp change in velocity due to drag

    }



    void compute_brownian_motion(Particle &particle){

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution (0.0,1.0);

        double x_rand = distribution(generator);
        double y_rand = distribution(generator);
        long double D = (kB * T)/(3 * M_PI * eta * (particle.radius*2)*std::pow(10,-6) ); //not sure
        long double k = std::sqrt(2 * D * time_step);
        k = noise_damper *k * std::pow(10,6); // because we are in micron scale /
        double dx = x_rand * k;
        double dy = y_rand * k;
        particle.tmp_noise = {dx,dy};

    }
    void compute_update_particle(Particle &particle) {
        Vector2d old_vel = particle.velocity;

        particle.velocity += particle.delta_v_circle(particle);
        force_damper[0] = 0.1; //1.0 / (10.0 * log(std::abs(particle.velocity[0])));//0.001; //(1.0/0.5 *std::abs(particle.velocity[0]));
        force_damper[1] =0.1;// 1.0 / (10.0 * log(std::abs(particle.velocity[1]))); //0.001; //(1.0/0.5 *std::abs(particle.velocity[1]));

        particle.velocity[0] = particle.velocity[0] * force_damper[0];
        particle.velocity[1] = particle.velocity[1] * force_damper[1];

        particle.tmp_update = (particle.velocity - old_vel);

    }

    //END OF COMPUTE AND NOT UPDATE FUNCTIONS:

};