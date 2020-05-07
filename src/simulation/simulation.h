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

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::Matrix2d;
using Eigen::Matrix;
using Eigen::MatrixXd;

class Simulation{
public:
    int num_steps = 1;
    double time_step = Particle::time_step;
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
    double stiffnes_constant = 20.0;
    //double rest_length = 2 * Particle::radius_for_spring;
    double damping_constant = 0.0;







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
        F_Newton_A[0] = delta_v_init_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_x_dx_A) - delta_v_init_A[1] * h2 * (d_f_circle_x_dy + d_f_spring_x_dy_A)
                        - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * d_f_spring_x_dx_A + v_k_minus_1_A[1] * (d_f_circle_x_dy + d_f_spring_x_dy_A)) )  );
        F_Newton_A[1] = delta_v_init_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_y_dy_A) - delta_v_init_A[0] * h2 * (d_f_circle_y_dx + d_f_spring_y_dx_A)
                        - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * d_f_spring_y_dy_A + v_k_minus_1_A[0] * (d_f_circle_y_dx + d_f_spring_y_dx_A)) )  );
        return F_Newton_A;
    }

    //The latest version
    void run_simulation_for_connected_Particles(std::vector<std::tuple <Particle&,Particle&,double> > &connected_particles){
        for(int i = 0; i<num_steps;i++) {
            for (auto& particle_pair : connected_particles) {
                UPDATE(particle_pair);
            }
        }
    }


    void UPDATE(std::tuple<Particle&,Particle&,double> &particle_pair){
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




        Matrix2d A_init_A = Identity - time_step * d_f_dv_A - std::pow(time_step,2) * d_f_dx_A;
        std::cout<<" A init A = "<<A_init_A<<std::endl;
        Matrix2d A_init_B = Identity - time_step * d_f_dv_B - std::pow(time_step,2) * d_f_dx_B;
        std::cout<<" A init  B= "<<A_init_B<<std::endl;

        Vector2d b_init_A = time_step * ( F_A  + time_step * d_f_dx_A * A.velocity);
        std::cout<<"b init A = "<<b_init_A<<std::endl;
        Vector2d b_init_B = time_step * ( F_B   + time_step * d_f_dx_B * B.velocity);
        std::cout<<"b init B = "<<b_init_B<<std::endl;

        Vector2d delta_v_init_A = A_init_A.colPivHouseholderQr().solve(b_init_A);
        Vector2d delta_v_init_B = A_init_B.colPivHouseholderQr().solve(b_init_B);
        std::cout<<" delta v A = "<<delta_v_init_A<<std::endl;
        std::cout<<" delta v B = "<<delta_v_init_B<<std::endl;




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



            Matrix2d Jacobian_A;

            Jacobian_A << 1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_x_dx_A + h* six_Pi_mu_r_A, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_A),
                    -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_A)  ,  1.0 + h * six_Pi_mu_r_A - h2 * d_f_spring_y_dy_A + h* six_Pi_mu_r_A;



            Vector2d old_eval_F = eval_F(delta_v_old_A,  h, six_Pi_mu_r_A, h2, d_f_spring_x_dx_A, d_f_spring_y_dx_A,  d_f_spring_x_dy_A,  d_f_spring_y_dy_A, d_f_circle_y_dx, d_f_circle_x_dy, F_A, v_k_minus_1_A);
            Matrix2d A_MAT = Identity - time_step * d_f_dv_A - std::pow(time_step,2) * d_f_dx_A;
            Vector2d b_vec = time_step * ( F_A  + time_step * d_f_dx_A * A.velocity);
            Vector2d delta_v_init_A = A_MAT.colPivHouseholderQr().solve(b_vec);
            //double old_eval = (A_MAT * delta_v_old_A - b_vec).norm();
            double old_eval = F_Newton_A.norm();

            Vector2d F_Newton_A_new;




            Vector2d delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;

            F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                            - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
            F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                            - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );

            //double new_eval = (A_MAT * delta_v_new_A - b_vec).norm();
            double new_eval = F_Newton_A_new.norm();

            while(new_eval > old_eval + alpha * t * old_eval_F.transpose() *(-Jacobian_A.inverse() * F_Newton_A)){
                t = beta * t;
                delta_v_new_A = delta_v_old_A -  t *Jacobian_A.inverse() * F_Newton_A;
                F_Newton_A_new[0] = delta_v_new_A[0] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0))) - delta_v_new_A[1] * h2 * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)
                                    - ( h * (  F_A[0] + h * (v_k_minus_1_A[0] * (d_f_spring_x_dx_A + d_fcircle_dxy_A(0,0)) + v_k_minus_1_A[1] * (d_fcircle_dxy_A(0,1) + d_f_spring_x_dy_A)) )  );
                F_Newton_A_new[1] = delta_v_new_A[1] * (1.0 + h * six_Pi_mu_r_A - h2 * (d_f_spring_y_dy_A+d_fcircle_dxy_A(1,1))) - delta_v_new_A[0] * h2 * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)
                                    - ( h * (  F_A[1] + h * (v_k_minus_1_A[1] * (d_f_spring_y_dy_A + d_fcircle_dxy_A(1,1)) + v_k_minus_1_A[0] * (d_fcircle_dxy_A(1,0) + d_f_spring_y_dx_A)) )  );


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

            std::cout<<" F Newton = "<<F_Newton_B<<std::endl;

            Matrix2d Jacobian_B;
            Jacobian_B << 1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_x_dx_B + h* six_Pi_mu_r_B, -h2 * (d_f_circle_x_dy + d_f_spring_x_dy_B), -h2 * (d_f_circle_y_dx + d_f_spring_y_dx_B)  ,  1.0 + h * six_Pi_mu_r_B - h2 * d_f_spring_y_dy_B + h* six_Pi_mu_r_B;



            Vector2d old_eval_F = eval_F(delta_v_old_B,  h, six_Pi_mu_r_B, h2, d_f_spring_x_dx_B, d_f_spring_y_dx_B,  d_f_spring_x_dy_B,  d_f_spring_y_dy_B, d_f_circle_y_dx, d_f_circle_x_dy, F_B, v_k_minus_1_B);
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

           // double new_eval = (A_MAT * delta_v_new_B - b_vec).norm();

           double new_eval = F_Newton_B_new.norm();

                int while_counter = 0;
            while(new_eval > old_eval + alpha * t * old_eval_F.transpose() *(-Jacobian_B.inverse() * F_Newton_B)){
                t = beta * t;
                delta_v_new_B = delta_v_old_B -  t *Jacobian_B.inverse() * F_Newton_B;

                F_Newton_B_new[0] = delta_v_new_B[0] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0))) - delta_v_new_B[1] * h2 * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)
                                    - ( h * (  F_B[0] + h * (v_k_minus_1_B[0] * (d_f_spring_x_dx_B + d_fcircle_dxy_B(0,0)) + v_k_minus_1_B[1] * (d_fcircle_dxy_B(0,1) + d_f_spring_x_dy_B)) )  );
                F_Newton_B_new[1] = delta_v_new_B[1] * (1.0 + h * six_Pi_mu_r_B - h2 * (d_f_spring_y_dy_B+d_fcircle_dxy_B(1,1))) - delta_v_new_B[0] * h2 * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)
                                    - ( h * (  F_B[1] + h * (v_k_minus_1_B[1] * (d_f_spring_y_dy_B + d_fcircle_dxy_B(1,1)) + v_k_minus_1_B[0] * (d_fcircle_dxy_B(1,0) + d_f_spring_y_dx_B)) )  );

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
        Vector2d delta_v_B = delta_v_old_B;
        std::cout<<" delta v A 2= "<<delta_v_A<<std::endl;
        std::cout<<" delta v B 2= "<<delta_v_B<<std::endl;

        //END OF NEWTONSMETHOD

        std::cout<<" A vel pre = "<<A.velocity<<std::endl;
        std::cout<<" B vel pre = "<<B.velocity<<std::endl;
        A.velocity += delta_v_A;
        std::cout<<" A vel post = "<<A.velocity<<std::endl;

        B.velocity += delta_v_B;
        std::cout<<" B vel post = "<<B.velocity<<std::endl;
        std::cout<<" A pos pre = "<<A.position<<std::endl;
        std::cout<<" B pos pre = "<<B.position<<std::endl;
        //A.velocity *= 0.1;
       // B.velocity *= 0.1;
        A.position += time_step * A.velocity;
        B.position += time_step * B.velocity;
        std::cout<<" A pos post = "<<A.position<<std::endl;
        std::cout<<" B pos post = "<<B.position<<std::endl;



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