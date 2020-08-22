#include<iostream>
#include<array>
#include<vector>
#include<cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#pragma once

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::Matrix2d;
using Eigen::Matrix;


class Particle{
public:
    Vector2d velocity, position;
    Vector2d tmp_drag;
    Vector2d tmp_noise;
    Vector2d tmp_update;
    double mass,radius;
    double six_Pi_mu_r;
    int index = -1;

    static double time_step;
    static double radius_for_spring;
    static double eta;
    bool is_drawn = false;
    bool visited  = false;
    bool visited_ekin = false;
    //everything in micro
    Particle(double mass = 10 ,double radius = 10, Vector2d position = {200,200}): mass(mass), radius(radius), position(position) {
        position = initial_position;//{350.0,300.0};
        velocity = initial_velocity;
        tmp_drag = {0,0};
        tmp_noise = {0,0};
        tmp_update = {0,0};
        six_Pi_mu_r = ( 6.0 * M_PI * eta * (radius*std::pow(10,-6))  * std::pow(10,6) );
        //changed initial radius to much bigger, previously 10^-6

    }
    Vector2d initial_position = {200,200};
    Vector2d initial_velocity = {0.0,0.0};


    Vector2d force(Particle &particle) {
        Vector2d force_;
        force_[0] = std::pow(particle.position[0], 2) + particle.velocity[0];
        force_[1] = 2.0 * particle.position[1] + particle.velocity[1];
        return force_;
    }


    Matrix2d dforce_dx(Particle &particle){
        Matrix2d dforce_dx_;
        double tmp = particle.position[0] * 2.0;
        dforce_dx_<<tmp,0,0,2.0;
        return dforce_dx_;
    }

    Matrix2d dforce_dv(Particle &particle){
        Matrix2d dforce_dv_;
        dforce_dv_<<1.0,0,0,1.0;
        return dforce_dv_;
    }


    Vector2d delta_v(Particle &particle) {
        Eigen::SparseMatrix<double> A(2, 2);
        Eigen::SparseMatrix<double> Identity(2, 2);
        Identity.setIdentity();
        A = Identity - time_step * dforce_dv(particle) - std::pow(time_step, 2) * dforce_dx(particle);
        Vector2d b;
        b = time_step * (force(particle) +time_step * dforce_dx(particle) * particle.velocity);
        Vector2d delta_v_;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
        solver.compute(A);
        delta_v_ = solver.solve(b);
        return delta_v_;
    }

    Vector2d push_force(Particle &particle){
        Vector2d force_;
        force_ = -particle.position; //push inwards
        return force_;
    }

    Matrix2d dforce_push_dx(Particle &particle){
        Matrix2d dforce_push_dx_;
        dforce_push_dx_<<-1.0,0,0,-1.0;
        return dforce_push_dx_;
    }
    Matrix2d dforce_push_dv(Particle &particle){
        Matrix2d dforce_push_dv_;
        dforce_push_dv_<<0,0,0,0;
        return dforce_push_dv_;
    }

    Vector2d delta_v_push(Particle &particle) {
        Eigen::SparseMatrix<double> A(2, 2);
        Eigen::SparseMatrix<double> Identity(2, 2);
        Identity.setIdentity();
        A = Identity - time_step * dforce_push_dv(particle) - std::pow(time_step, 2) * dforce_push_dx(particle);
        Vector2d b;
        b = time_step * (push_force(particle) +time_step * dforce_push_dx(particle) * particle.velocity);
        Vector2d delta_v_;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
        solver.compute(A);
        delta_v_ = solver.solve(b);
        return delta_v_;
    }




    //circle force field

    Vector2d circle_force2(Particle &particle){
        Vector2d force_;
        force_[0] = -particle.position[1];
        force_[1] = particle.position[0];
        Vector2d normforce = force_.normalized() * 100.0;
       // std::cout<<"FORCE = "<<normforce<<std::endl;
        return normforce;
    }

    Matrix2d circle_dforce_dx2(Particle &particle){
        Matrix2d dforce_circle_dx_;
        Matrix2d Identity;
        Identity <<1,0,0,1;
        double x = -particle.position[1];
        double y = particle.position[0];
        double x_dx = x*y/(std::pow( std::pow(x,2) + std::pow(y,2),(3.0/2.0)    )   );
        double x_dy = -1.0 * std::pow(x,2) / (std::pow( std::pow(x,2) + std::pow(y,2),(3.0/2.0)    )   );
        double y_dy = std::pow(y,2) / (std::pow( std::pow(x,2) + std::pow(y,2),(3.0/2.0)    )   );
        double y_dx = -1.0 * x * y /(std::pow( std::pow(x,2) + std::pow(y,2),(3.0/2.0)    )   );

       Vector2d force = particle.circle_force2(particle);
       double t0 = force.norm();
       dforce_circle_dx_ = 1.0/t0 * Identity - (1.0/std::pow(t0,3)) * (force * force.transpose());


        //dforce_circle_dx_<<x_dx,x_dy,y_dx,y_dy;
        //std::cout<<" D FORCE Dx = "<<dforce_circle_dx_<<std::endl;
        return dforce_circle_dx_ * 100.0;
    }
    Matrix2d circle_dforce_dv2(Particle &particle){
        Matrix2d dforce_circle_dv_;
        dforce_circle_dv_<<0,0,0,0;
        return dforce_circle_dv_;
    }

    Vector2d delta_v_circle2(Particle &particle) {
        Eigen::Matrix2d A;
        Eigen::Matrix2d Identity;
        Identity<<1,0,0,1;
        A = Identity - time_step * circle_dforce_dv2(particle) - std::pow(time_step, 2) * circle_dforce_dx2(particle);
        Vector2d b;
        b = time_step * (circle_force2(particle) +time_step * circle_dforce_dx2(particle) * particle.velocity);
        Vector2d delta_v_;

        delta_v_ = A.colPivHouseholderQr().solve(b);
       //std::cout<<"delta v = "<<delta_v_<<std::endl;
        return delta_v_;

    }





    Vector2d circle_force(Particle &particle){
        Vector2d force_;
        force_[0] = -particle.position[1];
        force_[1] = particle.position[0];

        return force_;
    }

    Matrix2d circle_dforce_dx(Particle &particle){
        Matrix2d dforce_circle_dx_;




        dforce_circle_dx_<<0,-1.0,1.0,0;

        return dforce_circle_dx_;
    }
    Matrix2d circle_dforce_dv(Particle &particle){
        Matrix2d dforce_circle_dv_;
        dforce_circle_dv_<<0,0,0,0;
        return dforce_circle_dv_;
    }

    Vector2d delta_v_circle(Particle &particle) {
        Eigen::Matrix2d A;
        Eigen::Matrix2d Identity;
        Identity<<1,0,0,1;
        A = Identity - time_step * circle_dforce_dv(particle) - std::pow(time_step, 2) * circle_dforce_dx(particle);
        Vector2d b;
        b = time_step * (circle_force(particle) +time_step * circle_dforce_dx(particle) * particle.velocity);
        Vector2d delta_v_;

        delta_v_ = A.colPivHouseholderQr().solve(b);
        //std::cout<<"delta v = "<<delta_v_<<std::endl;
        return delta_v_;

    }

    //end circle force field

    /* Vector2d new_x(Particle particle_) {
         Vector2d new_x_;
         new_x_ = particle_.position + time_step * (particle_.velocity + delta_v(particle_));
         return new_x_;
     }*/






};

double Particle::time_step =0.001; // seconds
double Particle::radius_for_spring = 10;
double Particle::eta = 1.0 * std::pow(10,-3);
