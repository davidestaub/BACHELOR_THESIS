#include<iostream>
#include<array>
#include<vector>
#include<cmath>
#include <Eigen/Core>
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
    double mass,radius;
    static double time_step;
    Particle(double mass = 1.0 * std::pow(10,-2) ,double radius = 1.0 * std::pow(10,-6)): mass(mass), radius(radius) {
        position = {1.0,1.0};
        velocity = {0.0,0.0};

    }


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

    /* Vector2d new_x(Particle particle_) {
         Vector2d new_x_;
         new_x_ = particle_.position + time_step * (particle_.velocity + delta_v(particle_));
         return new_x_;
     }*/



     void update_particle(Particle &particle) {
         std::cout<<"x-vel "<<particle.velocity[0]<<std::endl;
         std::cout<<"y-vel "<<particle.velocity[1]<<std::endl;
        std::cout<<"x-pos "<<particle.position[0]<<std::endl;
        std::cout<<"y-pos "<<particle.position[1]<<std::endl;

         particle.velocity += delta_v(particle);
         particle.position += (time_step * particle.velocity);


    }



public:
    Vector2d initial_position {1.0,1.0};
    Vector2d initial_velocity {0.0,0.0};
    double default_mass = 1.0 * std::pow(10,-4);
    double default_radius = 1.0 * std::pow(10,-5);

};

double Particle::time_step =0.01;