#include <iostream>
#include <array>
#include <vector>
#include"particles.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <chrono>
#include <random>

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::Matrix2d;
using Eigen::Matrix;

class Simulation {
public:
    int num_steps = 200;
    double time_step = Particle::time_step;
    std::vector<Particle> particles;
    double eta = 1.0 * std::pow(10,-3); //viscosity of water
    double T = 293; //room temperatur
    double kB = 1.38 * std::pow(10,-23); //Boltzmann constant




    Simulation(std::vector<Particle> particles): particles(particles)  {}

    void add_viscous_drag(Particle &particle){
        double tmp = 6.0 * M_PI * eta * particle.radius * time_step;
        Vector2d delta_v_drag = particle.velocity * (tmp/(1.0 - tmp));
        particle.velocity += delta_v_drag;
    }

    void add_brownian_motion(Particle &particle){

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution (0.0,1.0);

        double x_rand = distribution(generator);
        double y_rand = distribution(generator);
        double dx = x_rand * std::sqrt((kB * T * time_step)/(3 * M_PI * eta)) * std::sqrt(particle.radius);
        double dy = y_rand * std::sqrt((kB * T * time_step)/(3 * M_PI * eta)) * std::sqrt(particle.radius);
        particle.position[0] += dx;
        particle.position[1] += dy;
    }

    void run_simulation(std::vector<Particle> particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                particles[j].update_particle(particles[j]);
            }
        }
    }

    void run_simulation_with_drag(std::vector<Particle> particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                add_viscous_drag(particles[j]); // assumption add the new velocities together
                particles[j].update_particle(particles[j]);
            }
        }
    }

    void run_simulation_with_brownian_motion(std::vector<Particle> particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                particles[j].update_particle(particles[j]);
                add_brownian_motion(particles[j]);
            }
        }
    }

    void run_simulation_with_brownian_motion_only(std::vector<Particle> particles) {
        for (int i = 0; i < num_steps; i++) {
            for (int j = 0; j<particles.size(); j++) {
                add_brownian_motion(particles[j]);
               /* std::cout<<"xpos = "<<particles[j].position[0]<<std::endl;
                std::cout<<"ypos = "<<particles[j].position[1]<<std::endl;
                std::cout<<"xvel = "<<particles[j].velocity[0]<<std::endl;
                std::cout<<"yvel = "<<particles[j].velocity[1]<<std::endl;*/
            }
        }
    }

};