#include "application.h"
#include <imgui.h>
#include "simulation.h"

#include <iostream>
#include <math.h>
#include <deque>
#include <chrono>
#include <vector>



#include <fstream>
#include <string>

#include <Eigen/Core>
using Eigen::Vector2d;
using Eigen::Vector2f;

std::ofstream outputFile;
std::ofstream fs;

std::string filename = "simulation_data.csv";

class TestApp : public Application
{
#define COLOR_OUT    nvgRGBA(220,50,50,255)
#define COLOR_IN     nvgRGBA(50,50,220,255)
#define COLOR_SOLVED nvgRGBA(50,220,50,255)



public:
    TestApp(int w, int h, const char * title) : Application(title, w, h) {
        lastFrame = std::chrono::high_resolution_clock::now();

        //flat triangle from laura
       // simulation.connect(A,B,(A.radius+B.radius),connected_particles);
       // simulation.connect(B,C,(B.radius + C.radius), connected_particles);
       // simulation.connect(A,C,(A.radius + 2.0 * B.radius + C.radius), connected_particles);


        //frist triangle from laura
       // simulation.connect(A_t,B_t,(A_t.radius+B_t.radius),connected_particles);
       // simulation.connect(B_t,C_t,(B_t.radius + C_t.radius), connected_particles);
       // simulation.connect(A_t,C_t,(std::sqrt(std::pow(B_t.radius + A_t.radius,2) + std::pow(B_t.radius+C_t.radius,2))), connected_particles);

        //dumbel
        simulation.connect(A,B,(A.radius+B.radius),connected_particles);
        simulation.connect(A_2,B_2,(A_2.radius+B_2.radius),connected_particles);
        simulation.connect(A_3,B_3,(A_3.radius+B_3.radius),connected_particles);
        /*simulation.connect(A_4,B_4,(A_4.radius+B_4.radius),connected_particles);
        simulation.connect(A_5,B_5,(A_5.radius+B_5.radius),connected_particles);
        simulation.connect(A_6,B_6,(A_6.radius+B_6.radius),connected_particles);
        simulation.connect(A_7,B_7,(A_7.radius+B_7.radius),connected_particles);
        simulation.connect(A_8,B_8,(A_8.radius+B_8.radius),connected_particles);
        simulation.connect(A_9,B_9,(A_9.radius+B_9.radius),connected_particles);
        simulation.connect(A_10,B_10,(A_10.radius+B_10.radius),connected_particles);*/

        //rotator
       // simulation.connect(A,B,(A.radius+B.radius),connected_particles);
      //  simulation.connect(B,C,(B.radius+C.radius),connected_particles);
      //  simulation.connect(C,D,(C.radius+D.radius),connected_particles);
      //  simulation.connect(A,C,(std::sqrt(std::pow((A.radius + B.radius),2) + std::pow((B.radius + C.radius),2))),connected_particles);
      //  simulation.connect(B,D,(std::sqrt(std::pow((B.radius + C.radius),2) + std::pow((C.radius + D.radius),2))),connected_particles);

        //spiral triangle that does not do spirals
        // simulation.connect(A,B,(A.radius+B.radius),connected_particles);
       //    double fix_rest = std::sqrt(std::pow(mscale *(100.0 - (104.0 + 5.0/std::sqrt(2))),2) + std::pow(mscale *(100.0 - (100.0 + 5.0/std::sqrt(2))),2));
        //   simulation.connect(B,C,(fix_rest), connected_particles);
        //  simulation.connect(A,C,A.radius + C.radius, connected_particles);


        //spiral triangle
       // simulation.connect(A,B,(A.radius+B.radius),connected_particles);
       // double fix_rest_2 = std::sqrt(std::pow(mscale *(100.0 - (105.0 + 4.0/std::sqrt(2))),2) + std::pow(mscale *(100.0 - (100.0 + 4.0/std::sqrt(2))),2));
        //simulation.connect(B,C,(fix_rest_2), connected_particles);
        //simulation.connect(A,C,A.radius + C.radius, connected_particles);



        // L shaped particle system
       //simulation.connect(A,B,(A.radius+B.radius),connected_particles);
       // simulation.connect(B,C,(B.radius + 2.0 *A.radius + C.radius),connected_particles);
       // simulation.connect(A,C,(A.radius+C.radius),connected_particles);
       //simulation.connect(C,D,(C.radius+D.radius),connected_particles);
       // simulation.connect(Up,B,(Up.radius+B.radius),connected_particles);
       // simulation.connect(D,B,(2.0 *A.radius+B.radius + 2.0 * C.radius + D.radius),connected_particles);
       // simulation.connect(A,D,(A.radius+ 2.0 *C.radius + D.radius),connected_particles);
       // simulation.connect(Up,A,(std::sqrt(std::pow(B.radius + A.radius,2) + std::pow(B.radius+Up.radius,2))), connected_particles);



        rect = Box{{simulation.Boxcenter[0], simulation.Boxcenter[1]}, Vector2d(simulation.Boxsize[0], simulation.Boxsize[1]), nvgRGBA(50, 50, 50, 0), nvgRGBA(50, 50, 50, 200)};

        simulation.assign_index(connected_particles);
        simulation.initialize_position_vector_minus_1(connected_particles);
        simulation.count_neighbours(connected_particles);

    }





    void process() override {



        if(runSim) {
           //flat triangle from laura
          //  simulation.connect_new(A,B,(A.radius+B.radius),connected_particles,0);
          //  simulation.connect_new(B,C,(B.radius + C.radius), connected_particles,1);
          //  simulation.connect_new(A,C,(A.radius + 2.0 * B.radius + C.radius), connected_particles,2);


            //first triangle from laura
          //  simulation.connect_new(A_t,B_t,(A_t.radius+B_t.radius),connected_particles,0);
          //  simulation.connect_new(B_t,C_t,(B_t.radius + C_t.radius), connected_particles,1);
         //   simulation.connect_new(A_t,C_t,(std::sqrt(std::pow(B_t.radius + A_t.radius,2) + std::pow(B_t.radius+C_t.radius,2))), connected_particles,2);



            //dumbel
            simulation.connect_new(A,B,(A.radius + B.radius), connected_particles,0);
            simulation.connect_new(A_2,B_2,(A_2.radius+B_2.radius),connected_particles,1);
            simulation.connect_new(A_3,B_3,(A_3.radius+B_3.radius),connected_particles,2);
           /* simulation.connect_new(A_4,B_4,(A_4.radius+B_4.radius),connected_particles,3);
            simulation.connect_new(A_5,B_5,(A_5.radius+B_5.radius),connected_particles,4);
            simulation.connect_new(A_6,B_6,(A_6.radius+B_6.radius),connected_particles,5);
            simulation.connect_new(A_7,B_7,(A_7.radius+B_7.radius),connected_particles,6);
            simulation.connect_new(A_8,B_8,(A_8.radius+B_8.radius),connected_particles,7);
            simulation.connect_new(A_9,B_9,(A_9.radius+B_9.radius),connected_particles,8);
            simulation.connect_new(A_10,B_10,(A_10.radius+B_10.radius),connected_particles,9);*/


            //rotator
         //   simulation.connect_new(A,B,(A.radius+B.radius),connected_particles,0);
         //   simulation.connect_new(B,C,(B.radius+C.radius),connected_particles,1);
          //  simulation.connect_new(C,D,(C.radius+D.radius),connected_particles,2);
         //   simulation.connect_new(A,C,(std::sqrt(std::pow((A.radius + B.radius),2) + std::pow((B.radius + C.radius),2))),connected_particles,3);
         //   simulation.connect_new(B,D,(std::sqrt(std::pow((B.radius + C.radius),2) + std::pow((C.radius + D.radius),2))),connected_particles,4);



            //spiral triangle that does not do spirals
          //  simulation.connect_new(A,B,(A.radius+B.radius),connected_particles,0);
          //  double fix_rest = std::sqrt(std::pow(mscale *(100.0 - (104.0 + 5.0/std::sqrt(2))),2) + std::pow(mscale *(100.0 - (100.0 + 5.0/std::sqrt(2))),2));
          //  simulation.connect_new(B,C,(fix_rest), connected_particles,1);
          //  simulation.connect_new(A,C,A.radius + C.radius, connected_particles,2);

            //spiral triangle
           // simulation.connect_new(A,B,(A.radius+B.radius),connected_particles,0);
            //double fix_rest_2 = std::sqrt(std::pow(mscale *(100.0 - (105.0 + 4.0/std::sqrt(2))),2) + std::pow(mscale *(100.0 - (100.0 + 4.0/std::sqrt(2))),2));
            //simulation.connect_new(B,C,(fix_rest_2), connected_particles,1);
            //simulation.connect_new(A,C,A.radius + C.radius, connected_particles,2);


            // L shaped particle system
         // simulation.connect_new(A,B,(A.radius+B.radius),connected_particles,0);
         /*   simulation.connect_new(B,C,(B.radius + 2.0 *A.radius + C.radius),connected_particles,1);
            simulation.connect_new(A,C,(A.radius+C.radius),connected_particles,2);

            simulation.connect_new(C,D,(C.radius+D.radius),connected_particles,3);
            simulation.connect_new(Up,B,(Up.radius+B.radius),connected_particles,4);
            simulation.connect_new(D,B,(2.0 *A.radius+B.radius + 2.0 * C.radius + D.radius),connected_particles,5);
            simulation.connect_new(A,D,(A.radius+ 2.0 *C.radius + D.radius),connected_particles,6);
            simulation.connect_new(Up,A,(std::sqrt(std::pow(B.radius + A.radius,2) + std::pow(B.radius+Up.radius,2))), connected_particles,7);
*/


    // move image if right mouse button is pressed
        if(mouseState.rButtonPressed){
            auto dw = (int)(mouseState.lastMouseX - cursorPosDown[0]);
            auto dh = (int)(mouseState.lastMouseY - cursorPosDown[1]);
            translation[0] += dw/(double)base;
            translation[1] -= dh/(double)base;
            cursorPosDown[0] = mouseState.lastMouseX;
            cursorPosDown[1] = mouseState.lastMouseY;
        }



       //simulation.run_simulation_with_forcefield_only(particles); //no changes made to vel and pos, works now with damping and smaller timestep problme was to fast convergence to zero due to huge time step
       //simulation.run_simulation_with_drag_only(particles); //pos increases drastically, better now with smaller timestep
        //simulation.run_simulation_with_brownian_motion_only(particles); // works maybe a bit too much movement
        //simulation.run_simulation_with_drag(particles);

          //  simulation.run_simulation(particles);
        //simulation.run_simulation_with_spring_connected_Particles(A,B,connect_vector);
           // simulation.run_simulation_with_spring_connected_Particles_averaged_movement(A,B,connect_vector);


           simulation.run_simulation_for_connected_Particles(connected_particles);

            simulation.time_index++;

        }

    }

    void drawImGui() override {

        using namespace ImGui;

       /* Begin("Hello World!");
        TextWrapped("Use the arrow keys to move the first circle.");
        TextWrapped("Drag the other circle with the mouse.");
        TextWrapped("Put them both in the Box!");
        End();*/

        BeginMainMenuBar();
        if(BeginMenu("Debug")){
            Checkbox("draw cursor", &drawCursor);
            Checkbox("draw circles", &drawCircles);
            Checkbox(" Switch view", &switch_view);
            Checkbox("Safety off",&simulation.safety);
            Checkbox("Dimers Only",&simulation.dimers_only);
            SliderScalar("Move Camera along x-axis",ImGuiDataType_Double, &center_of_frame[0], &min_cof_x,&max_cof_x);
            SliderScalar("Move Camera along y-axis",ImGuiDataType_Double, &center_of_frame[1], &min_cof_y,&max_cof_y);
            SliderScalar("Move Camera along z-axis",ImGuiDataType_Double, &center_of_frame[2], &min_cof_z,&max_cof_z);
            SliderScalar("Zoom", ImGuiDataType_Double, &mscale, &min_my_zoom, &max_my_zoom);
            if(Checkbox("increase frequency by 100",&hundred_step)){
                simulation.lower_electrode.frequency +=100;
            }
            if(Checkbox("increase frequency by 500",&five_hundred_step)){
                simulation.lower_electrode.frequency +=500;
            }
            if(Checkbox("increase frequency by 1000",&thousand_step)){
                simulation.lower_electrode.frequency +=1000;
            }
            if(Checkbox("increase frequency by 10000",&ten_thousand_step)){
                simulation.lower_electrode.frequency +=10000;
            }
            if(Checkbox("increase frequency by 100000",&hundred_thousand_step)){
                simulation.lower_electrode.frequency +=100000;
            }
            if(Checkbox("Tracer on", &simulation.tracer_on)){
                simulation.dimer_tracers = false;
            }
            if(Checkbox("Dimer Tracers on", &simulation.dimer_tracers)){
                simulation.tracer_on = false;
            }
            if(Checkbox("Erase Dimer Tracers", &simulation.erase_dimer_tracers)){
                simulation.dimer_tracers = false;
            }
            if(Checkbox("Follow cam",&follow_cam)) {
                Vector3d com = Eigen::VectorXd::Zero(3);
                double count = 0.0;
                for (auto &particle_pair : connected_particles) {
                    if (!std::get<0>(particle_pair).visited) {
                        com += (1.0 / mscale) * std::get<0>(particle_pair).position;
                        std::get<0>(particle_pair).visited = true;
                        count = count + 1.0;
                    }
                    if (!std::get<1>(particle_pair).visited) {
                        com += (1.0 / mscale) * std::get<1>(particle_pair).position;
                        std::get<1>(particle_pair).visited = true;
                        count = count + 1.0;
                    }
                }
                if (count != 0.0) {
                    com = com / count;
                    center_of_frame = com;
                }
                simulation.reset_flags(connected_particles);
            }
            ImGui::EndMenu();
        }
        EndMainMenuBar();

        BeginMainMenuBar();
        if(BeginMenu("Particles")){
            int i = 0;
            for(auto& particle_pair: connected_particles){
                Particle& A = std::get<0>(particle_pair);
                Particle& B = std::get<1>(particle_pair);
                if(!A.visited){
                    std::string label = "Particle ";
                    label+= std::to_string(i);
                    const char* label_ = label.c_str();
                    if(CollapsingHeader(label_)){
                        std::string mass_label = label;
                        std::string mass = " mass";
                        mass_label+= mass;
                        const char* mass_label_ = mass_label.c_str();
                        if(SliderScalar(mass_label_,ImGuiDataType_Double,&A.mass,&min_mass,&max_mass)){

                        }

                        std::string charge_label = label;
                        std::string charge = " charge";
                        charge_label+= charge;
                        const char* charge_label_ = charge_label.c_str();
                        if(SliderScalar(charge_label_,ImGuiDataType_Double,&A.charge,&min_charge,&max_charge)){

                        }
                        std::string radius_label = label;
                        std::string radius = " radius";
                        radius_label += radius;
                        const char* radius_label_ = radius_label.c_str();
                        if(SliderScalar(radius_label_,ImGuiDataType_Double,&A.radius,&min_radius,&max_radius)){

                        }
                        std::string set_radius2_label = label;
                        std::string set_radius2 = " set radius to 2 micro meters";
                        set_radius2_label += set_radius2;
                        const char* set_radius2_label_ = set_radius2_label.c_str();
                        if(Checkbox(set_radius2_label_, &A.two_m)){
                            A.radius = 2.0 * std::pow(10,-6);
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                        }
                        std::string set_radius3_label = label;
                        std::string set_radius3 = " set radius to 3 micro meters";
                        set_radius3_label += set_radius3;
                        const char* set_radius3_label_ = set_radius3_label.c_str();
                        if(Checkbox(set_radius3_label_, &A.three_m)){
                            A.radius = 3.0 * std::pow(10,-6);
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                        }

                        std::string set_radius0_7_label = label;
                        std::string set_radius0_7 = " set radius to 0.7 micro meters";
                        set_radius0_7_label += set_radius0_7;
                        const char* set_radius0_7_label_ = set_radius0_7_label.c_str();
                        if(Checkbox(set_radius0_7_label_, &A.zero_seven_m)){
                            A.radius = 0.7 * std::pow(10,-6);
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                        }

                        std::string set_radius1_label = label;
                        std::string set_radius1 = " set radius to 1.0 micro meters";
                        set_radius1_label += set_radius1;
                        const char* set_radius1_label_ = set_radius1_label.c_str();
                        if(Checkbox(set_radius1_label_, &A.one_m)){
                            A.radius = 1.0 * std::pow(10,-6);
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                        }



                        std::string permittivity_label = label;
                        std::string permittivity = " permittivity";
                        permittivity_label += permittivity;
                        const char* permittivity_label_ = permittivity_label.c_str();
                        if(SliderScalar(permittivity_label_,ImGuiDataType_Double,&A.permittivity,&min_permittivity,&max_permittivity)){

                        }
                        std::string conductivity_label = label;
                        std::string conductivity = " stern layer conductance";
                        conductivity_label += conductivity;
                        const char* conductivity_label_ = conductivity_label.c_str();
                        if(SliderScalar(conductivity_label_,ImGuiDataType_Double,&A.stern_layer_conductance,&min_stern_layer_conductance,&max_stern_layer_conductance)){

                        }

                        std::string zeta_label = label;
                        std::string zeta = " zeta potential";
                        zeta_label += zeta;
                        const char* zeta_label_ = zeta_label.c_str();
                        if(SliderScalar(zeta_label_,ImGuiDataType_Double,&A.zeta_potential,&min_zeta,&max_zeta)){

                        }

                        std::string dissolve_label = label;
                        std::string dissolve = " dissolve";
                        dissolve_label+=dissolve;
                        const char* dissolve_label_ = dissolve_label.c_str();
                        Checkbox(dissolve_label_,&A.dissolve);

                        std::string silica_label = label;
                        std::string silica = " Silica";
                        silica_label+=silica;
                        const char* silica_label_ = silica_label.c_str();
                        if(Checkbox(silica_label_,&A.is_silica)){
                            A.permittivity = 3.9;
                            A.density = 2650.0;
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3); //new
                            A.stern_layer_conductance = 0.05 * std::pow(10,-9);
                            A.zeta_potential = -1.0 * 43 * std::pow(10,-3);
                        }
                        std::string PS_label = label;
                        std::string PS = " PS";
                        PS_label+=PS;
                        const char* PS_label_ = PS_label.c_str();
                        if(Checkbox(PS_label_, &A.is_PS)){
                            A.permittivity = 2.5;
                            A.density = 1050;
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                            A.stern_layer_conductance = 5 * std::pow(10,-9);
                           A.zeta_potential = -56.0 * std::pow(10,-3); //check if mv or V

                        }

                        std::string PNiPAM_microgel_label = label;
                        std::string PNiPAM_microgel = " PNiPAM microgel";
                        PNiPAM_microgel_label+=PNiPAM_microgel;
                        const char* PNiPAM_microgel_label_ = PNiPAM_microgel_label.c_str();
                        if(Checkbox(PNiPAM_microgel_label_, &A.is_PNiPAM_microgel)){
                            A.permittivity = 2000;
                            A.density = 1100;
                            A.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                            A.stern_layer_conductance = 5 * std::pow(10,-9);
                            A.zeta_potential = -10.0 * std::pow(10,-3); //check if mv or V

                        }


                    }
                    A.visited=true;
                    i++;
                }
                if(!B.visited){
                    std::string label = "Particle ";
                    label+= std::to_string(i);
                    const char* label_ = label.c_str();
                    if(CollapsingHeader(label_)){
                        std::string mass_label = label;
                        std::string mass = " mass";
                        mass_label+= mass;
                        const char* mass_label_ = mass_label.c_str();
                        if(SliderScalar(mass_label_,ImGuiDataType_Double,&B.mass,&min_mass,&max_mass)){

                        }
                        std::string charge_label = label;
                        std::string charge = " charge";
                        charge_label+= charge;
                        const char* charge_label_ = charge_label.c_str();
                        if(SliderScalar(charge_label_,ImGuiDataType_Double,&B.charge,&min_charge,&max_charge)){

                        }
                        std::string radius_label = label;
                        std::string radius = " radius";
                        radius_label += radius;
                        const char* radius_label_ = radius_label.c_str();
                        if(SliderScalar(radius_label_,ImGuiDataType_Double,&B.radius,&min_radius,&max_radius)){

                        }

                        std::string set_radius2_label = label;
                        std::string set_radius2 = " set radius to 2 micro meters";
                        set_radius2_label += set_radius2;
                        const char* set_radius2_label_ = set_radius2_label.c_str();
                        if(Checkbox(set_radius2_label_, &B.two_m)){
                            B.radius = 2.0 * std::pow(10,-6);
                            B.mass = B.density * (4.0/3.0) * M_PI * std::pow(B.radius,3);
                        }
                        std::string set_radius3_label = label;
                        std::string set_radius3 = " set radius to 3 micro meters";
                        set_radius3_label += set_radius3;
                        const char* set_radius3_label_ = set_radius3_label.c_str();
                        if(Checkbox(set_radius3_label_, &B.three_m)){
                            B.radius = 3.0 * std::pow(10,-6);
                            B.mass = B.density * (4.0/3.0) * M_PI * std::pow(B.radius,3);
                        }

                        std::string set_radius0_7_label = label;
                        std::string set_radius0_7 = " set radius to 0.7 micro meters";
                        set_radius0_7_label += set_radius0_7;
                        const char* set_radius0_7_label_ = set_radius0_7_label.c_str();
                        if(Checkbox(set_radius0_7_label_, &B.zero_seven_m)){
                            B.radius = 0.7 * std::pow(10,-6);
                            B.mass = B.density * (4.0/3.0) * M_PI * std::pow(B.radius,3);
                        }

                        std::string set_radius1_label = label;
                        std::string set_radius1 = " set radius to 1.0 micro meters";
                        set_radius1_label += set_radius1;
                        const char* set_radius1_label_ = set_radius1_label.c_str();
                        if(Checkbox(set_radius1_label_, &B.one_m)){
                            B.radius = 1.0 * std::pow(10,-6);
                            B.mass = B.density * (4.0/3.0) * M_PI * std::pow(B.radius,3);
                        }


                        std::string permittivity_label = label;
                        std::string permittivity = " permittivity";
                        permittivity_label += permittivity;
                        const char* permittivity_label_ = permittivity_label.c_str();
                        if(SliderScalar(permittivity_label_,ImGuiDataType_Double,&B.permittivity,&min_permittivity,&max_permittivity)){

                        }
                        std::string conductivity_label = label;
                        std::string conductivity = " stern layer conductance";
                        conductivity_label += conductivity;
                        const char* conductivity_label_ = conductivity_label.c_str();
                        if(SliderScalar(conductivity_label_,ImGuiDataType_Double,&B.stern_layer_conductance,&min_stern_layer_conductance,&max_stern_layer_conductance)){

                        }

                        std::string zeta_label = label;
                        std::string zeta = " zeta potential";
                        zeta_label += zeta;
                        const char* zeta_label_ = zeta_label.c_str();
                        if(SliderScalar(zeta_label_,ImGuiDataType_Double,&B.zeta_potential,&min_zeta,&max_zeta)){

                        }


                        std::string dissolve_label = label;
                        std::string dissolve = " dissolve";
                        dissolve_label+=dissolve;
                        const char* dissolve_label_ = dissolve_label.c_str();
                        Checkbox(dissolve_label_,&B.dissolve);

                        std::string silica_label = label;
                        std::string silica = " Silica";
                        silica_label+=silica;
                        const char* silica_label_ = silica_label.c_str();
                        if(Checkbox(silica_label_,&B.is_silica)){
                            B.permittivity = 3.9;
                            B.density = 2650.0;

                            B.mass = B.density * (4.0/3.0) * M_PI * std::pow(B.radius,3);
                            B.stern_layer_conductance = 0.05 * std::pow(10,-9);
                            A.zeta_potential = -1.0 * 43 * std::pow(10,-3);
                        }
                        std::string PS_label = label;
                        std::string PS = " PS";
                        PS_label+=PS;
                        const char* PS_label_ = PS_label.c_str();


                        if(Checkbox(PS_label_, &B.is_PS)){
                            B.permittivity = 2.5;
                            B.density = 1050;
                           // B.conductivity = 5 * std::pow(10,-9);
                           //B.conductivity = 0.001;
                     //      B.conductivity = 2.0 *  ((5 * std::pow(10,-9)) / B.radius);
                            B.mass = B.density * (4.0/3.0) * M_PI * std::pow(B.radius,3);
                            B.zeta_potential = -56.0 * std::pow(10,-3); //check if mv or V
                            B.stern_layer_conductance = 5 * std::pow(10,-9);
                        }

                        std::string PNiPAM_microgel_label = label;
                        std::string PNiPAM_microgel = " PNiPAM microgel";
                        PNiPAM_microgel_label+=PNiPAM_microgel;
                        const char* PNiPAM_microgel_label_ = PNiPAM_microgel_label.c_str();
                        if(Checkbox(PNiPAM_microgel_label_, &B.is_PNiPAM_microgel)){
                            B.permittivity = 2000;
                            B.density = 1100;
                            B.mass = A.density * (4.0/3.0) * M_PI * std::pow(A.radius,3);
                            B.stern_layer_conductance = 5 * std::pow(10,-9);
                            B.zeta_potential = -10.0 * std::pow(10,-3); //check if mv or V

                        }

                    }


                    B.visited=true;
                    i++;
                }

            }
            simulation.reset_flags(connected_particles);
            ImGui::EndMenu();
        }
        EndMainMenuBar();

       /* while(std::abs(mouseState.lastMouseX - std::get<0>(connected_particles[0]).position[0]  ) < 0.01 && std::abs(mouseState.lastMouseY - std::get<0>(connected_particles[0]).position[1]  ) < 0.01){
            std::string label = "Particle ";
            int i = std::get<0>(connected_particles[0]).index/3;
            label+= std::to_string(i);
            const char* label_ = label.c_str();
            Begin(label_);
            End();
        }*/

        BeginMainMenuBar();
        if(BeginMenu("Control")) {
            if (CollapsingHeader("Springs")) {
                if (SliderScalar("Stiffnes Constant", ImGuiDataType_Double, &simulation.stiffnes_constant,
                                 &min_stiffnes_constant, &max_stiffnes_constant))
                    std::cout << simulation.stiffnes_constant;

                if (SliderScalar("Damping Constant", ImGuiDataType_Double, &simulation.damping_constant,
                                 &min_damping_constant, &max_damping_constant)){}
            }

            if (CollapsingHeader("Drag")) {
                if (SliderScalar("Viscosity Multiplier", ImGuiDataType_Double, &simulation.viscosity_multiplier,
                                 &min_viscosity_multiplier, &max_viscosity_multiplier))
                    std::cout << simulation.viscosity_multiplier;
            }
            if (CollapsingHeader("EHD")) {
                if (SliderScalar("Permittivity of medium", ImGuiDataType_Double, &simulation.relative_permittivity,&min_epsilon_m, &max_epsilon_m)){}
                if (SliderScalar("Conductivity of medium", ImGuiDataType_Double, &simulation.conductivity,&min_sigma_m, &max_sigma_m)){}
            }

            if(CollapsingHeader("Experiments")){
                if(Checkbox(" Silica - Silica Experiment", &first_experiment)){
                    simulation.beta_ehd = 0.02;
                    simulation.lower_electrode.peak_voltage = 6.5;
                    simulation.upper_electrode.position[2] = simulation.lower_electrode.position[2] + 100.0 * mscale;
                }
                if(Checkbox(" SETUP PS - PS Experiment", &second_experiment)){
                    simulation.beta_ehd = 0.08;
                    //simulation.lower_electrode.frequency = 500.0;
                    simulation.lower_electrode.peak_voltage = 6.5;
                    simulation.lower_electrode.frequency = 300.0;
                    simulation.upper_electrode.position[2] = simulation.lower_electrode.position[2] + 100.0 * mscale;

                }
                if(Checkbox(" SETUP Microgel - PS Experiment", &second_experiment)){
                    simulation.beta_ehd = 0.1478;
                    //simulation.lower_electrode.frequency = 500.0;
                    simulation.lower_electrode.peak_voltage = 5.0;
                    simulation.lower_electrode.frequency = 300.0;
                    simulation.upper_electrode.position[2] = simulation.lower_electrode.position[2] + 100.0 * mscale;

                }
                /*if(Checkbox(" RUN FREQUENCY EXPERIMENT", &run_w_exp)){
                    std::cout<<" inside"<<std::endl;
                    if (simulation.safety){
                        std::cout<<" running"<<std::endl;
                        for (int w = 300; w < 11000; w += 5) {
                            simulation.lower_electrode.frequency = w;
                        }
                    }
                }*/
            }

            if (CollapsingHeader("Brownian Motion")) {
                if (SliderScalar("Brownian Multiplier", ImGuiDataType_Double, &simulation.brownian_motion_multiplier,
                                 &min_brownian_motion_multiplier, &max_brownian_motion_multiplier))
                    std::cout << simulation.brownian_motion_multiplier;
                if (SliderScalar("Temperatur [K]", ImGuiDataType_Double, &simulation.T,
                                 &min_T, &max_T)){}
            }
            if (CollapsingHeader("Simulation")) {
                if (SliderScalar("Time Step", ImGuiDataType_Double, &simulation.time_step, &min_time_step,
                                 &max_time_step)) {
                    std::cout << simulation.time_step;
                }
                if(Checkbox("Set Time Step to 1ms",&set_t_to_1ms)){
                    simulation.time_step = 1.0 * std::pow(10.0,-3.0);
                }
                if(Checkbox("Set Time Step to 2ms",&set_t_to_2ms)){
                    simulation.time_step = 2.0 * std::pow(10.0,-3.0);
                }
                if(Checkbox("Set Time Step to 5ms",&set_t_to_5ms)){
                    simulation.time_step = 5.0 * std::pow(10.0,-3.0);
                }
                if(Checkbox("Set Time Step to 10ms",&set_t_to_10ms)){
                    simulation.time_step = 10.0 * std::pow(10.0,-3.0);
                }
                if(Checkbox("Set Time Step to 20ms",&set_t_to_20ms)){
                    simulation.time_step = 20.0 * std::pow(10.0,-3.0);
                }
                if (SliderScalar("Max itteration", ImGuiDataType_S32, &simulation.max_iterations, &min_max_iterations,
                                 &max_max_iterations)) {
                    std::cout << simulation.max_iterations;
                }
            }
            if(CollapsingHeader("Electrodes")){
                if(CollapsingHeader("Lower Electrode")){
                    if(SliderScalar(" X position Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.position[0], &min_lower_electrode_x_position, &max_lower_electrode_x_position)){

                    }
                    if(SliderScalar(" Y position Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.position[1], &min_lower_electrode_y_position, &max_lower_electrode_y_position)){

                    }
                    if(SliderScalar(" Z position Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.position[2], &min_lower_electrode_y_position, &max_lower_electrode_y_position)){

                    }
                    if(SliderScalar(" Charge Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.charge, &min_lower_electrode_charge, &max_lower_electrode_charge)){

                    }
                    if(SliderScalar(" Voltage Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.voltage, &min_lower_electrode_voltage, &max_lower_electrode_voltage)){

                    }
                    if(SliderScalar(" Frequency Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.frequency, &min_lower_electrode_frequency, &max_lower_electrode_frequency)){

                    }

                    if(SliderScalar(" Peak Voltage Lower Electrode",ImGuiDataType_Double, &simulation.lower_electrode.peak_voltage, &min_lower_electrode_peak_voltage, &max_lower_electrode_peak_voltage)){

                    }
                }
                if(CollapsingHeader("Upper Electrode")){
                    if(SliderScalar(" X position Upper Electrode",ImGuiDataType_Double, &simulation.upper_electrode.position[0], &min_upper_electrode_x_position, &max_upper_electrode_x_position)){

                    }
                    if(SliderScalar(" Y position Upper Electrode",ImGuiDataType_Double, &simulation.upper_electrode.position[1], &min_upper_electrode_y_position, &max_upper_electrode_y_position)){

                    }
                    if(SliderScalar(" Charge upper Electrode",ImGuiDataType_Double, &simulation.upper_electrode.charge, &min_upper_electrode_charge, &max_upper_electrode_charge)){

                    }
                    if(SliderScalar(" Voltage Upper Electrode",ImGuiDataType_Double, &simulation.upper_electrode.voltage, &min_upper_electrode_voltage, &max_upper_electrode_voltage)){

                    }

                }
            }
            ImGui::EndMenu();
        }
        EndMainMenuBar();
    }

    void drawNanoVG() override {
        // draw Box
        nvgBeginPath(vg);
        nvgRect(vg, rect.center[0]-rect.size[0]/2.f + center_of_frame[0], rect.center[1]-rect.size[1]/2.f + center_of_frame[1], rect.size[0], rect.size[1]);
        nvgFillColor(vg, rect.colorFillBox);
        nvgFill(vg);
        nvgStrokeColor(vg, rect.colorStrokeBox);
        nvgStrokeWidth(vg, 10.0f);
        nvgStroke(vg);

        if(drawRectangle){
            auto drawRectangle = [this](const Rectangle &rectangle, bool is_lower_electrode,double scale, Vector3d distance){
                if(switch_view){
                    nvgBeginPath(vg);

                    nvgRect(vg,rectangle.position[0]  * scale + center_of_frame[0], rectangle.position[1] * scale + center_of_frame[1], rectangle.width,rectangle.depth);

                    if(is_lower_electrode){
                        nvgFillColor(vg,nvgRGBA(250,0,0,200));
                    }
                    else{
                        nvgFillColor(vg,nvgRGBA(0,0,250,200));
                    }
                    nvgFill(vg);
                    nvgStrokeColor(vg, rectangle.colorStroke);
                    nvgStrokeWidth(vg, 10.0f);
                    nvgStroke(vg);
                }else{
                    nvgBeginPath(vg);

                    nvgRect(vg,rectangle.position[0]  * scale + center_of_frame[0], rectangle.position[2] * scale + center_of_frame[2], rectangle.width,rectangle.height);

                    if(is_lower_electrode){
                        nvgFillColor(vg,nvgRGBA(250,0,0,200));
                    }
                    else{
                        nvgFillColor(vg,nvgRGBA(0,0,250,200));
                    }
                    nvgFill(vg);
                    nvgStrokeColor(vg, rectangle.colorStroke);
                    nvgStrokeWidth(vg, 10.0f);
                    nvgStroke(vg);
                }



            };
            Vector3d formation_center= (simulation.lower_electrode.position + simulation.upper_electrode.position) * 0.5;
            Vector3d distance_1 = formation_center - simulation.lower_electrode.position;
            Vector3d distance_2 = formation_center - simulation.upper_electrode.position;
            drawRectangle(Rectangle(simulation.lower_electrode.position, simulation.lower_electrode.width,simulation.lower_electrode.length,simulation.lower_electrode.depth),true, (1.0/mscale),distance_1);
            drawRectangle(Rectangle(simulation.upper_electrode.position, simulation.upper_electrode.width,simulation.upper_electrode.length,simulation.upper_electrode.depth),false,(1.0/mscale),distance_2);
        }

        if(drawCircles)
        {
            auto drawCircle = [this](const Circle &circle,int tmpfortest, double scale, bool is_tracer){
                if(switch_view){
                    nvgBeginPath(vg);
                    nvgCircle(vg, (circle.pos[0]    * scale) + center_of_frame[0], (circle.pos[1]   * scale) + center_of_frame[1], circle.radius *scale);
                    if(is_tracer){
                        std::ranlux48 gen;
                        std::uniform_int_distribution<int>  uniform_0_255(0, 255);
                        nvgFillColor(vg, nvgRGBA(102, 217, 239,100));
                    }else {

                        if (tmpfortest == 1) {
                            nvgFillColor(vg, nvgRGBA(150, 0, 0, 100));
                        } else {
                            nvgFillColor(vg, nvgRGBA(0, 0, 150, 100));
                        }
                    }
                    nvgFill(vg);
                    nvgStrokeColor(vg, circle.colorStroke);
                    nvgStrokeWidth(vg, 10.0f);
                    nvgStroke(vg);
                }else{
                    nvgBeginPath(vg);
                    nvgCircle(vg, (circle.pos[0]   * scale) + center_of_frame[0], (circle.pos[2]   * scale) + center_of_frame[2], circle.radius *scale);
                    if(is_tracer){
                        nvgFillColor(vg, nvgRGBA(150, 0, 0, 100));
                       // std::cout<<"true"<<std::endl;
                    }else {
                       // std::cout<<" false"<<std::endl;
                        if (tmpfortest == 1) {
                            nvgFillColor(vg, nvgRGBA(150, 0, 0, 100));
                        } else {
                            nvgFillColor(vg, nvgRGBA(0, 0, 150, 100));
                        }
                    }

                    nvgFill(vg);
                    nvgStrokeColor(vg, circle.colorStroke);
                    nvgStrokeWidth(vg, 10.0f);
                    nvgStroke(vg);
                }

            };

            /*for(const auto &p : connect_vector) //change to particles for std sim
                drawCircle(Circle(p.position,p.radius));*/

           /* Vector3d formation_center;
            formation_center.setZero();
            double objects_total = simulation.size/3;

            for(auto &pair : connected_particles){

                if(!std::get<0>(pair).visited) {
                    formation_center = formation_center + std::get<0>(pair).position;
                    std::get<0>(pair).visited = true;
                }
                if(!std::get<1>(pair).visited) {
                    formation_center = formation_center + std::get<1>(pair).position;
                    std::get<1>(pair).visited = true;
                }

            }
            formation_center[0] = formation_center[0]/(objects_total);
            formation_center[1] = formation_center[1]/(objects_total);
            formation_center[2] = formation_center[2]/(objects_total);

            Vector3d formation_center_2 = formation_center * 3;
            //std::cout<<" formation center 1 = "<<formation_center<<std::endl;
            //std::cout<<" formation center 2 = "<<formation_center_2<<std::endl;
            simulation.reset_flags(connected_particles);*/



            for(auto &pair : connected_particles) { //change to particles for std sim

                Particle& A = std::get<0>(pair);
                Particle& B = std::get<1>(pair);



               if(!A.visited) {
                  // Vector3d distance = formation_center - std::get<0>(pair).position;
                   //std::cout<<"distance = "<<distance<<std::endl;
                   drawCircle(Circle(A.position, A.radius), 1, (1.0/mscale),false);
                   A.visited = true;
               }
               if(!B.visited) {
                   //Vector3d distance = formation_center - std::get<1>(pair).position;
                   drawCircle(Circle(B.position, B.radius), 2, (1.0/mscale),false);
                   B.visited = true;
               }

            }
            simulation.reset_flags(connected_particles);

            if(simulation.tracer_on) {
                for (int i = 0; i < simulation.tracer.size(); i++) {
                    drawCircle(Circle(simulation.tracer[i], mscale* 0.6), 1, (1.0 / mscale),true);
                }
            }
            if(simulation.dimer_tracers){
                for(int i = 0; i<simulation.dimer_tracers_vector.size(); i++){
                    for(int j = 0; j<simulation.dimer_tracers_vector[i].size();j++){
                        drawCircle(Circle(simulation.dimer_tracers_vector[i][j], mscale * 0.6), 1, (1.0/mscale), true);
                    }
                }
            }
            /*for(auto &particle : particles){
                drawCircle(Circle(particle.position, particle.radius));
            }*/ //uncomment this if you want also single particles

          /*  for(auto &pair : connected_particles) { //change to particles for std sim
                    pair.first.is_drawn = false;
                    pair.second.is_drawn = false;
            }*/





            // deleted struct keyword and changed everthing to vector2d

        }

        if(drawCursor){
            nvgBeginPath(vg);
            nvgCircle(vg, mouseState.lastMouseX, mouseState.lastMouseY, 10.f);
            nvgFillColor(vg, rect.colorStrokeBox);
            nvgFill(vg);
        }

       /* if(DRAW_DRAG){
            for(auto& particle_pair: connected_particles) {
                Particle& A = std::get<0>(particle_pair);
                Particle& B = std::get<1>(particle_pair);
                if(!A.visited){

                    Vector3d current_drag = simulation.get_stokes_friction_force(A);
                    //drawVelocity({A.position[0] * (1.0/mscale), A.position[1]* (1.0/mscale)}, {current_drag[0] * (1.0/mscale), current_drag[1] * (1.0/mscale)}, nvgRGB(102, 217, 239));
                    A.visited = true;
                    nvgBeginPath(vg);
                    nvgMoveTo(vg,A.position[0] * (1.0/mscale),A.position[1] * (1.0/mscale));
                    nvgLineTo(vg,current_drag[0] * (1.0/mscale) + A.position[0] * (1.0/mscale), current_drag[1] * (1.0/mscale) + A.position[1] * (1.0/mscale));

                }
                if(!B.visited){
                    Vector3d current_drag = simulation.get_stokes_friction_force(B);
                    //drawVelocity({B.position[0] * (1.0/mscale), B.position[1] * (1.0/mscale)}, {current_drag[0] * (1.0/mscale), current_drag[1] * (1.0/mscale)}, nvgRGB(102, 217, 239));
                    B.visited = true;
                    nvgBeginPath(vg);
                    nvgMoveTo(vg,B.position[0] * (1.0/mscale),B.position[1] * (1.0/mscale));
                    nvgLineTo(vg,current_drag[0] * (1.0/mscale), current_drag[1] * (1.0/mscale));

                }
            }
            simulation.reset_flags(connected_particles);
        }*/
    }

protected:

    void resizeWindow(int w, int h) override {
        Application::resizeWindow(w, h);

    }

    void scrollWheel(double  xoffset, double yoffset) override {
        double zoomOld = zoom;
        zoom *= std::pow(1.10, yoffset);
        double cursorPos[2] = { mouseState.lastMouseX, mouseState.lastMouseY };
        for (int dim = 0; dim < 2; ++dim) {
            double c = cursorPos[dim]/(double) ((dim == 0) ? base : -base);
            translation[dim] = c - zoomOld/zoom * (c-translation[dim]);
        }
    }

    void mouseButtonPressed(int button, int mods) override {
        cursorPosDown[0] = mouseState.lastMouseX;
        cursorPosDown[1] = mouseState.lastMouseY;

        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            Vector2d cursor = fromScreen(mouseState.lastMouseX, mouseState.lastMouseY);
        }

    }

    void keyPressed(int key, int  /*mods*/) override {
        // play / pause with space bar
        if(key == GLFW_KEY_SPACE)
            runSim = !runSim;

        if (key == GLFW_KEY_R){
            runSim =!runSim;
            simulation.reset_simulation(connected_particles);
        }
        if(key == GLFW_KEY_LEFT){
            center_of_frame[0] = center_of_frame[0] +50.0;
        }
        if(key == GLFW_KEY_RIGHT){
            center_of_frame[0] = center_of_frame[0] -50.0;
        }

        if(key == GLFW_KEY_UP){
            center_of_frame[1] = center_of_frame[1] -50.0;
        }
        if(key == GLFW_KEY_DOWN){
            center_of_frame[1] = center_of_frame[1] +50.0;
        }

    }



    /*void mouseButtonPressed(int button, int mods) override {
        Vector2f x = Vector2f(mouseState.lastMouseX, mouseState.lastMouseY);
        if(button == GLFW_MOUSE_BUTTON_LEFT && circleMouse.isInside(x)) {
            draggingCircle = true;
            draggingCircleOffset = x - circleMouse.pos;
        }
    }*/ //commented out

    void mouseButtonReleased(int button, int mods) override {
        draggingCircle = false;
    }

private:
   // bool DRAW_THRUSTER = true;
  //  bool DRAW_VELOCITY = false;
   // bool DRAW_DRAG = true;
    int loadFonts(NVGcontext* vg)
    {
        int font;
        font = nvgCreateFont(vg, "sans", "../example/Roboto-Regular.ttf");
        if (font == -1) {
            printf("Could not add font regular.\n");
            return -1;
        }
        font = nvgCreateFont(vg, "sans-bold", "../example/Roboto-Bold.ttf");
        if (font == -1) {
            printf("Could not add font bold.\n");
            return -1;
        }
        return 0;
    }

  //  void drawVelocity(const Vector2d &s_, const Vector2d &v_, const NVGcolor &COLOR) { drawVector(s_, .1*v_, COLOR); }

  /*  double eps = .25; // * otherwise stuff doesn't show up
    double get_L_() { return double(std::min(height, width)); }
    double get_f_() { return get_L_() / (2. + 2 * eps); }
    double ZOOM_ = 1.f;
    double ZOOM_MIN_ = .1f;
    double ZOOM_MAX_ = 2.f;
    double ZOOM() { return 1. / ZOOM_; }
    Vector2f _2nvg(Vector2d xy) {
        // zoom*(-1. - eps, 1. + eps) -x-> (0, L)
        // ""                         -y-> (L, 0)
        xy /= ZOOM();
        Vector2f ret = (get_f_() * (xy + (1. + eps)*Vector2d::Ones())).cast<float>();
        ret.y() = get_L_() - ret.y();
        return ret;
    }

    float CIRCLE_RADIUS = 4.f;

    void drawVector(const Vector2d &s_, const Vector2d &F_, const NVGcolor &COLOR) {
        Vector2f s = _2nvg(s_);
        Vector2f t = _2nvg(s_ + F_);
        Vector2f st = t - s;
        Vector2f e = CIRCLE_RADIUS * Vector2f(-st.y(), st.x()).normalized();
        Vector2f sP = s + e;
        Vector2f sM = s - e;
        // --
        nvgReset(vg);
        nvgBeginPath(vg);
        nvgMoveTo(vg, t.x(), t.y());
        nvgLineTo(vg, sP.x(), sP.y());
        nvgLineTo(vg, sM.x(), sM.y());
        nvgLineTo(vg, t.x(), t.y());
        nvgFillColor(vg, COLOR);
        nvgFill(vg);
    }

*/

private:
    VectorXd fromScreen(int i, int j, int w, int h) const {
        VectorXd x(2);
        x[0] = ((double)i/(double)w - translation[0])*zoom/pixelRatio;
        x[1] = (-(double)j/(double)h - translation[1])*zoom/pixelRatio;
        return x;
    }

    template<class S>
    VectorXd fromScreen(S i, S j) const {
        return fromScreen((double)i, (double)j, base, base);
    }

    double toScreen(double s, int dim) const {
        return (s/zoom*pixelRatio + translation[dim]) * (double)((dim == 0) ? base : -base);
    }

    double toScreen(double s) const {
        return s/zoom*pixelRatio * base;
    }

public:
    Vector3d center_of_frame = {1000,1000,1000};
    bool runSim = false;
    double cursorPosDown[2]{};
    double translation[2] = {0.75*pixelRatio, -0.25*pixelRatio};
    double zoom = 10;
    int base;

    double min_my_zoom = std::pow(10,-6);
    double max_my_zoom = std::pow(10,-2);

    //Particle BIG_A = Particle(100,40, {520,-100});
    //Particle SMALL_B = Particle(10,15,{500,-90});

    double mscale = std::pow(10,-6);

    //flat triangle from laura
   //  Particle A = Particle(1050 ,3.0 * mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
   //  Particle B = Particle(1050,2.0* mscale,-1.0,{105.0 * mscale,100.0 * mscale, 100.0 * mscale});
   //  Particle C = Particle(1050,3.0* mscale,-1.0,{110.0 * mscale,100.0 * mscale, 100.0 * mscale});

    //first triangle from laura
  //  Particle A_t = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
  //Particle B_t = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
  //  Particle C_t = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,102.0 * mscale, 100.0 * mscale});



    //dumbel
   Particle A = Particle(1050 ,1.0 * mscale,-1.0,{101.7 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B = Particle(1050,0.7* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

    Particle A_2 = Particle(1050 ,1.0 * mscale,-1.0,{201.7 * mscale,200.0 * mscale, 100.0 * mscale});
    Particle B_2 = Particle(1050,0.7* mscale,-1.0,{200.0 * mscale,200.0 * mscale, 100.0 * mscale});

    Particle A_3 = Particle(1050 ,1.0 * mscale,-1.0,{1.7 * mscale,0.0 * mscale, 100.0 * mscale});
    Particle B_3 = Particle(1050,0.7* mscale,-1.0,{0.0 * mscale,0.0 * mscale, 100.0 * mscale});

    Particle A_4 = Particle(1050 ,1.0 * mscale,-1.0,{101.7 * mscale,-100.0 * mscale, 100.0 * mscale});
    Particle B_4 = Particle(1050,0.7* mscale,-1.0,{100.0 * mscale,-100.0 * mscale, 100.0 * mscale});

    Particle A_5 = Particle(1050 ,1.0 * mscale,-1.0,{301.7 * mscale,-100.0 * mscale, 100.0 * mscale});
    Particle B_5 = Particle(1050,0.7* mscale,-1.0,{300.0 * mscale,-100.0 * mscale, 100.0 * mscale});

    Particle A_6 = Particle(1050 ,1.0 * mscale,-1.0,{101.7 * mscale,-200.0 * mscale, 100.0 * mscale});
    Particle B_6 = Particle(1050,0.7* mscale,-1.0,{100.0 * mscale,-200.0 * mscale, 100.0 * mscale});

    Particle A_7 = Particle(1050 ,1.0 * mscale,-1.0,{301.7 * mscale,200.0 * mscale, 100.0 * mscale});
    Particle B_7 = Particle(1050,0.7* mscale,-1.0,{300.0 * mscale,200.0 * mscale, 100.0 * mscale});

    Particle A_8 = Particle(1050 ,1.0 * mscale,-1.0,{-151.7 * mscale,-300.0 * mscale, 100.0 * mscale});
    Particle B_8 = Particle(1050,0.7* mscale,-1.0,{-150.0 * mscale,-300.0 * mscale, 100.0 * mscale});

    Particle A_9 = Particle(1050 ,1.0 * mscale,-1.0,{1.7 * mscale,250.0 * mscale, 100.0 * mscale});
    Particle B_9 = Particle(1050,0.7* mscale,-1.0,{0.0 * mscale,250.0 * mscale, 100.0 * mscale});

    Particle A_10 = Particle(1050 ,1.0 * mscale,-1.0,{-101.7 * mscale,300.0 * mscale, 100.0 * mscale});
    Particle B_10 = Particle(1050,0.7* mscale,-1.0,{-100.0 * mscale,300.0 * mscale, 100.0 * mscale});





/*    Particle A_5 = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B_5 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

    Particle A_6 = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B_6 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

    Particle A_7 = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B_7 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

    Particle A_8 = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B_8 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

    Particle A_9 = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B_9 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

    Particle A_10 = Particle(1050 ,1.5 * mscale,-1.0,{102.5 * mscale,100.0 * mscale, 100.0 * mscale});
    Particle B_10 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});*/

   // Particle A_2 = Particle(1050 ,1.5 * mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
   // Particle B_2 = Particle(1050,1.0* mscale,-1.0,{/*(100.0 + 5/std::sqrt(2))*/ 100 * mscale,/*(100.0 + 5/std::sqrt(2))*/102.5 * mscale, 100.0 * mscale});

  //  Particle A_3 = Particle(1050 ,1.5 * mscale,-1.0,{97.5 * mscale,100.0 * mscale, 100.0 * mscale});
   // Particle B_3 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

   // Particle A_4 = Particle(1050 ,1.5 * mscale,-1.0,{100.0 * mscale,97.5 * mscale, 100.0 * mscale});
   // Particle B_4 = Particle(1050,1.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});

   //spiral triangle that doesn not do spirals
   //  Particle A = Particle(1050 ,2.0 * mscale,-1.0,{104.0 * mscale,100.0 * mscale, 100.0 * mscale});
   //  Particle B = Particle(1050,2.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
    // Particle C = Particle(1050,3.0* mscale,-1.0,{(104.0  + 5.0/std::sqrt(2.0))* mscale,(100.0 + 5.0/std::sqrt(2.0)) * mscale, 100.0 * mscale});

    //spiral triangle that doesn not do spirals
    //Particle A = Particle(1050 ,3.0 * mscale,-1.0,{105.0 * mscale,100.0 * mscale, 100.0 * mscale});
    //Particle B = Particle(1050,2.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
    //Particle C = Particle(1050,2.0* mscale,-1.0,{(105.0  + 4.0/std::sqrt(2.0))* mscale,(100.0 + 4.0/std::sqrt(2.0)) * mscale, 100.0 * mscale});


    //rotator
 // Particle A = Particle(1050 ,3.0 * mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
 // Particle B = Particle(1050,2.0* mscale,-1.0,{105.0 * mscale,100.0 * mscale, 100.0 * mscale});
 // Particle C = Particle(1050,2.0* mscale,-1.0,{105.0 * mscale,104.0 * mscale, 100.0 * mscale});
 // Particle D = Particle(1050,3.0* mscale,-1.0,{110.0 * mscale,104.0 * mscale, 100.0 * mscale});


    //l shaped
   // Particle A = Particle(1050 ,2.0 * mscale,-1.0,{104.0 * mscale,100.0 * mscale, 100.0 * mscale});
   // Particle B = Particle(1050,2.0* mscale,-1.0,{100.0 * mscale,100.0 * mscale, 100.0 * mscale});
   // Particle Up = Particle(1050,3.0* mscale,-1.0,{100.0 * mscale,105.0 * mscale, 100.0 * mscale});
   // Particle C = Particle(1050,2.0* mscale,-1.0,{108.0 * mscale,100.0 * mscale, 100.0 * mscale});
   // Particle D = Particle(1050,2.0* mscale,-1.0,{112.0 * mscale,100.0 * mscale, 100.0 * mscale});



    std::vector<Particle> particles;
    std::vector<std::tuple <Particle&,Particle&,double> > connected_particles;
    Simulation simulation;







    bool drawCursor = false;
    bool drawCircles = true;
    bool drawRectangle = true;

    double min_stiffnes_constant = 1.0;
    double max_stiffnes_constant = 4000.0;
    double min_damping_constant = 0.0;
    double max_damping_constant = 100.0;
    double min_viscosity_multiplier = 0.01;
    double max_viscosity_multiplier = 1000.0;
    double min_brownian_motion_multiplier = 0.0;
    double max_brownian_motion_multiplier = 10.0;
    double min_time_step = 0.000001;
    double max_time_step = 1.0;
    int min_max_iterations = 100;
    int max_max_iterations = 1000000;
    double min_lower_electrode_x_position = -1000.0 * mscale;
    double max_lower_electrode_x_position = 1000.0* mscale;
    double min_lower_electrode_y_position = -1000.0* mscale;
    double max_lower_electrode_y_position = 1000.0* mscale;
    double min_upper_electrode_x_position = -1000.0* mscale;
    double max_upper_electrode_x_position = 1000.0* mscale;
    double min_upper_electrode_y_position = -1000.0* mscale;
    double max_upper_electrode_y_position = 1000.0* mscale;
    double min_lower_electrode_charge = -10.0;
    double max_lower_electrode_charge = 10.0;
    double min_upper_electrode_charge = -10.0;
    double max_upper_electrode_charge = 10.0;

    double min_sigma_m = 0.5 * std::pow(10,-4);
    double max_sigma_m = 1.0 * std::pow(10,-2);
    double min_epsilon_m = 1.0;
    double max_epsilon_m = 100.0;

    double min_cof_x = -1000.0;
    double max_cof_x = 2000.0;
    double min_cof_y = -1000.0;
    double max_cof_y = 2000.0;
    double min_cof_z = -1000.0;
    double max_cof_z = 2000.0;

    double min_zeta = -100 * std::pow(10,-3);
    double max_zeta = 100 * std::pow(10,-3);


    double min_lower_electrode_voltage = 0.0;
    double max_lower_electrode_voltage = 300.0;
    double min_upper_electrode_voltage = -300.0;
    double max_upper_electrode_voltage = 0.0;
    double min_mass = std::pow(10.0,-15.0);
    double max_mass = 100 * mscale;
    double min_charge = -30;
    double max_charge = 30.0;
    double min_radius = 1.0 * mscale;
    double max_radius = 200.0 * mscale;
    double min_lower_electrode_frequency = std::pow(10.0,2.0);
    double max_lower_electrode_frequency = std::pow(10.0,9.0);
    double min_lower_electrode_peak_voltage = 0.0;
    double max_lower_electrode_peak_voltage = 1000.0;
    double min_permittivity = 0.0;
    double max_permittivity = 200.0;
    double min_conductivity = 0.05 * std::pow(10,-9);
    double max_conductivity = 10 * std::pow(10,-9);
    double min_T = 1.0;
    double max_T = 373.0;

    double min_stern_layer_conductance = 0.05 * std::pow(10,-9);
    double max_stern_layer_conductance = 10 * std::pow(10,-9);

    bool is_silica_1 = false;
    bool is_silica_2 = false;
    bool first_experiment = false;
    bool second_experiment = true;
    bool is_PS_1 = false;
    bool is_PS_2 = false;
    bool run_w_exp = false;
    bool hundred_step = false;
    bool five_hundred_step = false;
    bool thousand_step = false;
    bool ten_thousand_step = false;
    bool hundred_thousand_step = false;

    bool set_t_to_1ms = false;
    bool set_t_to_2ms = false;
    bool set_t_to_5ms = false;
    bool set_t_to_10ms = false;
    bool set_t_to_20ms = false;



    bool switch_view = false;
    bool follow_cam = false;







    struct Rectangle
    {
        Rectangle(Vector3d p, double w, double h, double d){
            this->position = p;
            this->width = w;
            this->height = h;
            this->depth = d;
        }
        Vector3d position;
        double height;
        double width;
        double depth;
        NVGcolor colorFill = nvgRGBA(50, 50, 50, 100) , colorStroke = nvgRGBA(50, 50, 50, 100);
    };

    struct Circle
    {
        Circle(Vector3d p, double r=10){ //took away struct
            this->pos = p;
            this->radius = r;
        }
        bool isInside(const Vector3d &x){
            return (x-pos).squaredNorm() <= radius*radius;
        }

        Vector3d pos;
        double radius;
        NVGcolor colorFill= nvgRGBA(50, 50, 50, 100) , colorStroke = nvgRGBA(50, 50, 50, 100);

    }; //deleted circleKey circle MOuse
    struct Box {
        bool isInside(const Vector2d &x, float slack = 0.f) {
            for (int i = 0; i < 2; ++i)
                if(x[i]-slack < center[i]-size[i]/2 || x[i]+slack > center[i]+size[i]/2)
                    return false;
            return true;
        }

        Vector2d center;
        Vector2d size;
        NVGcolor colorFillBox = nvgRGBA(0,0,0,255) , colorStrokeBox = nvgRGBA(50, 50, 50, 100);
    } rect;

    bool draggingCircle = false;
    Vector2d draggingCircleOffset;
    std::deque<std::pair<Vector2d, Circle>> circles;
    std::chrono::high_resolution_clock::time_point lastFrame;
};

int main(int, char**)
{





    //simulation.run_simulation_with_brownian_motion(particles);
    //simulation.run_simulation_with_brownian_motion_only(particles);
    std::cout<<"end of sim 1" <<std::endl<<std::endl;

    //simulation.run_simulation_with_drag(particles);


    TestApp app(1200, 800, "Assignment 0");
    app.run();

  /*  std::ofstream myfile;
    myfile.open ("../../../data.csv");
    myfile << "Total Energy"<<","<<"Kinetic Energy"<<","<<"Potential Energy"<<","<<"E"<<","<<"SPRING FORCE X"<<","<<"SPRING FORCE Y"<<","<<"DERIVATIVE"<<","<<"X"<<","<<"Y"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.total_energy.size(); i++){
        std::cout<<i<<std::endl;
        //std::cout<<app.simulation.spring_force_vec_over_time_x[i]<<std::endl;
        //std::cout<<app.simulation.spring_force_vec_over_time_x[i]/app.simulation.stiffnes_constant<<std::endl;
        myfile << app.simulation.total_energy[i]<<","<<app.simulation.kinetic_energy[i]<<","<<app.simulation.potential_energy[i]<<","<<app.simulation.e_vec[i]<<","<<app.simulation.spring_force_vec_over_time_x[i]/app.simulation.stiffnes_constant<<","<<app.simulation.spring_force_vec_over_time_y[i]/app.simulation.stiffnes_constant<<","<<app.simulation.spring_force_derivative_x_in_x[i]/app.simulation.stiffnes_constant<<","<<app.simulation.position_vec_over_time_in_x[i]<<","<<app.simulation.position_vec_over_time_in_y[i]<<","<<i<<std::endl;
    std::cout<<"after i"<<std::endl;
    }
    myfile.close();*/

    std::ofstream curvature;
    curvature.open("../../../curvature.csv");
    curvature<<"x"<<","<<"y"<<","<<"dumbbell_index"<<","<<"time_passed"<<std::endl;
    for(int i = 0; i<app.simulation.com_x_vector.size(); i++){
        curvature<<app.simulation.com_x_vector[i]<<","<<app.simulation.com_y_vector[i]<<","<<app.simulation.dumbbell_index_vector[i]<<","<<app.simulation.time_passed_vector[i]<<std::endl;
        std::cout<<" curvature x"<<app.simulation.com_x_vector[i]<<std::endl;
    }
    curvature.close();

    std::cout<<"ttest"<<std::endl;
    std::ofstream FMA;
    FMA.open ("../../../FMA.csv");
    FMA <<"F-M*a"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.F_ma.size(); i++){
        FMA<<app.simulation.F_ma[i]<<","<<i<<std::endl;
    }
    FMA.close();

    std::ofstream brownian;
    brownian.open("../../../brownian.csv");
    brownian <<"Brownian Force X"<<","<<"Brownian Force Y"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.b_force_x.size(); i++){
        brownian<<app.simulation.b_force_x[i]<<","<<app.simulation.b_force_y[i]<<","<<i<<std::endl;
    }
    brownian.close();

    std::ofstream ehd;
    ehd.open("../../../ehd.csv");
    ehd <<"EHD Force X"<<","<<"EHD Force Y"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.ehd_force_x.size(); i++){
        ehd<<app.simulation.ehd_force_x[i]<<","<<app.simulation.ehd_force_y[i]<<","<<i<<std::endl;
    }
    ehd.close();

    std::ofstream spring;
    spring.open("../../../spring.csv");
    spring <<"Spring Force X"<<","<<"Spring Force Y"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.spring_force_x.size(); i++){
        spring<<app.simulation.spring_force_x[i]<<","<<app.simulation.spring_force_y[i]<<","<<i<<std::endl;
        std::cout<<" spring "<<app.simulation.spring_force_x[i]<<std::endl;
    }
    spring.close();

    std::ofstream amount;
    amount.open("../../../amount.csv");
    amount <<"Amount X"<<","<<"Amount Y"<<","<<"Amount T"<<","<<"Norm"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.amount_x.size(); i++){
        amount<<app.simulation.amount_x[i]<<","<<app.simulation.amount_y[i]<<","<<app.simulation.amount_z[i]<<","<<app.simulation.direction_norm[i]<<","<<i<<std::endl;
    }
    amount.close();

    std::ofstream drag;
    drag.open("../../../drag.csv");
    drag <<"Drag X"<<","<<"Drag Y"<<","<<"Drag Px"<<","<<"Drag Py"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.drag_x.size(); i++){
        drag<<app.simulation.drag_x[i]<<","<<app.simulation.drag_y[i]<<","<<app.simulation.drag_px[i]<<","<<app.simulation.drag_py[i]<<","<<i<<std::endl;
    }
    drag.close();

    std::ofstream rotator_position;
    rotator_position.open("../../../rotator_position.csv");
    rotator_position<<"Position in X"<<","<<"Position in Y"<<","<<"Position in X 2"<<","<<"Position in Y 2"<<","<<"Time Step"<<std::endl;
    /*for(int i = 0; i<app.simulation.rotator_position_vec_over_time_in_x.size(); i++){
        rotator_position<<app.simulation.rotator_position_vec_over_time_in_x[i]<<","<<app.simulation.rotator_position_vec_over_time_in_y[i]<<","<<app.simulation.rotator_position_vec_over_time_in_x_2[i]<<","<<app.simulation.rotator_position_vec_over_time_in_y_2[i]<<","<<i<<std::endl;
    std::cout<<"rot"<<app.simulation.rotator_position_vec_over_time_in_y_2[i]<<std::endl;
    }*/
    rotator_position.close();

    std::ofstream position;
    position.open("../../../position.csv");
   position<<"Position in X"<<","<<"Position in Y"<<","<<"Position in X 2"<<","<<"Position in Y 2"<<","<<"Position in X 3"<<","<<"Position in Y 3"<<","<<"Position in X 4"<<","<<"Position in Y 4"<<","<<"Time Step"<<std::endl;
    for(int i = 0; i<app.simulation.position_vec_over_time_in_x.size(); i++){
        position<<app.simulation.position_vec_over_time_in_x[i]<<","<<app.simulation.position_vec_over_time_in_y[i]<<","<<app.simulation.position_vec_over_time_in_x_2[i]<<","<<app.simulation.position_vec_over_time_in_y_2[i]<<","<<app.simulation.position_vec_over_time_in_x_3[i]<<","<<app.simulation.position_vec_over_time_in_y_3[i]<<","<<app.simulation.position_vec_over_time_in_x_4[i]<<","<<app.simulation.position_vec_over_time_in_y_4[i]<<","<<i<<std::endl;
        std::cout<<"pos 4"<<app.simulation.position_vec_over_time_in_x_4[i]<<","<<app.simulation.position_vec_over_time_in_y_4[i]<<std::endl;
    }
    position.close();

    std::ofstream triangle_position;
    triangle_position.open("../../../triangle_position.csv");
    triangle_position<<"Triangle Position in X"<<","<<"Triangle Position in Y"<<","<<"Triangle Position in X 2"<<","<<"Triangle Position in Y 2"<<","<<"Triangle Position in X 3"<<","<<"Triangle Position in Y 3"<<","<<"COM Position in X 4"<<","<<"COM Position in Y 4"<<","<<"Time Step"<<std::endl;
  /* for(int i = 0; i<app.simulation.triangle_position_vec_over_time_in_x.size(); i++){
        double com_x = app.simulation.triangle_position_vec_over_time_in_x[i] + app.simulation.triangle_position_vec_over_time_in_x_2[i] +app.simulation.triangle_position_vec_over_time_in_x_3[i];
        double com_y = app.simulation.triangle_position_vec_over_time_in_y[i] + app.simulation.triangle_position_vec_over_time_in_y_2[i] +app.simulation.triangle_position_vec_over_time_in_y_3[i];
        com_x = com_x/3.0;
        com_y = com_y/3.0;
        triangle_position<<app.simulation.triangle_position_vec_over_time_in_x[i]<<","<<app.simulation.triangle_position_vec_over_time_in_y[i]<<","<<app.simulation.triangle_position_vec_over_time_in_x_2[i]<<","<<app.simulation.triangle_position_vec_over_time_in_y_2[i]<<","<<app.simulation.triangle_position_vec_over_time_in_x_3[i]<<","<<app.simulation.triangle_position_vec_over_time_in_y_3[i]<<","<<com_x<<","<<com_y<<","<<i<<std::endl;
        std::cout<<"com x"<<com_x<<","<<com_y<<std::endl;
    }*/
    triangle_position.close();



    std::ofstream v_w;
    v_w.open("../../../v_w.csv");
    v_w <<" Velocity"<<","<<"Frequency"<<std::endl;
    for(int i = 0; i<app.simulation.frequencies.size();i++){
        v_w<<app.simulation.velocities[i]<<","<<app.simulation.frequencies[i]<<std::endl;
        std::cout<<app.simulation.frequencies[i]<<std::endl;
    }
    v_w.close();

    std::ofstream v_w_2;
    v_w_2.open("../../../v_w_2.csv");
    v_w_2<<"U"<<","<<"frequency"<<","<<"U A"<<","<<"U B"<<","<<"RE CM"<<","<<"IM CM"<<std::endl;
    for(int i=0; i<app.simulation.U_vec.size(); i++){
        v_w_2<<app.simulation.U_vec[i]<<","<<app.simulation.frequencies2[i]<<","<<app.simulation.UA_vec[i]<<","<<app.simulation.UB_vec[i]<<","<<app.simulation.Re_CM[i]<<","<<app.simulation.Im_CM[i]<<std::endl;
    }
    v_w_2.close();

    std::ofstream v_E;
    v_E.open("../../../v_E.csv");
    v_E <<" Velocity"<<","<<"E Field"<<std::endl;
    for(int i = 0; i<app.simulation.Efields.size();i++){
        v_E<<app.simulation.velocities_2[i]<<","<<app.simulation.Efields[i]<<std::endl;
    }
    v_E.close();

    std::ofstream velocities;
    velocities.open ("../../../velocities.csv");
    std::cout<<"test2"<<std::endl;

    velocities << "VELX"<<","<<"VELY"<<","<<"friction"<<","<<"time_step"<<std::endl;
    /*for(int i = 0; i<app.simulation.velocities_over_time1_in_x.size();i++){
        velocities<<app.simulation.velocities_over_time1_in_x[i]<<","<<app.simulation.velocities_over_time1_in_y[i]<<","<<app.simulation.friction_force_over_time_x[i]<<","<<i<<std::endl;
        std::cout<<app.simulation.velocities_over_time1_in_x[i]<<std::endl;
    }
    velocities.close();*/
    std::cout<<" test test 1"<<std::endl;

    std::ofstream myfile2;
    myfile2.open("../../../data2.csv");
    myfile2<<"x"<<","<<"y"<<","<<"z1"<<","<<"z2"<<std::endl;
    /*for(int i = 0; i<899;i++){
        myfile2<<app.simulation.x_values[i]<<","<<app.simulation.y_values[i]<<","<<app.simulation.z_values1[i]<<","<<app.simulation.z_values2[i]<<std::endl;
        //std::cout<<app.simulation.x_values[i]<<std::endl;
       // std::cout<<app.simulation.y_values[i]<<std::endl;
        //std::cout<<app.simulation.z_values1[i]<<std::endl;
    }*/
    myfile2.close();
    std::cout<<" test test 1"<<std::endl;

    std::ofstream myfile3;
    myfile3.open("../../../data3.csv");
    /*for(int i = 0; i<app.simulation.x_values.size();i++){
        myfile3<<app.simulation.x_values[i]<<","<<app.simulation.y_values[i]<<","<<app.simulation.z_values[i]<<std::endl;
        //std::cout<<app.simulation.x_values[i]<<std::endl;
        //std::cout<<app.simulation.y_values[i]<<std::endl;
        //std::cout<<app.simulation.z_values1[i]<<std::endl;
    }*/
    myfile3.close();
    std::cout<<" test test 1"<<std::endl;

    std::ofstream best;
    best.open("../../../best.csv");
    best<<"x"<<","<<"y"<<","<<"z"<<std::endl;
   // best<<app.simulation.best_found[0]<<","<<app.simulation.best_found[1]<<","<<app.simulation.best_found[2]<<std::endl;
    best.close();

    std::cout<<" test test 1"<<std::endl;

    std::ofstream FD;
    FD.open("../../../FD.csv");
    /*FD<<"F_X"<<","<<"DEDX"<<","<<"Difference"<<","<<"FORCE_JACOBIAN"<<","<<"DFDX"<<","<<"DIFFERENCE2"<<","<<"i"<<std::endl;
    for(int i = 0; i<app.simulation.forces_FD.size();i++){
        FD<<app.simulation.forces_FD[i]<<","<<app.simulation.energy_FD[i]<<","<<std::abs(app.simulation.forces_FD[i] - app.simulation.energy_FD[i])<<","<<app.simulation.force_jacobian_FD[i]<<","<<app.simulation.dFx_dx_FD[i]<<","<<std::abs(app.simulation.force_jacobian_FD[i] - app.simulation.dFx_dx_FD[i])<<","<<i<<std::endl;
    }*/
    FD.close();
    std::cout<<" test test 1"<<std::endl;

    std::ofstream test_e;
    test_e.open("../../../test_e.csv");
    test_e<<"ex1"<<","<<"ex2"<<","<<"value"<<std::endl;
    /*for(int i = 0; i<app.simulation.e_z.size(); i++){
        test_e<<app.simulation.e_x[i]<<","<<app.simulation.e_y[i]<<","<<app.simulation.e_z[i]<<std::endl;
        //std::cout<<app.simulation.e_z[i];
    }*/
    test_e.close();
    std::cout<<" test test 1"<<std::endl;

    std::ofstream sigma_m;
    sigma_m.open("../../../sigma_m.csv");
    sigma_m<<"Medium conductivity"<<","<<"Velocity"<<std::endl;
    for(int i = 0; i<app.simulation.sigma_m.size(); i++){
        sigma_m<<app.simulation.sigma_m[i]<<","<<app.simulation.U_sm[i]<<std::endl;
    }
    sigma_m.close();
    std::cout<<" test test 1"<<std::endl;

    std::ofstream sigma_p;
    sigma_p.open("../../../sigma_p.csv");
    sigma_p<<"Particle conductivity"<<","<<"Velocity"<<std::endl;
    for(int i = 0; i<app.simulation.sigma_p.size(); i++){
        sigma_p<<app.simulation.sigma_p[i]<<","<<app.simulation.U_sp[i]<<std::endl;
    }
    sigma_p.close();
    std::cout<<" test test 1"<<std::endl;

    std::ofstream epsilon_m;
    epsilon_m.open("../../../epsilon_m.csv");
    epsilon_m<<"Relative medium permittivity"<<","<<"Velocity"<<std::endl;
    for(int i = 0; i<app.simulation.epsilon_m.size(); i++){
        epsilon_m<<app.simulation.epsilon_m[i]<<","<<app.simulation.U_em[i]<<std::endl;
    }
    epsilon_m.close();
    std::cout<<" test test epsilon m"<<std::endl;

    std::ofstream epsilon_p;
    epsilon_p.open("../../../epsilon_p.csv");
    epsilon_p<<"Relative particle permittivity"<<","<<"Velocity"<<std::endl;
    for(int i = 0; i<app.simulation.epsilon_p.size(); i++){
        epsilon_p<<app.simulation.epsilon_p[i]<<","<<app.simulation.U_ep[i]<<std::endl;
    }
    epsilon_p.close();
    std::cout<<" test test epsilon p"<<std::endl;

    std::ofstream R_A;
    R_A.open("../../../R_A.csv");
    R_A<<"Particle radius"<<","<<"Velocity"<<std::endl;
    for(int i = 0; i<app.simulation.R_A.size(); i++){
        R_A<<app.simulation.R_A[i]<<","<<app.simulation.U_ra[i]<<std::endl;
    }
    R_A.close();
    std::cout<<" test test final"<<std::endl;




    return 0;
}
