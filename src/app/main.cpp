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


       // particles.push_back(particle_1);
       // particles.push_back(particle_2);
       /* particles.push_back(particle_3);
        particles.push_back(particle_4);*/
       /* particles.push_back(particle_5);
        particles.push_back(particle_6);
        particles.push_back(particle_7);
        particles.push_back(particle_8);
        particles.push_back(particle_9);
        particles.push_back(particle_10);*/
       //connect_vector.push_back(A);
       //connect_vector.push_back(B);
      // particles.push_back(loner);
       simulation.connect(A,B,(A.radius+B.radius),connected_particles);
//       simulation.connect(B,C,(C.radius+B.radius),connected_particles);
     //  simulation.connect(A,C, (A.radius+C.radius),connected_particles);
/*
        simulation.connect(A2,B2,(A2.radius+B2.radius),connected_particles);
        simulation.connect(B2,C2,(C2.radius+B2.radius),connected_particles);
        simulation.connect(A2,C2,(A2.radius+C2.radius+2*B.radius),connected_particles);

        simulation.connect(A3,B3,(A3.radius+B3.radius),connected_particles);
        simulation.connect(B3,C3,(C3.radius+B3.radius),connected_particles);
        simulation.connect(A3,C3,(A3.radius+C3.radius),connected_particles);
        simulation.connect(C3,D3,C3.radius + D3.radius, connected_particles);
        simulation.connect(A3,D3,std::sqrt(std::pow(A3.radius,2)/2 + std::pow((A3.radius/std::sqrt(2) + 2*D3.radius + C.radius), 2)), connected_particles);
        simulation.connect(B3,D3,std::sqrt(std::pow(B3.radius,2)/2 + std::pow((B3.radius/std::sqrt(2) + 2*D3.radius + C.radius), 2)), connected_particles);


        simulation.connect(A6,B6,A6.radius + B6.radius, connected_particles);
        simulation.connect(A6,C6,A6.radius + C6.radius, connected_particles);
        simulation.connect(C6,B6,C6.radius + B6.radius, connected_particles);
        simulation.connect(D6,B6,D6.radius + B6.radius, connected_particles);
        simulation.connect(D6,A6,D6.radius + A6.radius, connected_particles);


        simulation.connect(A4,B4,A4.radius + B4.radius, connected_particles);
        simulation.connect(A4,C4,A4.radius + C4.radius, connected_particles);
       // simulation.connect(C4,B4,C4.radius + B4.radius, connected_particles);
        simulation.connect(B4,D4,D4.radius + B4.radius, connected_particles);
       // simulation.connect(D4,A4,D4.radius + A4.radius, connected_particles);
       simulation.connect(D4,C4,D4.radius + C4.radius, connected_particles);
       simulation.connect(B4,C4,std::sqrt(std::pow(B4.radius+D4.radius,2) + std::pow(C4.radius + D4.radius,2)),connected_particles);
        simulation.connect(A4,D4,std::sqrt(std::pow(C4.radius+D4.radius,2) + std::pow(C4.radius + A4.radius,2)),connected_particles);


        simulation.connect(BIG_A,SMALL_B,BIG_A.radius +SMALL_B.radius,connected_particles);



        simulation.connect(B_blitz,C_blitz,C_blitz.radius+B_blitz.radius,connected_particles);
        simulation.connect(A_blitz,C_blitz,std::sqrt(std::pow(A_blitz.radius,2) + std::pow(A_blitz.radius + 2 * B_blitz.radius + C_blitz.radius,2)),connected_particles);
        simulation.connect(C_blitz,D_blitz,std::sqrt(std::pow(C_blitz.radius,2) + std::pow(C_blitz.radius  + D_blitz.radius,2)),connected_particles);
        simulation.connect(A_blitz,B_blitz,std::sqrt(std::pow(B_blitz.radius,2) + std::pow(A_blitz.radius +   B_blitz.radius ,2)),connected_particles);
        simulation.connect(D_blitz,B_blitz,std::sqrt(std::pow(D_blitz.radius,2) + std::pow(D_blitz.radius + 2 * C_blitz.radius + B_blitz.radius,2)),connected_particles);
        simulation.connect(A_blitz,D_blitz,std::sqrt( std::pow(D_blitz.radius+A_blitz.radius,2) + std::pow(D_blitz.radius + 2*C_blitz.radius + 2* B_blitz.radius + A_blitz.radius,2) ), connected_particles);

*/

        rect = Box{{simulation.Boxcenter[0], simulation.Boxcenter[1]}, Vector2d(simulation.Boxsize[0], simulation.Boxsize[1]), nvgRGBA(50, 50, 50, 0), nvgRGBA(50, 50, 50, 200)};






        simulation.assign_index(connected_particles);
        simulation.initialize_position_vector_minus_1(connected_particles);

    }





    void process() override {
        my_data.open("data.csv");


        if(runSim) {


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
        my_data.close();
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
                        std::string permittivity_label = label;
                        std::string permittivity = " permittivity";
                        permittivity_label += permittivity;
                        const char* permittivity_label_ = permittivity_label.c_str();
                        if(SliderScalar(permittivity_label_,ImGuiDataType_Double,&A.permittivity,&min_permittivity,&max_permittivity)){

                        }
                        std::string conductivity_label = label;
                        std::string conductivity = " conductivity";
                        conductivity_label += conductivity;
                        const char* conductivity_label_ = conductivity_label.c_str();
                        if(SliderScalar(conductivity_label_,ImGuiDataType_Double,&A.conductivity,&min_conductivity,&max_conductivity)){

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
                        if(Checkbox(silica_label_,&is_silica_1)){
                            A.permittivity = 3.9;
                            A.density = 2650.0;
                            A.radius = 3.0 * mscale;
                            A.mass = M_PI * std::pow(A.radius,2) * A.density;
                            A.conductivity = 0.0;
                        }
                        std::string PS_label = label;
                        std::string PS = " PS";
                        PS_label+=PS;
                        const char* PS_label_ = PS_label.c_str();
                        if(Checkbox(PS_label_, &is_PS_1)){
                            A.permittivity = 2.5;
                            A.density = 1050;
                            A.radius = 3.0 * mscale;
                            A.mass = M_PI * std::pow(A.radius,2) * A.density;
                            A.conductivity = 5 * std::pow(10,-9);

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
                        std::string permittivity_label = label;
                        std::string permittivity = " permittivity";
                        permittivity_label += permittivity;
                        const char* permittivity_label_ = permittivity_label.c_str();
                        if(SliderScalar(permittivity_label_,ImGuiDataType_Double,&B.permittivity,&min_permittivity,&max_permittivity)){

                        }
                        std::string conductivity_label = label;
                        std::string conductivity = " conductivity";
                        conductivity_label += conductivity;
                        const char* conductivity_label_ = conductivity_label.c_str();
                        if(SliderScalar(conductivity_label_,ImGuiDataType_Double,&B.conductivity,&min_conductivity,&max_conductivity)){

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
                        if(Checkbox(silica_label_,&is_silica_2)){
                            B.permittivity = 3.9;
                            B.density = 2650.0;
                            B.radius = 2.0 * mscale;
                            B.mass = M_PI * std::pow(B.radius,2) * B.density;
                            B.conductivity = 0.0;
                        }
                        std::string PS_label = label;
                        std::string PS = " PS";
                        PS_label+=PS;
                        const char* PS_label_ = PS_label.c_str();


                        if(Checkbox(PS_label_, &is_PS_2)){
                            B.permittivity = 2.5;
                            B.density = 1050;
                            B.radius = 2.0 * mscale;
                            B.mass = M_PI * std::pow(A.radius,2) * A.density;
                            B.conductivity = 5 * std::pow(10,-9);

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

        BeginMainMenuBar();
        if(BeginMenu("Control")) {
            if (CollapsingHeader("Springs")) {
                if (SliderScalar("Stiffnes Constant", ImGuiDataType_Double, &simulation.stiffnes_constant,
                                 &min_stiffnes_constant, &max_stiffnes_constant))
                    std::cout << simulation.stiffnes_constant;
            }

            if (CollapsingHeader("Drag")) {
                if (SliderScalar("Viscosity Multiplier", ImGuiDataType_Double, &simulation.viscosity_multiplier,
                                 &min_viscosity_multiplier, &max_viscosity_multiplier))
                    std::cout << simulation.viscosity_multiplier;
            }

            if(CollapsingHeader("Experiments")){
                if(Checkbox(" Silica - Silica Experiment", &first_experiment)){
                    simulation.beta = 1.0;
                    simulation.lower_electrode.peak_voltage = 6.5;
                    simulation.upper_electrode.position[1] = simulation.lower_electrode.position[1] + 100.0 * mscale;
                }
                if(Checkbox(" PS - PS Experiment", &second_experiment)){
                    simulation.beta = 1.0;
                    simulation.lower_electrode.frequency = 500.0;
                    simulation.upper_electrode.position[1] = simulation.lower_electrode.position[1] + 100.0 * mscale;
                }
            }

            if (CollapsingHeader("Brownian Motion")) {
                if (SliderScalar("Brownian Multiplier", ImGuiDataType_Double, &simulation.brownian_motion_multiplier,
                                 &min_brownian_motion_multiplier, &max_brownian_motion_multiplier))
                    std::cout << simulation.brownian_motion_multiplier;
            }
            if (CollapsingHeader("Simulation")) {
                if (SliderScalar("Time Step", ImGuiDataType_Double, &simulation.time_step, &min_time_step,
                                 &max_time_step)) {
                    std::cout << simulation.time_step;
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
            auto drawRectangle = [this](const Rectangle &rectangle, bool is_lower_electrode,double scale, Vector2d distance){

                nvgBeginPath(vg);

                nvgRect(vg,rectangle.position[0] + distance[0] * scale + center_of_frame[0], rectangle.position[1]+ distance[1] * scale + center_of_frame[1], rectangle.width,rectangle.height);

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

            };
            Vector2d formation_center= (simulation.lower_electrode.position + simulation.upper_electrode.position) * 0.5;
            Vector2d distance_1 = formation_center - simulation.lower_electrode.position;
            Vector2d distance_2 = formation_center - simulation.upper_electrode.position;
            drawRectangle(Rectangle(simulation.lower_electrode.position, simulation.lower_electrode.width,simulation.lower_electrode.length),true, (1.0/mscale),distance_1);
            drawRectangle(Rectangle(simulation.upper_electrode.position, simulation.upper_electrode.width,simulation.upper_electrode.length),false,(1.0/mscale),distance_2);
        }

        if(drawCircles)
        {
            auto drawCircle = [this](const Circle &circle,int tmpfortest, double scale, Vector2d distance){
                nvgBeginPath(vg);
                nvgCircle(vg, (circle.pos[0]  + distance[0] * scale) + center_of_frame[0], (circle.pos[1] + distance[1] * scale) + center_of_frame[1], circle.radius *scale);
                if(tmpfortest ==1) {
                    nvgFillColor(vg, nvgRGBA(150, 0, 0, 100));
                }else{
                    nvgFillColor(vg, nvgRGBA(0, 0, 150, 100));
                }
                nvgFill(vg);
                nvgStrokeColor(vg, circle.colorStroke);
                nvgStrokeWidth(vg, 10.0f);
                nvgStroke(vg);
            };

            /*for(const auto &p : connect_vector) //change to particles for std sim
                drawCircle(Circle(p.position,p.radius));*/

            Vector2d formation_center;
            formation_center.setZero();
            double objects_total = simulation.size;
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
            formation_center[0] = formation_center[0]/objects_total;
            formation_center[1] = formation_center[1]/objects_total;
            simulation.reset_flags(connected_particles);



            for(auto &pair : connected_particles) { //change to particles for std sim





               if(!std::get<0>(pair).visited) {
                   Vector2d distance = formation_center - std::get<0>(pair).position;
                   drawCircle(Circle(std::get<0>(pair).position, std::get<0>(pair).radius), 1, (1.0/mscale), distance);
                   std::get<0>(pair).visited = true;
               }
               if(!std::get<1>(pair).visited) {
                   Vector2d distance = formation_center - std::get<1>(pair).position;
                   drawCircle(Circle(std::get<1>(pair).position, std::get<1>(pair).radius), 2, (1.0/mscale), distance);
                   std::get<1>(pair).visited = true;
               }

            }
            simulation.reset_flags(connected_particles);
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
    Vector2d center_of_frame = {1000,1000};
    bool runSim = false;
    double cursorPosDown[2]{};
    double translation[2] = {0.75*pixelRatio, -0.25*pixelRatio};
    double zoom = 10;
    int base;

    //Particle BIG_A = Particle(100,40, {520,-100});
    //Particle SMALL_B = Particle(10,15,{500,-90});

    double mscale = std::pow(10,-6);


    Particle A = Particle(50.0 *mscale ,20.0 * mscale,1.0,{100.0 * mscale,100.0 * mscale});
    Particle B = Particle(50.0 * mscale,20.0* mscale,1.0,{140.0 * mscale,100.0 * mscale});
    //Particle C = Particle(50.0,20.0,1.0,{150.0,0.0});

  /*  Particle A2 = Particle(20,20,{100,100});
    Particle B2 = Particle(20,20,{0,0});
    Particle C2 = Particle(20,20,{-100,-100});

    Particle A3 = Particle(20,20,{-200,-200});
    Particle B3 = Particle(20,20,{-300,-300});
    Particle C3 = Particle(20,20,{-400,-400});
    Particle D3 = Particle(20,20,{-400,-400});


    Particle A4 = Particle(20,20,{200,-200});
    Particle B4 = Particle(20,20,{300,-300});
    Particle C4 = Particle(20,20,{400,-400});
    Particle D4 = Particle(20,20,{400,-400});


    Particle A_blitz = Particle(20,20,{100,-200});
    Particle B_blitz = Particle(20,20,{200,-300});
    Particle C_blitz = Particle(20,20,{300,-400});
    Particle D_blitz = Particle(20,20,{300,-400});


    Particle A6 = Particle(20,20,{-200,200});
    Particle B6 = Particle(20,20,{-300,300});
    Particle C6 = Particle(20,20,{-400,400});
    Particle D6 = Particle(20,20,{-400,400});
    Particle E6 = Particle(20,20,{-500,400});
    Particle F6 = Particle(20,20,{-400,200});*/

    std::vector<Particle> connect_vector;

    /*Particle particle_1;
    Particle particle_2 = Particle(10,10,{400,-400});
    Particle particle_3 = Particle(10,20,{-200,-200});
    Particle particle_4 = Particle(10,10,{500,200});
    Particle particle_5 = Particle(10,10,{200,500});
    Particle particle_6 = Particle(10,10,{10,20});
    Particle particle_7 = Particle(10,10,{0,-300});
    Particle particle_8 = Particle(10,10,{50,50});
    Particle particle_9 = Particle(10,10,{470,-100});
    Particle particle_10 = Particle(10,10,{-600,300});

    Particle loner = Particle(10,30,{-600,-600});*/
    std::vector<Particle> particles;
    std::vector<std::tuple <Particle&,Particle&,double> > connected_particles;
    Simulation simulation;
    std::ofstream my_data;






    bool drawCursor = false;
    bool drawCircles = true;
    bool drawRectangle = true;

    double min_stiffnes_constant = 0.0;
    double max_stiffnes_constant = 100000.0;
    double min_viscosity_multiplier = 0.0;
    double max_viscosity_multiplier = 1000.0;
    double min_brownian_motion_multiplier = 0.0;
    double max_brownian_motion_multiplier = 10.0;
    double min_time_step = 0.0* mscale;
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

    double min_lower_electrode_voltage = 0.0;
    double max_lower_electrode_voltage = 300.0;
    double min_upper_electrode_voltage = -300.0;
    double max_upper_electrode_voltage = 0.0;
    double min_mass = 0.1 * mscale;
    double max_mass = 100 * mscale;
    double min_charge = -30;
    double max_charge = 30.0;
    double min_radius = 1.0 * mscale;
    double max_radius = 200.0 * mscale;
    double min_lower_electrode_frequency = 0.0;
    double max_lower_electrode_frequency = 12000.0;
    double min_lower_electrode_peak_voltage = 0.0;
    double max_lower_electrode_peak_voltage = 1000.0;
    double min_permittivity = 0.0;
    double max_permittivity = 200.0;
    double min_conductivity = 0.0;
    double max_conductivity = 1.0 * std::pow(10,-3);

    bool is_silica_1 = false;
    bool is_silica_2 = false;
    bool first_experiment = false;
    bool second_experiment = false;
    bool is_PS_1 = false;
    bool is_PS_2 = false;





    struct Rectangle
    {
        Rectangle(Vector2d p, double w, double h){
            this->position = p;
            this->width = w;
            this->height = h;
        }
        Vector2d position;
        double height;
        double width;
        NVGcolor colorFill = nvgRGBA(50, 50, 50, 100) , colorStroke = nvgRGBA(50, 50, 50, 100);
    };

    struct Circle
    {
        Circle(Vector2d p, double r=10){ //took away struct
            this->pos = p;
            this->radius = r;
        }
        bool isInside(const Vector2d &x){
            return (x-pos).squaredNorm() <= radius*radius;
        }

        Vector2d pos;
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

    std::ofstream myfile;
    myfile.open ("../../../data.csv");
    myfile << "Total Energy"<<","<<"Kinetic Energy"<<","<<"Potential Energy"<<","<<"E"<<","<<"SPRING FORCE X"<<","<<"SPRING FORCE Y"<<","<<"DISTANCE"<<","<<"DERIVATIVE"<<","<<"X"<<","<<"Y"<<","<<"time step"<<std::endl;
    std::cout<<"test"<<std::endl;
    for(int i = 0; i<app.simulation.total_energy.size(); i++){
        //std::cout<<app.simulation.spring_force_vec_over_time_x[i]<<std::endl;
        //std::cout<<app.simulation.spring_force_vec_over_time_x[i]/app.simulation.stiffnes_constant<<std::endl;
        myfile << app.simulation.total_energy[i]<<","<<app.simulation.kinetic_energy[i]<<","<<app.simulation.potential_energy[i]<<","<<app.simulation.e_vec[i]<<","<<app.simulation.spring_force_vec_over_time_x[i]/app.simulation.stiffnes_constant<<","<<app.simulation.spring_force_vec_over_time_y[i]/app.simulation.stiffnes_constant<<","<<app.simulation.A_B_DISTANCE[i]<<","<<app.simulation.spring_force_derivative_x_in_x[i]/app.simulation.stiffnes_constant<<","<<app.simulation.position_vec_over_time_in_x[i]<<","<<app.simulation.position_vec_over_time_in_y[i]<<","<<i<<std::endl;
    }

    std::ofstream v_w;
    v_w.open("../../../v_w.csv");
    v_w <<" Velocity"<<","<<"Frequency"<<std::endl;
    for(int i = 0; i<app.simulation.frequencies.size();i++){
        v_w<<app.simulation.velocities[i]<<","<<app.simulation.frequencies[i]<<std::endl;
    }

    std::ofstream v_E;
    v_E.open("../../../v_E.csv");
    v_E <<" Velocity"<<","<<"E Field"<<std::endl;
    for(int i = 0; i<app.simulation.Efields.size();i++){
        v_E<<app.simulation.velocities_2[i]<<","<<app.simulation.Efields[i]<<std::endl;
    }

    std::ofstream velocities;
    velocities.open ("../../../velocities.csv");
    std::cout<<"test2"<<std::endl;

    velocities << "VELX"<<","<<"VELY"<<","<<"friction"<<","<<"time_step"<<std::endl;
    for(int i = 0; i<app.simulation.velocities_over_time1_in_x.size();i++){
        velocities<<app.simulation.velocities_over_time1_in_x[i]<<","<<app.simulation.velocities_over_time1_in_y[i]<<","<<app.simulation.friction_force_over_time_x[i]<<","<<i<<std::endl;
        std::cout<<app.simulation.velocities_over_time1_in_x[i]<<std::endl;
    }
    velocities.close();

    std::ofstream myfile2;
    myfile2.open("../../../data2.csv");
    myfile2<<"x"<<","<<"y"<<","<<"z1"<<","<<"z2"<<std::endl;
    for(int i = 0; i<899;i++){
        myfile2<<app.simulation.x_values[i]<<","<<app.simulation.y_values[i]<<","<<app.simulation.z_values1[i]<<","<<app.simulation.z_values2[i]<<std::endl;
        //std::cout<<app.simulation.x_values[i]<<std::endl;
       // std::cout<<app.simulation.y_values[i]<<std::endl;
        //std::cout<<app.simulation.z_values1[i]<<std::endl;
    }

    std::ofstream myfile3;
    myfile3.open("../../../data3.csv");
    for(int i = 0; i<app.simulation.x_values.size();i++){
        myfile3<<app.simulation.x_values[i]<<","<<app.simulation.y_values[i]<<","<<app.simulation.z_values[i]<<std::endl;
        //std::cout<<app.simulation.x_values[i]<<std::endl;
        //std::cout<<app.simulation.y_values[i]<<std::endl;
        //std::cout<<app.simulation.z_values1[i]<<std::endl;
    }

    std::ofstream best;
    best.open("../../../best.csv");
    best<<"x"<<","<<"y"<<","<<"z"<<std::endl;
    best<<app.simulation.best_found[0]<<","<<app.simulation.best_found[1]<<","<<app.simulation.best_found[2]<<std::endl;
    myfile.close();
    v_w.close();
    v_E.close();
    myfile2.close();
    myfile3.close();
    best.close();

    std::ofstream FD;
    FD.open("../../../FD.csv");
    FD<<"F_X"<<","<<"DEDX"<<","<<"Difference"<<","<<"FORCE_JACOBIAN"<<","<<"DFDX"<<","<<"DIFFERENCE2"<<","<<"i"<<std::endl;
    for(int i = 0; i<app.simulation.forces_FD.size();i++){
        FD<<app.simulation.forces_FD[i]<<","<<app.simulation.energy_FD[i]<<","<<std::abs(app.simulation.forces_FD[i] - app.simulation.energy_FD[i])<<","<<app.simulation.force_jacobian_FD[i]<<","<<app.simulation.dFx_dx_FD[i]<<","<<std::abs(app.simulation.force_jacobian_FD[i] - app.simulation.dFx_dx_FD[i])<<","<<i<<std::endl;
    }

    std::ofstream test_e;
    test_e.open("../../../test_e.csv");
    test_e<<"ex1"<<","<<"ex2"<<","<<"value"<<std::endl;
    for(int i = 0; i<app.simulation.e_z.size(); i++){
        test_e<<app.simulation.e_x[i]<<","<<app.simulation.e_y[i]<<","<<app.simulation.e_z[i]<<std::endl;
        //std::cout<<app.simulation.e_z[i];
    }



    test_e.close();
    return 0;
}
