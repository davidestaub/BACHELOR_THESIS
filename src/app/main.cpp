#include "application.h"
#include <imgui.h>
#include "simulation.h"

#include <iostream>
#include <math.h>
#include <deque>
#include <chrono>
#include <vector>

#include <Eigen/Core>
using Eigen::Vector2d;

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
       simulation.connect(B,C,(C.radius+B.radius),connected_particles);
        simulation.connect(A,C,(A.radius+C.radius),connected_particles);

        /*simulation.connect(A2,B2,(A2.radius+B2.radius),connected_particles);
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








    }



    void process() override {

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
        if(BeginMenu("debug")){
            Checkbox("draw cursor", &drawCursor);
            Checkbox("draw circles", &drawCircles);
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

        if(drawCircles)
        {
            auto drawCircle = [this](const Circle &circle){
                nvgBeginPath(vg);
                nvgCircle(vg, circle.pos[0] + center_of_frame[0], circle.pos[1]+ center_of_frame[1], circle.radius);
                nvgFillColor(vg, circle.colorFill);
                nvgFill(vg);
                nvgStrokeColor(vg, circle.colorStroke);
                nvgStrokeWidth(vg, 10.0f);
                nvgStroke(vg);
            };

            /*for(const auto &p : connect_vector) //change to particles for std sim
                drawCircle(Circle(p.position,p.radius));*/
            for(auto &pair : connected_particles) { //change to particles for std sim
               // if(!pair.first.is_drawn) {
                    drawCircle(Circle(std::get<0>(pair).position, std::get<0>(pair).radius));
                 //   pair.first.is_drawn = true;
               // }
               // if(!pair.second.is_drawn) {
                    drawCircle(Circle(std::get<1>(pair).position, std::get<1>(pair).radius));
                //    pair.second.is_drawn = true;
               // }

            }
            for(auto &particle : particles){
                drawCircle(Circle(particle.position, particle.radius));
            }

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

    void scrollWheel(double  /*xoffset*/, double yoffset) override {
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
            simulation.reset_simulation(connect_vector);
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

private:
    Vector2d center_of_frame = {1000,1000};
    bool runSim = false;
    double cursorPosDown[2]{};
    double translation[2] = {0.75*pixelRatio, -0.25*pixelRatio};
    double zoom = 24;
    int base;

    Particle BIG_A = Particle(100,40, {520,-100});
    Particle SMALL_B = Particle(10,15,{500,-90});


    Particle A = Particle(20,20,{200,200});
    Particle B = Particle(20,20,{300,300});
    Particle C = Particle(20,20,{400,400});

    Particle A2 = Particle(20,20,{100,100});
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
    Particle F6 = Particle(20,20,{-400,200});

    std::vector<Particle> connect_vector;

    Particle particle_1;
    Particle particle_2 = Particle(10,10,{400,-400});
    Particle particle_3 = Particle(10,20,{-200,-200});
    Particle particle_4 = Particle(10,10,{500,200});
    Particle particle_5 = Particle(10,10,{200,500});
    Particle particle_6 = Particle(10,10,{10,20});
    Particle particle_7 = Particle(10,10,{0,-300});
    Particle particle_8 = Particle(10,10,{50,50});
    Particle particle_9 = Particle(10,10,{470,-100});
    Particle particle_10 = Particle(10,10,{-600,300});

    Particle loner = Particle(10,30,{-600,-600});
    std::vector<Particle> particles;
    std::vector<std::tuple <Particle&,Particle&,double> > connected_particles;
    Simulation simulation;






    bool drawCursor = false;
    bool drawCircles = true;

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

    return 0;
}
