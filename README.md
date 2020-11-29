
1. Clone this repository and load submodules:
    ```
    git clone --recurse-submodules YOUR_GIT_URL
    ```
2. Run CMake to create project files
    ```
    cd NAME_OF_YOUR_REPO
    mkdir build && cd build
    cmake ..
    ```
3. Compile code and run executable
    - for MacOS and Linux:
        ```
        make
        ./src/app/main.cpp
        ```

INSTRUCTIONS:
    A default situation is given in the main.cpp file (a PNiPAM microgel - Polystyrene dumbbell).
    If you want to change this to simulate something else (e.g the L-shaped robot or the rotator etc.) you need to uncomment the correct lines in the code!
    These lines start at line 37, line 94 and line 1045. It works the following. Say you want to simulate not 1 but 2 dumbbells.:
    1. uncomment the correct line (they are labeles with comments) around line 1045. In this case uncomment a line to create another dumbbell.
    Next go to line 37 and uncomment the corresponding line(s). In this case uncomment one line corresponding to the dumbell that you want to run (again they are labeled with comments).1
    Finally uncomment the corresponding line around line 94. If anything is unclear do not hesitate to write me an email (PS. all the provided scenarios work so if something does not work it probably means that a wrong line was uncommented)


    Now run the simulation as described above. At first nothing will move, you have to do several things:

    1. Select the type of experiment you want to run under the tab Experiments (I suggest to take the third one)
    2. Select the materials of your particles under the tab particles (I suggest PS for the first and microgel for the second)
    3. Hit space to run the simulation.
    4. Under the tab 'debug' select the box 'safety off'(I implemented this for more complex robots to be 100% sure that all the springs are at rest length before starting the experiment to nut falsify the data)
    5. Enjoy :) (you can play around with all the checkboxes and sliders and).

    To restart the simulation press R.

    Note there is one thing you should not do and that is to set the box 'dimer tracer on' to true if you run an experiment with a robots that is not a dimer (dumbbell) this will crash the app. Just use the normal tracer checkbox.

    There are a lot of checkboxes and sliders, I think they should be self explanatory. If any question or problem appears feel free to contact me under dastaub@student.ethz.ch

 
