# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug

# Include any dependencies generated for this target.
include src/simulation/CMakeFiles/simulation.dir/depend.make

# Include the progress variables for this target.
include src/simulation/CMakeFiles/simulation.dir/progress.make

# Include the compile flags for this target's objects.
include src/simulation/CMakeFiles/simulation.dir/flags.make

src/simulation/CMakeFiles/simulation.dir/particles.cpp.o: src/simulation/CMakeFiles/simulation.dir/flags.make
src/simulation/CMakeFiles/simulation.dir/particles.cpp.o: ../src/simulation/particles.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/simulation/CMakeFiles/simulation.dir/particles.cpp.o"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simulation.dir/particles.cpp.o -c /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation/particles.cpp

src/simulation/CMakeFiles/simulation.dir/particles.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulation.dir/particles.cpp.i"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation/particles.cpp > CMakeFiles/simulation.dir/particles.cpp.i

src/simulation/CMakeFiles/simulation.dir/particles.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulation.dir/particles.cpp.s"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation/particles.cpp -o CMakeFiles/simulation.dir/particles.cpp.s

src/simulation/CMakeFiles/simulation.dir/simulation.cpp.o: src/simulation/CMakeFiles/simulation.dir/flags.make
src/simulation/CMakeFiles/simulation.dir/simulation.cpp.o: ../src/simulation/simulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/simulation/CMakeFiles/simulation.dir/simulation.cpp.o"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simulation.dir/simulation.cpp.o -c /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation/simulation.cpp

src/simulation/CMakeFiles/simulation.dir/simulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulation.dir/simulation.cpp.i"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation/simulation.cpp > CMakeFiles/simulation.dir/simulation.cpp.i

src/simulation/CMakeFiles/simulation.dir/simulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulation.dir/simulation.cpp.s"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation/simulation.cpp -o CMakeFiles/simulation.dir/simulation.cpp.s

# Object files for target simulation
simulation_OBJECTS = \
"CMakeFiles/simulation.dir/particles.cpp.o" \
"CMakeFiles/simulation.dir/simulation.cpp.o"

# External object files for target simulation
simulation_EXTERNAL_OBJECTS =

src/simulation/libsimulation.a: src/simulation/CMakeFiles/simulation.dir/particles.cpp.o
src/simulation/libsimulation.a: src/simulation/CMakeFiles/simulation.dir/simulation.cpp.o
src/simulation/libsimulation.a: src/simulation/CMakeFiles/simulation.dir/build.make
src/simulation/libsimulation.a: src/simulation/CMakeFiles/simulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libsimulation.a"
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && $(CMAKE_COMMAND) -P CMakeFiles/simulation.dir/cmake_clean_target.cmake
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/simulation/CMakeFiles/simulation.dir/build: src/simulation/libsimulation.a

.PHONY : src/simulation/CMakeFiles/simulation.dir/build

src/simulation/CMakeFiles/simulation.dir/clean:
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation && $(CMAKE_COMMAND) -P CMakeFiles/simulation.dir/cmake_clean.cmake
.PHONY : src/simulation/CMakeFiles/simulation.dir/clean

src/simulation/CMakeFiles/simulation.dir/depend:
	cd /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/src/simulation /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation /Users/davidestaub/Desktop/BACHELOR_THESIS/BACHELOR_THESIS/cmake-build-debug/src/simulation/CMakeFiles/simulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/simulation/CMakeFiles/simulation.dir/depend

