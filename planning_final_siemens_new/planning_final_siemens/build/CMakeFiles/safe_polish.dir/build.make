# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jiansong/planning_final_siemens/planning_final_siemens

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jiansong/planning_final_siemens/planning_final_siemens/build

# Include any dependencies generated for this target.
include CMakeFiles/safe_polish.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/safe_polish.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/safe_polish.dir/flags.make

CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.o: CMakeFiles/safe_polish.dir/flags.make
CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.o: ../src/safe_polish_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jiansong/planning_final_siemens/planning_final_siemens/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.o -c /home/jiansong/planning_final_siemens/planning_final_siemens/src/safe_polish_main.cpp

CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiansong/planning_final_siemens/planning_final_siemens/src/safe_polish_main.cpp > CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.i

CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiansong/planning_final_siemens/planning_final_siemens/src/safe_polish_main.cpp -o CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.s

# Object files for target safe_polish
safe_polish_OBJECTS = \
"CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.o"

# External object files for target safe_polish
safe_polish_EXTERNAL_OBJECTS =

safe_polish: CMakeFiles/safe_polish.dir/src/safe_polish_main.cpp.o
safe_polish: CMakeFiles/safe_polish.dir/build.make
safe_polish: src/libquadprog.a
safe_polish: src/libyaml-cpp.a
safe_polish: CMakeFiles/safe_polish.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jiansong/planning_final_siemens/planning_final_siemens/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable safe_polish"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/safe_polish.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/safe_polish.dir/build: safe_polish

.PHONY : CMakeFiles/safe_polish.dir/build

CMakeFiles/safe_polish.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/safe_polish.dir/cmake_clean.cmake
.PHONY : CMakeFiles/safe_polish.dir/clean

CMakeFiles/safe_polish.dir/depend:
	cd /home/jiansong/planning_final_siemens/planning_final_siemens/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jiansong/planning_final_siemens/planning_final_siemens /home/jiansong/planning_final_siemens/planning_final_siemens /home/jiansong/planning_final_siemens/planning_final_siemens/build /home/jiansong/planning_final_siemens/planning_final_siemens/build /home/jiansong/planning_final_siemens/planning_final_siemens/build/CMakeFiles/safe_polish.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/safe_polish.dir/depend

