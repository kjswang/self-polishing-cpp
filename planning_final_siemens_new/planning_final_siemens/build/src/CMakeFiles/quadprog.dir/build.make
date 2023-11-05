# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build

# Include any dependencies generated for this target.
include src/CMakeFiles/quadprog.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/quadprog.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/quadprog.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/quadprog.dir/flags.make

src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o: src/CMakeFiles/quadprog.dir/flags.make
src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o: /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/src/lib/QP_lib/QuadProg++.cc
src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o: src/CMakeFiles/quadprog.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o"
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o -MF CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o.d -o CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o -c /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/src/lib/QP_lib/QuadProg++.cc

src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.i"
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/src/lib/QP_lib/QuadProg++.cc > CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.i

src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.s"
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/src/lib/QP_lib/QuadProg++.cc -o CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.s

# Object files for target quadprog
quadprog_OBJECTS = \
"CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o"

# External object files for target quadprog
quadprog_EXTERNAL_OBJECTS =

src/libquadprog.a: src/CMakeFiles/quadprog.dir/lib/QP_lib/QuadProg++.cc.o
src/libquadprog.a: src/CMakeFiles/quadprog.dir/build.make
src/libquadprog.a: src/CMakeFiles/quadprog.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libquadprog.a"
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src && $(CMAKE_COMMAND) -P CMakeFiles/quadprog.dir/cmake_clean_target.cmake
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quadprog.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/quadprog.dir/build: src/libquadprog.a
.PHONY : src/CMakeFiles/quadprog.dir/build

src/CMakeFiles/quadprog.dir/clean:
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src && $(CMAKE_COMMAND) -P CMakeFiles/quadprog.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/quadprog.dir/clean

src/CMakeFiles/quadprog.dir/depend:
	cd /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/src /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src /home/jiansong/HTPC-1/planning_final_siemens_new/planning_final_siemens/build/src/CMakeFiles/quadprog.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/CMakeFiles/quadprog.dir/depend

