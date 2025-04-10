# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/hellbats/Code/Phase-Transition-Model

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/hellbats/Code/Phase-Transition-Model/build

# Include any dependencies generated for this target.
include CMakeFiles/build.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/build.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/build.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/build.dir/flags.make

CMakeFiles/build.dir/src/Lattice.c.o: CMakeFiles/build.dir/flags.make
CMakeFiles/build.dir/src/Lattice.c.o: /home/hellbats/Code/Phase-Transition-Model/src/Lattice.c
CMakeFiles/build.dir/src/Lattice.c.o: CMakeFiles/build.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/build.dir/src/Lattice.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/build.dir/src/Lattice.c.o -MF CMakeFiles/build.dir/src/Lattice.c.o.d -o CMakeFiles/build.dir/src/Lattice.c.o -c /home/hellbats/Code/Phase-Transition-Model/src/Lattice.c

CMakeFiles/build.dir/src/Lattice.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/build.dir/src/Lattice.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/hellbats/Code/Phase-Transition-Model/src/Lattice.c > CMakeFiles/build.dir/src/Lattice.c.i

CMakeFiles/build.dir/src/Lattice.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/build.dir/src/Lattice.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/hellbats/Code/Phase-Transition-Model/src/Lattice.c -o CMakeFiles/build.dir/src/Lattice.c.s

CMakeFiles/build.dir/src/Parameters.c.o: CMakeFiles/build.dir/flags.make
CMakeFiles/build.dir/src/Parameters.c.o: /home/hellbats/Code/Phase-Transition-Model/src/Parameters.c
CMakeFiles/build.dir/src/Parameters.c.o: CMakeFiles/build.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/build.dir/src/Parameters.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/build.dir/src/Parameters.c.o -MF CMakeFiles/build.dir/src/Parameters.c.o.d -o CMakeFiles/build.dir/src/Parameters.c.o -c /home/hellbats/Code/Phase-Transition-Model/src/Parameters.c

CMakeFiles/build.dir/src/Parameters.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/build.dir/src/Parameters.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/hellbats/Code/Phase-Transition-Model/src/Parameters.c > CMakeFiles/build.dir/src/Parameters.c.i

CMakeFiles/build.dir/src/Parameters.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/build.dir/src/Parameters.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/hellbats/Code/Phase-Transition-Model/src/Parameters.c -o CMakeFiles/build.dir/src/Parameters.c.s

CMakeFiles/build.dir/src/Random.c.o: CMakeFiles/build.dir/flags.make
CMakeFiles/build.dir/src/Random.c.o: /home/hellbats/Code/Phase-Transition-Model/src/Random.c
CMakeFiles/build.dir/src/Random.c.o: CMakeFiles/build.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/build.dir/src/Random.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/build.dir/src/Random.c.o -MF CMakeFiles/build.dir/src/Random.c.o.d -o CMakeFiles/build.dir/src/Random.c.o -c /home/hellbats/Code/Phase-Transition-Model/src/Random.c

CMakeFiles/build.dir/src/Random.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/build.dir/src/Random.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/hellbats/Code/Phase-Transition-Model/src/Random.c > CMakeFiles/build.dir/src/Random.c.i

CMakeFiles/build.dir/src/Random.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/build.dir/src/Random.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/hellbats/Code/Phase-Transition-Model/src/Random.c -o CMakeFiles/build.dir/src/Random.c.s

CMakeFiles/build.dir/src/Sim.c.o: CMakeFiles/build.dir/flags.make
CMakeFiles/build.dir/src/Sim.c.o: /home/hellbats/Code/Phase-Transition-Model/src/Sim.c
CMakeFiles/build.dir/src/Sim.c.o: CMakeFiles/build.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/build.dir/src/Sim.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/build.dir/src/Sim.c.o -MF CMakeFiles/build.dir/src/Sim.c.o.d -o CMakeFiles/build.dir/src/Sim.c.o -c /home/hellbats/Code/Phase-Transition-Model/src/Sim.c

CMakeFiles/build.dir/src/Sim.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/build.dir/src/Sim.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/hellbats/Code/Phase-Transition-Model/src/Sim.c > CMakeFiles/build.dir/src/Sim.c.i

CMakeFiles/build.dir/src/Sim.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/build.dir/src/Sim.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/hellbats/Code/Phase-Transition-Model/src/Sim.c -o CMakeFiles/build.dir/src/Sim.c.s

CMakeFiles/build.dir/src/UnitCellTemplate.c.o: CMakeFiles/build.dir/flags.make
CMakeFiles/build.dir/src/UnitCellTemplate.c.o: /home/hellbats/Code/Phase-Transition-Model/src/UnitCellTemplate.c
CMakeFiles/build.dir/src/UnitCellTemplate.c.o: CMakeFiles/build.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/build.dir/src/UnitCellTemplate.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/build.dir/src/UnitCellTemplate.c.o -MF CMakeFiles/build.dir/src/UnitCellTemplate.c.o.d -o CMakeFiles/build.dir/src/UnitCellTemplate.c.o -c /home/hellbats/Code/Phase-Transition-Model/src/UnitCellTemplate.c

CMakeFiles/build.dir/src/UnitCellTemplate.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/build.dir/src/UnitCellTemplate.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/hellbats/Code/Phase-Transition-Model/src/UnitCellTemplate.c > CMakeFiles/build.dir/src/UnitCellTemplate.c.i

CMakeFiles/build.dir/src/UnitCellTemplate.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/build.dir/src/UnitCellTemplate.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/hellbats/Code/Phase-Transition-Model/src/UnitCellTemplate.c -o CMakeFiles/build.dir/src/UnitCellTemplate.c.s

# Object files for target build
build_OBJECTS = \
"CMakeFiles/build.dir/src/Lattice.c.o" \
"CMakeFiles/build.dir/src/Parameters.c.o" \
"CMakeFiles/build.dir/src/Random.c.o" \
"CMakeFiles/build.dir/src/Sim.c.o" \
"CMakeFiles/build.dir/src/UnitCellTemplate.c.o"

# External object files for target build
build_EXTERNAL_OBJECTS =

build: CMakeFiles/build.dir/src/Lattice.c.o
build: CMakeFiles/build.dir/src/Parameters.c.o
build: CMakeFiles/build.dir/src/Random.c.o
build: CMakeFiles/build.dir/src/Sim.c.o
build: CMakeFiles/build.dir/src/UnitCellTemplate.c.o
build: CMakeFiles/build.dir/build.make
build: /usr/lib/gcc/x86_64-linux-gnu/14/libgomp.so
build: /usr/lib/x86_64-linux-gnu/libpthread.a
build: CMakeFiles/build.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C executable build"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/build.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/build.dir/build: build
.PHONY : CMakeFiles/build.dir/build

CMakeFiles/build.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/build.dir/cmake_clean.cmake
.PHONY : CMakeFiles/build.dir/clean

CMakeFiles/build.dir/depend:
	cd /home/hellbats/Code/Phase-Transition-Model/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/hellbats/Code/Phase-Transition-Model /home/hellbats/Code/Phase-Transition-Model /home/hellbats/Code/Phase-Transition-Model/build /home/hellbats/Code/Phase-Transition-Model/build /home/hellbats/Code/Phase-Transition-Model/build/CMakeFiles/build.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/build.dir/depend

