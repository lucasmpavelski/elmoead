# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/lucasmp/projects/elmoead_gpu

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lucasmp/projects/elmoead_gpu

# Include any dependencies generated for this target.
include lib/hv/CMakeFiles/hv.dir/depend.make

# Include the progress variables for this target.
include lib/hv/CMakeFiles/hv.dir/progress.make

# Include the compile flags for this target's objects.
include lib/hv/CMakeFiles/hv.dir/flags.make

lib/hv/CMakeFiles/hv.dir/hv.c.o: lib/hv/CMakeFiles/hv.dir/flags.make
lib/hv/CMakeFiles/hv.dir/hv.c.o: lib/hv/hv.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucasmp/projects/elmoead_gpu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lib/hv/CMakeFiles/hv.dir/hv.c.o"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/hv.dir/hv.c.o   -c /home/lucasmp/projects/elmoead_gpu/lib/hv/hv.c

lib/hv/CMakeFiles/hv.dir/hv.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/hv.dir/hv.c.i"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lucasmp/projects/elmoead_gpu/lib/hv/hv.c > CMakeFiles/hv.dir/hv.c.i

lib/hv/CMakeFiles/hv.dir/hv.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/hv.dir/hv.c.s"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lucasmp/projects/elmoead_gpu/lib/hv/hv.c -o CMakeFiles/hv.dir/hv.c.s

lib/hv/CMakeFiles/hv.dir/hv.c.o.requires:

.PHONY : lib/hv/CMakeFiles/hv.dir/hv.c.o.requires

lib/hv/CMakeFiles/hv.dir/hv.c.o.provides: lib/hv/CMakeFiles/hv.dir/hv.c.o.requires
	$(MAKE) -f lib/hv/CMakeFiles/hv.dir/build.make lib/hv/CMakeFiles/hv.dir/hv.c.o.provides.build
.PHONY : lib/hv/CMakeFiles/hv.dir/hv.c.o.provides

lib/hv/CMakeFiles/hv.dir/hv.c.o.provides.build: lib/hv/CMakeFiles/hv.dir/hv.c.o


lib/hv/CMakeFiles/hv.dir/avl.c.o: lib/hv/CMakeFiles/hv.dir/flags.make
lib/hv/CMakeFiles/hv.dir/avl.c.o: lib/hv/avl.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucasmp/projects/elmoead_gpu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object lib/hv/CMakeFiles/hv.dir/avl.c.o"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/hv.dir/avl.c.o   -c /home/lucasmp/projects/elmoead_gpu/lib/hv/avl.c

lib/hv/CMakeFiles/hv.dir/avl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/hv.dir/avl.c.i"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lucasmp/projects/elmoead_gpu/lib/hv/avl.c > CMakeFiles/hv.dir/avl.c.i

lib/hv/CMakeFiles/hv.dir/avl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/hv.dir/avl.c.s"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lucasmp/projects/elmoead_gpu/lib/hv/avl.c -o CMakeFiles/hv.dir/avl.c.s

lib/hv/CMakeFiles/hv.dir/avl.c.o.requires:

.PHONY : lib/hv/CMakeFiles/hv.dir/avl.c.o.requires

lib/hv/CMakeFiles/hv.dir/avl.c.o.provides: lib/hv/CMakeFiles/hv.dir/avl.c.o.requires
	$(MAKE) -f lib/hv/CMakeFiles/hv.dir/build.make lib/hv/CMakeFiles/hv.dir/avl.c.o.provides.build
.PHONY : lib/hv/CMakeFiles/hv.dir/avl.c.o.provides

lib/hv/CMakeFiles/hv.dir/avl.c.o.provides.build: lib/hv/CMakeFiles/hv.dir/avl.c.o


# Object files for target hv
hv_OBJECTS = \
"CMakeFiles/hv.dir/hv.c.o" \
"CMakeFiles/hv.dir/avl.c.o"

# External object files for target hv
hv_EXTERNAL_OBJECTS =

lib/hv/libhv.a: lib/hv/CMakeFiles/hv.dir/hv.c.o
lib/hv/libhv.a: lib/hv/CMakeFiles/hv.dir/avl.c.o
lib/hv/libhv.a: lib/hv/CMakeFiles/hv.dir/build.make
lib/hv/libhv.a: lib/hv/CMakeFiles/hv.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucasmp/projects/elmoead_gpu/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libhv.a"
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && $(CMAKE_COMMAND) -P CMakeFiles/hv.dir/cmake_clean_target.cmake
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hv.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/hv/CMakeFiles/hv.dir/build: lib/hv/libhv.a

.PHONY : lib/hv/CMakeFiles/hv.dir/build

lib/hv/CMakeFiles/hv.dir/requires: lib/hv/CMakeFiles/hv.dir/hv.c.o.requires
lib/hv/CMakeFiles/hv.dir/requires: lib/hv/CMakeFiles/hv.dir/avl.c.o.requires

.PHONY : lib/hv/CMakeFiles/hv.dir/requires

lib/hv/CMakeFiles/hv.dir/clean:
	cd /home/lucasmp/projects/elmoead_gpu/lib/hv && $(CMAKE_COMMAND) -P CMakeFiles/hv.dir/cmake_clean.cmake
.PHONY : lib/hv/CMakeFiles/hv.dir/clean

lib/hv/CMakeFiles/hv.dir/depend:
	cd /home/lucasmp/projects/elmoead_gpu && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucasmp/projects/elmoead_gpu /home/lucasmp/projects/elmoead_gpu/lib/hv /home/lucasmp/projects/elmoead_gpu /home/lucasmp/projects/elmoead_gpu/lib/hv /home/lucasmp/projects/elmoead_gpu/lib/hv/CMakeFiles/hv.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/hv/CMakeFiles/hv.dir/depend

