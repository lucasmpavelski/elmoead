################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Toolkit/ExampleProblems.cpp \
../Toolkit/ExampleShapes.cpp \
../Toolkit/ExampleTransitions.cpp \
../Toolkit/FrameworkFunctions.cpp \
../Toolkit/Misc.cpp \
../Toolkit/ShapeFunctions.cpp \
../Toolkit/TransFunctions.cpp 

OBJS += \
./Toolkit/ExampleProblems.o \
./Toolkit/ExampleShapes.o \
./Toolkit/ExampleTransitions.o \
./Toolkit/FrameworkFunctions.o \
./Toolkit/Misc.o \
./Toolkit/ShapeFunctions.o \
./Toolkit/TransFunctions.o 

CPP_DEPS += \
./Toolkit/ExampleProblems.d \
./Toolkit/ExampleShapes.d \
./Toolkit/ExampleTransitions.d \
./Toolkit/FrameworkFunctions.d \
./Toolkit/Misc.d \
./Toolkit/ShapeFunctions.d \
./Toolkit/TransFunctions.d 


# Each subdirectory must supply rules for building sources it contributes
Toolkit/%.o: ../Toolkit/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


