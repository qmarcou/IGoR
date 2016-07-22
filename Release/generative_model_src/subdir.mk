################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../generative_model_src/Aligner.cpp \
../generative_model_src/Deletion.cpp \
../generative_model_src/Dinuclmarkov.cpp \
../generative_model_src/Errorrate.cpp \
../generative_model_src/GenModel.cpp \
../generative_model_src/Genechoice.cpp \
../generative_model_src/Insertion.cpp \
../generative_model_src/IntStr.cpp \
../generative_model_src/Model_Parms.cpp \
../generative_model_src/Model_marginals.cpp \
../generative_model_src/Rec_Event.cpp \
../generative_model_src/Singleerrorrate.cpp \
../generative_model_src/Utils.cpp \
../generative_model_src/main.cpp 

OBJS += \
./generative_model_src/Aligner.o \
./generative_model_src/Deletion.o \
./generative_model_src/Dinuclmarkov.o \
./generative_model_src/Errorrate.o \
./generative_model_src/GenModel.o \
./generative_model_src/Genechoice.o \
./generative_model_src/Insertion.o \
./generative_model_src/IntStr.o \
./generative_model_src/Model_Parms.o \
./generative_model_src/Model_marginals.o \
./generative_model_src/Rec_Event.o \
./generative_model_src/Singleerrorrate.o \
./generative_model_src/Utils.o \
./generative_model_src/main.o 

CPP_DEPS += \
./generative_model_src/Aligner.d \
./generative_model_src/Deletion.d \
./generative_model_src/Dinuclmarkov.d \
./generative_model_src/Errorrate.d \
./generative_model_src/GenModel.d \
./generative_model_src/Genechoice.d \
./generative_model_src/Insertion.d \
./generative_model_src/IntStr.d \
./generative_model_src/Model_Parms.d \
./generative_model_src/Model_marginals.d \
./generative_model_src/Rec_Event.d \
./generative_model_src/Singleerrorrate.d \
./generative_model_src/Utils.d \
./generative_model_src/main.d 


# Each subdirectory must supply rules for building sources it contributes
generative_model_src/%.o: ../generative_model_src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -ljemalloc -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


