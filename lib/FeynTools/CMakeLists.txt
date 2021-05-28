add_library(FeynTools Config.cpp halfint.cpp HelicityAmplitudes.cpp Histogram.cpp MultiArray.cpp 
    RandomNumberGenerator.cpp Spinors.cpp udouble.cpp units.cpp utils.cpp Vectors.cpp wavefunc.cpp
    Config.hpp halfint.hpp HelicityAmplitudes.hpp Histogram.hpp MultiArray.hpp RandomNumberGenerator.hpp 
    Spinors.hpp udouble.hpp units.hpp utils.hpp Vectors.hpp wavefunc.hpp
    )

target_include_directories(FeynTools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})