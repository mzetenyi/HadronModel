set(SOURCES
    Config.cpp
    Config.hpp
    halfint.cpp
    halfint.hpp
    HelicityAmplitudes.cpp
    HelicityAmplitudes.hpp
    Histogram.cpp
    Histogram.hpp
    MultiArray.cpp
    MultiArray.hpp
    RandomNumberGenerator.cpp
    RandomNumberGenerator.hpp
    Spinors.cpp
    Spinors.hpp
    udouble.cpp
    udouble.hpp
    units.cpp
    units.hpp
    utils.cpp
    utils.hpp
    Vectors.cpp
    Vectors.hpp
    wavefunc.cpp
    wavefunc.hpp
)

add_compile_options(-O3)
add_library(FeynTools ${SOURCES})

target_include_directories(FeynTools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_property(TARGET FeynTools PROPERTY CXX_STANDARD 11)