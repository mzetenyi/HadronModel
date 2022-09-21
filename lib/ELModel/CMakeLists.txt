set(SOURCES
    FormFactors.cpp
    FormFactors.hpp
    Gamma.cpp
    Gamma.hpp
    Vrancx.cpp
    Vrancx.hpp)

add_compile_options(-O3)
add_library(ELModel ${SOURCES})

target_include_directories(ELModel PUBLIC ../FeynTools)
set_property(TARGET ELModel PROPERTY CXX_STANDARD 17)