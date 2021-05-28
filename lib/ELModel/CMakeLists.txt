set(SOURCES 
    Gamma.cpp 
    Gamma.hpp 
    Vrancx.cpp 
    Vrancx.hpp)


add_library(ELModel ${SOURCES})

target_include_directories(ELModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_property(TARGET ELModel PROPERTY CXX_STANDARD 11)