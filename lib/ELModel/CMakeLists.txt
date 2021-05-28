add_library(ELModel Gamma.cpp Vrancx.cpp Gamma.hpp Vrancx.hpp)

target_include_directories(ELModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})