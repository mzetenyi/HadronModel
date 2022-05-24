add_compile_options(-O3)

add_executable(NN_NDelta NN_NDelta.cpp FeynmanRules.cpp NN_NDelta.hpp FeynmanRules.hpp)

target_include_directories(NN_NDelta PUBLIC ../../lib/FeynTools)

target_link_libraries(NN_NDelta FeynTools)