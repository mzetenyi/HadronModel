add_compile_options(-O3)

add_executable(NN_NDelta NN_NDelta_main.cpp Amplitude_NN_NDelta.cpp FeynmanRules.cpp Amplitude_NN_NDelta.hpp FeynmanRules.hpp)
target_include_directories(NN_NDelta PUBLIC ../../lib/FeynTools)

target_link_libraries(NN_NDelta FeynTools)