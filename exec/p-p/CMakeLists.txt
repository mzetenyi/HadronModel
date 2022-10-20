add_compile_options(-O3)

add_executable(NN_NDelta NN_NDelta_main.cpp Amplitude_NN_NDelta.cpp FeynmanRules.cpp Amplitude_NN_NDelta.hpp FeynmanRules.hpp)
target_compile_features(NN_NDelta PRIVATE cxx_std_17)
target_include_directories(NN_NDelta PUBLIC ../../lib/FeynTools)

target_link_libraries(NN_NDelta FeynTools)