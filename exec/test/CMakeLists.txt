add_compile_options(-O3)

add_executable(test test.cpp)
target_compile_features(test PRIVATE cxx_std_17)
target_include_directories(NN_NDelta PUBLIC ../../lib/FeynTools)

target_link_libraries(test FeynTools)