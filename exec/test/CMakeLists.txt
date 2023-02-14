add_compile_options(-O3)

add_executable(test test.cpp)
target_compile_features(test PRIVATE cxx_std_17)
target_include_directories(test PUBLIC ../../lib/FeynTools)

target_link_libraries(test FeynTools)

add_executable(test32projector test32projector.cpp)
target_compile_features(test32projector PRIVATE cxx_std_17)
target_include_directories(test32projector PUBLIC ../../lib/FeynTools)

target_link_libraries(test32projector FeynTools)

add_executable(testTransformations testTransformations.cpp)
target_compile_features(testTransformations PRIVATE cxx_std_17)
target_include_directories(testTransformations PUBLIC ../../lib/FeynTools)

target_link_libraries(testTransformations FeynTools)