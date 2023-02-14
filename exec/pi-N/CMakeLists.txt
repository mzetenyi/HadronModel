add_compile_options(-O3)

add_executable(piN_elastic Amplitude_piN_elastic.cpp piN_elastic.cpp piN_elastic_main.cpp Amplitude_piN_elastic.hpp 
piN_elastic.hpp)

target_include_directories(piN_elastic PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_compile_features(piN_elastic PRIVATE cxx_std_17)

target_link_libraries(piN_elastic ELModel)
target_link_libraries(piN_elastic FeynTools)

add_executable(R_piN_decay Amplitude_R_piN.cpp R_piN_decay.cpp R_piN_decay_main.cpp Amplitude_R_piN.hpp 
R_piN_decay.hpp)

target_include_directories(R_piN_decay PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_compile_features(R_piN_decay PRIVATE cxx_std_17)

target_link_libraries(R_piN_decay ELModel)
target_link_libraries(R_piN_decay FeynTools)

add_executable(piN_Ndilep BornTerms.cpp cutoffs.cpp observables.cpp piN_Ndilep.cpp piN_Ndilep_main.cpp piN_Ngammastar.cpp
    BornTerms.hpp cutoffs.hpp observables.hpp piN_Ndilep.hpp piN_Ngammastar.hpp)

target_include_directories(piN_Ndilep PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Ndilep ELModel)
target_link_libraries(piN_Ndilep FeynTools)

add_executable(piN_Npipi BornTerms.cpp cutoffs.cpp observables.cpp piN_Npipi.cpp piN_Npipi_main.cpp piN_Ngammastar.cpp
    BornTerms.hpp cutoffs.hpp observables.hpp piN_Npipi.hpp piN_Ngammastar.hpp)

target_include_directories(piN_Npipi PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Npipi ELModel)
target_link_libraries(piN_Npipi FeynTools)

add_executable(piN_Npipi_generator BornTerms.cpp cutoffs.cpp observables.cpp piN_Npipi.cpp piN_Npipi_generator.cpp piN_Ngammastar.cpp
    BornTerms.hpp cutoffs.hpp observables.hpp piN_Npipi.hpp piN_Ngammastar.hpp)

target_include_directories(piN_Npipi_generator PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Npipi_generator ELModel)
target_link_libraries(piN_Npipi_generator FeynTools)

add_executable(pionPhotoprod BornTerms.cpp cutoffs.cpp observables.cpp pionPhotoprod.cpp pionPhotoprod_main.cpp piN_Ngammastar.cpp
    BornTerms.hpp cutoffs.hpp observables.hpp pionPhotoprod.hpp piN_Ngammastar.hpp)

target_include_directories(pionPhotoprod PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(pionPhotoprod ELModel)
target_link_libraries(pionPhotoprod FeynTools)
