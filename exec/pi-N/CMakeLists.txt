add_compile_options(-O3)
add_executable(piN_Ndilep BornTerms.cpp  cutoffs.cpp  observables.cpp  piN_Ndilep.cpp  piN_Ndilep_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  piN_Ndilep.hpp  piN_Ngammastar.hpp)

target_include_directories(piN_Ndilep PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Ndilep FeynTools)
target_link_libraries(piN_Ndilep ELModel)

add_executable(piN_Npipi BornTerms.cpp  cutoffs.cpp  observables.cpp  piN_Npipi.cpp  piN_Npipi_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  piN_Npipi.hpp  piN_Ngammastar.hpp)

target_include_directories(piN_Npipi PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Npipi FeynTools)
target_link_libraries(piN_Npipi ELModel)

add_executable(pionPhotoprod BornTerms.cpp  cutoffs.cpp  observables.cpp  pionPhotoprod.cpp  pionPhotoprod_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  pionPhotoprod.hpp  piN_Ngammastar.hpp)

target_include_directories(pionPhotoprod PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(pionPhotoprod FeynTools)
target_link_libraries(pionPhotoprod ELModel)

