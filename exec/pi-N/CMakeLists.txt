add_compile_options(-O3)
add_executable(piN_Ndilep BornTerms.cpp  cutoffs.cpp  observables.cpp  piN_Ndilep.cpp  piN_Ndilep_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  piN_Ndilep.hpp  piN_Ngammastar.hpp)

target_include_directories(piN_Ndilep PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Ndilep FeynTools)
target_link_libraries(piN_Ndilep ELModel)