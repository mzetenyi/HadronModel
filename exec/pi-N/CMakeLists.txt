add_compile_options(-O3)

add_executable(piN_Ndilep BornTerms.cpp  cutoffs.cpp  observables.cpp  piN_Ndilep.cpp  piN_Ndilep_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  piN_Ndilep.hpp  piN_Ngammastar.hpp)

target_include_directories(piN_Ndilep PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Ndilep ELModel)
target_link_libraries(piN_Ndilep FeynTools)

add_executable(piN_Npipi BornTerms.cpp  cutoffs.cpp  observables.cpp  piN_Npipi.cpp  piN_Npipi_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  piN_Npipi.hpp  piN_Ngammastar.hpp)

target_include_directories(piN_Npipi PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(piN_Npipi ELModel)
target_link_libraries(piN_Npipi FeynTools)

add_executable(pionPhotoprod BornTerms.cpp  cutoffs.cpp  observables.cpp  pionPhotoprod.cpp  pionPhotoprod_main.cpp  piN_Ngammastar.cpp 
    BornTerms.hpp  cutoffs.hpp  observables.hpp  pionPhotoprod.hpp  piN_Ngammastar.hpp)

target_include_directories(pionPhotoprod PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(pionPhotoprod ELModel)
target_link_libraries(pionPhotoprod FeynTools)

add_executable(pionPhotoprodAnalytic pionPhotoprodAnalytic.cpp pionPhotoprodAnalytic_main.cpp)

target_include_directories(pionPhotoprodAnalytic PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(pionPhotoprodAnalytic ELModel)
target_link_libraries(pionPhotoprodAnalytic FeynTools)

add_executable(pionPhotoprodTest pionPhotoprodTest.cpp BornTerms.cpp BornTerms.hpp)

target_include_directories(pionPhotoprodTest PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(pionPhotoprodTest ELModel)
target_link_libraries(pionPhotoprodTest FeynTools)

add_executable(rhoPhase rhoPhase.cpp BornTerms.cpp  cutoffs.cpp  observables.cpp  pionPhotoprod.cpp piN_Ngammastar.cpp 
BornTerms.hpp  cutoffs.hpp  observables.hpp  pionPhotoprod.hpp  piN_Ngammastar.hpp)

target_include_directories(rhoPhase PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(rhoPhase ELModel)
target_link_libraries(rhoPhase FeynTools)

add_executable(fit_resonance_couplings fit_resonance_couplings.cpp cutoffs.cpp cutoffs.hpp resdec.cpp resdec.hpp)

target_include_directories(fit_resonance_couplings PUBLIC ../../lib/FeynTools ../../lib/ELModel)

target_link_libraries(fit_resonance_couplings ELModel)
target_link_libraries(fit_resonance_couplings FeynTools)
