cmake_minimum_required(VERSION 3.17)
project(CKKS)

set(CMAKE_CXX_STANDARD 17)
# seal attack
add_executable(sealattack
        my_code/examples.h
        my_code/sealmain.cpp my_code/sealattack.cpp my_code/sealattack.h my_code/utils.cpp my_code/utils.h)
find_package(SEAL 3.6 REQUIRED)
target_link_libraries(sealattack ntl gmp)
target_link_libraries(sealattack SEAL::seal)
target_include_directories(sealattack PRIVATE SEAL/native/src)

# helib attack
add_executable(helibattack my_code/helibmain.cpp my_code/helibattack.cpp my_code/helibattack.h my_code/utils.h my_code/utils.cpp)
find_package(helib)
target_link_libraries(helibattack helib)
target_include_directories(helibattack PRIVATE HElib/src)
target_include_directories(helibattack PRIVATE HElib/include)

#heaan attack
add_executable(heaanattack my_code/heaanmain.cpp my_code/utils.cpp my_code/utils.h)
target_link_libraries(heaanattack /home/msh/CKKS/HEAAN/HEAAN/lib/libHEAAN.a)
target_link_libraries(heaanattack ntl gmp pthread gf2x)
target_include_directories(heaanattack PRIVATE HEAAN/HEAAN/src)

#gf2 inv test
add_executable(gf2invtest my_code/gf2invtest.cpp)
target_link_libraries(gf2invtest ntl gmp pthread gf2x)

#ntru tr9 test
add_executable(tr9 my_code/tr9valid.cpp)
target_link_libraries(tr9 ntl gmp pthread)

#heaan-rns attack
add_executable(heaanrnsattack my_code/heaanrnsmain.cpp my_code/utils.cpp my_code/utils.h)
target_link_libraries(heaanrnsattack ntl pthread gmp)
target_link_libraries(heaanrnsattack /home/msh/CKKS/FullRNS-HEAAN/lib/libFRNSHEAAN.a)
target_include_directories(heaanrnsattack PRIVATE FullRNS-HEAAN/src)