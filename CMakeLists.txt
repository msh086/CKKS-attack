cmake_minimum_required(VERSION 3.17)
project(CKKS)

set(CMAKE_CXX_STANDARD 17)
# seal attack
add_executable(sealattack
        my_code/examples.h
        my_code/sealmain.cpp my_code/sealattack.cpp my_code/sealattack.h my_code/utils.cpp my_code/utils.h my_code/helibmain.cpp)
find_package(SEAL 3.6 REQUIRED)
target_link_libraries(sealattack ntl gmp)
target_link_libraries(sealattack SEAL::seal)
include_directories(my_code)
include_directories(my_code/tclap)
add_subdirectory(SEAL EXCLUDE_FROM_ALL)

# helib attack
add_executable(helibattack my_code/helibmain.cpp)
find_package(helib)
target_link_libraries(helibattack helib)
include_directories(HElib/src)
include_directories(HElib/include)