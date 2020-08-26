cmake_minimum_required(VERSION 2.8)
# You should not specify the compiler in the CMakeLists.txt file !
#SET(CMAKE_CXX_COMPILER /usr/bin/clang++-6.0)
#SET(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries/linux/bin/intel64/icc)
# You should not have specified the compiler in CMakeLists.txt

# Boost
find_package(Boost 1.50.0 REQUIRED)
list(APPEND LIBRARIES ${Boost_LIBRARIES})
list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
# a gcc special 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fno-finite-math-only -mrecip")
# and an icc special
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --funroll-all-loops")
# and a clang special
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fno-finite-math-only -mrecip")

# Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND INCLUDES Aboria/src)
list(APPEND INCLUDES Aboria/third-party)

# Tinyply
list(APPEND INCLUDES tinyply/source)

# yaml-cpp
list(APPEND INCLUDES yaml-cpp/include)
list(APPEND INCLUDES yaml-cpp/src)

list(APPEND INCLUDES src)

include_directories(src ${INCLUDES})

file(GLOB all_SRC
    "src/*.cpp"
	"tinyply/source/tinyply.cpp"
    "yaml-cpp/src/*.cpp"
    "yaml-cpp/include/*.cpp"
    )

set(SOURCE ${all_SRC})

add_executable(meshless ${SOURCE})

target_link_libraries(meshless ${LIBRARIES})
