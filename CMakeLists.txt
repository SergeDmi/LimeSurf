cmake_minimum_required(VERSION 2.8)

# Boost
find_package(Boost 1.50.0 REQUIRED)
list(APPEND LIBRARIES ${Boost_LIBRARIES})
list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

# VTK
find_package(VTK REQUIRED)
if (VTK_FOUND)
    add_definitions(-DHAVE_VTK)
endif(VTK_FOUND)
list(APPEND LIBRARIES ${VTK_LIBRARIES})
list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11")

# Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND INCLUDES Aboria/src)
list(APPEND INCLUDES Aboria/third-party)

# Cytosim
list(APPEND INCLUDES Cytosim)

# LimeSurf
list(APPEND INCLUDES .)

include_directories(src ${INCLUDES})

file(GLOB all_SRC
    "src/*.cpp"
    "Cytosim/*.cc"
)

#MESSAGE( STATUS ${all_SRC})

set(SOURCE ${all_SRC}
)

add_executable(meshless ${SOURCE})
target_link_libraries(meshless ${LIBRARIES})
