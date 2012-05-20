# ---------------------------------------------------------------------------
#           Maxwell-Boltzmann distribution demonstration (2D/3D)
#   Copyright (C) 2012  Jose Luis Blanco
# ---------------------------------------------------------------------------
PROJECT(maxwell_boltzmann)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)  # Required by CMake 2.7+
endif(COMMAND cmake_policy)

FIND_PACKAGE(MRPT REQUIRED base gui)

# Declare the target (an executable)
ADD_EXECUTABLE(boltz-3d maxwell-boltzmann.cpp)
TARGET_LINK_LIBRARIES(boltz-3d ${MRPT_LIBS})

# Declare the target (an executable)
ADD_EXECUTABLE(boltz-2d maxwell-boltzmann.cpp)
TARGET_LINK_LIBRARIES(boltz-2d ${MRPT_LIBS})
SET_TARGET_PROPERTIES(boltz-2d PROPERTIES COMPILE_FLAGS "-DSIMUL_2D")

# Set optimized building:
IF(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native")
ENDIF(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")

