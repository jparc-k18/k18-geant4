#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(S2S)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
#set(ROOT_INCLUDE_DIR /home/gogami/root/root-5.34.36-install/include/root)
#set(ROOT_LIBRARIES -L/home/gogami/root/root-5.34.36-install/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic)
set(ROOT_INCLUDE_DIR /sw/packages/root/6.20.08/include)
set(ROOT_LIBRARIES -L/sw/packages/root/6.20.08/lib/ -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
add_executable(S2S S2S.cc ${sources} ${headers})
target_link_libraries(S2S ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX /home/had/yudai/private/E27_ichikawa_new1_1/S2S/k18geant4)
install(TARGETS S2S RUNTIME DESTINATION bin)



