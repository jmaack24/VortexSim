
include_directories(.)

set(VSIM_SRC
  AdMoulPC.cpp
  CenterVort.cpp
  Constants.cpp
  FlowField.cpp
  Hamiltonian.cpp
  Logger.cpp
  MomentInertia.cpp
  Moments.cpp
  OutputHQ.cpp
  PointVortex.cpp
  RuKu4.cpp
  ScreenedFlow.cpp
  ScreenedHamiltonian.cpp
  SphereFlow.cpp
  Storage.cpp
  Timer.cpp
  TwoLayerCenter.cpp
  TwoLayerFlow.cpp
  TwoLayerHamiltonian.cpp
  TwoLayerMoments.cpp
  UnboundedFlow.cpp
  VortexSim.cpp
  main.cpp)

add_executable(vsim ${VSIM_SRC})

find_package(GSL REQUIRED)
target_include_directories(vsim PUBLIC "${GSL_INCLUDE_DIR}")
# target_include_directories(vsim PUBLIC "${GSL_LIBRARIES}")
target_link_libraries(vsim PRIVATE "${GSL_LIBRARIES}")
