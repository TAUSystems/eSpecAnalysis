# CMakeList.txt : CMake project for eSpecAnalysisCMake, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (eSpecAnalysis "src/main.cpp" )
# add_library(eSpecAnalysis SHARED "src/main.cpp" )

target_link_libraries (eSpecAnalysis ${OpenCV_LIBS})
target_link_libraries (eSpecAnalysis ${PYTHON_LIBRARIES})
target_link_libraries (eSpecAnalysis OpenMP::OpenMP_CXX)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property (TARGET eSpecAnalysis PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
