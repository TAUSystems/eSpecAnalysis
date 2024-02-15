# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)
# set (CMAKE_CXX_STANDARD 11)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy (SET CMP0141 NEW)
  set (CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("eSpecAnalysis")
find_package (OpenCV REQUIRED core calib3d highgui imgproc features2d)
find_package (PythonLibs REQUIRED)
find_package (OpenMP REQUIRED)
#find_package (Python3 REQUIRED COMPONENTS Interpreter NumPy)

include_directories(${OpenMP_LIBRARY_DIR})
include_directories (${OpenCV_INCLUDE_DIR})
include_directories (${Python_INCLUDE_DIRS})
include_directories ("/usr/local/include/python3.9")
include_directories ("/usr/local/lib/python3.9/site-packages/numpy/core/include")
# The matplotlibcpp.h file from the matplotlib-cpp repository should be put into 
# <project_root>/eSpecAnalysis/include
include_directories ("eSpecAnalysis/include")

# Include sub-projects.
add_subdirectory ("eSpecAnalysis")
