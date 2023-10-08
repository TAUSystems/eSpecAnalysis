# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)


# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy (SET CMP0141 NEW)
  set (CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("eSpecAnalysis")
find_package (OpenCV REQUIRED)
find_package (PythonLibs REQUIRED)

include_directories (${OpenCV_INCLUDE_DIR})
include_directories (${Python_INCLUDE_DIRS})
include_directories ("..//matplotlib-cpp-master")


# Include sub-projects.
add_subdirectory ("eSpecAnalysis")

