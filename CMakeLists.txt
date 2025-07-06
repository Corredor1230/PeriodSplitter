# Specify the minimum CMake version required
cmake_minimum_required(VERSION 3.29)

# Define the project name and set the C++ standard to C++17
project(WavParser VERSION 1.0 LANGUAGES CXX)

# Set the required C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#NEW
include_directories("C:/Libraries/libsndfile/include")
include_directories("C:/Libraries/LibPyin/Debug")

# Specify the source files and the executable name
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/source")
set(EXECUTABLE_NAME PeriodicSplitter)

# Manually specify the paths for libsndfile include and lib directories
# Adjust these paths to where you installed libsndfile
set(LIBSNDFILE_INCLUDE_DIR "C:/Libraries/libsndfile/include")
set(LIBSNDFILE_LIBRARY_DIR "C:/Libraries/libsndfile/lib")
set(LIBPYIN_INCLUDE_DIR "C:/Libraries/LibPyin/source")
set(LIBPYIN_LIBRARY_DIR "C:/Libraries/LibPyin/Debug")
set(FFTW_INCLUDE_DIR "C:/Libraries/fftw")
set(FFTW_LIBRARY_DIR "C:/Libraries/fftw")


# Add the source file(s)
add_executable(${EXECUTABLE_NAME} ${SOURCE_DIR}/main.cpp 
				${SOURCE_DIR}/Splitter.cpp 
				${SOURCE_DIR}/Splitter.h 
				${SOURCE_DIR}/Correlator.cpp 
				${SOURCE_DIR}/Correlator.h
				${SOURCE_DIR}/csv.h)

#target_link_libraries(${EXECUTABLE_NAME} PRIVATE "C:/Libraries/libsndfile/lib/sndfile.lib")

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")

# Include the libsndfile headers
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${LIBSNDFILE_INCLUDE_DIR})
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${FFTW_INCLUDE_DIR})
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${LIBPYIN_INCLUDE_DIR})

# Link to the libsndfile library
target_link_libraries(${EXECUTABLE_NAME} PRIVATE "${LIBSNDFILE_LIBRARY_DIR}/sndfile.lib")
target_link_libraries(${EXECUTABLE_NAME} PRIVATE "${FFTW_LIBRARY_DIR}/libfftw3-3.lib")
target_link_libraries(${EXECUTABLE_NAME} PRIVATE "${LIBPYIN_LIBRARY_DIR}/LibPyin.lib")

# If you also have a dynamic version (.dll) of libsndfile, make sure to copy the DLL to your executable directory.