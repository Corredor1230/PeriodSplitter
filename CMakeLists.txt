# Specify the minimum CMake version required
cmake_minimum_required(VERSION 3.29)

# Define the project name and set the C++ standard to C++17
project(WavParser VERSION 1.0 LANGUAGES CXX)

# Set the required C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Specify the source files and the executable name
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/source")
set(EXECUTABLE_NAME PeriodicSplitter)

# add_library(glad STATIC ${SOURCE_DIR}/src/glad.c)
# target_link_directories(glad PUBLIC ${SOURCE_DIR}/include)

find_package(OpenGL REQUIRED)

# Manually specify the paths for libsndfile include and lib directories
# Adjust these paths to where you installed libsndfile
set(LIBSNDFILE_INCLUDE_DIR "C:/Libraries/libsndfile/include")
set(LIBSNDFILE_LIBRARY_DIR "C:/Libraries/libsndfile/lib")
set(LIBPYIN_INCLUDE_DIR "C:/Libraries/LibPyin/source")
set(LIBPYIN_LIBRARY_DIR "C:/Libraries/LibPyin/Debug")
set(FFTW_INCLUDE_DIR "C:/Libraries/fftw")
set(FFTW_LIBRARY_DIR "C:/Libraries/fftw")
set(GLFW3_INCLUDE_DIR "C:/Libraries/vcpkg/packages/glfw3_x64-windows/include/GLFW")
set(GLFW3_LIBRARY_DIR "C:/Libraries/vcpkg/packages/glfw3_x64-windows/lib")
set(IMGUI_LIBRARY_DIR "C:/Libraries/vcpkg/packages/imgui_x64-windows/lib")
set(IMGUI_INCLUDE_DIR "C:/Libraries/vcpkg/packages/imgui_x64-windows/include")

# --- Define your logical components ---

set(ANALYSIS_SOURCES
    ${SOURCE_DIR}/analysis/HarmonicTracker.cpp
    ${SOURCE_DIR}/analysis/HarmonicTracker.h
    ${SOURCE_DIR}/analysis/NoiseTracker.cpp
    ${SOURCE_DIR}/analysis/NoiseTracker.h
    ${SOURCE_DIR}/analysis/SitranoAnalysis.cpp
    ${SOURCE_DIR}/analysis/SitranoAnalysis.h
	${SOURCE_DIR}/analysis/OvertoneFinder.h
	${SOURCE_DIR}/analysis/OvertoneFinder.cpp
)

set(DSP_SOURCES
    ${SOURCE_DIR}/dsp/PeriodCutter.cpp 
    ${SOURCE_DIR}/dsp/PeriodCutter.h
    ${SOURCE_DIR}/dsp/FIRFilter.cpp
    ${SOURCE_DIR}/dsp/FIRFilter.h
    ${SOURCE_DIR}/dsp/PitchFinder.h
    ${SOURCE_DIR}/dsp/PitchFinder.cpp
)

set(GUI_SOURCES
    ${SOURCE_DIR}/gui/GUI.h
    ${SOURCE_DIR}/gui/GUI.cpp
    ${SOURCE_DIR}/gui/SpectrumView.cpp
    ${SOURCE_DIR}/gui/SpectrumView.h
)

set(SUPPORT_SOURCES
    ${SOURCE_DIR}/support/Splitter.cpp 
    ${SOURCE_DIR}/support/Splitter.h
	${SOURCE_DIR}/support/AudioChunkStrategy.cpp
	${SOURCE_DIR}/support/AudioChunkStrategy.h
	${SOURCE_DIR}/support/PeriodLoopStrategy.cpp
	${SOURCE_DIR}/support/PeriodLoopStrategy.h
	${SOURCE_DIR}/support/SinglePeriodStrategy.cpp
	${SOURCE_DIR}/support/SinglePeriodStrategy.h
	${SOURCE_DIR}/support/IWindowStrategy.h
)

set(INCLUDE_FILES
    ${SOURCE_DIR}/include/csv.h
    ${SOURCE_DIR}/include/SitranoHeader.h
)

# --- Create the executable ---

# Start by creating the target with just the main file
add_executable(${EXECUTABLE_NAME} ${SOURCE_DIR}/main.cpp)

# --- Add all other sources to the target ---
target_sources(${EXECUTABLE_NAME} PRIVATE
    ${ANALYSIS_SOURCES}
    ${DSP_SOURCES}
    ${GUI_SOURCES}
    ${SUPPORT_SOURCES}
    ${INCLUDE_FILES}
)

# Add the source file(s)
# add_executable(${EXECUTABLE_NAME} ${SOURCE_DIR}/main.cpp 
# 				${SOURCE_DIR}/include/csv.h
# 				${SOURCE_DIR}/include/SitranoHeader.h
# 				${SOURCE_DIR}/analysis/HarmonicTracker.h
# 				${SOURCE_DIR}/analysis/HarmonicTracker.cpp
# 				${SOURCE_DIR}/analysis/NoiseTracker.cpp
# 				${SOURCE_DIR}/analysis/NoiseTracker.h
# 				${SOURCE_DIR}/analysis/SitranoAnalysis.cpp
# 				${SOURCE_DIR}/analysis/SitranoAnalysis.h
# 				${SOURCE_DIR}/dsp/PeriodCutter.cpp 
# 				${SOURCE_DIR}/dsp/PeriodCutter.h
# 				${SOURCE_DIR}/dsp/FIRFilter.cpp
# 				${SOURCE_DIR}/dsp/FIRFilter.h
# 				${SOURCE_DIR}/dsp/PitchFinder.h
# 				${SOURCE_DIR}/dsp/PitchFinder.cpp
# 				${SOURCE_DIR}/support/Splitter.cpp 
# 				${SOURCE_DIR}/support/Splitter.h 
# 				${SOURCE_DIR}/gui/GUI.h
# 				${SOURCE_DIR}/gui/GUI.cpp
# 				${SOURCE_DIR}/gui/SpectrumView.cpp
# 				${SOURCE_DIR}/gui/SpectrumView.h
# 				)

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")

# Include the libsndfile headers
target_include_directories(${EXECUTABLE_NAME} PRIVATE 
    ${LIBSNDFILE_INCLUDE_DIR}
    ${FFTW_INCLUDE_DIR}
    ${LIBPYIN_INCLUDE_DIR}
    ${GLFW3_INCLUDE_DIR}
    ${IMGUI_INCLUDE_DIR}
	${OPENGL_LIBRARIES}
	${SOURCE_DIR}
)

target_link_libraries(${EXECUTABLE_NAME} PRIVATE 
    "${LIBSNDFILE_LIBRARY_DIR}/sndfile.lib"
    "${FFTW_LIBRARY_DIR}/libfftw3-3.lib"
	"${FFTW_LIBRARY_DIR}/libfftw3f-3.lib"
	"${FFTW_LIBRARY_DIR}/libfftw3l-3.lib"
    "${LIBPYIN_LIBRARY_DIR}/LibPyin.lib"
    "${GLFW3_LIBRARY_DIR}/glfw3dll.lib"
    "${IMGUI_LIBRARY_DIR}/imgui.lib"
    opengl32
)