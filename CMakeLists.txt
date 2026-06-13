# Versión mínima de CMake
cmake_minimum_required(VERSION 3.29)

# Nombre del proyecto y C++ 17
project(Sihat VERSION 1.0 LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- ¡LA MAGIA DE VCPKG EMPIEZA AQUÍ! ---
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(fftw3 CONFIG REQUIRED)
find_package(fftw3f CONFIG REQUIRED)
find_package(SndFile CONFIG REQUIRED)
find_package(RtAudio CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

# --- RUTAS A BIBLIOTECAS (Solo para las que no están en vcpkg) ---
set(LIBPYIN_INCLUDE_DIR "C:/Libraries/LibPyin/source")
set(LIBPYIN_LIBRARY_DIR "C:/Libraries/LibPyin/Debug")

# --- Rutas del Proyecto (Esto está bien) ---
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/source")
set(EXECUTABLE_NAME PeriodicSplitter)
set(PROJECT_INCLUDES ${CMAKE_SOURCE_DIR}/include)

# --- Fuentes del Proyecto (Sin cambios) ---
set(ANALYSIS_SOURCES
    ${SOURCE_DIR}/analysis/HarmonicTracker.cpp
    ${SOURCE_DIR}/analysis/HarmonicTracker.h
    ${SOURCE_DIR}/analysis/NoiseTracker.cpp
    ${SOURCE_DIR}/analysis/NoiseTracker.h
    ${SOURCE_DIR}/analysis/SitranoAnalysis.cpp
    ${SOURCE_DIR}/analysis/SitranoAnalysis.h
    ${SOURCE_DIR}/analysis/OvertoneFinder.h
    ${SOURCE_DIR}/analysis/OvertoneFinder.cpp
    ${SOURCE_DIR}/analysis/TransientRecognition.cpp
    ${SOURCE_DIR}/analysis/TransientRecognition.h
    ${SOURCE_DIR}/analysis/TransientAnalysis.cpp
    ${SOURCE_DIR}/analysis/TransientAnalysis.h
)

set(DSP_SOURCES
    ${SOURCE_DIR}/dsp/PeriodCutter.cpp 
    ${SOURCE_DIR}/dsp/PeriodCutter.h
    ${SOURCE_DIR}/dsp/FIRFilter.cpp
    ${SOURCE_DIR}/dsp/FIRFilter.h
    ${SOURCE_DIR}/dsp/PitchFinder.h
    ${SOURCE_DIR}/dsp/PitchFinder.cpp
    ${SOURCE_DIR}/dsp/ZPFilter.h
    ${SOURCE_DIR}/dsp/WTransform.h
    ${SOURCE_DIR}/dsp/WTransform.cpp
    ${SOURCE_DIR}/dsp/MorletWavelet.h
    ${SOURCE_DIR}/dsp/MedianFilter.cpp
    ${SOURCE_DIR}/dsp/MedianFilter.h
)

set(GUI_SOURCES
    ${SOURCE_DIR}/gui/GUI.h
    ${SOURCE_DIR}/gui/GUI.cpp
    ${SOURCE_DIR}/gui/ConfigEditor.h
    ${SOURCE_DIR}/gui/SihatEditor.cpp
    ${SOURCE_DIR}/gui/SihatEditor.h
	${SOURCE_DIR}/gui/EditorViewModel.h
    ${SOURCE_DIR}/gui/glad.c
)

set(SUPPORT_SOURCES
    ${SOURCE_DIR}/support/Splitter.cpp 
    ${SOURCE_DIR}/support/Splitter.h
    ${SOURCE_DIR}/support/AudioChunkStrategy.h
    ${SOURCE_DIR}/support/PeriodLoopStrategy.h
    ${SOURCE_DIR}/support/SinglePeriodStrategy.h
    ${SOURCE_DIR}/support/IWindowStrategy.h
    ${SOURCE_DIR}/support/SihatLogger.h
)

set(INCLUDE_FILES
    ${SOURCE_DIR}/include/csv.h
    ${SOURCE_DIR}/include/SitranoHeader.h
    ${SOURCE_DIR}/include/SitranoHeader.cpp
    ${SOURCE_DIR}/include/ResynthHeader.h
    ${SOURCE_DIR}/include/json.hpp
    ${SOURCE_DIR}/include/ResynthConfig.h
    ${SOURCE_DIR}/include/InterpConfig.h
    ${SOURCE_DIR}/include/InterpHeader.h
)

set(FILE_FILES
    ${SOURCE_DIR}/file/Sihatcsv.h
    ${SOURCE_DIR}/file/SihatFile.h
    ${SOURCE_DIR}/file/SihatJson.h
    ${SOURCE_DIR}/file/JsonAdapters.h
)

set(RESYNTH_FILES
    ${SOURCE_DIR}/resynthesis/SihatResynth.h
    ${SOURCE_DIR}/resynthesis/SihatResynth.cpp
)

set(INTERP_FILES
    ${SOURCE_DIR}/interpolation/Interpolator.h
    ${SOURCE_DIR}/interpolation/Interpolator.cpp
    ${SOURCE_DIR}/interpolation/AudioPlayer.h
    ${SOURCE_DIR}/interpolation/AudioPlayer.cpp
)

set(APPLICATION_SOURCES
    ${SOURCE_DIR}/SihatApplication.cpp
    ${SOURCE_DIR}/SihatApplication.h
)

# --- Crear el ejecutable ---
# ¡¡CAMBIO IMPORTANTE!!
# Ponemos TODAS las fuentes en add_executable.
# Esto es más robusto que target_sources.
add_executable(${EXECUTABLE_NAME}
    ${SOURCE_DIR}/main.cpp
    ${APPLICATION_SOURCES}
    ${ANALYSIS_SOURCES}
    ${DSP_SOURCES}
    ${GUI_SOURCES}
    ${SUPPORT_SOURCES}
    ${INCLUDE_FILES}
    ${FILE_FILES}
    ${RESYNTH_FILES}
    ${INTERP_FILES}
)

target_compile_definitions(${EXECUTABLE_NAME} PRIVATE NOMINMAX)

if(MSVC)
    target_compile_options(${EXECUTABLE_NAME} PRIVATE /bigobj)
endif()

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")

# --- Incluir directorios (Híbrido) ---
target_include_directories(${EXECUTABLE_NAME} PRIVATE 
    # vcpkg añade automáticamente las rutas de include
    ${LIBPYIN_INCLUDE_DIR}
    ${SOURCE_DIR}
    ${PROJECT_INCLUDES}
)

# --- Enlazar bibliotecas (Híbrido) ---
# ¡¡CAMBIOS IMPORTANTES AQUÍ!!
target_link_libraries(${EXECUTABLE_NAME} PRIVATE 
    # Bibliotecas encontradas por vcpkg
    imgui::imgui
    glfw
    SndFile::sndfile
    # FFTW3::fftw3        
    FFTW3::fftw3f       
    RtAudio::rtaudio
    # Bibliotecas del sistema
    OpenGL::GL
    Eigen3::Eigen
    
    # Tu biblioteca manual
    "${LIBPYIN_LIBRARY_DIR}/LibPyin.lib"
)

add_custom_command(TARGET ${EXECUTABLE_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${SOURCE_DIR}/configFiles/default.json"
    $<TARGET_FILE_DIR:${EXECUTABLE_NAME}>/default.json
)

if(WIN32)
    add_custom_command(TARGET ${EXECUTABLE_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${LIBPYIN_LIBRARY_DIR}/LibPyin.dll"  # <--- Make sure this filename matches exactly!
        $<TARGET_FILE_DIR:${EXECUTABLE_NAME}>
    )
endif()