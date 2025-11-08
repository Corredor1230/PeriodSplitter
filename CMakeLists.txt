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
find_package(SndFile CONFIG REQUIRED) # Corregido (antes libsndfile)

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
    ${SOURCE_DIR}/gui/ConfigEditor.h
    ${SOURCE_DIR}/gui/SihatEditor.cpp
    ${SOURCE_DIR}/gui/SihatEditor.h
    ${SOURCE_DIR}/gui/glad.c
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
    ${SOURCE_DIR}/include/SitranoHeader.cpp
    ${SOURCE_DIR}/include/json.hpp
)

set(FILE_FILES
    ${SOURCE_DIR}/file/Sihatcsv.h
    ${SOURCE_DIR}/file/SihatFile.cpp
    ${SOURCE_DIR}/file/SihatFile.h
    ${SOURCE_DIR}/file/SihatJson.cpp
    ${SOURCE_DIR}/file/SihatJson.h
    ${SOURCE_DIR}/file/JsonAdapters.h
)

# --- Crear el ejecutable ---
# ¡¡CAMBIO IMPORTANTE!!
# Ponemos TODAS las fuentes en add_executable.
# Esto es más robusto que target_sources.
add_executable(${EXECUTABLE_NAME}
    ${SOURCE_DIR}/main.cpp
    ${ANALYSIS_SOURCES}
    ${DSP_SOURCES}
    ${GUI_SOURCES}
    ${SUPPORT_SOURCES}
    ${INCLUDE_FILES}
    ${FILE_FILES}
)

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
    
    # Bibliotecas del sistema
    OpenGL::GL
    
    # Tu biblioteca manual
    "${LIBPYIN_LIBRARY_DIR}/LibPyin.lib"
)