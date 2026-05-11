set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSandHypoplasticApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories(${KRATOS_SOURCE_DIR}/kratos)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${KRATOS_SOURCE_DIR}/applications/MPMApplication)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/MPMApplication)

set(KRATOS_SAND_HYPOPLASTIC_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/sand_hypoplastic_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sand_hypoplastic_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/path_a_fortran_dll/sand_hypoplastic_fortran_dll_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/path_b_cpp_port/sand_hypoplastic_kernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/path_b_cpp_port/sand_hypoplastic_cpp_law.cpp
)

set(KRATOS_SAND_HYPOPLASTIC_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/sand_hypoplastic_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

add_library(KratosSandHypoplasticCore SHARED ${KRATOS_SAND_HYPOPLASTIC_APPLICATION_CORE})
target_link_libraries(KratosSandHypoplasticCore PUBLIC KratosCore KratosMPMCore)
set_target_properties(KratosSandHypoplasticCore PROPERTIES COMPILE_DEFINITIONS "SAND_HYPOPLASTIC_APPLICATION=EXPORT,API")

pybind11_add_module(KratosSandHypoplasticApplication MODULE THIN_LTO ${KRATOS_SAND_HYPOPLASTIC_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosSandHypoplasticApplication PRIVATE KratosSandHypoplasticCore)
set_target_properties(KratosSandHypoplasticApplication PROPERTIES PREFIX "")

IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosSandHypoplasticCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosSandHypoplasticApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosSandHypoplasticApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosSandHypoplasticApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/SandHypoplasticApplication.py KratosMultiphysics/SandHypoplasticApplication/__init__.py)

get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME})

install(TARGETS KratosSandHypoplasticCore DESTINATION libs)
install(TARGETS KratosSandHypoplasticApplication DESTINATION libs)

set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosSandHypoplasticCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosSandHypoplasticApplication" PARENT_SCOPE)
