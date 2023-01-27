cmake_minimum_required(VERSION 3.7...3.25)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
 
set (mason2_args ${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS})
list (APPEND mason2_args "-DSEQAN_BUILD_SYSTEM=APP:mason2")

include (ExternalProject)
ExternalProject_Add (
    mason2
    PREFIX mason2
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/seqan"
    CMAKE_ARGS "${mason2_args}"
    INSTALL_COMMAND ""
)
unset (mason2_args)

ExternalProject_Get_property(mason2 BINARY_DIR)
install(DIRECTORY "${BINARY_DIR}/bin/" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" USE_SOURCE_PERMISSIONS)

# add seqan3 as interface
add_library (common_interface INTERFACE)
target_link_libraries ("common_interface" INTERFACE seqan3::seqan3)
target_include_directories ("common_interface" INTERFACE ../include)
 
# build app generate_query_matches with seqan3
add_executable ("generate_query_matches" generate_query_matches.cpp)
target_link_libraries ("generate_query_matches" "common_interface")

# build app translate_ref with seqan3
add_executable ("translate_ref" translate_ref.cpp)
target_link_libraries ("translate_ref" "common_interface")