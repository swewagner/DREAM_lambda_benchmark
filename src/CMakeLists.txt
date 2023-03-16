cmake_minimum_required(VERSION 3.7...3.25)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# add seqan3 as interface
add_library (common_interface INTERFACE)
target_link_libraries ("common_interface" INTERFACE seqan3::seqan3)
target_include_directories ("common_interface" INTERFACE ../include)
 
# build app generate_query_matches with seqan3
add_executable ("generate_query_matches" generate_query_matches.cpp)
target_link_libraries ("generate_query_matches" "common_interface")

# build app insert_query_matches with seqan3
add_executable ("insert_query_matches" insert_query_matches.cpp)
target_link_libraries ("insert_query_matches" "common_interface")

# build app translate_ref with seqan3
add_executable ("translate_ref" translate_ref.cpp)
target_link_libraries ("translate_ref" "common_interface")