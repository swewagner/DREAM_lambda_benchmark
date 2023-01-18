cmake_minimum_required (VERSION 3.4)
project (generate_query_matches CXX)
 
# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
 
# build app with seqan3
add_executable (generate_query_matches generate_query_matches.cpp)
target_link_libraries (generate_query_matches seqan3::seqan3)

# build app with seqan3
add_executable (translate_ref translate_ref.cpp)
target_link_libraries (translate_ref seqan3::seqan3)