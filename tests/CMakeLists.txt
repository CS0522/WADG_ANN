set(CMAKE_CXX_STANDARD 11)

add_executable(test_nsg_index test_nsg_index.cpp)
target_link_libraries(test_nsg_index ${PROJECT_NAME} -ltcmalloc)

add_executable(test_nsg_search test_nsg_search.cpp)
target_link_libraries(test_nsg_search ${PROJECT_NAME})

add_executable(test_nsg_optimized_search test_nsg_optimized_search.cpp)
target_link_libraries(test_nsg_optimized_search ${PROJECT_NAME} -ltcmalloc)

add_executable(test_wadg_search test_wadg_search.cpp)
target_link_libraries(test_wadg_search ${PROJECT_NAME})

add_executable(test_cal_top_K_precision test_cal_top_K_precision.cpp)
add_executable(test_cal_top_K_precision_bm test_cal_top_K_precision_bm.cpp)

# mimalloc
add_executable(test_nsg_search_bm test_nsg_search_bm.cpp)
target_link_libraries(test_nsg_search_bm ${PROJECT_NAME})
target_link_libraries(test_nsg_search_bm PUBLIC mimalloc)

add_executable(test_nsg_index_bm test_nsg_index_bm.cpp)
target_link_libraries(test_nsg_index_bm ${PROJECT_NAME} -ltcmalloc)
target_link_libraries(test_nsg_index_bm PUBLIC mimalloc)