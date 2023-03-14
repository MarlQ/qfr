# command line transcription app
add_executable(${PROJECT_NAME}_app ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)

# link against main project library
target_link_libraries(${PROJECT_NAME}_app PRIVATE ${PROJECT_NAME})

# add OpenMP flags
target_link_libraries(${PROJECT_NAME}_app PRIVATE "C:/Program Files/LLVM/lib/libomp.lib")
#set(OpenMP_CXX_FLAGS "-fopenmp")
#target_link_libraries(${PROJECT_NAME}_app PRIVATE ${OpenMP_CXX_FLAGS})
#target_compile_options(${PROJECT_NAME}_app PUBLIC ${OpenMP_CXX_FLAGS})

# enable interprocedural optimization if it is supported
enable_lto(${PROJECT_NAME}_app)

# create symlinks
add_custom_command(TARGET ${PROJECT_NAME}_app
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_app ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_app
                   COMMENT "Creating symlink for ${PROJECT_NAME}_app"
                   VERBATIM)
