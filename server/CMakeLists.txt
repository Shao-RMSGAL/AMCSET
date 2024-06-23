add_executable(amcset_server
    ${PROJECT_SOURCE_DIR}/server/src/main.cpp
    ${PROJECT_SOURCE_DIR}/server/src/amcset_server.cpp
    )

target_include_directories(amcset_server PRIVATE
    ${PROJECT_SOURCE_DIR}/server/include/
    )

target_link_libraries(amcset_server PUBLIC compiler_flags)
