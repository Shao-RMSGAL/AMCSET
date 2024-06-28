add_executable(amcset-server
    ${PROJECT_SOURCE_DIR}/server/src/main.cpp
    ${PROJECT_SOURCE_DIR}/server/src/amcset_server.cpp
    ${PROJECT_SOURCE_DIR}/server/include/amcset_server.h
    )

target_include_directories(amcset-server PRIVATE
    ${PROJECT_SOURCE_DIR}/server/include/
    )

target_link_libraries(amcset-server PRIVATE compiler_flags)
