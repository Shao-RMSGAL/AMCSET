add_library(amcset_common
    ${PROJECT_SOURCE_DIR}/common/src/amcset_common.cpp
    )

target_include_directories(amcset_common PUBLIC
    ${PROJECT_SOURCE_DIR}/common/include/
    )

target_link_libraries(amcset_common PUBLIC compiler_flags)
