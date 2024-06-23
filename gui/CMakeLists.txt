add_executable(amcset_gui
    ${PROJECT_SOURCE_DIR}/gui/src/main.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/amcset_gui.cpp
    )

target_include_directories(amcset_gui PRIVATE
    ${PROJECT_SOURCE_DIR}/gui/include/
    )

target_link_libraries(amcset_gui PUBLIC compiler_flags)
