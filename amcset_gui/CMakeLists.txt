# cmake_minimum_required(VERSION 3.15)

# project(acmset_gui VERSION 0.1 LANGUAGES CXX)


set(TS_FILES amcset_gui_en_US.ts)

set(PROJECT_SOURCES
        main.cpp
        amcset_main_window.cpp
        amcset_main_window.h
        amcset_main_window.ui
)

qt_add_executable(amcset_gui
    ${PROJECT_SOURCES}
)

target_link_libraries(amcset_gui PRIVATE Qt6::Widgets compiler_flags)

set_target_properties(amcset_gui PROPERTIES
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)
