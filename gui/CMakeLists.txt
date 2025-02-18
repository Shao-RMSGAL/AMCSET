set(CMAKE_AUTORCC ON)  # Automatically handle .qrc files

add_executable(amcset-gui
    ${PROJECT_SOURCE_DIR}/gui/src/main.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/amcset_gui.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/application.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/mainwindow.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/layerwindow.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/canvas.cpp
    ${PROJECT_SOURCE_DIR}/gui/include/amcset_gui.h
    ${PROJECT_SOURCE_DIR}/gui/include/application.h
    ${PROJECT_SOURCE_DIR}/gui/include/mainwindow.h
    ${PROJECT_SOURCE_DIR}/gui/include/layerwindow.h
    ${PROJECT_SOURCE_DIR}/gui/include/canvas.h
    ${PROJECT_SOURCE_DIR}/gui/resources/resources.qrc
    )

target_include_directories(amcset-gui PRIVATE
    ${PROJECT_SOURCE_DIR}/gui/include/
    )

target_link_libraries(amcset-gui PRIVATE compiler_flags)
