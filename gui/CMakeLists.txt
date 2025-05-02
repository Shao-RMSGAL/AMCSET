set(CMAKE_AUTORCC ON)  # Automatically handle .qrc files

add_executable(amcset-gui
    ${PROJECT_SOURCE_DIR}/gui/src/main.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/amcset_gui.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/application.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/mainwindow.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/layerwindow.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/layer_delegate.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/canvas.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/3dview.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/counter.cpp
    ${PROJECT_SOURCE_DIR}/gui/src/simulation.cpp
    ${PROJECT_SOURCE_DIR}/gui/include/amcset_gui.h
    ${PROJECT_SOURCE_DIR}/gui/include/application.h
    ${PROJECT_SOURCE_DIR}/gui/include/mainwindow.h
    ${PROJECT_SOURCE_DIR}/gui/include/layerwindow.h
    ${PROJECT_SOURCE_DIR}/gui/include/layer_delegate.h
    ${PROJECT_SOURCE_DIR}/gui/include/canvas.h
    ${PROJECT_SOURCE_DIR}/gui/include/3dview.h
    ${PROJECT_SOURCE_DIR}/gui/include/counter.h
    ${PROJECT_SOURCE_DIR}/gui/include/simulation.h
    ${PROJECT_SOURCE_DIR}/gui/resources/resources.qrc
    )

target_include_directories(amcset-gui PRIVATE
    ${PROJECT_SOURCE_DIR}/gui/include/
    )

target_link_libraries(amcset-gui PRIVATE compiler_flags)

# In your CMakeLists.txt, after your target is defined:
if(WIN32)
  # Find the windeployqt executable
  get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
  get_filename_component(QT_BIN_DIR "${QT_QMAKE_EXECUTABLE}" DIRECTORY)
  find_program(WINDEPLOYQT_EXECUTABLE windeployqt HINTS "${QT_BIN_DIR}")

  # Run windeployqt after build
  add_custom_command(TARGET amcset-gui POST_BUILD
    COMMAND "${WINDEPLOYQT_EXECUTABLE}"
            --debug
            --no-translations
            "$<TARGET_FILE:amcset-gui>"
    COMMENT "Running windeployqt..."
  )
endif()
