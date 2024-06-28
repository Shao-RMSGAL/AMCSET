add_library(amcset-common STATIC
    ${PROJECT_SOURCE_DIR}/common/src/amcset_common.cpp
    ${PROJECT_SOURCE_DIR}/common/include/amcset_common.h
    )

# state that common need PIC when the default is shared libraries
set_target_properties(amcset-common PROPERTIES
                        POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
                        )

target_include_directories(amcset-common PUBLIC
    ${PROJECT_SOURCE_DIR}/common/include/
    )

target_link_libraries(amcset-common PUBLIC compiler_flags)

# define the symbol stating we are using the declspec(dllexport) when
# building on windows
target_compile_definitions(amcset-common PRIVATE "EXPORTING_AMCSET")
