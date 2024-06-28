add_library(amcset_common STATIC
    ${PROJECT_SOURCE_DIR}/common/src/amcset_common.cpp
    ${PROJECT_SOURCE_DIR}/common/include/amcset_common.h
    )

# state that common need PIC when the default is shared libraries
set_target_properties(amcset_common PROPERTIES
                        POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
                        )

target_include_directories(amcset_common PUBLIC
    ${PROJECT_SOURCE_DIR}/common/include/
    )

target_link_libraries(amcset_common PUBLIC compiler_flags)

# define the symbol stating we are using the declspec(dllexport) when
# building on windows
target_compile_definitions(amcset_common PRIVATE "EXPORTING_MYMATH")
