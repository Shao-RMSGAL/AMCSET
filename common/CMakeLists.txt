add_library(amcset-common STATIC
    ${PROJECT_SOURCE_DIR}/common/src/amcset_common.cpp
    ${PROJECT_SOURCE_DIR}/common/include/amcset_common.h
    ${PROJECT_SOURCE_DIR}/common/include/isotope_data.h
    ${PROJECT_SOURCE_DIR}/common/include/amcset_utilities.h
    ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    ${PROJECT_SOURCE_DIR}/common/assets/electron_screening_parameters.csv
    ${PROJECT_SOURCE_DIR}/common/assets/mott_scattering_parameters.csv
    )

# state that common need PIC when the default is shared libraries
set_target_properties(amcset-common PROPERTIES
                        POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
                        )


target_link_libraries(amcset-common PUBLIC compiler_flags)

# define the symbol stating we are using the declspec(dllexport) when
# building on windows
target_compile_definitions(amcset-common PRIVATE "EXPORTING_AMCSET")

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Add the custom command to generate the isotope database
add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    COMMAND
    ${Python3_EXECUTABLE}
    ${PROJECT_SOURCE_DIR}/common/scripts/generate_isotope_data.py
            ${PROJECT_SOURCE_DIR}/common/assets/isotope_data.txt
            ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    DEPENDS ${PROJECT_SOURCE_DIR}/common/assets/isotope_data.txt
            ${PROJECT_SOURCE_DIR}/common/scripts/generate_isotope_data.py
    COMMENT "Generating isotope database"
)


target_include_directories(amcset-common PUBLIC
    ${PROJECT_SOURCE_DIR}/common/include/
    ${PROJECT_SOURCE_DIR}/common/generated/
    )
