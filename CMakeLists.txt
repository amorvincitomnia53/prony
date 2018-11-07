cmake_minimum_required(VERSION 3.12)
project(enshu)

set(CMAKE_CXX_STANDARD 17)


add_compile_options(-Wall -Wextra)



find_package(PkgConfig)
pkg_check_modules(gtkmm IMPORTED_TARGET gtkmm-3.0 )
pkg_check_modules(epoxy IMPORTED_TARGET epoxy)

find_package(Boost)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/oglplus)

#add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/prony.glade
#        COMMAND  cp ${CMAKE_CURRENT_LIST_DIR}/prony.glade ${CMAKE_CURRENT_BINARY_DIR}/
#        DEPENDS ${CMAKE_CURRENT_LIST_DIR}/prony.glade
#        )


add_executable(prony main.cpp ) #${CMAKE_CURRENT_BINARY_DIR}/prony.glade
target_link_libraries(prony PkgConfig::gtkmm PkgConfig::epoxy Boost::boost oglplus)
