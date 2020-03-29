cmake_minimum_required(VERSION 3.5)

project(app)

find_package(Json REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(${PROJECT_NAME}
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    guiLib simulation
)


