cmake_minimum_required(VERSION 3.10)
project(reconstructionArmCurrentAnalysis)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(PythonLibs REQUIRED)

# 指定头文件搜索路径
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PYTHON_INCLUDE_DIRS}
)

# 设置可执行文件输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# 指定源文件
set(PROCESSDATA src/armCurrentDataAnalysis.cpp src/main.cpp)

# 添加测试程序的可执行文件
add_executable(processData ${PROCESSDATA})

# 如果您的项目依赖于Eigen库，并且Eigen库位于include目录下
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
    target_link_libraries(processData 
                        ${EIGEN3_LIBRARIES}                    
                        ${PYTHON_LIBRARIES}
                        python3.10.so
                        pthread
                        dl
    )
else()
    message(WARNING "Eigen3 not found. Please specify the path to Eigen3.")
endif()
