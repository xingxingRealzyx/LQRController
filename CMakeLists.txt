cmake_minimum_required(VERSION 3.10)
project(lqr_controller VERSION 1.0)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 添加include目录
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/3rdparty/NanoMatrix
)

# 构建测试程序
add_executable(test_lqr test_lqr.cpp)

# 如果有第三方依赖（比如线性代数库），需要在这里添加
# 例如：
# find_package(Eigen3 REQUIRED)
# target_link_libraries(test_lqr Eigen3::Eigen)

# 设置编译选项
target_compile_options(test_lqr PRIVATE 
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# 安装规则（可选）
install(TARGETS test_lqr
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
) 