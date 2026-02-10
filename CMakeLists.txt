# 指定CMake的最低版本要求
cmake_minimum_required(VERSION 3.10)

# 项目名称和语言设置
project(Quad2ManifoldSpline)

# 包含头文件目录
include_directories(
    ${PROJECT_SOURCE_DIR}/3rdParty/eigen-3.4.0
    ${PROJECT_SOURCE_DIR}/3rdParty/freeglut/include
    ${PROJECT_SOURCE_DIR}/3rdParty
    ${PROJECT_SOURCE_DIR}/include
)

# 平台特定的设置
if(WIN32)
    # Windows平台设置
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)  # 自动导出所有符号
    add_definitions(-DInitMeshLib_EXPORTS)
    
    # 指定静态库路径
    link_directories(${PROJECT_SOURCE_DIR}/3rdParty/freeglut/lib/x64)
    
    # MSVC编译器特定选项
    if(MSVC)
        add_compile_options(/bigobj)
    endif()
else()
    # Linux/macOS平台设置
    add_compile_options(-fPIC)  # 位置无关代码
    
    # 对于macOS可能需要指定freeglut路径
    if(APPLE)
        find_package(GLUT REQUIRED)
    endif()
endif()

# 收集源文件
set(SRCS
    include/NUBE_quadMesh2NURBS_API.h
    include/CCG_QMS_headerTool.h
    include/CCG_QMS_meshOperation.h
    include/NUBE_printInformation.h
    include/CCG_QMS_meshOperation.h
   include/CCG_QMS_meshTool.h
    include/CCG_QMS_quad2ManifoldSpline.h
    include/CCG_QMS_meshView.h
    include/BsplingCuvApprox.h
    include/CCG_QMS_bsplineTool.h
    src/main.cpp
    src/NUBE_quadMesh2NURBS_API.cpp
    src/CCG_QMS_meshOperation.cpp
    src/NUBE_printInformation.cpp
    src/CCG_QMS_meshOperation.cpp
    src/CCG_QMS_meshTool.cpp
    src/CCG_QMS_quad2ManifoldSpline.cpp
    src/CCG_QMS_meshView.cpp
    src/CCG_QMS_bsplineTool.cpp
    3rdParty/MeshLib/core/bmp/RgbImage.cpp
    src/BsplingCuvApprox.cpp
)

# 添加可执行文件
add_executable(${PROJECT_NAME} ${SRCS})

# 添加动态库
add_library(${PROJECT_NAME}_lib SHARED ${SRCS})

# 设置目标属性
set_target_properties(${PROJECT_NAME}_lib PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# 链接库
if(WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE freeglut)
    target_link_libraries(${PROJECT_NAME}_lib PRIVATE freeglut)
else()
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL GLUT::GLUT -lGLU)
    target_link_libraries(${PROJECT_NAME}_lib PRIVATE OpenGL::GL GLUT::GLUT -lGLU)
endif()

install(TARGETS ${PROJECT_NAME}_lib
    DESTINATION ${CMAKE_SOURCE_DIR}/bin
)
