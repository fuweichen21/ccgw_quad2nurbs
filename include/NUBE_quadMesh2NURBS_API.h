/*
* (C++) FuweiChen
* 2025/6/21
* API:simple mesh lib test to DLL
*/

#ifndef _SIMPLE_QUADMESHTONURBS_API_H_
#define _SIMPLE_QUADMESHTONURBS_API_H_

#if defined(_WIN32) || defined(_WIN64)
#define CCGL_API __declspec(dllexport)
#else
#define CCGL_API __attribute__((visibility("default")))
#endif

#include <iostream>
extern "C"
{
	//CCGL_API bool simpleMeshLib(const std::string prefix_FilePath_Mesh);
	//CCGL_API bool QuadMeshToNURBS_API(const std::string prefix_FilePath_TriangleMesh, const std::string prefix_FilePath_QuadMesh, double _approxi_thread, double _approxi_smooth_weight, const bool printToscreen_Log, FILE* filePointer_Log);
	CCGL_API bool QuadMeshToNURBS_API(const std::string prefix_FilePath_TriangleMesh, const std::string prefix_FilePath_QuadMesh, double _approxi_thread, double _approxi_smooth_weight, const bool printToscreen_Log, FILE* filePointer_Log,const bool spline_surface_fitting);
}


#endif // !_SIMPLE_QUADMESHTONURBS_API_H_





