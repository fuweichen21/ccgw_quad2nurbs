/*
* (C++) FuweiChen
* 2025/6/21
* A header for mesh to view
*/

#ifndef _CCG_QMS_VIEWER_H_
#define _CCG_QMS_VIEWER_H_

#include"CCG_QMS_meshTool.h"
#include"CCG_QMS_quad2ManifoldSpline.h"
#include"CCG_QMS_bsplineTool.h"

namespace CCG_QMSLib
{
	typedef CCG_QMS_Vertex	V;
	typedef CCG_QMS_Edge	E;
	typedef CCG_QMS_Face	F;
	typedef CCG_QMS_HalfEdge H;
	typedef CCG_QMS_Mesh<V, E, F, H> M;
	typedef typename M::It			It;

	/*view mesh, texture*/
	void viewMesh(M* pMesh, int argc, char* argv[]);

	/*
	* view mesh
	*/
	void viewMesh(M* pMesh);

	/*view tri mesh and quad mesh*/
	void viewMesh(M* triMesh, M* quadMesh, int argc, char* argv[]);

	/*view some scalar on quad Mesh*/
	void viewScalarOnMeshFace(M* pMesh, int argc, char* argv[]);

	/*view sampling on quad*/
	void viewSampling(M* pMesh, std::vector<QB_sampling> sts);
}


#endif // !_CCG_QMS_VIEWER_H_


