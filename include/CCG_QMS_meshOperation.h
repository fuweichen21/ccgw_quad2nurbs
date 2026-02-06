/*****************************************************************//**
 * \file   CCG_QMS_meshOperation.h
 * \brief  The header file for some mesh operation
 * 
 * \author FuweiChen
 * \date   2025\6\21
 *********************************************************************/

#ifndef _CCG_QMS_MESHOPERATION_HEADER_H_
#define _CCG_QMS_MESHOPERATION_HEADER_H_

#include"CCG_QMS_meshTool.h"
#include"CCG_QMS_quad2ManifoldSpline.h"

namespace CCG_QMSLib
{
	typedef CCG_QMS_Vertex	V;
	typedef CCG_QMS_Edge	E;
	typedef CCG_QMS_Face	F;
	typedef CCG_QMS_HalfEdge H;
	typedef CCG_QMS_Mesh<V, E, F, H> M;
	typedef typename M::It			It;

	//Attribute to quad mesh to manifold spline
	/*-----------------------------------------------------------------------------------------------*/
	/*
	* Mark the feature edges by sharp edges
	*/
	void markFeatureEdgesBySharpEdges(M* pMesh);

	void arrangeMeshVertexOrder(M* pMesh);

	/*
	* Enlarge the coordinates by
	* a certain order of magnitude
	*/
	void magnifyMeshPoint(M* pMesh, double scale);

	/*
	* Marking T-Layout edges on mesh
	* T-layout patch information is required
	* on input file
	*/
	void markTLayout_quadMesh(M* quadMesh);

	/*
	* statistic the number of TLayout patches
	* on quad mesh
	*/
	int statisticTLayoutPatchNum(M* quadMesh);
	/*
	* statistic the number of quad mesh faces on each TLayout
	*/
	void statisticTLayoutPatch_quadNum(M* quadMesh);

	void compute_mesh_normal(M* pMesh);

	/*
	* Marking extraordinary points
	* especially, the boundary (degree is 1) is not 
	* regarded as the singularity
	*/
	void markExtraordinaryPts_cw(M* pMesh);

	/*
	* Mark all singularities in the mesh
	*/
	void markExtraordinaryPts(M* pMesh);

	/*Marking halfedges' local id, making the mapping between control points from the bezier surface and local id*/
	/*
	* ///////////////////////////////////////////
		*     |                     |
		*     |                     |
		*     |        [2]          |
		* ---(4)--------------------(3)---
		*     |  10       11        |
		*  [3]|                     |[1]
		*     |  6        7         |
		* ---(1)--------------------(2)---
		*     |       [0]           |
		*     |                     |
		*     |                     |
		* ///////////////////////////////////////////
		* (1)(2)(3)(4) the order of quad mesh face's vertex, sequence from the input file
		* [[0][1][2][3]halfedge local id
		* [0]'s source is£¨1£©
		* 6£¬7£¬10£¬11 are face control points
		*/
	void markHalfedgeLocalId(M* pMesh);

	/*Marking corner (degree is 1) at quad mesh, */
	void markCorner_quadMesh(M* pMesh);

	/*Marking corner (degree is 1) and feature corner at quad mesh, */
	void markCorner_feature_quadMesh(M* pMesh);

	/*Avoid feature vertex to be the singular vertex*/
	void remarkSingularPts(M* pMesh);

	/*Computing the degree of vertex on the mesh*/
	void computeV_degree_mesh(M* pMesh);

	/*Enlarge the parameters 'uv' in the mesh by a certain magnitude*/
	void magnifyMeshUV(M* pMesh, double scale);

	/*
	* Establish the relationship between the boundary points of the triangular mesh 
	* and the boundary edges of the quadrilateral mesh based on the parameterized results.
	* Establish the relationship directly based on the "trajectory" on the "halfedges".
	* The considered object is the boundary halfedge on the triangular mesh
	*/
	void createTriQuadBoundaryVertexRealationByParameterization4(M* triMesh, M* quadMesh);
	

	/*
		* Obtaining the T layout based on the quadrilateral mesh
		* Using the idea of obtaining the motorcycle graph
		* To minimize the number of surfaces obtained
		* Not considering using the feature lines as the boundaries of the patches
		*/
	void obtainTLayout_motorcycleGraph_featureFree(M* pMesh);

	/*
		* Obtaining the T layout based on the quadrilateral mesh
		* Using the idea of obtaining the motorcycle graph
		* To minimize the number of surfaces obtained
		* Considering using the feature lines as the boundaries of the patches
		*/
	void obtainTLayout_motorcycleGraph(M* pMesh);

	/*
		* Obtaining the T layout based on the quadrilateral mesh
		* Using the idea of obtaining the motorcycle graph
		* To minimize the number of surfaces obtained
		* Considering using the feature lines as the boundaries of the patches
		* suppose that the feature lines constructs a quad layout
		*/
	void obtainTLayout_motorcycleGraph_featureAware(M* pMesh);

	/*
		* Obtain the layout based on quadrilateral mesh
		* Make the number of patches consistent with the number of quadrilateral meshes.
		*/
	void obtainLayout_consistencyQuadMeshFaceNum(M* pMesh);
	/*-----------------------------------------------------------------------------------------------*/


	//Attribute to boundary curve fitting
	/*-----------------------------------------------------------------------------------------------*/
	void copyInitialMeshVPos(M* pMesh);
	/*Compute the cos of every angles on the mesh*/
	void computeCos(M* pMesh);
	/*Compute the Gauss curvature for every vertex*/
	void computeGaussCur(M* pMesh);
	/*Restore the quadrilateral mesh verrtex coordinates of the flat areas at the boundary to their original positions*/
	void recoverFlatBoundaryVertexPoint(M* pMesh);
	/*-----------------------------------------------------------------------------------------------*/

	/*-----------------------------------------------------------------------------------------------*/
	/*
	* By using the parametric results, the triangular grid points are located on the quadrilateral grid, 
	* and the corresponding sampling point information is obtained.
	*/
	void obtainSamplingFromTriMeshByParameterization(M* triMesh, M* quadMesh, std::vector<QB_sampling>& sts);

	/*Assign the fitted results to the initial quadrilateral mesh points*/
	void updateMeshPosByControlMesh(M* pMesh, CCG_QMS_model& model);

	/*
	* Based on the parametric results of quadrilateral and triangular meshes,
	* identify the quadrilaterals whose four vertices' parameters are within the parametric range of the same triangle.
	*/
	std::vector<F*> FindQuadsInsideTris_Brute_WithBuckets(M* quadMesh, M* triMesh, double eps);
	/*-----------------------------------------------------------------------------------------------*/

	/*------------------------------------------Quad mesh qulity calculation------------------------------------------*/
	/*
	* brief Calculates the Scaled Jacobian for quadrilateral face.
	*/
	void calculateScaledJacobian_quadMesh(M* pMesh);
	/*
	* brief Calculates the Edge Ratio for quadrilateral face.
	*/
	void calculateEdgeRatio_quadMesh(M* pMesh);
	/*--------------------------------------------------------------------------------------------------------------*/
}


#endif // !_CCG_QMS_MESHOPERATION_HEADER_H_

