/*
* (C++) FuweiChen
* The mesh header for mesh basic operation
* 2025/6/21
*/

#ifndef _CCG_QMS_MESHTOOL_H_
#define _CCG_QMS_MESHTOOL_H_

#include"CCG_QMS_headerTool.h"
#include <unordered_map>

using namespace MeshLib;

namespace CCG_QMSLib
{
	class CCG_QMS_Vertex;
	class CCG_QMS_dge;
	class CCG_QMS_HalfEdge;
	class CCG_QMS_Face;
	
	class CCG_QMS_Vertex :public MeshLib::CVertex
	{
	public:
		CCG_QMS_Vertex() {
			//Attribute to quad mesh to manifold spline
			/*----------------------------------------*/
			m_feature = false;
			m_degree = 0;
			m_ifSingular = false;
			m_ifCorner = false;
			m_originalId = -1;
			m_ifSharp = false;
			m_ifVisit = false;
			m_boundaryFaceId = -1;
			m_boundaryTargetVId = -1;
			/*----------------------------------------*/

			//Attribute to boundary curve fitting
			/*-----------------------------------------------------------------------------------------------*/
			m_guassCur = 0.0;
			/*-----------------------------------------------------------------------------------------------*/

			//Attribute to fitting
			/*-----------------------------------------------------------------------------------------------*/
			m_ifFixed = false;
			/*-----------------------------------------------------------------------------------------------*/
		};
		~CCG_QMS_Vertex() {};

		void _from_string();
		void _to_string();

		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		bool& feature() { return m_feature; }
		int& degree() { return m_degree; }
		bool& ifSingular() { return m_ifSingular; }
		bool& ifCorner() { return m_ifCorner; }
		int& originalId() { return m_originalId; }
		bool& ifSharp() { return m_ifSharp; }
		bool& ifVisit() { return m_ifVisit; }
		int& boundaryFaceId() { return m_boundaryFaceId; }
		int& boundaryTargetVId() { return m_boundaryTargetVId; }
		/*----------------------------------------*/

		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		double& gaussCur() { return m_guassCur; }
		CPoint& initialPos() { return m_initialPos; }
		/*-----------------------------------------------------------------------------------------------*/

		//Attribute to fitting
		/*-----------------------------------------------------------------------------------------------*/
		bool& fixed() { return m_ifFixed; }
		/*-----------------------------------------------------------------------------------------------*/

	private:
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		/*Marking feature*/
		bool m_feature;
		/*The number of linking faces*/
		int m_degree;
		/*Is it a singularity? ture or false*/
		bool m_ifSingular;
		/*Is it a corner (degree is 1)? ture or false*/
		bool m_ifCorner;
		/*Record the IDs of the original mesh points corresponding to the simplified triangular mesh points*/
		int m_originalId;
		bool m_ifSharp;
		bool m_ifVisit;
		/*Record triangular mesh bouondary vertex Id mapping to the face id from quad mesh*/
		int m_boundaryFaceId; 
		/*Record triangular mesh bouondary vertex Id mapping to the halfedge target id from quad mesh*/
		int m_boundaryTargetVId;
		/*----------------------------------------*/

		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		double m_guassCur;
		CPoint m_initialPos;//copy initial position
		/*-----------------------------------------------------------------------------------------------*/

		//Attribute to fitting
		/*-----------------------------------------------------------------------------------------------*/
		bool m_ifFixed;
		/*-----------------------------------------------------------------------------------------------*/

	};

	class CCG_QMS_Edge :public MeshLib::CEdge
	{
	public:
		CCG_QMS_Edge() {
			//Attribute to quad mesh to manifold spline
			/*----------------------------------------*/
			m_sharp = false;
			m_constrinedBoundary = false;
			m_feature = false;
			m_tLayout = false;
			m_knotInterval = 1.0;
			/*----------------------------------------*/
			//Attribute to boundary curve fitting
			/*-----------------------------------------------------------------------------------------------*/
			m_visit = false;
			/*-----------------------------------------------------------------------------------------------*/
		};
		~CCG_QMS_Edge() {};
		void _from_string();
		void _to_string();
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		bool& sharp() { return m_sharp; }
		bool& constrinedBoundary() { return m_constrinedBoundary; }
		bool& feature() { return m_feature; }
		bool& tLayout() { return m_tLayout; }
		double& knotInterval() { return m_knotInterval; }
		/*----------------------------------------*/
		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		bool& visit() { return m_visit; };
		/*-----------------------------------------------------------------------------------------------*/
	protected:
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		/*Marking edges*/
		bool m_sharp;
		/*constrianed patch boundary edge*/
		bool m_constrinedBoundary;
		/*Marking feature edges*/
		bool m_feature;
		/*Marking T-layout edges*/
		bool m_tLayout;
		/*
		* Define the knot intervals along the edges of the mesh,
		* and it is required that the node interval values for the opposite sides of each quadrilateral be equal
		*/
		double m_knotInterval;
		/*----------------------------------------*/

		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		bool m_visit;
		/*-----------------------------------------------------------------------------------------------*/
	};

	class CCG_QMS_Face :public MeshLib::CFace
	{
	public:
		CCG_QMS_Face() {
			//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
			m_patchIndex = -1;
			m_subPatchIndex = -1;
			m_f_spline_patchId = 0;
			m_localId = -1;
			m_tempMark = false;
			m_f_patchId = 0;
		/*----------------------------------------*/
		/*--------------Quality-------------------------*/
			m_quality_scaledJocbian = 1.0;
			m_quality_edgeRatio = 1.0;
			m_scalar = 1.0;
		/*----------------------------------------*/
		};
		~CCG_QMS_Face() {};
		CPoint& normal() { return m_normal; }
		void _from_string();
		void _to_string();
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		int& patchIndex() { return m_patchIndex; }
		int& subPatchIndex() { return m_subPatchIndex; }
		int& f_spline_patchId() { return m_f_spline_patchId; }
		int& localId() { return m_localId; }
		bool& tempMark() { return m_tempMark; }
		int& f_patchId() { return m_f_patchId; }
		/*----------------------------------------*/
		/*--------------Quality-------------------------*/
		double& quality_scaledJocbian() { return m_quality_scaledJocbian; }
		double& quality_edgeRatio() { return m_quality_edgeRatio; }
		double& scalar() { return m_scalar; }
		/*----------------------------------------*/

	protected:
		CPoint m_normal;
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		/*The piecewise index of the initial triangular mesh on quad mesh, starting from 0*/
		int m_patchIndex;
		/*The piecewise index of the initial triangular mesh patch, starting from 0*/
		int m_subPatchIndex;
		/*Indicates the patch of the spline where the face is located*/
		int m_f_spline_patchId;
		/*Establish the association between the Bezier surface patch and the B-spline surface*/
		int m_localId;
		/*temp marking*/
		bool m_tempMark;
		int m_f_patchId;//表示当前面所在的patch
		/*----------------------------------------*/
		/*--------------Quality-------------------------*/
		double m_quality_scaledJocbian;
		double m_quality_edgeRatio;
		double m_scalar;
		/*----------------------------------------*/
	};

	class CCG_QMS_HalfEdge :public MeshLib::CHalfEdge
	{
	public:
		CCG_QMS_HalfEdge() {
			//Attribute to quad mesh to manifold spline
			/*----------------------------------------*/
			m_localId = 0;
			m_boundaryFaceId = -1;
			m_boundaryTargetVId = -1;
			m_id_FeatureTrajectory = -1;
			m_NURBS_patchId = -1;
			m_isoParameterLineType = -1;
			/*----------------------------------------*/

			//Attribute to boundary curve fitting
			/*-----------------------------------------------------------------------------------------------*/
			m_angleCos = 0.0;
			/*-----------------------------------------------------------------------------------------------*/

			//Attribute to quantify the surface quality
		    /*-----------------------------------------------------------------------------------------------*/
			m_visit = false;
		    /*-----------------------------------------------------------------------------------------------*/
		};
		~CCG_QMS_HalfEdge() {};
		void _from_string();
		void _to_string();
		CPoint2& uv() { return m_uv; }
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		int& localId() { return m_localId; }
		int& boundaryFaceId() { return m_boundaryFaceId; }
		int& boundaryTargetVId() { return m_boundaryTargetVId; }
		int& id_FeatureTrajectory() { return m_id_FeatureTrajectory; }
		int& NURBS_patchId() { return m_NURBS_patchId; }
		int& isoParameterLineType() { return m_isoParameterLineType; }
		std::vector<int>& originalIds() { return m_originalIds; }
		/*----------------------------------------*/

		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		double& angleCos() { return m_angleCos; }
		/*-----------------------------------------------------------------------------------------------*/

		//Attribute to quantify the surface quality
		/*-----------------------------------------------------------------------------------------------*/
		bool& visit() { return m_visit; };
		/*-----------------------------------------------------------------------------------------------*/

	private:
		CPoint2 m_uv; //uv
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		/*
		* Local id, 1,2,3,4
		* Mappint 1 to the first control point
		* 2 to the fourth control point
		* 3 to the sixteenth control point
		* 4 to the thirteenth control point
		*/
		int m_localId;
		/*Record the face ID on the quadrilateral mesh corresponding to the halfedge of the triangular mesh boundary*/
		int m_boundaryFaceId;
		/*Record the target ID corresponding to the half edge of the triangular mesh boundary and the half edge of the quadrilateral mesh boundary*/
		int m_boundaryTargetVId;
		/*The trajectory's id of feature on the triangular mesh*/
		int  m_id_FeatureTrajectory;
		/*Record the NURBS patch ID where the halfedge of the quadrilateral mesh boundary is located*/
		int m_NURBS_patchId;
		/*
		* Record the parametric lines of the NURBS patch 
		* (represented in integer form: 1 indicates u = u_min; 2 indicates v = v_min; 3 indicates u = u_max; 4 indicates v = v_max)
		*/
		int m_isoParameterLineType;
		/*Store the originalIds on the halfedge of the boundary of the quadrilateral mesh*/
		std::vector<int> m_originalIds;
		/*----------------------------------------*/

		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		double m_angleCos;
		/*-----------------------------------------------------------------------------------------------*/

		//Attribute to quantify the surface quality
		/*-----------------------------------------------------------------------------------------------*/
		bool m_visit;
		/*-----------------------------------------------------------------------------------------------*/
	};

	template<typename VV, typename EE, typename FF, typename HH>
	class CCG_QMS_Mesh :public MeshLib::CDynamicMesh<VV, EE, FF, HH>
	{
	public:
		typedef CCG_QMS_Mesh<VV, EE, FF, HH>			M;
		typedef MyLib::My_MeshIterator<M>	It;
		typedef VV							V;
		typedef EE							E;
		typedef FF							F;
		typedef HH							HE;

	public:
		CCG_QMS_Mesh() {
			//Attribute to quad mesh to manifold spline
			/*----------------------------------------*/
			m_tLayout_num = 0;
			/*----------------------------------------*/
		};
		~CCG_QMS_Mesh() {};

		/**
		* @brief Build Index: Classify triangle faces according to Patch and SubPatch indices.
		* @return A nested hash map containing the classified faces.
		*/
		void classifyFacesByPatch();

	public:
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		int& tLayout_num() { return m_tLayout_num; }
		std::map<int, int>& tLayoutPatches_quadNum() { return m_tLayoutPatches_quadNum; }
		std::unordered_map<int, std::unordered_map<int, std::vector<F*>>>& patchMap_mesh() { return m_patchMap_mesh; }
		void updateMapVertex();/*update m_map_vert*/
		/*----------------------------------------*/
	protected:
		//Attribute to quad mesh to manifold spline
		/*----------------------------------------*/
		//The number of TLayout patches
		int m_tLayout_num;
		//The number of quad mesh faces on each TLayout patch
		std::map<int, int> m_tLayoutPatches_quadNum;	
		//PatchMap: Patch Index -> SubPatch Index -> List of mesh faces
		std::unordered_map<int, std::unordered_map<int, std::vector<F*>>> m_patchMap_mesh;
		/*----------------------------------------*/
	};
	template<typename VV, typename EE, typename FF, typename HH>
	inline void CCG_QMS_Mesh<VV, EE, FF, HH>::classifyFacesByPatch()
	{
		m_patchMap_mesh.clear();
		for (auto f : It::MFIterator(this))
		{
			int pIdx = f->patchIndex();
			int sIdx = f->subPatchIndex();
			m_patchMap_mesh[pIdx][sIdx].push_back(f);
		}
	}
	template<typename VV, typename EE, typename FF, typename HH>
	inline void CCG_QMS_Mesh<VV, EE, FF, HH>::updateMapVertex()
	{
		m_map_vert.clear();
		for (auto v : this->vertices())
		{
			m_map_vert.insert(std::pair<int, V*>(v->id(), v));
		}
	}
}

#endif // !_CCG_QMS_MESHTOOL_H_


