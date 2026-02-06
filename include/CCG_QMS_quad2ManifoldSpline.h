/*****************************************************************
 * \file   CCG_QMS_quad2ManifoldSpline.h
 * \brief  
 * A header for quad to manifold spline
 * \author FuweiChen
 * \date   2025/6/28
 *********************************************************************/

#ifndef _CCG_QMS_QUAD2MANIFOLDSPLINE_H_
#define _CCG_QMS_QUAD2MANIFOLDSPLINE_H_

#include"CCG_QMS_meshTool.h"
#include<memory>

namespace CCG_QMSLib{

	/*
	* Rename some classes of mesh, easy to use
	*/
	typedef CCG_QMS_Vertex	V;
	typedef CCG_QMS_Edge	E;
	typedef CCG_QMS_Face	F;
	typedef CCG_QMS_HalfEdge H;
	typedef CCG_QMS_Mesh<V, E, F, H> M;
	typedef typename M::It			It;

	/*Declare some classes*/
	/*Integral model class*/
	class  CCG_QMS_model;
	/*Bezier surface class*/
	class CCG_QMS_bezierSurf;
	/*B spline surface class*/
	class CCG_QMS_bsplineSurf;
	/*Boundary of one B spline surface*/
	class NUBE_bsurfBoudnary;
	/*Segment of one B spline surface boundary*/
	class NUBE_bsurfSegment;
	/*Sampling point*/
	class QB_sampling;

	/*Bezier surface type list*/
	enum BezierSurfType
	{
		/*All vertex on face is regular*/
		RegularSurf = 0,
		/*There is at least one vertex  on face is sigular*/
		SingularSurf = 1,
		/*All vertex on face is regular, but linking singularity*/
		LinkingSingularSurf = 2
	};

	struct BE_linkingCtlPoints
	{
		std::map<int, double> linkingCtlPId_weights;

		/*Operator'-'*/
		BE_linkingCtlPoints operator - (BE_linkingCtlPoints& bts)
		{
			BE_linkingCtlPoints tempCW;

			for (auto cwIt : linkingCtlPId_weights)
			{
				std::pair<int, double> tempV = cwIt;
				tempV.second *= 1.0;

				if (tempCW.linkingCtlPId_weights.find(tempV.first) == tempCW.linkingCtlPId_weights.end())
				{
					tempCW.linkingCtlPId_weights.insert(tempV);
				}
				else
				{
					std::map<int, double>::iterator tIt = tempCW.linkingCtlPId_weights.find(tempV.first);
					tIt->second += (tempV.second);
				}
			}
			
			for (auto cwIt : bts.linkingCtlPId_weights)
			{
				std::pair<int, double> tempV = cwIt;
				tempV.second *= -1.0;

				if (tempCW.linkingCtlPId_weights.find(tempV.first) == tempCW.linkingCtlPId_weights.end())
				{
					tempCW.linkingCtlPId_weights.insert(tempV);
				}
				else
				{
					std::map<int, double>::iterator tIt = tempCW.linkingCtlPId_weights.find(tempV.first);
					tIt->second += (tempV.second);
				}
			}
			return tempCW;
		}

		/*Operator'+'*/
		BE_linkingCtlPoints operator + (BE_linkingCtlPoints& bts)
		{
			BE_linkingCtlPoints tempCW;

			for (auto cwIt : linkingCtlPId_weights)
			{
				std::pair<int, double> tempV = cwIt;
				tempV.second *= 1.0;

				if (tempCW.linkingCtlPId_weights.find(tempV.first) == tempCW.linkingCtlPId_weights.end())
				{
					tempCW.linkingCtlPId_weights.insert(tempV);
				}
				else
				{
					std::map<int, double>::iterator tIt = tempCW.linkingCtlPId_weights.find(tempV.first);
					tIt->second += (tempV.second);
				}
			}
			
			for (auto cwIt : bts.linkingCtlPId_weights)
			{
				std::pair<int, double> tempV = cwIt;
				tempV.second *= 1.0;

				if (tempCW.linkingCtlPId_weights.find(tempV.first) == tempCW.linkingCtlPId_weights.end())
				{
					tempCW.linkingCtlPId_weights.insert(tempV);
				}
				else
				{
					std::map<int, double>::iterator tIt = tempCW.linkingCtlPId_weights.find(tempV.first);
					tIt->second += (tempV.second);
				}
			}
			return tempCW;
		}

		/*Operator'*'scalar*/
		BE_linkingCtlPoints operator * (double r)
		{
			BE_linkingCtlPoints tempCW;

			for (auto cwIt : linkingCtlPId_weights)
			{
				std::pair<int, double> tempV = cwIt;
				tempV.second *= r;

				if (tempCW.linkingCtlPId_weights.find(tempV.first) == tempCW.linkingCtlPId_weights.end())
				{
					tempCW.linkingCtlPId_weights.insert(tempV);
				}
				else
				{
					std::map<int, double>::iterator tIt = tempCW.linkingCtlPId_weights.find(tempV.first);
					tIt->second += (tempV.second);
				}
			}
			return tempCW;
		}

	};

	struct NURBS_boundaryCurve
	{
		int patchId;//nurbs patch id
		int curveType;//( 1: u=u_min; 2: v = v_min; 3: u = u_max; 4:  v=v_max )
		std::vector<int> originalTriMeshVIds;//record the original tri mesh vetex id
		std::vector<int> triMeshVIds;//record the tri mesh vetex id
	};

	// Sampling points of boundary curve
	struct BoundaryCurveSample {
		CPoint2 uvParam;
		CPoint point;
		int curveIndex;
	};

	/*Create the class for the declaration*/
	class CCG_QMS_model
	{
	public:
		CCG_QMS_model(){
			m_diffG0_sampling = NULL;
			m_diffG1_sampling = NULL;
			m_diff_duu_sampling = NULL;
			m_diff_duv_sampling = NULL;
			m_diff_dvv_sampling = NULL;
			m_boundaryEdgeRatio = 0.0;
			m_minCtrlPtDist = 0.0;
			m_maxHausdorff = 0.0;
		};
		~CCG_QMS_model(){
			// Release memory for bezier surfaces
			/*std::cout << "0:" << m_bfs.capacity() << std::endl;
			std::cout << "0:" << m_bfs.size() << std::endl;*/
			/*for (auto &ptr : m_bfs) {
				delete ptr;
			}*/
			/*std::cout << "1:" << m_bfs.capacity() << std::endl;
			std::cout << "1:" << m_bfs.size() << std::endl;*/
			/*m_bfs.clear();
			m_bfs.shrink_to_fit();*/
			/*std::cout << "2:" << m_bfs.capacity() << std::endl;
			std::cout << "2:" << m_bfs.size() << std::endl;*/
			
			//// Release memory for bspline surfaces
			//for (auto ptr : m_bsplineSurfs) {
			//	delete []ptr;
			//}
			//m_bsplineSurfs.clear();

			/*for (typename std::list<CCG_QMS_bezierSurf*>::iterator biter = m_bfs_test.begin(); biter != m_bfs_test.end(); ++biter)
			{
				CCG_QMS_bezierSurf* bf = *biter;
				delete bf;
			}
			m_bfs_test.clear();*/
		};
		/*
		* obtain the number of control points from quad mesh,
		* equal to the number of vertex on quad mesh
		*/
		void obtainCtlPts(M* pMesh);
		
		/*
		* build bezier surface by leveraging the topological and geometric information of the mesh
		* reference to 《Isogeometric analysis using G-spline surfaces with arbitrary
		* unstructured quadrilateral layout》
		*/
		void buildBezierSurfs_cw(M* pMesh);

		/*
		* Raise the degree of all Bezier surfaces in the model
		*/
		void raiseBezierSurfDegree3To5_all(M* pMesh);

		/*Sort all bezier surface according to their ids*/
		void sortBezierSurfaces();

		/*
		* Adjust to the sequence of halfedges to adapt the construction of basis function using
		* 《Isogeometric analysis using G-spline surfaces with arbitrary
		* unstructured quadrilateral layout》
		*/
		void sortBezierCtlPts_cw(M* pMesh);

		/*
		* G1 constriants,update control points，equation（38）on
		* 《Isogeometric analysis using G-spline surfaces with arbitrary
		* unstructured quadrilateral layout》
		*/
		void G1Constraints_optimization_check(M* pMesh);

		/*
		* G1 constriants,update control points，equation（38）on
		* Isogeometric analysis using G-spline surfaces with arbitrary
		* unstructured quadrilateral layout》
		* keep the relationship between the Bezier control points and the quadrilateral vertex points
		*/
		void G1Constraints_optimization_cws(M* pMesh);

		/*Obtain the knot vector of the Bezier surface*/
		void computeBezierSurfacesKnots(M* pMesh);

		/*
		* Sort the control points in the Bezier surface in a counterclockwise direction based on the first point in the quadrilateral mesh file.
		*/
		void sortBezierCtlPts(M* pMesh);

		/*
		* By using the control point ID and weights, the control points on each Bezier surface are calculated.
		* and scale down to the initial size of the model
		*/
		void obtainBezierCptsBy_cpts(double scale);

		/**
		 * Obtain B-spline surface patches using the topology of a quadrilateral mesh.
		 * It can handle special cases where there are no singular points in the mesh.
		 * The direction of the patch is consistent with the direction of the quadrilateral mesh.
		 * Provide the IDs of the four vertices of the B-spline surface corresponding to the four grid points of the quadrilateral.
		 * In order to create the quadrilateral mesh corresponding to the B-spline surface patch, 
		 * it is necessary to ensure that the id of each half edge on the surface is in the correct order:（0,0）--2-->(1,0) --3-->(1,1) --4-->(0,1)--1-->(0,0)
		 * Obtain the boundary edge information of each B-spline surface
		 */
		void obtainBSplineSurfacePatches_general_UVConsistency_topology_boundary(M* pMesh);

		/**
		 * Obtain B-spline surface patches using the topology of a quadrilateral mesh.
		 * It can handle special cases where there are no singular points in the mesh.
		 * The direction of the patch is consistent with the direction of the quadrilateral mesh.
		 * Record the NURBS patch ID where each halfedge of the quadrilateral mesh boundary is located and
		 * isoparametric lines of the NURBS patch（express in integer form：1---u=u_min; 2---v = v_min; 3---u = u_max; 4---v=v_max）
		 */
		void obtainBSplineSurfacePatches_general_UVConsistency_recordQuadBoundaryNURBSBoundaryCurveRelation(M* pMesh);

		/*Remove the duplicate knots inside the B-spline surface*/
		void removeMultiKnots_BSplines();

		/*Calculate the values of the non-repeating knots in the B-spline surface and the corresponding multiplicities.*/
		void computeBSplineSurfaceUniqueKnotMultiNum();

		/**
		 * Export the B-spline surface in the.iges format file.
		 */
		void outputBSpline_iges(const char* output);

		/**
		* Export the B-spline surface in the.iges format file.
		* version 5.3
		*/
		void outputBSpline_iges_5_3(const char* output);

		/*
		* Output the B-spline surface to the following fixed format
		* ============================================
		* Number_of_unique_knots_u (int),
		* Number_of_unique_knots_v (int),
		* Number_of_ctrl_pts_u (int),
		* Number_of_ctrl_pts_v (int),
		* Unique_knots_u vector (doubles),
		* Unique_knots_v vector (doubles),
		* Multiplicity_knots_u vector (ints),
		* Multiplicity_knots_v vector (ints),
		* Control points grid (doubles), varying faster along U.
		* Control points’ weights. (Optional if all being 1)
		* ============================================
		* Now let me spell it out, the STANDARD knots vector of
		* {0,0,0,0,0.5,1,1,1,1} as in text books will be represented
		* by a UNIQUE knot vector {0,0.5,1} and an accompanying knot
		* MULTIPLICITY vector of the same size: {4,1,4}.
		* The degrees will be inferred from the size of the standard
		* knot vector and number of control points along both directions.
		* Look at an example below:
		* ================================
		* 4;                      // Number_of_unique_knots_u (int),
		* 3;                      // Number_of_unique_knots_v (int),
		* 6;                      // Number_of_ctrl_pts_u (int),
		* 4;                      // Number_of_ctrl_pts_v (int),
		* 0.0, 0.33, 0.67, 1.0;   // Unique_knots_u vector (doubles),
		* 0.0, 0.5, 1.0;          // Unique_knots_v vector (doubles),
		* 4,1,1,4;                // Multiplicity_knots_u vector (ints),
		* 3,1,3;                  // Multiplicity_knots_v vector (ints),
		* 0.0,  0.0, 0.0;   // row 0 along u. Control points grid (doubles)
		* 5.0,  0.0, 1.0;
		* 20.0,  0.0, 2.0;
		* 30.0,  0.0, 2.0;
		* 45.0,  0.0, 1.0;
		* 50.0,  0.0, 0.0;
		* 0.0,  5.0, 1.0;   // row 1 along u
		* 5.0,  5.0, 2.0;
		* 20.0,  5.0, 3.0;
		* 30.0,  5.0, 3.0;
		* 45.0,  5.0, 2.0;
		* 50.0,  5.0, 1.0;
		* 0.0, 25.0, 1.0;   // row 2 along u
		* 5.0, 25.0, 2.0;
		* 20.0, 25.0, 3.0;
		* 30.0, 25.0, 3.0;
		* 45.0, 25.0, 2.0;
		* 50.0, 25.0, 1.0;
		* 0.0, 30.0, 0.0;   // row 3 along u
		* 5.0, 30.0, 1.0;
		* 20.0, 30.0, 2.0;
		* 30.0, 30.0, 2.0;
		* 45.0, 30.0, 1.0;
		* 50.0, 30.0, 0.0;
		* ================================
		*/
		void outputBSplineSurf_fixedFormat3(const char* output);//Regarding the updated output after the removal of internal duplicate knots, and update the periodic and closed property

		/*
		* Read simplifiedId and initialId
		* Obatain NURBS patch bouondary curve's triangular mesh vertex, and output to .txt file
		* .txt format：
		* Total_NURBS_BoundaryCurve_number 37. (BCV: BoundaryCurveVertex; IMBV: Initial Input Mesh Boundary Vertex)
		* BoundaryCurve[0]: Patch 1, Side 2. (BCV# = 2)
		* BCV[0] IMBV 25 (Collapsed# = 2):
		*  25, 17;
		* BCV[1] IMBV 36 (Collapsed# = 1):
		*  36;
		* BoundaryCurve[1]: Patch 1, Side 3. (BCV# = 4)
		* .
		* .
		* .
		* Taking into account that the common points on the boundary line will be used multiple times
		*/
		void obtainNURBS_boundaryCurve_output3(M* triMesh, M* quadMesh, const char* input, const char* output);

		/*
		* Obtaining nurbs boundary curve information by quad mesh boundary edge
		*/
		void obtainNURBS_boundaryCurve(M* quadMesh);

		/*
		*计算未简化三角网格边界点到整个NURBS边界曲线的Hausdorff距离
		*
		* @param model BE_model对象，包含所有NURBS曲面信息
		* @param objTriMesh 原始三角网格
		* @param sampleStep 在参数域上的采样步长
		* @param sortedDistances 输出参数，存储排序后的点与其最小距离
		* @return double Hausdorff距离值
		*/
		double computeHausdorffDistance_allNURBSBoundaryCurve(M* objTriMesh, double sampleStep, std::vector<std::pair<int, double>>& sortedDistances);

		/*
		*Calculate the Hausdorff distance from the boundary points of the un-simplified triangular mesh to the entire NURBS boundary curve.
		*
		* @param model CCG_QMS_model
		* @param objTriMesh
		* @param sampleStep
		* @param sortedDistances
		* @return double Hausdorff diatance
		* add some points on boundary edges on triangular mesh
		*/
		double computeHausdorffDistance_allNURBSBoundaryCurve2(M* objTriMesh, double sampleStep, std::vector<std::pair<int, double>>& sortedDistances);

		//Output the sorted distance values to a .txt file, with each line formatted as: Point ID, Distance Value
		void outputSortedDistances(const std::vector<std::pair<int, double>>& sortedDistances, const std::string& filename);


		std::vector<CPoint>& controlPoints() { return m_controlPoints; }
		std::vector<std::shared_ptr<CCG_QMS_bezierSurf>>& bfs() { return m_bfs; }
		std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>>& bsplineSurfs() { return m_bsplineSurfs; }
		std::vector<NURBS_boundaryCurve>& nbcs() { return m_nbcs; }
		std::map<int, std::vector<int>>& vertexMap() { return m_vertexMap; }
		std::vector<std::vector<BoundaryCurveSample>>& getBoundarySamples() {return m_boundarySamples;}

		//Attribute to boundary curve fitting
		/*-----------------------------------------------------------------------------------------------*/
		/*
		* Calculate the representation of the sampling points on the model
		*/
		void computeRealSamplingPos(int faceId, int firstVId, CPoint2 uv, std::vector<int>& controlIndexs, std::vector<double>& controlWeights);

		/*
		* 直接缩放每个bezier曲面上的控制点
		* 到模型初始的尺寸
		*/
		void obtainBezierCptsBy_CWs_scale(M* pMesh, double scale);
		/*-----------------------------------------------------------------------------------------------*/

		//Attribute to quantify the surface quality
		/*-----------------------------------------------------------------------------------------------*/
		/**
		 * create all boundary curves for model.
		 * take into account periodic spline
		 */
		void create_bsurBoundarys_general(M* pMesh);		
		void createTopology_bspline_TMesh(M* pMesh);
		void outputTopology_bspline_TMesh_obj(M* pMesh, const char* output);
		/**
		 * Quantize the bouondary fitting property.
		 * add vertex on tri mesh boundary edges
		 * using sampling points from boundary segments of model's boundary
		 * computing vertex's minimum distance to all sampling points
		 * computing boundary's hausdorff distance by computing their sampling points to all vertex from tri mesh boundary vertex and edeges
		 */
		double computeHausdorffDistance_allNURBSBoundaryCurve3(M* objTriMesh, int _numSamplesPerEdge, std::vector<std::pair<int, double>>& sortedDistances);

		/**
		 * Traverse all QB_bsplineSurf，calculate：
		 *  1. Maximal boundaryEdgeRatio
		 *  2. Minimal minCtrlPtDist
		 *  3. Maximal maxHausdorff
		 */
		void computeGlobalSurfaceMetrics();
		/**
		 * Write the following five pieces of data into a text file in a clear and readable format:
		 *   - qbs_diffG0_sampling
		 *   - qbs_diffG1_sampling
		 *   - qbs_boundaryEdgeRatio
		 *   - qbs_minCtrlPtDist
		 *   - qbs_maxHausdorff
		 *
		 * @param filename The path of the output text file (including the file name)
		 * @return If the file is successfully opened and written to, return true; otherwise, return false.
		 */
		void outputMetricsToFile(const std::string& filename);

		M& bsplineSurfacesTopology_quad() { return m_bsplineSurfacesTopology_quad; }
		QB_sampling*& diffG0_sampling() { return m_diffG0_sampling; };
		QB_sampling*& diffG1_sampling() { return m_diffG1_sampling; };
		QB_sampling*& diff_duu_sampling() { return m_diff_duu_sampling; };
		QB_sampling*& diff_duv_sampling() { return m_diff_duv_sampling; };
		QB_sampling*& diff_dvv_sampling() { return m_diff_dvv_sampling; };

		double& boundaryEdgeRatio() { return m_boundaryEdgeRatio; }
		double& minCtrlPtDist() { return m_minCtrlPtDist; }
		double& maxHausdorff() { return m_maxHausdorff; }
		/*-----------------------------------------------------------------------------------------------*/

		/*-----------------------------------------------------------------------------------------------*/
		/*
		* Fitting
		*/
		std::vector< std::shared_ptr<QB_sampling>>& samplings() { return m_samplings; }
		void computeSamplingPtUV(M* pMesh);
		/**
		 * Calculate the integral on the Bezier surface through the formula derived by oneself.
		 */
		void computeBezierSufIntegral_mine(std::shared_ptr<CCG_QMS_bezierSurf> bf, std::map<int, BE_linkingCtlPoints>& ctl_ws);
		/*
		* Approximate based on the sampling points, optimize to obtain the corresponding control point set.
		* The optimization conditions include thin plate energy and keeping the boundary points of the control mesh (initial quadrilateral mesh) unchanged.
		* Check the fitting algorithm.
		*/
		void approximateSurface_optimize_check(M* pMesh, double approxi_thread, double approxapproxi_smooth_weight);
		/*
		* Approximate based on the sampling points,
		* optimize to obtain the corresponding control point set,
		* add thin plate energy to the optimization conditions, and check the fitting algorithm.
		*/
		void approximateSurface_optimize_check2(M* pMesh, double approxi_thread, double approxapproxi_smooth_weight);
		/*-----------------------------------------------------------------------------------------------*/

	private:
		/**
		 * control points.
		 */
		std::vector<CPoint> m_controlPoints;

		/*Bezier surface list*/
		std::vector<std::shared_ptr<CCG_QMS_bezierSurf>> m_bfs;
		/**
		 * B spline surface list
		 */
		std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>> m_bsplineSurfs;

		/*store all nurbs boundary curve*/
		std::vector<NURBS_boundaryCurve> m_nbcs;

		std::map<int, std::vector<int>> m_vertexMap;

		std::vector<std::vector<BoundaryCurveSample>> m_boundarySamples;

		//Attribute to quantify the surface quality
		/*-----------------------------------------------------------------------------------------------*/
		/**
		 * Store the topology of the B-spline surface in the form of a quadrilateral (polygon) mesh.
		 */
		M m_bsplineSurfacesTopology_quad;
		QB_sampling* m_diffG0_sampling;
		QB_sampling* m_diffG1_sampling;
		QB_sampling* m_diff_duu_sampling;
		QB_sampling* m_diff_duv_sampling;
		QB_sampling* m_diff_dvv_sampling;

		double m_boundaryEdgeRatio;  // The maximum value of the ratio of the longest edge to the shortest edge among all the polygons.
		double m_minCtrlPtDist;  // The minimum distance between any two points in the control mesh
		double m_maxHausdorff;  // The hausdorff distance
		/*-----------------------------------------------------------------------------------------------*/

		/*-----------------------------------------------------------------------------------------------*/
		/*
		* Fitting
		*/		
		std::vector<std::shared_ptr<QB_sampling>> m_samplings;
		/*-----------------------------------------------------------------------------------------------*/

	};

	class CCG_QMS_bezierSurf
	{
	public:
		CCG_QMS_bezierSurf(){
			m_type = RegularSurf;
			m_degree[0] = 3, m_degree[1] = 3;
			m_id = -1;
			m_maxU = 1.0;
			m_maxV = 1.0;
			m_cpts.resize((m_degree[0] + 1) * (m_degree[1] + 1));
			m_linkCtlPWs.resize((m_degree[0] + 1) * (m_degree[1] + 1));
		};
		~CCG_QMS_bezierSurf() {};

		BezierSurfType& type() { return m_type; }
		int& id() { return m_id; }
		int& degree(int i) { return m_degree[i]; }
		std::vector<MeshLib::CPoint>& cpts() { return m_cpts; }
		std::vector<BE_linkingCtlPoints>& linkCptPWs() { return m_linkCtlPWs; }
		std::vector<BE_linkingCtlPoints>& singular_linkCptPWs() { return m_singular_linkCtlPWs; }
		std::vector<MeshLib::CPoint>& cpts_singular() { return m_cpts_singular; }
		std::vector<double>& knotsU() { return m_knotsU; }
		std::vector<double>& knotsV() { return m_knotsV; }
		double& maxU() { return m_maxU; }
		double& maxV() { return m_maxV; }

		/*Raise the degree of bezier surface, 3 to 5*/
		void raiseDegree3to5(M* pMesh);

	private:
		BezierSurfType m_type;
		/*
		* The degree of two directions, 
		* default m_degree[0] is at U
		* m_degree[1] is at V
		*/
		int m_degree[2] = { 3,3 };
		int m_id;
		/*
		* Control points
		* The index of control ponts can be noted as 1 ~ (be_degree[0] + 1) * (be_degree[1] + 1)
		* Halfedge local id is 0，1，2，3, starting from the halfedge at the quad face, counting CCW
		* Bi-3 bezier：
		* 1---the first corner, computing from the halfedge 0 and its symmetric halfedge
		* 2，3---two edge control points at the first edge,computing from the halfedge 0 and its symmetric halfedge
		* 4---the second corner, computing from the halfedge 1 and its symmetric halfedge
		* 5，9---two edge control points at the fourth edge,computing from the halfedge 3 and its symmetric halfedge
		* 6，7---the first and second face control points,computing from the halfedge 0
		* 8，12---two edge control points at the second edge,computing from the halfedge 1 and its symmetric halfedge
		* 10，11---the third and fourth face control points,computing from the halfedge 2
		* 16---the third corner, computing from the halfedge 2
		* 13---the fourth corner, computing from the halfedge 3
		* 14，15---two edge control points at the third edge,computing from the halfedge 2 and its symmetric halfedge
		*/
		std::vector<MeshLib::CPoint> m_cpts;
		/*Mapping of control points and quad mesh vertex*/
		std::vector<BE_linkingCtlPoints> m_linkCtlPWs;
		/*Mapping of control points and quad mesh vertex at singular surface*/
		std::vector<BE_linkingCtlPoints> m_singular_linkCtlPWs;
		/*
		* For the singular surface patches, the control points of the surface are reordered
		* the index is noted by 1~(be_degree[0] + 1) * (be_degree[1] + 1)
		* starting from singularity, the local id of halfedges respectively are 0，1，2，3;(logically, not change local id)
		* an example bi-5 bezier surface：
		* 0,1,2,3,4,5 (local id 0 halfedge)
		* 5,11,17,23,29,35 (local id 1 halfedge)
		* 35,34,33,32,31,30 (local id 2 halfedge)
		* 30,24,18,12,6,0 (local id 3 halfedge)
		*/
		std::vector<MeshLib::CPoint> m_cpts_singular;

		/*kont vector in U*/
		std::vector<double> m_knotsU;
		/*kont vector in V*/
		std::vector<double> m_knotsV;
		/*the max knot value in U*/
		double m_maxU;
		/*the max knot value in V*/
		double m_maxV;
		
	};

	class CCG_QMS_bsplineSurf
	{
	public:
		CCG_QMS_bsplineSurf() {
			m_id = 0;
			m_periodicAlongU = false;
			m_periodicAlongV = false;
			m_closedAlongU = false;
			m_closedAlongV = false;
			m_quadMeshFaceId = -1;
			m_quadLayoutIndex = -1;
			m_quadIds.resize(4);
			//Attribute to quantify the surface quality
			/*-----------------------------------------------------------------------------------------------*/
			m_diffG0_sampling = NULL;
			m_diffG1_sampling = NULL;
			m_diff_duu_sampling = NULL;
			m_diff_duv_sampling = NULL;
			m_diff_dvv_sampling = NULL;
			m_boundaryEdgeRatio = 0.0; // Initialize the boundary edge ratio to 0.0.
			m_minCtrlPtDist = std::numeric_limits<double>::infinity(); // Initialize the minimum control point distance to positive infinity.
			m_maxHausdorff = 0.0;
			/*-----------------------------------------------------------------------------------------------*/

		};
		~CCG_QMS_bsplineSurf() {};
		int& id() { return m_id; }
		bool& periodicAlongU() { return m_periodicAlongU; }
		bool& periodicAlongV() { return m_periodicAlongV; }
		bool& closedAlongU() { return m_closedAlongU; }
		bool& closedAlongV() { return m_closedAlongV; }
		std::vector<std::vector<CPoint>>& cpts() { return m_cpts; }
		std::vector<int>& degree() { return m_degree; }
		int& quadMeshFaceId() { return m_quadMeshFaceId; }
		std::vector<double>& knotsU() { return m_knotsU; }
		std::vector<double>& knotsV() { return m_knotsV; }
		std::vector<int>& quadIds() { return m_quadIds; }
		std::map<double, int>& knotsU_unique_num() { return m_knotsU_unique_num; }
		std::map<double, int>& knotsV_unique_num() { return m_knotsV_unique_num; }
		int& quadLayoutIndex() { return m_quadLayoutIndex; }
		std::vector<NUBE_bsurfBoudnary*>& boundarys() { return m_boundarys; }

		/*Jingpeng Yin 2025/5/19, knot removal algorithm under given tolerance*/
		void RemoveCurveKnot_UVConsistency2(double error_threshold);
		/*Calculate the values of the non-redundant knots and their corresponding redundancy counts based on the knot vectors.*/
		void computeUniqueKnotNum();
		/**
		* Using the De Boor algorithm, calculate the B-spline surface points at the given parameters u and v
		*/
		CPoint DeBoorAlgorithmSurface(const std::vector<std::vector<CPoint>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, double u, double v);

		//Attribute to quantify the surface quality
		/*-----------------------------------------------------------------------------------------------*/
		QB_sampling*& diffG0_sampling() { return m_diffG0_sampling; };
		QB_sampling*& diffG1_sampling() { return m_diffG1_sampling; };
		QB_sampling*& diff_duu_sampling() { return m_diff_duu_sampling; };
		QB_sampling*& diff_duv_sampling() { return m_diff_duv_sampling; };
		QB_sampling*& diff_dvv_sampling() { return m_diff_dvv_sampling; };
		double& boundaryEdgeRatio() { return m_boundaryEdgeRatio; } 
		double& minCtrlPtDist() { return m_minCtrlPtDist; }
		double& maxHausdorff() { return m_maxHausdorff; }
		/*-----------------------------------------------------------------------------------------------*/

	private:
		//The id and id-1 of the B-spline surface are consistent with the position of the corresponding vector.
		int m_id;
		/*periodic along U*/
		bool m_periodicAlongU;
		/*periodic along V*/
		bool m_periodicAlongV;
		/*closed along U*/
		bool m_closedAlongU;
		/*closed along V*/
		bool m_closedAlongV;
		/*control points*/
		std::vector<std::vector<CPoint>> m_cpts;
		/*The order of the two directions, with the default value being 3*/
		std::vector<int> m_degree = { 3,3 };
		/*Find the id of a quadrilateral mesh face that corresponds to the current B-spline surface.*/
		int m_quadMeshFaceId;
		//kont vector in U (qbb_cpts[0].size() + qbb_degree[0] + 1)
		std::vector<double> m_knotsU;
		//kont vector in V (qbb_cpts[1].size() + qbb_degree[1] + 1)
		std::vector<double> m_knotsV;
		/*
		* Store the IDs of the four vertices of the B-spline surface corresponding to the quadrilateral mesh vertex, 
		* ensuring that the IDs of each half-edge on the surface are in the correct order:
		* （0,0）[0]--2-->(1,0)[1] --3-->(1,1)[2] --4-->(0,1)[3]--1-->(0,0)[0]
		*/
		std::vector<int> m_quadIds;
		/*The values of non-repeating knots in the U direction and their corresponding multiplicities*/
		std::map<double, int> m_knotsU_unique_num;
		/*The values of non-repeating knots in the V direction and their corresponding multiplicities*/
		std::map<double, int> m_knotsV_unique_num;
		/*the current NURBS belongs to which T-layout patch*/
		int m_quadLayoutIndex;
		/**
		 * boudnary on this B spline surface.
		 * expected size is 4
		 */
		std::vector<NUBE_bsurfBoudnary*> m_boundarys;
		//Attribute to quantify the surface quality
		/*-----------------------------------------------------------------------------------------------*/
		QB_sampling* m_diffG0_sampling;
		QB_sampling* m_diffG1_sampling;
		QB_sampling* m_diff_duu_sampling;
		QB_sampling* m_diff_duv_sampling;
		QB_sampling* m_diff_dvv_sampling;
		double m_boundaryEdgeRatio;  //Ratio of the longest edge to the shortest edge
		double m_minCtrlPtDist;  // The minimum distance between any two points in the control point mesh
		double m_maxHausdorff;  // Store the maximum value of the Hausdorff dimension for each boundary on the current surface.
		/*-----------------------------------------------------------------------------------------------*/
	};

	class NUBE_bsurfBoudnary
	{
	public:
		NUBE_bsurfBoudnary() {
			m_boundary = false;
			m_localId = -1;
			m_bsurfId = -1;
			qbs_diffG0_sampling = NULL;
			qbs_diffG1_sampling = NULL;
			qbs_diff_duu_sampling = NULL;
			qbs_diff_duv_sampling = NULL;
			qbs_diff_dvv_sampling = NULL;
			m_hausdorffDist = 0.0;
		};
		~NUBE_bsurfBoudnary() {};
		bool& boundary() { return m_boundary; }
		std::vector<H*>& halfedges() { return m_halfedges; }
		int& localId() { return m_localId; }
		int& bsurfId() { return m_bsurfId; }
		std::vector<NUBE_bsurfSegment*>& segements() { return m_segments; }
		QB_sampling*& diffG0_sampling() { return qbs_diffG0_sampling; };
		QB_sampling*& diffG1_sampling() { return qbs_diffG1_sampling; };
		QB_sampling*& diff_duu_sampling() { return qbs_diff_duu_sampling; };
		QB_sampling*& diff_duv_sampling() { return qbs_diff_duv_sampling; };
		QB_sampling*& diff_dvv_sampling() { return qbs_diff_dvv_sampling; };
		double& hausdorffDist() { return m_hausdorffDist; }
	private:
		/**
		 * marking if is boundary at model.
		 */
		bool m_boundary;

		/**
		 * all halfedges on this boundary from a quadrilateral mesh.
		 * sequential arrangement in CCW
		 */
		std::vector<H*> m_halfedges;

		/**
		 * index and parameterization information on its b spline surface.
		 * index
		 * visited by its b spline surface
		 * parameterization information
		 * 1: u=u_min; 2: v=v_min; 3: u=u_max; 4: v=v_max;
		 */
		int m_localId;

		/**
		 * B spline surface id, counting from 1.
		 */
		int m_bsurfId;

		/**
		 * segments on this boundary.
		 */
		std::vector<NUBE_bsurfSegment*> m_segments;

		QB_sampling* qbs_diffG0_sampling;
		QB_sampling* qbs_diffG1_sampling;
		QB_sampling* qbs_diff_duu_sampling;
		QB_sampling* qbs_diff_duv_sampling;
		QB_sampling* qbs_diff_dvv_sampling;

		// Store the Hausdorff distance of this boundary
		double m_hausdorffDist;
	};

	/**
	 * class implementation
	 * Segment of one B spline surface boundary.
	 */
	class NUBE_bsurfSegment
	{
	public:
		NUBE_bsurfSegment() {
			m_bsurfId = -1;
			m_bsurfBoundaryLoacalId = -1;
			m_id = -1;
			m_firstHalfedgeIndex = -1;
			m_endHalfedgeIndex = -1;
			m_symBsurfSeg = NULL;
		};
		~NUBE_bsurfSegment() {};

		int& bsurfId() { return m_bsurfId; }
		int& bsurfBoundaryLoacalId() { return m_bsurfBoundaryLoacalId; }
		int& id() { return m_id; }
		NUBE_bsurfSegment*& symBsurfSeg() { return m_symBsurfSeg; }
		std::vector<QB_sampling*>& samplings() { return m_samplings; }
		int& firstHalfedgeIndex() { return m_firstHalfedgeIndex; }
		int& endHalfedgeIndex() { return m_endHalfedgeIndex; }
	private:

		/**
		 * B spline surface id, counting from 1.
		 */
		int m_bsurfId;

		/**
		 * index and parameterization information on its b spline surface.
		 * index
		 * visited by its b spline surface
		 * parameterization information
		 * 1: u=u_min; 2: v=v_min; 3: u=u_max; 4: v=v_max;
		 */
		int m_bsurfBoundaryLoacalId;

		/**
		 * index on its b spline surface boundary.
		 * visited by its b spline surface boundary
		 */
		int m_id;

		/**
		 * the index of starting halfedge on its boundary's halfedges
		 * counting from 1.
		 */
		int m_firstHalfedgeIndex;

		/**
		 * the index of ending halfedge on its boundary's halfedges.
		 */
		int m_endHalfedgeIndex;

		/**
		 * symmetry segment.
		 */
		NUBE_bsurfSegment* m_symBsurfSeg;

		/**
		 * sampling points from this segment.
		 * sequential arrangement in CCW.
		 */
		std::vector<QB_sampling*> m_samplings;

	};

	/**
	 *Definition of sampling point class.
	 */
	class QB_sampling
	{
	public:
		QB_sampling() {
			qbs_fId = -1;
			qbs_vId = -1;
			for (int i = 0; i < qbs_ws.size(); i++)
			{
				qbs_ws[i] = 0.0;
			}
			qbs_mark = true;
			qbs_boundary = false;
			qbs_patchId = -1;
			qbs_kappa1 = 0.0;
			qbs_kappa2 = 0.0;

			qbs_dualSampling = NULL;

			qbs_diffG0 = 0.0;
			qbs_diffG1 = 0.0;
			qbs_diff_duu = 0.0;
			qbs_diff_duv = 0.0;
			qbs_diff_dvv = 0.0;
		};
		~QB_sampling() {};
		int& fId() { return qbs_fId; }
		int& vId() { return qbs_vId; }
		CPoint& pos() { return qbs_pos; }
		CPoint& pos_approx() { return qbs_pos_approx; }
		CPoint& normal() { return qbs_n; }
		std::vector<double>& ws() { return qbs_ws; }
		CPoint2& uv() { return qbs_uv; }
		CPoint2& uv_init() { return qbs_uv_init; }
		std::vector<int>& index() { return qbs_index; }
		std::vector<double>& weight() { return qbs_weight; }
		bool& mark() { return qbs_mark; }
		bool& boundary() { return qbs_boundary; }
		int& patchId() { return qbs_patchId; }
		double& kappa1() { return qbs_kappa1; }
		double& kappa2() { return qbs_kappa2; }
		CPoint& t1() { return qbs_t1; }
		CPoint& t2() { return qbs_t2; }
		CPoint& du() { return qbs_du; }       // The first-order derivative in the u direction
		CPoint& dv() { return qbs_dv; }       // The first-order derivative in the v direction
		CPoint& duu() { return qbs_duu; }     // The second-order derivative in the u direction
		CPoint& duv() { return qbs_duv; }     // Mixed second-order derivative
		CPoint& dvv() { return qbs_dvv; }     // The second-order derivative in the v direction

		QB_sampling*& dualSampling() { return qbs_dualSampling; }
		double& diffG0() { return qbs_diffG0; }
		double& diffG1() { return qbs_diffG1; }
		double& diff_duu() { return qbs_diff_duu; }
		double& diff_duv() { return qbs_diff_duv; }
		double& diff_dvv() { return qbs_diff_dvv; }
	private:
		/*
		* The id corresponding to the quadrilateral mesh face.
		*/
		int qbs_fId;
		/*
		* The ID corresponding to the patch.
		*/
		int qbs_patchId;
		/*
		* The ID of the first vertex in the corresponding surface.
		*/
		int qbs_vId;
		/*
		* The spatial location of the sampling points.
		*/
		CPoint qbs_pos;
		/*
		* The spatial position fitted out by the sampling points.
		*/
		CPoint qbs_pos_approx;
		/*The parameter position weights of the sampling points*/
		std::vector<double> qbs_ws;
		/*
		Parameters of the sampling point.
		*/
		CPoint2	qbs_uv;
		/*
		The initial parameters of the sampling points, corresponding to the results given by the global parameterization
		*/
		CPoint2	qbs_uv_init;
		/*
		* The normal vector corresponding to the spatial coordinates of the sampling point
		*/
		CPoint qbs_n;
		/*Corresponding control point index*/
		std::vector<int> qbs_index;

		/*Corresponding control point weight*/
		std::vector<double> qbs_weight;

		/*Indicate whether the current sampling point can be used for approximation*/
		bool qbs_mark;

		/*Mark whether the current sampling point is a boundary sampling point*/
		bool qbs_boundary;

		/*The principal curvature and principal direction of the sampling point*/
		double qbs_kappa1, qbs_kappa2;         // principal curvature
		CPoint qbs_t1, qbs_t2;        // principal direction

		CPoint qbs_du;    // The first-order derivative in the u direction
		CPoint qbs_dv;    // The first-order derivative in the v direction
		CPoint qbs_duu;   // The second-order derivative in the u direction
		CPoint qbs_duv;   // Mixed second-order derivative
		CPoint qbs_dvv;   // The second-order derivative in the v direction

		QB_sampling* qbs_dualSampling;

		double qbs_diffG0;
		double qbs_diffG1;
		double qbs_diff_duu;
		double qbs_diff_duv;
		double qbs_diff_dvv;
	};

}
#endif // !_CCG_QMS_QUAD2MANIFOLDSPLINE_H_
