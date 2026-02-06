#ifndef _BSPLINECURVAPPROX_H_
#define _BSPLINECURVAPPROX_H_

#include"../3rdParty/MeshLib/core/Mesh/BaseMesh.h"
#include"CCG_QMS_meshTool.h"
#include <vector>
#include <map>
#include "CCG_QMS_quad2ManifoldSpline.h"

namespace CCG_QMSLib
{
	typedef CCG_QMS_Vertex	V;
	typedef CCG_QMS_Edge	E;
	typedef CCG_QMS_Face	F;
	typedef CCG_QMS_HalfEdge H;
	typedef CCG_QMS_Mesh<V, E, F, H> M;
	typedef typename M::It			It;

	class CLoopSegment
	{
	public:
		/*!
			Constructor of the CLoopSegment
			\param pMesh  pointer to the current mesh
			\param pH halfedge on the boundary loop
		*/
		CLoopSegment(M* pMesh, H* pH);
		/*!
		   Destructor of CLoop.
		*/
		~CLoopSegment()
		{
			m_halfedges.clear();
		};

		/*!
		The vector of haledges on the current boundary loop segment.
		*/
		std::vector<H*>& halfedges()
		{
			return m_halfedges;
		}
		/*!
			Starting vertex
		*/
		V* start()
		{
			if (m_halfedges.size() == 0) return NULL;
			H* he = m_halfedges[0];
			return m_pMesh->halfedgeSource(he);
		}
		/*!
			ending vertex
		*/
		V* end()
		{
			if (m_halfedges.size() == 0) return NULL;
			size_t n = m_halfedges.size();
			H* he = m_halfedges[n - 1];
			return m_pMesh->halfedgeTarget(he);
		}
		//triangular mesh arc id
		int& getTriArcId()
		{
			return 	 triArc_id;
		}
		//quadrilateral loop id
		int& getQuadArcId()
		{
			return quadArc_id;
		}
		std::list<V*>& vertex()
		{
			return m_vertexs;
		}

		/*!
			The length of the current boundary loop.
		*/

	protected:
		/*!
			The mesh pointer
		*/
		M* m_pMesh;

		/*!
			The vector of consecutive halfedges along the boundary loop.
		*/
		std::vector<H*>							  m_halfedges;
		int triArc_id;
		int quadArc_id;
		V* m_vertex;
		std::list<V*>							  m_vertexs;
		double m_vertexLength;
	};
	class CLoop
	{
	public:
		/*!
			Constructor of the CLoop
			\param pMesh  pointer to the current mesh
			\param pH halfedge on the boundary loop
		*/
		CLoop(M* pMesh, H* pH);
		/*!
			Constructor of the CLoop
			\param pMesh  pointer to the current mesh
		*/
		CLoop(M* pMesh) {
			m_pMesh = pMesh;
			m_length = 0;
			m_pHalfedge = NULL;
			m_vertexLength = 0;
			triLoop_id = -1;
			quadLoop_id = -1;
			m_vertex = NULL;
		};
		/*!
		   Destructor of CLoop.
		*/
		~CLoop();

		/*!
		The list of haledges on the current boundary loop.
		*/
		std::list<H*>& halfedges()
		{
			return m_halfedges;
		}
		std::list<V*>& vertex()
		{
			return m_vertexs;
		}

		/*!
			The length of the current boundary loop.
		*/
		double length()
		{
			return m_length;
		}
		/*!
			The vector of segments on this loop
		*/
		std::vector<CLoop*>& segments()
		{
			return m_segments;
		}
		/*!
			divide the loop to segments
			\param markers, the array of markers, the first marker is the starting one
		*/
		void divide(std::vector<V*>& markers);

		//Obtain the loop ID of the triangular mesh
		int& getTriLoopId()
		{
			return 	 triLoop_id;
		}
		//Obtain the loop ID of the quadrilateral mesh
		int& getQuadLoopId()
		{
			return quadLoop_id;
		}

		double vertexLength()
		{
			return m_vertexLength;
		}
	protected:
		/*!
			Pointer to the current mesh.
		*/
		M* m_pMesh;
		/*! The length of the current boundary loop.
		*/
		double											  m_length;
		/*!
			Starting halfedge of the current boundary loop.
		*/
		H* m_pHalfedge;
		/*!
			List of consecutive halfedges along the boundary loop.
		*/
		std::list<H*>							  m_halfedges;
		/*!
			Vector of segments
		*/
		std::vector<CLoop*>							  m_segments;


		int triLoop_id;
		int quadLoop_id;

		V* m_vertex;
		std::list<V*>							  m_vertexs;
		double m_vertexLength;

	};
	class CFeatureLoop
	{

	public:
		/*!
		CBoundary constructor
		\param pMesh pointer to the current mesh
		*/
		CFeatureLoop(M* pMesh);

		/*!
		CBoundary destructor
		*/
		~CFeatureLoop();
		/*!
		The list of boundary loops.
		*/
		std::vector<CLoop*>& loops()
		{
			return m_loops;
		}
		// Add arc segments
		std::vector<CLoopSegment*>& arcs() {
			return m_arcs;
		}

		M* getMesh() const { return m_pMesh; } // Provide read-only access

	protected:
		/*!
			Pointer to the current mesh.
		*/
		M* m_pMesh;
		/*!
			List of boundary loops.
		*/
		std::vector<CLoop*> m_loops;
		// Add arc segments
		std::vector<CLoopSegment*> m_arcs;


	};
	enum BezierCurveType
	{
		RegularCurve = 0,
		SingularCurve = 1,
		LinkingSingularCurve = 2
	};
	class BE_bezierCurve_approxi
	{
	public:
		BE_bezierCurve_approxi() {
			/* The degree of the curve is 3 */
			be_degree = 3;
			/* Initialize the type of the curve as a regular curve */
			be_type = RegularCurve;
			/* The id of the curve is initialized as -1, indicating that it is not associated with any specific data. */
			be_id = -1;
			/* Initialize the number of control points to (be_degree + 1) */
			be_cpts.resize(be_degree + 1);
			/* The default maximum value of the knot vector is 1.0 */
			be_maxU = 1.0;
			be_knotsU = { 0.0, 0.0 ,0.0,0.0,1.0,1.0,1.0,1.0 };
		};
		~BE_bezierCurve_approxi() {};

		/* Obtain the type of the curve and allow for modification */
		BezierCurveType& type() { return be_type; }
		/* Obtain the ID of the curve and allow modification */
		int& id() { return be_id; }
		/* Obtain the control points of the curve and allow for modification */
		std::vector<MeshLib::CPoint>& cpts() { return be_cpts; }
		
		/* Return and allow the modification of the degree of the curve */
		int& degree() { return be_degree; }
		/* Return and allow modification of the knot vector of the curve */
		std::vector<double>& knotsU() { return be_knotsU; }
		double& maxU() { return be_maxU; }
		/*!
	 * Calculate the point on the curve with parameter u
	 * @param u Curve parameter, range [0, 1]
	 * @return The coordinate of the curve point corresponding to parameter u
	 */
		CPoint evaluate(double u) const;

	private:

		/* Degree of curve */
		int be_degree = 3;
		/* Type of curve */
		BezierCurveType be_type;
		/* id */
		int be_id;
		/*
		 * Control points
		 * The positions of the control points can be recorded using labels ranging from 1 to (be_degree + 1).
		 * For example, a bicubic Bezier curve:
		 * 1---First point
		 * 2---Second point
		 * 3---Third point
		 * 4---Fourth point
		 */
		std::vector<MeshLib::CPoint> be_cpts;

		/* Knot vector */
		std::vector<double> be_knotsU;
		/* Maximum value of the knot */
		double be_maxU;
	};
	class QB_bsplineCurve_approxi
	{
	public:
		QB_bsplineCurve_approxi() {
			qbc_periodic = false;  // Initialize as a non-periodic curve
			qbc_closed = false;    // Initialize as a non-closed curve
		}
		~QB_bsplineCurve_approxi() {}

		// Getter and Setter function
		std::vector<CPoint>& cpts() { return qbc_cpts; }
		int& degree() { return qbc_degree; }
		std::vector<double>& knots() { return qbc_knots; }
		bool& periodic() { return qbc_periodic; }
		bool& closed() { return qbc_closed; }

	private:
		/*
		 * Control point coordinates
		 * Store the set of control points of the curve, with the storage structure being a 1D array
		 */
		std::vector<CPoint> qbc_cpts;

		/*
		 * The degree of the curve, with a default value of 3
		 */
		int qbc_degree = 3;

		/*
		 * Knot vector, with a size of qbc_cpts.size() + qbc_degree + 1
		 */
		std::vector<double> qbc_knots;

		/*
		 * Is it a periodic curve?
		 */
		bool qbc_periodic;

		/*
		 * Is it a closed curve?
		 */
		bool qbc_closed;
	};
	class QB_curve_sampling;//Declaration of the sampling point class

	struct BoundaryFittingArc {
		std::vector<H*> orderedHalfedges;  // 有序的半边序列
	};

	class BE_curve_model
	{
	public:
		BE_curve_model() {
			m_pMesh = NULL;
			be_triFeatureLoops = NULL;
			be_quadFeatureLoops = NULL;
		};
		~BE_curve_model() {};

		/*Access and modify all the Bezier curves in the model*/
		std::vector<BE_bezierCurve_approxi*>& bcs() { return be_bcs; }
		std::vector<QB_curve_sampling*>& samplings() { return qbm_samplings; }
		std::vector<CPoint>& controlPoints() { return qbm_controlPoints; }
		std::vector<QB_bsplineCurve_approxi*>& bsplineCurves() { return qbm_bsplineCurves; }

		/*Access and modify the feature line loop set of the triangular mesh*/
		CFeatureLoop*& triFeatureLoops() { return be_triFeatureLoops; }

		/*Access and modify the feature line loop set of the quadrilateral mesh*/
		CFeatureLoop*& quadFeatureLoops() { return be_quadFeatureLoops; }

		/*Normalize both the triangular mesh points and the quadrilateral mesh points together*/
		void normalize_triQuadMesh(M* triMesh, M* quadMesh);

		/*Establish the initial control points of the Bezier curve*/
		void setControlPoints(M* pMesh);

		/*
		* Marking feature points on quad mesh boundry
		* according to the input information
		*/
		void markFeaturePts_input(M* pMesh);

		//Calculate the angle function
		static double calculateAngle(const CPoint& a, const CPoint& b, const CPoint& c);
		//Triangle corner detection function
		bool checkTriCorner(const CPoint& quadPoint, V* triPoint) const;

		/*Rebuild the control points of the Bezier curve*/
		void resetControlPoints(M* pMesh);
		/*Determine the corresponding relationship between the feature rings of triangular mesh and the feature rings of quadrilateral mesh*/
		void mapTriLoopToQuadLoop(M* quadMesh, M* triMesh);

		/*Access and modify the corresponding relationship between the feature rings of triangular meshes and those of quadrilateral meshes*/
		std::map<int, int>& getTriToQuadLoopMap() {
			return triToQuadLoopMap;
		}
		/*Access and modify the corresponding relationship between the characteristic arcs of triangular meshes and those of quadrilateral meshes*/
		std::map<int, int>& getTriToQuadArcMap() {
			return triToQuadArcMap;
		}
		/*Determine the corresponding relationship between the triangular mesh points and the quadrilateral mesh halfedges.*/
		void triPointToQuadHalfedge(M* quadMesh, M* triMesh);
		/*
		* Specify boundary points on tri mesh mapping to the quad mesh boundary halfedge
		* but let not boundary points on tri mesh as sampling points
		*/
		void triPointToQuadHalfedge2(M* quadMesh, M* triMesh);

		/*
		* specify all boundary halfedges from triangle mesh
		* corresponding to a boundary halfedge from quad mesh
		*/
		void triBoundaryHalfedgesMapping2QuadBoundaryHalfedges();

		/**
	  * 定位贴体性较差的四边形网格边界边
	  * @param quadMesh 四边形网格
	  * @param minPointsThreshold 边界边上待拟合点数量的阈值，默认为4
	  * @param errorThreshold 拟合误差阈值，默认为5e-5
	  */
		void locatePoorFittingBoundaryEdges(M* quadMesh, int minPointsThreshold = 4, double errorThreshold = 5e-5);

		//Add interpolation points to the half-edges of the un-mapped quadrilateral
		void addInterpolatedPointsForUnmappedHalfedges(M* pMesh);
		/*
		* Add sampling points on quad boundry halfedges to satisfy known condition
		* according to "samplingPts2QuadHalfedges" mapping
		*/
		void addInterpolatedPointsForUnmappedHalfedges2();

		//Declaration of auxiliary function
		std::pair<V*, double> findClosestTriPoint(CLoop* triLoop, const MeshLib::CPoint& p0, const MeshLib::CPoint& p1);
		//std::pair<V*, double> findClosestTriPoint(CLoop* triLoop, const CPoint& p0, const CPoint& p1);
		std::pair<H*, H*> findAdjacentHalfedges(CLoop* triLoop, V* vert);
		std::vector<CPoint> interpolateHalfedges(H* he1, H* he2, int count);
		H* findClosestQuadHalfedge(CLoop* quadLoop, const CPoint& point);
		void addSamplingPoint(const CPoint& point, int loopId, H* quadHe);

		/*Access and modify the corresponding relationships between the triangular mesh points and the quadrilateral mesh edges.*/
		std::map<V*, H*>& getTriPointToQuadHalfedgeMap()
		{
			return triPointToQuadHalfedgeMap;
		}

		/*Access and modify the corresponding relationship between data points and the half-edges of the quadrilateral mesh*/
		std::map<QB_curve_sampling*, H*>& getDataPointToQuadHalfedgeMap()
		{
			return dataPointToQuadHalfedgeMap;
		}

		/*Reorder the IDs of the half-edges of the characteristic loops in the quadrilateral mesh*/
		void initializeQuadHalfedgeIdMap(M* quadMesh);
		void initializeQuadHalfedgeIdMap();
		/*Initialize tri boundary halfedges id according to feature halfedges' sequence*/
		void initializeTriHalfedgeIdMap();

		/*Access and modify the corresponding relationship between the half-edges of the feature rings in the quadrilateral mesh and their IDs*/
		std::map<H*, int>& getQuadHalfedgeToIdMap()
		{
			return quadHalfedgeToIdMap;
		}

		/*Calculate the distance from a point to a line*/
		double pointToSegmentDistance(const CPoint& point, const CPoint& segStart, const CPoint& segEnd);

		/*Reorder the IDs of the characteristic ring points in the quadrilateral mesh*/
		void initializeQuadVertexIdMap(M* quadMesh);

		/*Access and modify the corresponding relationship between the feature ring points and their IDs in the quadrilateral mesh*/
		std::map<V*, int>& getQuadLoopVertexIdToMap()
		{
			return quadLoopVertexIdToMap;
		}

		/*Access and modify the corresponding relationship between the id-> parameter of the triangular feature ring points*/
		std::map<int, double>& getTriPointIdToParamMap()
		{
			return triPointIdToParamMap;
		}
		void parametrize(M* quadMesh);
		void parametrize3(M* quadMesh);
		void parametrize2(M* quadMesh);
		void parametrize4(M* quadMesh);

		/**
		 * Parameterization of centrifugal force parameters
		 * @param startPoint The starting coordinate of one of the halfedges of the quadrilateral
		 * @param endPoint The targeting coordinate of one of the halfedges of the quadrilateral
		 * @param dataPointsSet of triangular mesh data points (already sorted by distance)
		 * @return The parameter value vector (within the range [0, 1]) corresponding to each data point
		 */
		std::vector<double> computeCentripetalParams(
			const CPoint& startPoint,
			const CPoint& endPoint,
			const std::vector<QB_curve_sampling*>& dataPoints);

		/*Set up equations, based on the relationship between the control points and the quadrilateral mesh points, to solve for the new quadrilateral mesh points*/
		void compute_quad_vertex(M* quadMesh);
		/*Set up equations, based on the relationship between the control points and the quadrilateral mesh points, to solve for the new quadrilateral mesh points*/
		void compute_quad_vertex1(M* quadMesh);
		void compute_quad_vertex2(M* quadMesh);
		/*
		* cfw add 20250613
		* change equations for some sampling points on halfedges linking feature poionts
		*/
		void compute_quad_vertex3(M* quadMesh);

		/*
		* cfw add 20250613
		* change equations for some sampling points on halfedges linking feature poionts
		* speed improvement
		*/
		void compute_quad_vertex4(M* quadMesh);

		//边界特征弧求解
		void compute_quad_vertex5(M* quadMesh);
		
		double computeTotalError();
		double computeTriTotalError();
		/* Maintain the sharp feature points of the quadrilateral grid */
		void preserveQuadCornerPoints(M* pMesh);

		/*
		* Check degenerated boudnary vertex which causes the quadrilateral mesh to fold
		* at fitting boundary curve
		* and update it to initial position
		*/
		void CheckDegeneratedBoudnaryVertexAndUpdate(M* pMesh);

		/*
		* Check orthogonality after boundary fitting qulity improvement
		* modifying the vertex where orthogonality is poor
		*/
		void postProcessing_checkAndModifyOrthogonality_boundaryFitting(M* pMesh);

		/*Access and modify the corresponding relationship between the id-> parameter of the triangular feature ring points*/
		std::map<int, double>& getDataPointIdToParamMap()
		{
			return dataPointIdToParamMap;
		}

		//Access and modify the mapping between the original points and the fitted points
		std::map<int, std::pair<CPoint, CPoint>>& getOriginalToApproxPointsMap() {
			return originalToApproxPoints;
		}
		void computeTriangleFeaturePointFitting();

		/*obatin and modify mapping between quadHalfedges and triHalfedges*/
		std::map<H*, std::pair<H*, H*>>& getMapQuadHalfedgesAndTriHalfedges()
		{
			return mapQuadHalfedgesAndTriHalfedges;
		}

		// 提供访问结果的接口
		const std::map<H*, std::pair<double, std::vector<QB_curve_sampling*>>>& getPoorFittingEdges() const {
			return poorFittingEdges;
		}

		// 获取端点统计信息
		const std::map<V*, int>& getEndpointCountMap() const {
			return endpointCountMap;
		}

		// 边界特征弧
		const std::vector<BoundaryFittingArc>& getBoundaryFittingArcs() const {
			return m_boundaryFittingArcs;
		}

		// 获取单条弧（按索引）
		const BoundaryFittingArc* getBoundaryFittingArc(size_t index) const {
			if (index < m_boundaryFittingArcs.size()) {
				return &m_boundaryFittingArcs[index];
			}
			return nullptr;  // 索引无效时返回空
		}

		// 返回顶点映射的常量引用
		const std::map<V*, int>& boundaryVertexIdMap() const {
			return m_boundaryVertexIdMap;
		}

		/*访问和修改四边形网格中边界特征环半边和其id的对应关系*/
		std::map<H*, int>& getQuadBoundaryArcHalfedgeToIdMap()
		{
			return quadBoundaryArcHalfedgeToIdMap;
		}

	private:
		M* m_pMesh;
		/*Store all the Bezier curves in the model*/
		std::vector<BE_bezierCurve_approxi*> be_bcs;

		/*Store all the B-spline curves in the model*/
		std::vector<QB_bsplineCurve_approxi*> qbm_bsplineCurves;

		/*Sampling point set*/
		std::vector<QB_curve_sampling*> qbm_samplings;

		/*Control point set*/
		std::vector<CPoint> qbm_controlPoints;

		/*The Loop set storing the feature lines of the triangular mesh*/
		CFeatureLoop* be_triFeatureLoops;

		/*The Loop set storing the feature lines of the quadrilateral mesh*/
		CFeatureLoop* be_quadFeatureLoops;

		/*Mapping relationship from triangular mesh LoopID to quadrilateral mesh LoopID*/
		std::map<int, int> triToQuadLoopMap;
		std::map<int, int> triToQuadArcMap;

		/*Mapping relationship: Triangle feature ring point -> Quadrilateral feature ring half-edge*/
		std::map<V*, H*> triPointToQuadHalfedgeMap;

		/*Mapping relationship: Data points -> Quadrilateral feature loop half-edge (including triangular mesh points and newly interpolated data points)*/
		std::map<QB_curve_sampling*, H*> dataPointToQuadHalfedgeMap;

		/*mapping relation: quad boundary halfedge to sampling on It*/
		std::map<H*, std::vector<QB_curve_sampling*>> samplingPts2QuadHalfedges;

		/*Mapping of halfedges of the quadrilateral mesh to IDs in a ring*/
		std::map<H*, int> quadHalfedgeToIdMap;
		// tri mesh halfedges id mapping to its loop
		std::map<H*, int> triHalfedgeToIdMap;

		//Mapping of quadrilateral points -> ID in the ring
		std::map<V*, int> quadLoopVertexIdToMap;

		//Mapping relationship: Triangle feature ring point ID -> Parameter
		std::map<int, double> triPointIdToParamMap;

		//Mapping relationship: Data point (sampling) ID -> Parameter
		std::map<int, double> dataPointIdToParamMap;

		//Store the mapping relationship between the halfedges of the quadrilateral mesh and the Bezier curve
		std::map<H*, BE_bezierCurve_approxi*> quadHalfedgeToBezierCurveMap;

		//Create the correspondence relationship between the IDs of the quadrilateral grids that store sharp feature points and their coordinates
		std::map<int, CPoint> sharpFeaturePoints;

		//Store the mapping between the original points and the fitted points
		std::map<int, std::pair<CPoint, CPoint>> originalToApproxPoints;

		/* recoard mapping between quadHalfedges and triHalfedges*/
		std::map<H*, std::pair<H*, H*>> mapQuadHalfedgesAndTriHalfedges;

		// 添加成员变量存储贴体性较差的边界边信息
		std::map<H*, std::pair<double, std::vector<QB_curve_sampling*>>> poorFittingEdges;

		// 添加端点统计map
		std::map<V*, int> endpointCountMap; // 记录顶点作为端点的出现次数

		std::vector<BoundaryFittingArc> m_boundaryFittingArcs;  // 边界特征弧集合

		// 使用unordered_map存储顶点映射（更高效）
		std::map<V*, int> m_boundaryVertexIdMap;

		// 四边形半边->环中ID的映射
		std::map<H*, int> quadBoundaryArcHalfedgeToIdMap;

		// 四边形边界特征弧的点->弧中ID的映射(以便列方程组时使用)
		std::map<V*, int> quadBoundaryArcVertexIdToMap;
	};

	/**
	 * The definition of the curve sampling point class.
	*/
	class QB_curve_sampling
	{
	public:
		QB_curve_sampling() {
			static int currentId = 0;//-1????
			qbc_vId = currentId++;  // Set the initial point ID to the current ID and increment the current ID by one.???
			qbc_mark = true;        // The sampling points are set as default and are available.
			for (int i = 0; i < qbc_ws.size(); i++) {
				qbc_ws[i] = 0.0;    // Initialize weights
			}
			qbc_quadLoopId = -1;
			qbc_quadArcId = -1;
			qbc_quadHalfedgesId = -1;
		}
		
		~QB_curve_sampling() {}

		//Access function
		int& vId() { return qbc_vId; }
		int& quadLoopId() { return qbc_quadLoopId; }
		int& quadArcId() { return qbc_quadArcId; }
		int& quadHalfedgesId() { return qbc_quadHalfedgesId; }
		CPoint& pos() { return qbc_pos; }
		CPoint& pos_approx() { return qbc_pos_approx; }
		CPoint& tangent() { return qbc_tangent; }
		std::vector<double>& ws() { return qbc_ws; }
		std::vector<int>& index() { return qbc_index; }
		std::vector<double>& weight() { return qbc_weight; }
		bool& mark() { return qbc_mark; }

	private:
		/*
		 * The unique ID of the corresponding point.
		 */
		int qbc_vId;

		/*
		 * The corresponding loop's ID of the point.
		 */
		int qbc_quadLoopId;

		/*
		 * The corresponding id of the arc for the point.
		 */
		int qbc_quadArcId;

		/*
		 * The corresponding halfedge ID of the loop for the point.
		 */
		int qbc_quadHalfedgesId;

		/*
		 * The location of the sampling points
		 */
		CPoint qbc_pos;

		/*
		 * The position obtained through sampling point fitting.
		 */
		CPoint qbc_pos_approx;

		/*
		 * The tangential vector of the sampling point (the direction of the curve)
		 */
		CPoint qbc_tangent;

		/*
		 * Weight parameter set
		 */
		std::vector<double> qbc_ws;

		/*
		 * Control point index set
		 */
		std::vector<int> qbc_index;

		/*
		 * The weight corresponding to the control point.
		 */
		std::vector<double> qbc_weight;

		/*
		 * Mark whether the current sampling point is valid (available).
		 */
		bool qbc_mark;
	};

	/**
	 * B spline curve.
	 */
	class BE_surface_model : public CCG_QMS_model
	{
	public:
		BE_surface_model() {};
		~BE_surface_model() {};
		void computeSamplingPosForTriangleFeatures(M* pMesh, CFeatureLoop* trifeatureLoop, std::map<int, double> ttp, std::map<V*, H*> ttq);

		//Obtain the UV parameter coordinates of the feature data points
		CPoint2 getUVFromSamplingPoint(V* v, std::map<int, double> ttp, std::map<V*, H*> ttq);
		// Access and modify the mapping between the original points and the fitted points
		std::map<int, std::pair<CPoint, CPoint>>& getOriginalToApproxPointsMap() {
			return originalToApproxPoints;
		}
	private:
		// Store the mapping of the original points and the fitted points.
		std::map<int, std::pair<CPoint, CPoint>> originalToApproxPoints;
	};
}
#endif // !_BSPLINECURVAPPROX_H_