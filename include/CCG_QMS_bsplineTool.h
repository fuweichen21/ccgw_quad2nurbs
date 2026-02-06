/*****************************************************************//**
 * \file   CCG_QMS_bsplineTool.h
 * \brief  some operation about b spline surface 
 * 
 * \author FuweiChen
 * \date   2025/6/30
 *********************************************************************/

#ifndef _CCG_QMS_BSPLINETOOL_H_
#define _CCG_QMS_BSPLINETOOL_H_

#include"CCG_QMS_quad2ManifoldSpline.h"
#include <unordered_map>

namespace CCG_QMSLib {

	/*Using the subPatchIndex attribute on the quadrilateral mesh as a bridge, establish the relationship between NURBS patches and T-layout*/
	void BPatchesTLayoutIndex(CCG_QMS_model& model, M* quadMesh);

	/*
	* Format：
	*###############################################################
	Total Number of TLayout Regions: xxx
	Total Number of Quad Cells: xxx
	Total Number of NURBS Patches: xxx
	Total Number of NURBS Control Points: xxx

	TLayout Region 1
	=================================================
		 Quad Cell 1

	---
			   NURBS Patch 1: xxx Control Points (xxx X xxx grid)
			   NURBS Patch 2: xxx Control Points (xxx X xxx grid)
			   …
	=================================================
		 Quad Cell 2

	---
			   NURBS Patch 1: xxx Control Points (xxx X xxx grid)
			   NURBS Patch 2: xxx Control Points (xxx X xxx grid)
			   …

	TLayout Region 2
	=================================================
		 Quad Cell 1

	---
			   NURBS Patch 1: xxx Control Points (xxx X xxx grid)
			   NURBS Patch 2: xxx Control Points (xxx X xxx grid)
			   …
	=================================================
		 Quad Cell 2

	---
			   NURBS Patch 1: xxx Control Points (xxx X xxx grid)
			   NURBS Patch 2: xxx Control Points (xxx X xxx grid)
			   …
	… …
	###############################################################

	---
	*/
	void outputRelationBSpline_TLayout_fixedFormat(const char* output, CCG_QMS_model& model, M* quadMesh);

	/**
	* .Determine the index of the node interval where parameter u is located
	*
	* \param n u-direction # control points
	* \param p	u-direction degree
	* \param u	u-dirction para
	* \param U	u-direction knot vector
	* \return
	*/
	int FindSpan1(int degree, const std::vector<double>& knots, double u);

	//Attribute to quantify the surface quality
	/*-----------------------------------------------------------------------------------------------*/
	/**
	 * Obtaing sampling points on boundary curve of b spline surface
	 * according to given step length.
	 * parallel for accelerating
	 */
	void obtainSamplingPoints_bsurfBoundarys_parrallel(CCG_QMS_model& model, double samplingStepLen);

	/**
	 *  Obtaing sampling points on boundary curve segment of b spline surface
	 * according to given step length.
	 */
	void obtainSamplingPoints_bsurfBoundarySegments(NUBE_bsurfBoudnary* bBoundaryCurve, NUBE_bsurfSegment* bBoundaryCurveSeg, const std::vector<std::vector<CPoint>>& cp, const std::vector<double>& kontsu, const std::vector<double>& kontsv, int degree_u, int degree_v, double stepL, double u_max, double v_max);

	/**
	 * Computing sampling points information (normal,du,dv,duu,duv,dvv...) from boundary curve segments.
	 */
	void computingSamplingInformation_bsurfBoundarySegment(NUBE_bsurfSegment* bBoundaryCurveSeg, std::shared_ptr<CCG_QMS_bsplineSurf> bsurf);

	/**
	 * obtaining sampling points's dual sampling point for boundary curve segment.
	 */
	void obtainingDualSamplingPoint_samplingPt_segment(NUBE_bsurfSegment* bBoundaryCurveSeg);

	/*
	* Use binary search to accurately locate the parameter range and avoid boundary errors in linear search.
	*/
	int FindSpan(int degree, const std::vector<double>& knots, double u);


	/**
	* Using the de Boor algorithm, calculate the B-spline surface point at the given parameters u and v.
	*/
	CPoint DeBoorAlgorithmSurface1(const std::vector<std::vector<CPoint>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, double u, double v);

	/**
	 * Calculate the principal curvatures and normal vectors at a given point on the B-spline surface.
	 * bpsur_m1: principal curvature 1
	 * bpsur_m2: principal curvature 2
	 * bpser_n: normal vector
	 * n: Number of control points in the U direction
	 * p: Degree of U direction
	 * U: Knot vector of U direction
	 * m: Number of control points in the V direction
	 * q: Degree of V direction
	 * V: Knot vector of V direction
	 * P： control points coordinate
	 * u: parameterization in U direction
	 * v: parameterization in V direction
	 * d: The order of the highest partial derivative
	 * S: Save the partial derivative vectors from order 0 to order d
	 */
	void getMaincurvature_normal(double& k1, double& k2, CPoint& m1, CPoint& m2, CPoint& bpser_n, int n, int p, std::vector<double> U, int m, int q, std::vector<double> V, std::vector<std::vector<CPoint>> P, double u, double v, int d, std::vector<std::vector<CPoint>>& S);
	/**
	 * 曲面直到d次的偏导.
	 * n: U方向控制点数量
	 * p: U方向次数
	 * U: u方向节点向量
	 * m: V方向控制点数量
	 * q: V方向次数
	 * V: V方向节点向量
	 * P： 控制点坐标
	 * u: U方向参数值
	 * v: V方向参数值
	 * d: 最高偏导数的次数
	 * S: 保存从0到d阶的偏导矢
	 */
	void getsurfaceDerivat_allDegree(int n, int p, std::vector<double> U, int m, int q, std::vector<double> V, std::vector<std::vector<CPoint>> P, double u, double v, int d, std::vector<std::vector<CPoint>>& S);

	/**
	* .计算非零B样条基函数及其导数
	*
	* \param i 第几个B样条基函数
	* \param u u方向参数值
	* \param p B样条基函数次数
	* \param n u方向控制点数量
	* \param U U方向的节点向量
	* \param Der 存储基函数的导数
 */
	void DersBasisFuns(int i, double u, int p, int n, std::vector<double>U, std::vector<std::vector<double>>& Der);

	/**
	* .确定参数u所在的节点区间的下标
	*
	* \param n u方向的控制点数量
	* \param p	u方向的次数
	* \param u	u方向参数值
	* \param U	u方向节点向量
	* \return
	*/
	int FindSpan(int n, int p, double u, std::vector<double> U);

	// statistics max difference for all boundary curve of G0 with T node
	void findMaxPositionDiffForEachEdge_T(CCG_QMS_model& model, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxPosDiff);

	/**
	 * statistics max difference for all boundary curve of G1 with T node.
	 */
	void findMaxAngleForEachEdge_T(CCG_QMS_model& model, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxAngles);

	/**
	 * statistics max difference for all boundary curve of G2 with T node.
	 */
	void findMaxSecondDerivativeDiffForEachEdge_T(
		CCG_QMS_model& model,
		std::unordered_map<NUBE_bsurfBoudnary*, double>& maxDuuDiff,
		std::unordered_map<NUBE_bsurfBoudnary*, double>& maxDuvDiff,
		std::unordered_map<NUBE_bsurfBoudnary*, double>& maxDvvDiff);

	// Calculate the ratio of the longest to the shortest boundary edge length of the spline surface patch
	void computeBoundaryEdgeRatio(const std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>>& allSurfs);

	/**
	 * For each QB_bsplineSurf, calculate the minimum Euclidean distance between each pair of its control points
	 * hen, the result is stored in the minCtrlPtDist() of this surface patch
	 *
	 * allSurfs：A container containing all the QB_bsplineSurf* items in the project
	 */
	void computeControlPointMinDistance(const std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>>& allSurfs);

	// A template function for converting a map to a vector, sorting it by value in descending order, and writing to a file
	template<typename T>
	void dumpSorted(const std::unordered_map<CEdge*, T>& m, const std::string& filename) {
		// 1) 转成 vector<pair<CEdge*, T>>  
		std::vector<std::pair<CEdge*, T>> v(m.begin(), m.end());
		// 2) 按第二项（差值）降序排序  
		std::sort(v.begin(), v.end(),
			[](auto& a, auto& b) { return a.second > b.second; }
		);
		// 3) 写文件  
		std::ofstream ofs(filename);
		if (!ofs) {
			std::cerr << "Cannot open " << filename << " for writing\n";
			return;
		}
		for (auto& p : v) {
			ofs << p.first->id()  // Edge id  
				<< "\t"
				<< p.second       // 对应的最大差值或角度  
				<< "\n";
		}
		ofs.close();
	}
	template<typename T>
	void dumpSorted_T(const std::unordered_map<NUBE_bsurfBoudnary*, T>& m, const std::string& filename)
	{
		// 1) 转成 vector<pair<NUBE_bsurfBoudnary*, T>>  
		std::vector<std::pair<NUBE_bsurfBoudnary*, T>> v(m.begin(), m.end());
		// 2) 按第二项（差值）降序排序  
		std::sort(v.begin(), v.end(),
			[](auto& a, auto& b) { return a.second > b.second; }
		);
		// 3) 写文件  
		std::ofstream ofs(filename);
		if (!ofs) {
			std::cerr << "Cannot open " << filename << " for writing\n";
			return;
		}
		for (auto& p : v) {
			ofs << p.first->bsurfId()  // bsurf id  
				<< "\t"
				<< p.first->localId()  // local id  
				<< "\t"
				<< p.second       // 对应的最大差值或角度  
				<< "\n";
		}
		ofs.close();
	}

	// 将点投影到三角形上
	CPoint2 projectPointToTriangle(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c);

	// 将点投影到线段上
	CPoint2 projectPointToSegment(const CPoint2& p, const CPoint2& a, const CPoint2& b);

	// 计算点到线段的最短距离
	double distanceToSegment(const CPoint2& p, const CPoint2& a, const CPoint2& b);

	double distanceToTriangle(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c);

	/**
		* 计算点在三角形中的重心坐标
		*/
	void computeBarycentricCoordinates(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c,
		double& u, double& v, double& w);

	/**
		* 计算采样点在四边形四个顶点上的权重
		*/
	void computeQuadWeights(M* quadMesh, F* quadFace, const std::vector<H*>& fhs, const CPoint2& uv, QB_sampling& samplingPoint);


	/*
		* 计算三个点构成的三角形的有向面积
		*/
	double signedArea(CPoint2 p0, CPoint2 p1, CPoint2 p2);

	/**
		* 判断点是否在三角形内部
		*/
	bool isPointInTriangle(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c);

	/**
		* 通过UV参数坐标在三角网格中定位对应点并计算其空间坐标(如果在三角网格找不到的点使用距离最近的三角网格网格点来计算)
		* @param triMesh 三角网格
		* @param quadMesh 四边形网格
		* @param quadFace 四边形面
		* @param uv 需要定位的UV参数坐标
		* @param samplingPoint 用于存储结果的采样点结构
		* @return 是否成功找到对应的三角网格点
		*/
	bool locateTrianglePointByUV1(M* triMesh, M* quadMesh, F* quadFace, const CPoint2& uv, QB_sampling& samplingPoint);

	void sampleQuadFace(
		M* triMesh,      // 已经全局参数化过的三角网格
		M* quadMesh,     // 四边形网格
		F* quadFace,     // 待采样的四边形面指针
		const std::vector<CPoint2>& sampleUVs,  // 参数域上的若干 (u,v) 采样坐标
		std::vector<QB_sampling>& outSamplings // 输出的采样点列表
	);

	/*
		 * 对每个四边形面检查已有采样点数量，如果少于 minSamples，就在该面的参数域上补足采样。
		 *  - triMesh: 已全局参数化的三角网格
		 *  - quadMesh: 四边形网格
		 *  - minSamples: 每个四边形面至少要有的采样点数
		 *  - sts: 输入/输出的采样列表（QB_sampling）。函数会统计其中每个面已有的采样点数，
		 *         并在不足时调用 FitTool::sampleQuadFace 生成新的采样并追加到 sts。
		 */
	void ensureQuadFaceSampling(
		M* triMesh,
		M* quadMesh,
		int                           minSamples,
		std::vector<QB_sampling>& sts);

	// 定义一个结构体来存储控制点和相关数据
	struct ControlPointData {
		CPoint point1;
		CPoint point2;
		int he1_localId;
		int he2_localId;
		int f1_id;
		int f2_id;
		E* e;
	};

	/*-----------------------------------------------------------------------------------------------*/
}

#endif // !_CCG_QMS_BSPLINETOOL_H_

