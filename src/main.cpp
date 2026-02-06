#include<iostream>
#include"NUBE_printInformation.h"
#include"CCG_QMS_meshOperation.h"
#include"CCG_QMS_quad2ManifoldSpline.h"
#include"CCG_QMS_bsplineTool.h"
#include"BsplingCuvApprox.h"
#include"CCG_QMS_meshView.h"

#ifdef _WIN32

#include <windows.h>

#include <psapi.h>
#include "main.h"

#elif __linux__

#else

#endif

static size_t _Ask_OS_memory_usage()
{

#ifdef _WIN32

	DWORD currentProcessID = GetCurrentProcessId();

	HANDLE hProcess;

	PROCESS_MEMORY_COUNTERS pmc;



	// Open the process with necessary access rights

	hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, currentProcessID);

	if (NULL == hProcess)

		return 0;

	if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))

		return pmc.WorkingSetSize;

	else

		return 0;

#elif __linux__

	return 0;

#else

#endif

}

#define CRTDBG MAP ALLOC
using namespace CCG_QMSLib;
using namespace NUBELib;

// test quad mesh 2 manifold spline
// control print information
// bi 5
// output boundary curve relation about initial vertice on tri mesh
// spline surface quality improvment
// face 1-1 patch done 
// .txt 0 0 1 0 to 1 0 1 0 
// no boundary curve approximation and output hausdorff distance result of initial trimesh boundary vertex (sampling points on boundary edges)
// add surface fitting
// add feature line approximation
void test250407_594_3_3_12512192258(int argc, char* argv[])
{
	std::string prefix_FilePath_TriangleMesh = argv[1];
	std::string prefix_FilePath_QuadMesh = argv[2];
	bool printToscreen_Log = true;
	FILE* filePointer_Log = NULL;
	double _approxi_thread = 1e-5;
	double _approxi_smooth_weight = 1.0;
	bool ifApproximation = true;
	NUBE_printInformation myPrint;
	clock_t start, finish;
	double totaltime;
	start = clock();
	/*读取四边形网格 */
	M quadMesh;
	quadMesh.read_m(prefix_FilePath_QuadMesh.c_str());
	arrangeMeshVertexOrder(&quadMesh);
	std::string print1 = "# V on QuadMesh: " + std::to_string(quadMesh.numVertices()) + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	print1 = "# F on QuadMesh: " + std::to_string(quadMesh.numFaces()) + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	for (auto v : It::MVIterator(&quadMesh))
	{
		if (v->feature())
		{
			//std::cout << "feature: " << v->id() << std::endl;
			v->feature() = false;
		}
	}
	//std::cout << "# boundary V on quadMesh: " << statisticBoundaryVertexNum(&quadMesh) << std::endl;
	//isParrallelogram_uv(&quadMesh);
	M triMesh;
	triMesh.read_m(prefix_FilePath_TriangleMesh.c_str());
	//std::cout << "triMesh #V: " << triMesh.numVertices() << " #E: " << triMesh.numEdges() << " #F: " << triMesh.numFaces() << std::endl;
	/*-------------------------------------------特征曲线约束开始------------------------------------------------------------------*/
	//copyInitialMeshVPos(&quadMesh);//cfw add 202404040107,为了固定边界平坦地方的网格点坐标
	//// 创建 CFeatureLoop 对象
	//CFeatureLoop* quadfeatureLoop = new CFeatureLoop(&quadMesh);
	////可视化四边形网格和特征边
	//viewMesh_featureLoop(&quadMesh, *quadfeatureLoop);

	////三角网格特征提取
	//CFeatureLoop* trifeatureLoop = new CFeatureLoop(&triMesh);
	//viewMesh_featureLoop(&triMesh, *trifeatureLoop);
	//// 获取四边形特征边环的数量
	//std::cout << "=========================================================" << std::endl;
	//std::cout << "Number of quadfeature loops: " << quadfeatureLoop->loops().size() << std::endl;
	////手动遍历特征边环
	//for (auto loop : quadfeatureLoop->loops()) {
	//	std::cout << " Loop contains " << loop->halfedges().size() << " half-edges." << std::endl;
	//}
	//std::cout << "Number of quadfeature arcs: " << quadfeatureLoop->arcs().size() << std::endl;
	//for (auto arc : quadfeatureLoop->arcs()) {
	//	std::cout << " Arc contains " << arc->halfedges().size() << " half-edges." << std::endl;
	//}
	//// 获取三角形特征边环的数量
	//std::cout << "Number of trifeature loops: " << trifeatureLoop->loops().size() << std::endl;
	////手动遍历特征边环
	//for (auto loop : trifeatureLoop->loops()) {
	//	std::cout << " Loop contains " << loop->halfedges().size() << " half-edges." << std::endl;
	//}
	//std::cout << "Number of trifeature arcs: " << trifeatureLoop->arcs().size() << std::endl;
	//for (auto arc : trifeatureLoop->arcs()) {
	//	std::cout << " Arc contains " << arc->halfedges().size() << " half-edges." << std::endl;
	//}

	////创建一个bezier样条曲面模型
	//BE_curve_model model;
	//model.quadFeatureLoops() = quadfeatureLoop;
	//model.triFeatureLoops() = trifeatureLoop;
	////可视化三角形网格和特征边
	//viewMesh_featureLoop(&triMesh, *trifeatureLoop);
	////匹配三角网格特征线和四边形网格特征线
	//model.mapTriLoopToQuadLoop(&quadMesh, &triMesh);
	////可视化特征线匹配关系
	////viewMesh_mapFeatureLoop(&quadMesh, *quadfeatureLoop, *trifeatureLoop, model);

	////匹配三角形网格点到四边形特征环中半边
	//model.triPointToQuadHalfedge(&quadMesh, &triMesh);
	////可视化三角形网格点到四边形特征环中半边
	//std::vector<QB_curve_sampling*> qbm_samplings;
	//qbm_samplings = model.samplings();
	//std::map<V*, H*> pointToHalfedgeMap;
	//pointToHalfedgeMap = model.getTriPointToQuadHalfedgeMap();
	//std::map<H*, int> quadHalfedgeInLoopID;
	//quadHalfedgeInLoopID = model.getQuadHalfedgeToIdMap();
	//std::cout << "qbm_samplings.size()_before" << qbm_samplings.size() << std::endl;

	//// 在调用函数时
	////std::map<H*, V*> quadToTriMap;  // 创建映射容器
	////model.addInterpolatedPointsForUnmappedHalfedges(&quadMesh, quadToTriMap);
	////viewMesh(&triMesh, &quadMesh, quadToTriMap);

	//// 添加新功能：为未映射的半边添加插值点
	//model.addInterpolatedPointsForUnmappedHalfedges(&quadMesh);
	//model.addInterpolatedPointsForFeatureArcs1(&quadMesh);

	//qbm_samplings = model.samplings();
	//pointToHalfedgeMap = model.getTriPointToQuadHalfedgeMap();
	//quadHalfedgeInLoopID = model.getQuadHalfedgeToIdMap();
	//std::cout << "qbm_samplings.size()_after" << qbm_samplings.size() << std::endl;
	////viewMesh_triPointToQuadHalfedge(&quadMesh, qbm_samplings, pointToHalfedgeMap, quadHalfedgeInLoopID, *quadfeatureLoop);

	////建立特征线的初始bezier控制点(尖锐点保持）
	//model.setControlPoints(&quadMesh);
	//model.markFeaturePts_input(&quadMesh);
	//std::vector<BE_bezierCurve_approxi*> be_bcs;
	//be_bcs = model.bcs();
	//std::cout << "be_bcs.size():" << be_bcs.size() << std::endl;
	////可视化四边形网格特征边的bezier控制点
	////viewMesh_featureLoopCPts(&quadMesh, *quadfeatureLoop, be_bcs);

	////重新确定四边形网格点
	//model.initializeQuadVertexIdMap(&quadMesh);	//特征线上四边形网格点重新排序(四边形网格点到特征环中点的新id映射）
	//model.parametrize2(&quadMesh);//数据点参数化

	////viewFeatureBezierCurve(be_bcs, &triMesh);

	//// 创建存储排序距离的向量
	//std::vector<std::pair<int, double>> sortedCurveDistances;
	////计算三角形网格边到样条曲线的hausdroff距离
	//double hausdorffCurveDistBefore = model.computeCurveHausdorffDistance(&triMesh, 0.05, sortedCurveDistances);
	//std::cout << "hausdorffCurveDistBefore: " << hausdorffCurveDistBefore << endl;

	////model.compute_quad_vertex(&quadMesh);//更新四边形网格点
	////model.compute_quad_vertex1(&quadMesh);//更新四边形网格点,天然正则化
	////model.compute_quad_vertex2(&quadMesh);//更新四边形网格点，加速求解eigen库
	//model.compute_quad_vertex7(&quadMesh);//最初版本
	//model.preserveQuadCornerPoints(&quadMesh);//保持四边形网格尖锐点
	//recoverFlatBoundaryVertexPoint(&quadMesh);//cfw add 202504040124 将边界平坦点保持初始坐标
	//model.CheckDegeneratedBoudnaryVertexAndUpdate(&quadMesh);//cfw add 202507112138

	////viewMesh(&quadMesh);
	////可视化新四边形网格
	//qbm_samplings = model.samplings();
	//pointToHalfedgeMap = model.getTriPointToQuadHalfedgeMap();
	//quadHalfedgeInLoopID = model.getQuadHalfedgeToIdMap();
	////viewMesh_triPointToQuadHalfedge(&quadMesh, qbm_samplings, pointToHalfedgeMap, quadHalfedgeInLoopID, *quadfeatureLoop);

	////可视化三角形网格和对应的bezier样条曲线
	//model.quadFeatureLoops() = quadfeatureLoop;
	////更新bezier控制点
	//model.resetControlPoints(&quadMesh);

	//be_bcs = model.bcs();
	////可视化四边形网格特征边的bezier控制点
	////viewMesh_featureLoopCPts(&quadMesh, *quadfeatureLoop, be_bcs);
	//std::cout << "be_bcs.size():" << be_bcs.size() << std::endl;
	////model.computeTotalError();//计算采样点拟合误差

	////计算三角形网格边到样条曲线的hausdroff距离
	//double hausdorffCurveDist = model.computeCurveHausdorffDistance(&triMesh, 0.05, sortedCurveDistances);
	//std::cout << "hausdorffCurveDistAfter: " << hausdorffCurveDist << endl;

	////model.computeTriTotalError();//计算三角形数据点拟合误差

	////viewFeatureBezierCurve(be_bcs, &triMesh);

	////viewMesh(&quadMesh);
	////viewFeatureBezierCurve(be_bcs, &quadMesh2);



	////清理内存
	////delete quadfeatureLoop;
	////delete trifeatureLoop;

	/*-------------------------------------------特征曲线约束结束------------------------------------------------------------------*/
	markFeatureEdgesBySharpEdges(&quadMesh);
	//viewMesh(&quadMesh, argc, argv);

	markTLayout_quadMesh(&quadMesh);
	print1 = "# tLayout on quadMesh: " + std::to_string(statisticTLayoutPatchNum(&quadMesh)) + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	statisticTLayoutPatch_quadNum(&quadMesh);
	//viewMesh(&quadMesh);

	std::string filename = prefix_FilePath_QuadMesh;
	/*--------------------------------------------*/
	if (ifApproximation)
	{
		//preprocessing quad mesh
		/*------------------calculate some quality on quad mesh start-----------------------------------------------------*/
		calculateScaledJacobian_quadMesh(&quadMesh);
		calculateEdgeRatio_quadMesh(&quadMesh);
		for (auto f : It::MFIterator(&quadMesh))
		{
			double sj = f->quality_scaledJocbian();
			double er = f->quality_edgeRatio();
			if (sj < 0.5 || er > 10)
			{
				for(auto fv:It::FVIterator(&quadMesh,f))
				{
					fv->fixed() = true;
					//std::cout << "bad face vertex id: " << fv->id() << std::endl;
				}
			}
			//std::cout << "face id: " << f->id() << " scaledJacobian: " << sj << " edgeRatio: " << er << std::endl;
		}
		/*------------------calculate some quality on quad mesh end-----------------------------------------------------*/
		BE_surface_model model_surf_approx;
		clock_t start4, finish4;
		double totaltime4;
		start4 = clock();
		model_surf_approx.obtainCtlPts(&quadMesh);
		compute_mesh_normal(&quadMesh);
		//model_surf.buildBezierSurfs(&quadMesh);
		model_surf_approx.buildBezierSurfs_cw(&quadMesh);
		remarkSingularPts(&quadMesh);
		//model.raiseBezierSurfDegree3To5(&quadMesh);
		model_surf_approx.raiseBezierSurfDegree3To5_all(&quadMesh);
		model_surf_approx.sortBezierCtlPts_cw(&quadMesh);
		computeV_degree_mesh(&quadMesh);
		model_surf_approx.G1Constraints_optimization_cws(&quadMesh);
		finish4 = clock();
		totaltime4 = (double)(finish4 - start4) / CLOCKS_PER_SEC;
		//std::cout << "Obtain G1 NURBS is finished! It takes " << totaltime4 << "s！" << std::endl;
		model_surf_approx.computeBezierSurfacesKnots(&quadMesh);
		model_surf_approx.sortBezierCtlPts(&quadMesh);
		int if_fixedBoundary = 1;
		/*approximation*/
		clock_t start3, finish3;
		double totaltime3;
		start3 = clock();
		/*读取初始三角网格 */
		//M triMesh;
		//triMesh.read_m(prefix_FilePath_TriangleMesh.c_str());
		//std::cout << "# V on triMesh: " << triMesh.numVertices() << std::endl;
		//std::cout << "# boundary V on triMesh: " << statisticBoundaryVertexNum(&triMesh) << std::endl;
		/*对网格中的参数uv进行一定数量级的放大*/
		//magnifyMeshUV(&quadMesh, 1e5);
		//magnifyMeshUV(&triMesh, 1e5);
		/*通过参数化结果在四边形网格上定位出三角形网格点，获取对应的采样点信息*/
		std::vector<QB_sampling> sts;
		obtainSamplingFromTriMeshByParameterization(&triMesh, &quadMesh, sts);
		print1 = "# samplingPts from tri: " + std::to_string(sts.size()) + "\n";
		myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);

		if (triMesh.patchMap_mesh().size() == 0)
		{
			triMesh.classifyFacesByPatch();
		}
		ensureQuadFaceSampling(&triMesh, &quadMesh, 20, sts);  // 确保每一片区域采样点足够

		/*部分四边形网格点作为采样点*/
		//bool part_sampling_quad = 0;
		//obtainSampling_quad(&quadMesh, sts, part_sampling_quad);

		print1 = "# samplingPts: " + std::to_string(sts.size()) + "\n";
		myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
		//std::cout << "# boundary sampling: " << statisticBoundarySamplingNum(sts) << std::endl;
		//viewMesh(&triMesh, &quadMesh, argc, argv);
		//viewMesh(&triMesh,argc,argv);
		//viewSampling(&quadMesh, sts);
		for (int i = 0; i < sts.size(); i++)
		{
			model_surf_approx.samplings().push_back(std::make_shared<QB_sampling>(sts[i]));
		}
		model_surf_approx.computeSamplingPtUV(&quadMesh);
		//statisticSampleNum_model(&quadMesh, model.samplings());
		finish3 = clock();
		totaltime3 = (double)(finish3 - start3) / CLOCKS_PER_SEC;
		std::cout << "Computing sampling UV is finished! It takes " << totaltime3 << "s！" << std::endl;
		clock_t start5, finish5;
		double totaltime5;
		start5 = clock();
		//model.approximateSurface();
		//model.approximateSurface_optimize3(&quadMesh);
		/*读取用户输入拟合精度*/
		double approxi_thread = _approxi_thread;
		/*读取用户输入光顺权重*/
		double approxi_smooth_weight = _approxi_smooth_weight;
		if (if_fixedBoundary == 1)
		{
			model_surf_approx.approximateSurface_optimize_check(&quadMesh, approxi_thread, approxi_smooth_weight);
		}
		else
		{
			model_surf_approx.approximateSurface_optimize_check2(&quadMesh, approxi_thread, approxi_smooth_weight);
		}
		//model.approximateSurface_optimize_check(&quadMesh,approxi_thread,approxi_smooth_weight);
		//model.approximateSurface_optimize_check2(&quadMesh, approxi_thread, approxi_smooth_weight);
		//model.computeSamplingPos_approx();
		//model.approximateSurface_boundary(&quadMesh, approxi_thread, approxi_smooth_weight);
		finish5 = clock();
		totaltime5 = (double)(finish5 - start5) / CLOCKS_PER_SEC;
		print1 = "Approximation is finished! It takes "  + std::to_string(totaltime5) + "s！" + "\n";
		myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
		updateMeshPosByControlMesh(&quadMesh, model_surf_approx);
		std::string outputApproxQuadMeshFileName = filename.substr(0, filename.length() - 2) + "_approxQuadMesh.m";
		quadMesh.write_m(outputApproxQuadMeshFileName.c_str());
	}
	/*--------------------------------------------*/

	//创建一个bezier样条曲面模型
	BE_surface_model model_surf;
	clock_t start4, finish4;
	double totaltime4;
	start4 = clock();
	model_surf.obtainCtlPts(&quadMesh);
	compute_mesh_normal(&quadMesh);
	//model_surf.buildBezierSurfs(&quadMesh);
	model_surf.buildBezierSurfs_cw(&quadMesh);
	remarkSingularPts(&quadMesh);
	//model.raiseBezierSurfDegree3To5(&quadMesh);
	model_surf.raiseBezierSurfDegree3To5_all(&quadMesh);
	std::string outputBSplineFileName = filename.substr(0, filename.length() - 2) + "_G0_bspline_tri_q.igs";
	model_surf.sortBezierCtlPts_cw(&quadMesh);
	computeV_degree_mesh(&quadMesh);
	model_surf.G1Constraints_optimization_check(&quadMesh);
	outputBSplineFileName = filename.substr(0, filename.length() - 2) + "_G1_CCG_bspline_tri_q.igs";
	finish4 = clock();
	totaltime4 = (double)(finish4 - start4) / CLOCKS_PER_SEC;
	//std::cout << "Obtain G1 NURBS is finished! It takes " << totaltime4 << "s！" << std::endl;
	model_surf.computeBezierSurfacesKnots(&quadMesh);
	model_surf.sortBezierCtlPts(&quadMesh);
	model_surf.obtainBezierCptsBy_cpts(1.0);

	//obtainQuadLayout_extend2(&quadMesh);
	obtainTLayout_motorcycleGraph_featureAware(&quadMesh);
	//obtainLayout_consistencyQuadMeshFaceNum(&quadMesh);
	//viewMesh(&quadMesh);
	//model.obtainBSplineSurfacePatches_general(&quadMesh);
	//model.obtainBSplineSurfacePatches_general_UVConsistency(&quadMesh);
	//model_surf.obtainBSplineSurfacePatches_general_UVConsistency_recordQuadBoundaryNURBSBoundaryCurveRelation(&quadMesh);
	model_surf.obtainBSplineSurfacePatches_general_UVConsistency_topology_boundary(&quadMesh);
	model_surf.removeMultiKnots_BSplines();
	model_surf.computeBSplineSurfaceUniqueKnotMultiNum();

	BPatchesTLayoutIndex(model_surf, &quadMesh);

	//model_surf.outputBSpline_iges(outputBSplineFileName.c_str());
	model_surf.outputBSpline_iges_5_3(outputBSplineFileName.c_str());
	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	print1 = "Reconstructing NURBS is finished! It takes " + std::to_string(totaltime) + "s！" + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	//model.outputBezierSurfaces_iges(&quadMesh, outputBezierFileName.c_str(), 3, 3);
	std::string outputBSplineSurface_fixedFormatFileName = filename.substr(0, filename.length() - 2) + "_CCG_spline_tri_q.txt";
	outputBSplineSurface_fixedFormatFileName = filename.substr(0, filename.length() - 2) + "_CCG_spline_tri_q.txt";
	model_surf.outputBSplineSurf_fixedFormat3(outputBSplineSurface_fixedFormatFileName.c_str());
	//model.readBSplineSurf_fixedFormat(outputBSplineSurface_fixedFormatFileName.c_str());
	std::string outputRelationBSpline_TLayout_fixedFormatFileName = filename.substr(0, filename.length() - 2) + "_CCG_relation_tri_q.txt";
	outputRelationBSpline_TLayout_fixedFormatFileName = filename.substr(0, filename.length() - 2) + "_CCG_relation_tri_q.txt";
	outputRelationBSpline_TLayout_fixedFormat(outputRelationBSpline_TLayout_fixedFormatFileName.c_str(), model_surf, &quadMesh);

	/*-------------------------------------样条曲面质量分析开始------------------------------------------start*/
	//model_surf.create_bsurBoundarys_general(&quadMesh);

	////model.createTopology_bspline(&quadMesh);
	//std::string outputBsplineTopologyTMeshFileName = filename.substr(0, filename.length() - 2) + "_BsplineTopologyTMesh.obj";
	//model_surf.outputTopology_bspline_TMesh_obj(&quadMesh, outputBsplineTopologyTMeshFileName.c_str());
	//model_surf.createTopology_bspline_TMesh(&quadMesh);
	////std::cout << "111" << std::endl;
	//std::string outputBsplineTopologyQuadMeshFileName = filename.substr(0, filename.length() - 2) + "_BsplineTopologyQuad.m";
	///*for (auto v : It::MVIterator(&model_surf.bsplineSurfacesTopology_quad()))
	//{
	//	v->point() = v->point() / scale;
	//}*/
	//model_surf.bsplineSurfacesTopology_quad().write_m(outputBsplineTopologyQuadMeshFileName.c_str());

	//std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>> allSurfs = model_surf.bsplineSurfs();
	//// ===================== 控制点一致性检查 ======================
	////std::vector<ControlPointData> unequalControlPoints;
	////checkContronPoints(model_surf, unequalControlPoints);
	////std::cout << "-----unequalControlPoints: " << unequalControlPoints.size() << std::endl;

	////magnifyMeshPoint(&quadMesh, 1.0 / scale);
	////viewControlPoints(&quadMesh, unequalControlPoints);
	//// ============================================================

	//// ===================== 边界采样点连续性分析 ==================
	//double stepL = 0.5;
	//obtainSamplingPoints_bsurfBoundarys_parrallel(model_surf, stepL);
	//// ===================== G0连续性分析 ==========================
	//std::unordered_map<NUBE_bsurfBoudnary*, double> maxPosDiff;
	//findMaxPositionDiffForEachEdge_T(model_surf, maxPosDiff);

	//// ===================== G1连续性分析 ==========================
	//std::unordered_map<NUBE_bsurfBoudnary*, double> maxAngles;
	//findMaxAngleForEachEdge_T(model_surf, maxAngles);

	//// ===================== G2连续性分析 ==========================
	//std::unordered_map<NUBE_bsurfBoudnary*, double> maxDuuDiff, maxDuvDiff, maxDvvDiff;
	//findMaxSecondDerivativeDiffForEachEdge_T(model_surf, maxDuuDiff, maxDuvDiff, maxDvvDiff);
	//// ============================================================

	//// 1、计算并打印所有曲面片的边界边长比
	//computeBoundaryEdgeRatio(allSurfs);

	//// 2、计算控制点之间的最小距离 s
	//computeControlPointMinDistance(allSurfs);

	//// 4. 调用对每条 NURBS 边界计算 Hausdorff 距离
	////computeHausdorffPerBoundary(&triMesh, allSurfs, sampleStep);

	//// 创建存储排序距离的向量
	//std::vector<std::pair<int, double>> sortedDistances;

	//// 计算Hausdorff距离
	//double hausdorffDist = model_surf.computeHausdorffDistance_allNURBSBoundaryCurve3(&triMesh, 3, sortedDistances);

	//// 输出排序后的距离到文件
	//std::string outputDistanceFile = filename.substr(0, filename.length() - 2) + "_hausdorff_distances.txt";
	//model_surf.outputSortedDistances(sortedDistances, outputDistanceFile);

	//model_surf.computeGlobalSurfaceMetrics();

	//// =================== 分析数据输出到.txt ======================
	//std::string maxPosDiffFileName = filename.substr(0, filename.length() - 2) + "_maxPosDiff.txt";
	//std::string maxAnglesFileName = filename.substr(0, filename.length() - 2) + "_maxAngles.txt";
	////std::string maxDuuDiffFileName = filename.substr(0, filename.length() - 2) + "_maxDuuDiff.txt";
	////std::string maxDuvDiffFileName = filename.substr(0, filename.length() - 2) + "_maxDuvDiff.txt";
	////std::string maxDvvDiffFileName = filename.substr(0, filename.length() - 2) + "_maxDvvDiff.txt";
	//dumpSorted_T(maxPosDiff, maxPosDiffFileName.c_str());
	//dumpSorted_T(maxAngles, maxAnglesFileName.c_str());
	////dumpSorted_T(maxDuuDiff, maxDuuDiffFileName.c_str());
	////dumpSorted_T(maxDuvDiff, maxDuvDiffFileName.c_str());
	////dumpSorted_T(maxDvvDiff, maxDvvDiffFileName.c_str());
	//std::string outputBsurfMetricsName = filename.substr(0, filename.length() - 2) + "_report.txt";
	//model_surf.outputMetricsToFile(outputBsurfMetricsName.c_str());
	/*-------------------------------------样条曲面质量分析结束------------------------------------------start*/

}

// test quad mesh 2 manifold spline
// control print information
// bi 5
// output boundary curve relation about initial vertice on tri mesh
// spline surface quality improvment
// face 1-1 patch done 
// .txt 0 0 1 0 to 1 0 1 0 
// no boundary curve approximation and output hausdorff distance result of initial trimesh boundary vertex (sampling points on boundary edges)
// add surface fitting
// add feature line approximation
void test250407_594_3_3(int argc, char* argv[])
{
	std::string prefix_FilePath_TriangleMesh = argv[1];
	std::string prefix_FilePath_QuadMesh = argv[2];
	bool printToscreen_Log = true;
	FILE* filePointer_Log = NULL;
	double _approxi_thread = 1e-5;
	double _approxi_smooth_weight = 1.0;
	bool ifApproximation = true;
	NUBE_printInformation myPrint;
	clock_t start, finish;
	double totaltime;
	start = clock();
	/*读取四边形网格 */
	M quadMesh;
	quadMesh.read_m(prefix_FilePath_QuadMesh.c_str());
	std::string print1 = "# V on QuadMesh: " + std::to_string(quadMesh.numVertices()) + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	print1 = "# F on QuadMesh: " + std::to_string(quadMesh.numFaces()) + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	for (auto v : It::MVIterator(&quadMesh))
	{
		if (v->feature())
		{
			//std::cout << "feature: " << v->id() << std::endl;
			v->feature() = false;
		}
	}
	//std::cout << "# boundary V on quadMesh: " << statisticBoundaryVertexNum(&quadMesh) << std::endl;
	//isParrallelogram_uv(&quadMesh);
	//arrangeMeshVertexOrder(&quadMesh);
	M triMesh;
	triMesh.read_m(prefix_FilePath_TriangleMesh.c_str());
	//std::cout << "triMesh #V: " << triMesh.numVertices() << " #E: " << triMesh.numEdges() << " #F: " << triMesh.numFaces() << std::endl;
	/*-------------------------------------------特征曲线约束开始------------------------------------------------------------------*/
	//copyInitialMeshVPos(&quadMesh);//cfw add 202404040107,为了固定边界平坦地方的网格点坐标
	//// 创建 CFeatureLoop 对象
	//CFeatureLoop* quadfeatureLoop = new CFeatureLoop(&quadMesh);
	////可视化四边形网格和特征边
	//viewMesh_featureLoop(&quadMesh, *quadfeatureLoop);

	////三角网格特征提取
	//CFeatureLoop* trifeatureLoop = new CFeatureLoop(&triMesh);
	//viewMesh_featureLoop(&triMesh, *trifeatureLoop);
	//// 获取四边形特征边环的数量
	//std::cout << "=========================================================" << std::endl;
	//std::cout << "Number of quadfeature loops: " << quadfeatureLoop->loops().size() << std::endl;
	////手动遍历特征边环
	//for (auto loop : quadfeatureLoop->loops()) {
	//	std::cout << " Loop contains " << loop->halfedges().size() << " half-edges." << std::endl;
	//}
	//std::cout << "Number of quadfeature arcs: " << quadfeatureLoop->arcs().size() << std::endl;
	//for (auto arc : quadfeatureLoop->arcs()) {
	//	std::cout << " Arc contains " << arc->halfedges().size() << " half-edges." << std::endl;
	//}
	//// 获取三角形特征边环的数量
	//std::cout << "Number of trifeature loops: " << trifeatureLoop->loops().size() << std::endl;
	////手动遍历特征边环
	//for (auto loop : trifeatureLoop->loops()) {
	//	std::cout << " Loop contains " << loop->halfedges().size() << " half-edges." << std::endl;
	//}
	//std::cout << "Number of trifeature arcs: " << trifeatureLoop->arcs().size() << std::endl;
	//for (auto arc : trifeatureLoop->arcs()) {
	//	std::cout << " Arc contains " << arc->halfedges().size() << " half-edges." << std::endl;
	//}

	////创建一个bezier样条曲面模型
	//BE_curve_model model;
	//model.quadFeatureLoops() = quadfeatureLoop;
	//model.triFeatureLoops() = trifeatureLoop;
	////可视化三角形网格和特征边
	//viewMesh_featureLoop(&triMesh, *trifeatureLoop);
	////匹配三角网格特征线和四边形网格特征线
	//model.mapTriLoopToQuadLoop(&quadMesh, &triMesh);
	////可视化特征线匹配关系
	////viewMesh_mapFeatureLoop(&quadMesh, *quadfeatureLoop, *trifeatureLoop, model);

	////匹配三角形网格点到四边形特征环中半边
	//model.triPointToQuadHalfedge(&quadMesh, &triMesh);
	////可视化三角形网格点到四边形特征环中半边
	//std::vector<QB_curve_sampling*> qbm_samplings;
	//qbm_samplings = model.samplings();
	//std::map<V*, H*> pointToHalfedgeMap;
	//pointToHalfedgeMap = model.getTriPointToQuadHalfedgeMap();
	//std::map<H*, int> quadHalfedgeInLoopID;
	//quadHalfedgeInLoopID = model.getQuadHalfedgeToIdMap();
	//std::cout << "qbm_samplings.size()_before" << qbm_samplings.size() << std::endl;

	//// 在调用函数时
	////std::map<H*, V*> quadToTriMap;  // 创建映射容器
	////model.addInterpolatedPointsForUnmappedHalfedges(&quadMesh, quadToTriMap);
	////viewMesh(&triMesh, &quadMesh, quadToTriMap);

	//// 添加新功能：为未映射的半边添加插值点
	//model.addInterpolatedPointsForUnmappedHalfedges(&quadMesh);
	//model.addInterpolatedPointsForFeatureArcs1(&quadMesh);

	//qbm_samplings = model.samplings();
	//pointToHalfedgeMap = model.getTriPointToQuadHalfedgeMap();
	//quadHalfedgeInLoopID = model.getQuadHalfedgeToIdMap();
	//std::cout << "qbm_samplings.size()_after" << qbm_samplings.size() << std::endl;
	////viewMesh_triPointToQuadHalfedge(&quadMesh, qbm_samplings, pointToHalfedgeMap, quadHalfedgeInLoopID, *quadfeatureLoop);

	////建立特征线的初始bezier控制点(尖锐点保持）
	//model.setControlPoints(&quadMesh);
	//model.markFeaturePts_input(&quadMesh);
	//std::vector<BE_bezierCurve_approxi*> be_bcs;
	//be_bcs = model.bcs();
	//std::cout << "be_bcs.size():" << be_bcs.size() << std::endl;
	////可视化四边形网格特征边的bezier控制点
	////viewMesh_featureLoopCPts(&quadMesh, *quadfeatureLoop, be_bcs);

	////重新确定四边形网格点
	//model.initializeQuadVertexIdMap(&quadMesh);	//特征线上四边形网格点重新排序(四边形网格点到特征环中点的新id映射）
	//model.parametrize2(&quadMesh);//数据点参数化

	////viewFeatureBezierCurve(be_bcs, &triMesh);

	//// 创建存储排序距离的向量
	//std::vector<std::pair<int, double>> sortedCurveDistances;
	////计算三角形网格边到样条曲线的hausdroff距离
	//double hausdorffCurveDistBefore = model.computeCurveHausdorffDistance(&triMesh, 0.05, sortedCurveDistances);
	//std::cout << "hausdorffCurveDistBefore: " << hausdorffCurveDistBefore << endl;

	////model.compute_quad_vertex(&quadMesh);//更新四边形网格点
	////model.compute_quad_vertex1(&quadMesh);//更新四边形网格点,天然正则化
	////model.compute_quad_vertex2(&quadMesh);//更新四边形网格点，加速求解eigen库
	//model.compute_quad_vertex7(&quadMesh);//最初版本
	//model.preserveQuadCornerPoints(&quadMesh);//保持四边形网格尖锐点
	//recoverFlatBoundaryVertexPoint(&quadMesh);//cfw add 202504040124 将边界平坦点保持初始坐标
	//model.CheckDegeneratedBoudnaryVertexAndUpdate(&quadMesh);//cfw add 202507112138

	////viewMesh(&quadMesh);
	////可视化新四边形网格
	//qbm_samplings = model.samplings();
	//pointToHalfedgeMap = model.getTriPointToQuadHalfedgeMap();
	//quadHalfedgeInLoopID = model.getQuadHalfedgeToIdMap();
	////viewMesh_triPointToQuadHalfedge(&quadMesh, qbm_samplings, pointToHalfedgeMap, quadHalfedgeInLoopID, *quadfeatureLoop);

	////可视化三角形网格和对应的bezier样条曲线
	//model.quadFeatureLoops() = quadfeatureLoop;
	////更新bezier控制点
	//model.resetControlPoints(&quadMesh);

	//be_bcs = model.bcs();
	////可视化四边形网格特征边的bezier控制点
	////viewMesh_featureLoopCPts(&quadMesh, *quadfeatureLoop, be_bcs);
	//std::cout << "be_bcs.size():" << be_bcs.size() << std::endl;
	////model.computeTotalError();//计算采样点拟合误差

	////计算三角形网格边到样条曲线的hausdroff距离
	//double hausdorffCurveDist = model.computeCurveHausdorffDistance(&triMesh, 0.05, sortedCurveDistances);
	//std::cout << "hausdorffCurveDistAfter: " << hausdorffCurveDist << endl;

	////model.computeTriTotalError();//计算三角形数据点拟合误差

	////viewFeatureBezierCurve(be_bcs, &triMesh);

	////viewMesh(&quadMesh);
	////viewFeatureBezierCurve(be_bcs, &quadMesh2);



	////清理内存
	////delete quadfeatureLoop;
	////delete trifeatureLoop;

	/*-------------------------------------------特征曲线约束结束------------------------------------------------------------------*/
	//markFeatureEdgesBySharpEdges(&quadMesh);
	////viewMesh(&quadMesh, argc, argv);

	//markTLayout_quadMesh(&quadMesh);
	//print1 = "# tLayout on quadMesh: " + std::to_string(statisticTLayoutPatchNum(&quadMesh)) + "\n";
	//myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	//statisticTLayoutPatch_quadNum(&quadMesh);
	////viewMesh(&quadMesh);

	//std::string filename = prefix_FilePath_QuadMesh;
	///*--------------------------------------------*/
	//if (ifApproximation)
	//{
	//	//preprocessing quad mesh
	//	/*------------------calculate some quality on quad mesh start-----------------------------------------------------*/
	//	calculateScaledJacobian_quadMesh(&quadMesh);
	//	calculateEdgeRatio_quadMesh(&quadMesh);
	//	for (auto f : It::MFIterator(&quadMesh))
	//	{
	//		double sj = f->quality_scaledJocbian();
	//		double er = f->quality_edgeRatio();
	//		if (sj < 0.5 || er > 10)
	//		{
	//			for (auto fv : It::FVIterator(&quadMesh, f))
	//			{
	//				fv->fixed() = true;
	//				//std::cout << "bad face vertex id: " << fv->id() << std::endl;
	//			}
	//		}
	//		//std::cout << "face id: " << f->id() << " scaledJacobian: " << sj << " edgeRatio: " << er << std::endl;
	//	}
	//	/*------------------calculate some quality on quad mesh end-----------------------------------------------------*/
	//	BE_surface_model model_surf_approx;
	//	clock_t start4, finish4;
	//	double totaltime4;
	//	start4 = clock();
	//	model_surf_approx.obtainCtlPts(&quadMesh);
	//	compute_mesh_normal(&quadMesh);
	//	//model_surf.buildBezierSurfs(&quadMesh);
	//	model_surf_approx.buildBezierSurfs_cw(&quadMesh);
	//	remarkSingularPts(&quadMesh);
	//	//model.raiseBezierSurfDegree3To5(&quadMesh);
	//	model_surf_approx.raiseBezierSurfDegree3To5_all(&quadMesh);
	//	model_surf_approx.sortBezierCtlPts_cw(&quadMesh);
	//	computeV_degree_mesh(&quadMesh);
	//	model_surf_approx.G1Constraints_optimization_cws(&quadMesh);
	//	finish4 = clock();
	//	totaltime4 = (double)(finish4 - start4) / CLOCKS_PER_SEC;
	//	//std::cout << "Obtain G1 NURBS is finished! It takes " << totaltime4 << "s！" << std::endl;
	//	model_surf_approx.computeBezierSurfacesKnots(&quadMesh);
	//	model_surf_approx.sortBezierCtlPts(&quadMesh);
	//	int if_fixedBoundary = 1;
	//	/*approximation*/
	//	clock_t start3, finish3;
	//	double totaltime3;
	//	start3 = clock();
	//	/*读取初始三角网格 */
	//	//M triMesh;
	//	//triMesh.read_m(prefix_FilePath_TriangleMesh.c_str());
	//	//std::cout << "# V on triMesh: " << triMesh.numVertices() << std::endl;
	//	//std::cout << "# boundary V on triMesh: " << statisticBoundaryVertexNum(&triMesh) << std::endl;
	//	/*对网格中的参数uv进行一定数量级的放大*/
	//	magnifyMeshUV(&quadMesh, 1e5);
	//	magnifyMeshUV(&triMesh, 1e5);
	//	/*通过参数化结果在四边形网格上定位出三角形网格点，获取对应的采样点信息*/
	//	std::vector<QB_sampling> sts;
	//	obtainSamplingFromTriMeshByParameterization(&triMesh, &quadMesh, sts);
	//	print1 = "# samplingPts from tri: " + std::to_string(sts.size()) + "\n";
	//	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);

	//	if (triMesh.patchMap_mesh().size() == 0)
	//	{
	//		triMesh.classifyFacesByPatch();
	//	}
	//	ensureQuadFaceSampling(&triMesh, &quadMesh, 20, sts);  // 确保每一片区域采样点足够

	//	/*部分四边形网格点作为采样点*/
	//	//bool part_sampling_quad = 0;
	//	//obtainSampling_quad(&quadMesh, sts, part_sampling_quad);

	//	print1 = "# samplingPts: " + std::to_string(sts.size()) + "\n";
	//	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	//	//std::cout << "# boundary sampling: " << statisticBoundarySamplingNum(sts) << std::endl;
	//	//viewMesh(&triMesh, &quadMesh, argc, argv);
	//	//viewMesh(&triMesh,argc,argv);
	//	//viewSampling(&quadMesh, sts);
	//	for (int i = 0; i < sts.size(); i++)
	//	{
	//		model_surf_approx.samplings().push_back(std::make_shared<QB_sampling>(sts[i]));
	//	}
	//	model_surf_approx.computeSamplingPtUV(&quadMesh);
	//	//statisticSampleNum_model(&quadMesh, model.samplings());
	//	finish3 = clock();
	//	totaltime3 = (double)(finish3 - start3) / CLOCKS_PER_SEC;
	//	std::cout << "Computing sampling UV is finished! It takes " << totaltime3 << "s！" << std::endl;
	//	clock_t start5, finish5;
	//	double totaltime5;
	//	start5 = clock();
	//	//model.approximateSurface();
	//	//model.approximateSurface_optimize3(&quadMesh);
	//	/*读取用户输入拟合精度*/
	//	double approxi_thread = _approxi_thread;
	//	/*读取用户输入光顺权重*/
	//	double approxi_smooth_weight = _approxi_smooth_weight;
	//	if (if_fixedBoundary == 1)
	//	{
	//		model_surf_approx.approximateSurface_optimize_check(&quadMesh, approxi_thread, approxi_smooth_weight);
	//	}
	//	else
	//	{
	//		model_surf_approx.approximateSurface_optimize_check2(&quadMesh, approxi_thread, approxi_smooth_weight);
	//	}
	//	//model.approximateSurface_optimize_check(&quadMesh,approxi_thread,approxi_smooth_weight);
	//	//model.approximateSurface_optimize_check2(&quadMesh, approxi_thread, approxi_smooth_weight);
	//	//model.computeSamplingPos_approx();
	//	//model.approximateSurface_boundary(&quadMesh, approxi_thread, approxi_smooth_weight);
	//	finish5 = clock();
	//	totaltime5 = (double)(finish5 - start5) / CLOCKS_PER_SEC;
	//	print1 = "Approximation is finished! It takes " + std::to_string(totaltime5) + "s！" + "\n";
	//	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	//	updateMeshPosByControlMesh(&quadMesh, model_surf_approx);
	//	std::string outputApproxQuadMeshFileName = filename.substr(0, filename.length() - 2) + "_approxQuadMesh.m";
	//	quadMesh.write_m(outputApproxQuadMeshFileName.c_str());
	//}
	/*--------------------------------------------*/

	//创建一个bezier样条曲面模型
	BE_surface_model model_surf;
	clock_t start4, finish4;
	double totaltime4;
	start4 = clock();
	model_surf.obtainCtlPts(&quadMesh);
	compute_mesh_normal(&quadMesh);
	//model_surf.buildBezierSurfs(&quadMesh);
	model_surf.buildBezierSurfs_cw(&quadMesh);
	remarkSingularPts(&quadMesh);
	//model.raiseBezierSurfDegree3To5(&quadMesh);
	model_surf.raiseBezierSurfDegree3To5_all(&quadMesh);
	std::string filename = prefix_FilePath_QuadMesh;
	std::string outputBSplineFileName = filename.substr(0, filename.length() - 2) + "_G0_bspline_tri_q.igs";
	model_surf.sortBezierCtlPts_cw(&quadMesh);
	computeV_degree_mesh(&quadMesh);
	model_surf.G1Constraints_optimization_check(&quadMesh);
	outputBSplineFileName = filename.substr(0, filename.length() - 2) + "_G1_CCG_bspline_tri_q.igs";
	finish4 = clock();
	totaltime4 = (double)(finish4 - start4) / CLOCKS_PER_SEC;
	//std::cout << "Obtain G1 NURBS is finished! It takes " << totaltime4 << "s！" << std::endl;
	model_surf.computeBezierSurfacesKnots(&quadMesh);
	model_surf.sortBezierCtlPts(&quadMesh);
	model_surf.obtainBezierCptsBy_cpts(1.0);

	//obtainQuadLayout_extend2(&quadMesh);
	obtainTLayout_motorcycleGraph_featureAware(&quadMesh);
	//obtainLayout_consistencyQuadMeshFaceNum(&quadMesh);
	//viewMesh(&quadMesh);
	//model.obtainBSplineSurfacePatches_general(&quadMesh);
	//model.obtainBSplineSurfacePatches_general_UVConsistency(&quadMesh);
	//model_surf.obtainBSplineSurfacePatches_general_UVConsistency_recordQuadBoundaryNURBSBoundaryCurveRelation(&quadMesh);
	model_surf.obtainBSplineSurfacePatches_general_UVConsistency_topology_boundary(&quadMesh);
	model_surf.removeMultiKnots_BSplines();
	model_surf.computeBSplineSurfaceUniqueKnotMultiNum();

	//BPatchesTLayoutIndex(model_surf, &quadMesh);

	//model_surf.outputBSpline_iges(outputBSplineFileName.c_str());
	model_surf.outputBSpline_iges_5_3(outputBSplineFileName.c_str());
	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	print1 = "Reconstructing NURBS is finished! It takes " + std::to_string(totaltime) + "s！" + "\n";
	myPrint.PrintLogInformation(print1, printToscreen_Log, filePointer_Log);
	//model.outputBezierSurfaces_iges(&quadMesh, outputBezierFileName.c_str(), 3, 3);

	/*-------------------------------------样条曲面质量分析开始------------------------------------------start*/
	//model_surf.create_bsurBoundarys_general(&quadMesh);

	////model.createTopology_bspline(&quadMesh);
	//std::string outputBsplineTopologyTMeshFileName = filename.substr(0, filename.length() - 2) + "_BsplineTopologyTMesh.obj";
	//model_surf.outputTopology_bspline_TMesh_obj(&quadMesh, outputBsplineTopologyTMeshFileName.c_str());
	//model_surf.createTopology_bspline_TMesh(&quadMesh);
	////std::cout << "111" << std::endl;
	//std::string outputBsplineTopologyQuadMeshFileName = filename.substr(0, filename.length() - 2) + "_BsplineTopologyQuad.m";
	///*for (auto v : It::MVIterator(&model_surf.bsplineSurfacesTopology_quad()))
	//{
	//	v->point() = v->point() / scale;
	//}*/
	//model_surf.bsplineSurfacesTopology_quad().write_m(outputBsplineTopologyQuadMeshFileName.c_str());

	//std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>> allSurfs = model_surf.bsplineSurfs();
	//// ===================== 控制点一致性检查 ======================
	////std::vector<ControlPointData> unequalControlPoints;
	////checkContronPoints(model_surf, unequalControlPoints);
	////std::cout << "-----unequalControlPoints: " << unequalControlPoints.size() << std::endl;

	////magnifyMeshPoint(&quadMesh, 1.0 / scale);
	////viewControlPoints(&quadMesh, unequalControlPoints);
	//// ============================================================

	//// ===================== 边界采样点连续性分析 ==================
	//double stepL = 0.5;
	//obtainSamplingPoints_bsurfBoundarys_parrallel(model_surf, stepL);
	//// ===================== G0连续性分析 ==========================
	//std::unordered_map<NUBE_bsurfBoudnary*, double> maxPosDiff;
	//findMaxPositionDiffForEachEdge_T(model_surf, maxPosDiff);

	//// ===================== G1连续性分析 ==========================
	//std::unordered_map<NUBE_bsurfBoudnary*, double> maxAngles;
	//findMaxAngleForEachEdge_T(model_surf, maxAngles);

	//// ===================== G2连续性分析 ==========================
	//std::unordered_map<NUBE_bsurfBoudnary*, double> maxDuuDiff, maxDuvDiff, maxDvvDiff;
	//findMaxSecondDerivativeDiffForEachEdge_T(model_surf, maxDuuDiff, maxDuvDiff, maxDvvDiff);
	//// ============================================================

	//// 1、计算并打印所有曲面片的边界边长比
	//computeBoundaryEdgeRatio(allSurfs);

	//// 2、计算控制点之间的最小距离 s
	//computeControlPointMinDistance(allSurfs);

	//// 4. 调用对每条 NURBS 边界计算 Hausdorff 距离
	////computeHausdorffPerBoundary(&triMesh, allSurfs, sampleStep);

	//// 创建存储排序距离的向量
	//std::vector<std::pair<int, double>> sortedDistances;

	//// 计算Hausdorff距离
	//double hausdorffDist = model_surf.computeHausdorffDistance_allNURBSBoundaryCurve3(&triMesh, 3, sortedDistances);

	//// 输出排序后的距离到文件
	//std::string outputDistanceFile = filename.substr(0, filename.length() - 2) + "_hausdorff_distances.txt";
	//model_surf.outputSortedDistances(sortedDistances, outputDistanceFile);

	//model_surf.computeGlobalSurfaceMetrics();

	//// =================== 分析数据输出到.txt ======================
	//std::string maxPosDiffFileName = filename.substr(0, filename.length() - 2) + "_maxPosDiff.txt";
	//std::string maxAnglesFileName = filename.substr(0, filename.length() - 2) + "_maxAngles.txt";
	////std::string maxDuuDiffFileName = filename.substr(0, filename.length() - 2) + "_maxDuuDiff.txt";
	////std::string maxDuvDiffFileName = filename.substr(0, filename.length() - 2) + "_maxDuvDiff.txt";
	////std::string maxDvvDiffFileName = filename.substr(0, filename.length() - 2) + "_maxDvvDiff.txt";
	//dumpSorted_T(maxPosDiff, maxPosDiffFileName.c_str());
	//dumpSorted_T(maxAngles, maxAnglesFileName.c_str());
	////dumpSorted_T(maxDuuDiff, maxDuuDiffFileName.c_str());
	////dumpSorted_T(maxDuvDiff, maxDuvDiffFileName.c_str());
	////dumpSorted_T(maxDvvDiff, maxDvvDiffFileName.c_str());
	//std::string outputBsurfMetricsName = filename.substr(0, filename.length() - 2) + "_report.txt";
	//model_surf.outputMetricsToFile(outputBsurfMetricsName.c_str());
	/*-------------------------------------样条曲面质量分析结束------------------------------------------start*/

}

void test_find(int argc, char* argv[]) {
	clock_t start, finish;
	double totaltime;
	start = clock();
	/*读取四边形网格 */
	M quadMesh;
	/*提示用户输入四边形网格数据*/
	if (argc < 3)
	{
		std::cout << " sorry, input data missing!" << std::endl;
	}
	quadMesh.read_m(argv[2]);
	std::cout << "# V on QuadMesh: " << quadMesh.numVertices() << std::endl;
	std::cout << "# F on QuadMesh: " << quadMesh.numFaces() << std::endl;
	//std::cout << "# boundary V on quadMesh: " << statisticBoundaryVertexNum(&quadMesh) << std::endl;

	/*读取初始三角网格 */
	M triMesh;
	triMesh.read_m(argv[1]);
	std::cout << "# V on triMesh: " << triMesh.numVertices() << std::endl;
	//std::cout << "# boundary V on triMesh: " << statisticBoundaryVertexNum(&triMesh) << std::endl;

	/*---------------FindQuadsInsideTris start---------------------------------------*/
	//auto hit = FindQuadsInsideTris_Brute_WithBuckets(&quadMesh, &triMesh, 1e-12);

	//std::cout << "[Brute] hit quad faces = " << hit.size() << "\n";
	////std::cout << "  first few face ids: ";
	//for (size_t i = 0; i < hit.size(); ++i) {
	//	hit[i]->tempMark() = true;
	//	//std::cout << hit[i]->id() << " ";
	//}
	/*---------------FindQuadsInsideTris end---------------------------------------*/

	/*------------------calculate some quality on quad mesh start-----------------------------------------------------*/
	calculateScaledJacobian_quadMesh(&quadMesh);
	calculateEdgeRatio_quadMesh(&quadMesh);
	for(auto f : It::MFIterator(&quadMesh))
	{
		double sj = f->quality_scaledJocbian();
		double er = f->quality_edgeRatio();
		if (sj < 0.5 || er > 10)
		{
			f->tempMark() = true;
		}
		//std::cout << "face id: " << f->id() << " scaledJacobian: " << sj << " edgeRatio: " << er << std::endl;
	}
	/*------------------calculate some quality on quad mesh end-----------------------------------------------------*/

	//std::cout << std::endl;
	//viewMesh(&triMesh);
	//viewMesh(&quadMesh, &quadMesh, argc, argv);
	//viewScalarOnMeshFace(&quadMesh, argc, argv);
}

int main(int argc, char* argv[])
{
	float memuse_4 = _Ask_OS_memory_usage();
	std::cout << "memuse_4:" << memuse_4 << std::endl;
#ifdef MEMORY_LEAK_DETECTION
	/* Memory leak detection*/
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	test250407_594_3_3_12512192258(argc, argv);
	//test250407_594_3_3(argc, argv);
	//test_find(argc, argv);
	float memuse_5 = _Ask_OS_memory_usage();
	std::cout << "memuse_5:" << memuse_5 << std::endl;
	return 1;
}
