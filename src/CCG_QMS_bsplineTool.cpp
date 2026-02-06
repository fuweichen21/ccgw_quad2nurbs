#include"CCG_QMS_bsplineTool.h"
#include<omp.h>

void CCG_QMSLib::BPatchesTLayoutIndex(CCG_QMS_model& model, M* quadMesh)
{
	for (auto bp : model.bsplineSurfs())
	{
		F* bpf = quadMesh->idFace(bp->quadMeshFaceId());
		bp->quadLayoutIndex() = bpf->subPatchIndex();
	}

	/*check*/
	/*for (auto bp : model.bsplineSurfs())
	{
		std::cout <<"bp->quadLayoutIndex(): " << bp->quadLayoutIndex() << std::endl;
	}*/
}

void CCG_QMSLib::outputRelationBSpline_TLayout_fixedFormat(const char* output, CCG_QMS_model& model, M* quadMesh)
{
	std::fstream _os(output, std::fstream::out);
	if (_os.fail())
	{
		fprintf(stderr, "Error is opening file %s\n", output);
		return;
	}

	/*Output of basic information*/
	_os << "Total Number of TLayout Regions: " << quadMesh->tLayout_num() << std::endl;
	_os << "Total Number of Quad Cells: " << quadMesh->numFaces() << std::endl;
	_os << "Total Number of NURBS Patches:  " << model.bsplineSurfs().size() << std::endl;
	int totalCtlPtsNum = 0;
	for (auto bs : model.bsplineSurfs())
	{
		totalCtlPtsNum += (bs->cpts().size() * bs->cpts()[0].size());
	}
	_os << "Total Number of NURBS Control Points: " << totalCtlPtsNum << std::endl;

	/*output the relation between TLayout and  NURBS surface*/
	/*build relation*/
	/*
	* Create a two-dimensional array.
	* The first dimension represents the index of TLayout, 
	* and the second dimension represents the id of the B-spline surface patches contained on each TLayout patch.
	*/
	std::vector<std::vector<int>> relationTLayoutBSplines(quadMesh->tLayout_num());
	for (int i = 0; i < model.bsplineSurfs().size(); i++)
	{
		relationTLayoutBSplines[model.bsplineSurfs()[i]->quadLayoutIndex()].push_back(i);
	}
	/*output relationship*/
	for (int i = 0; i < quadMesh->tLayout_num(); i++)
	{
		_os << "TLayout Region " << i + 1 << std::endl;
		_os << "    Quad Cell " << quadMesh->tLayoutPatches_quadNum()[i] << std::endl;
		for (int j = 0; j < relationTLayoutBSplines[i].size(); j++)
		{
			_os << "        NURBS Patch " << j + 1 << ": " << model.bsplineSurfs()[relationTLayoutBSplines[i][j]]->cpts()[0].size() * model.bsplineSurfs()[relationTLayoutBSplines[i][j]]->cpts().size() <<
				"  ( " << (int)model.bsplineSurfs()[relationTLayoutBSplines[i][j]]->knotsU().back() << " X " << (int)model.bsplineSurfs()[relationTLayoutBSplines[i][j]]->knotsV().back() << " grid )" << std::endl;
		}
	}
	_os.close();
}

int CCG_QMSLib::FindSpan1(int degree, const std::vector<double>& knots, double u)
{
	int n = knots.size() - degree - 1;
	if (u >= knots[n]) return n - 1;
	if (u <= knots[degree]) return degree;
	int low = degree;
	int high = n;
	int mid = (low + high) / 2;
	while (u < knots[mid] || u >= knots[mid + 1]) {
		if (u < knots[mid]) high = mid;
		else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

//Attribute to quantify the surface quality
/*-----------------------------------------------------------------------------------------------*/
void CCG_QMSLib::obtainSamplingPoints_bsurfBoundarys_parrallel(CCG_QMS_model& model, double samplingStepLen)
{
	// First parallel section: sampling point processing
#pragma omp parallel for
	for (auto bsurf : model.bsplineSurfs())
	{
		for (auto bsurfBoundary : bsurf->boundarys())
		{
			for (auto bsurfBoundarySeg : bsurfBoundary->segements())
			{
				//sampling point parameterization & position
				std::vector<std::vector<CPoint>> cp = bsurf->cpts();
				int degree_u = bsurf->degree()[0]; // the degree of U
				int degree_v = bsurf->degree()[1]; // the degree of V
				std::vector<double> kontsu = bsurf->knotsU();//the knots of U
				std::vector<double> kontsv = bsurf->knotsV();//the knots of V
				double u_max = bsurf->knotsU().back();
				double v_max = bsurf->knotsV().back();
				obtainSamplingPoints_bsurfBoundarySegments(bsurfBoundary, bsurfBoundarySeg, cp, kontsu, kontsv, degree_u, degree_v, samplingStepLen, u_max, v_max);
				computingSamplingInformation_bsurfBoundarySegment(bsurfBoundarySeg, bsurf);
			}
		}
	}

	// Second parallel section: computing sampling points information
#pragma omp parallel for
	for (auto bsurf : model.bsplineSurfs())
	{
		for (auto bsurfBoundary : bsurf->boundarys())
		{
			for (auto bsurfBoundarySeg : bsurfBoundary->segements())
			{
				//G0、G1、duu、duv、dvv...
				obtainingDualSamplingPoint_samplingPt_segment(bsurfBoundarySeg);
			}
		}
	}

	// Third section: computing boundary information (needs reduction for max values)
	QB_sampling* qbs_diffG0_sampling_model = NULL;
	QB_sampling* qbs_diffG1_sampling_model = NULL;
	QB_sampling* qbs_diff_duu_sampling_model = NULL;
	QB_sampling* qbs_diff_duv_sampling_model = NULL;
	QB_sampling* qbs_diff_dvv_sampling_model = NULL;

	// We'll use parallel reduction for finding maximum values
#pragma omp parallel
	{
		QB_sampling* local_diffG0 = NULL;
		QB_sampling* local_diffG1 = NULL;
		QB_sampling* local_diff_duu = NULL;
		QB_sampling* local_diff_duv = NULL;
		QB_sampling* local_diff_dvv = NULL;

#pragma omp for nowait
		for (auto bsurf : model.bsplineSurfs())
		{
			QB_sampling* qbs_diffG0_sampling_bsurf = NULL;
			QB_sampling* qbs_diffG1_sampling_bsurf = NULL;
			QB_sampling* qbs_diff_duu_sampling_bsurf = NULL;
			QB_sampling* qbs_diff_duv_sampling_bsurf = NULL;
			QB_sampling* qbs_diff_dvv_sampling_bsurf = NULL;

			bool tempMark = false;//recoard bsurf have one non-boudnry curve at least

			for (auto bsurfBoundary : bsurf->boundarys())
			{
				QB_sampling* qbs_diffG0_sampling_boundary = NULL;
				QB_sampling* qbs_diffG1_sampling_boundary = NULL;
				QB_sampling* qbs_diff_duu_sampling_boundary = NULL;
				QB_sampling* qbs_diff_duv_sampling_boundary = NULL;
				QB_sampling* qbs_diff_dvv_sampling_boundary = NULL;

				for (auto bsurfBoundarySeg : bsurfBoundary->segements())
				{
					if (bsurfBoundarySeg->symBsurfSeg() != NULL)
					{
						tempMark = true;
						for (auto sp : bsurfBoundarySeg->samplings())
						{
							if (qbs_diffG0_sampling_boundary == NULL)
							{
								qbs_diffG0_sampling_boundary = sp;
							}
							else
							{
								qbs_diffG0_sampling_boundary = qbs_diffG0_sampling_boundary->diffG0() < sp->diffG0() ? sp : qbs_diffG0_sampling_boundary;
							}

							if (qbs_diffG1_sampling_boundary == NULL)
							{
								qbs_diffG1_sampling_boundary = sp;
							}
							else
							{
								qbs_diffG1_sampling_boundary = qbs_diffG1_sampling_boundary->diffG1() < sp->diffG1() ? sp : qbs_diffG1_sampling_boundary;
							}

							if (qbs_diff_duu_sampling_boundary == NULL)
							{
								qbs_diff_duu_sampling_boundary = sp;
							}
							else
							{
								qbs_diff_duu_sampling_boundary = qbs_diff_duu_sampling_boundary->diff_duu() < sp->diff_duu() ? sp : qbs_diff_duu_sampling_boundary;
							}

							if (qbs_diff_duv_sampling_boundary == NULL)
							{
								qbs_diff_duv_sampling_boundary = sp;
							}
							else
							{
								qbs_diff_duv_sampling_boundary = qbs_diff_duv_sampling_boundary->diff_duv() < sp->diff_duv() ? sp : qbs_diff_duv_sampling_boundary;
							}

							if (qbs_diff_dvv_sampling_boundary == NULL)
							{
								qbs_diff_dvv_sampling_boundary = sp;
							}
							else
							{
								qbs_diff_dvv_sampling_boundary = qbs_diff_dvv_sampling_boundary->diff_dvv() < sp->diff_dvv() ? sp : qbs_diff_dvv_sampling_boundary;
							}
						}
					}
					else
					{
						bsurfBoundary->boundary() = true;
					}
				}

				if (!bsurfBoundary->boundary())
				{
					bsurfBoundary->diffG0_sampling() = qbs_diffG0_sampling_boundary;
					bsurfBoundary->diffG1_sampling() = qbs_diffG1_sampling_boundary;
					bsurfBoundary->diff_duu_sampling() = qbs_diff_duu_sampling_boundary;
					bsurfBoundary->diff_duv_sampling() = qbs_diff_duv_sampling_boundary;
					bsurfBoundary->diff_dvv_sampling() = qbs_diff_dvv_sampling_boundary;

					if (qbs_diffG0_sampling_bsurf == NULL)
					{
						qbs_diffG0_sampling_bsurf = bsurfBoundary->diffG0_sampling();
					}
					else
					{
						qbs_diffG0_sampling_bsurf = qbs_diffG0_sampling_bsurf->diffG0() < bsurfBoundary->diffG0_sampling()->diffG0() ? bsurfBoundary->diffG0_sampling() : qbs_diffG0_sampling_bsurf;
					}

					if (qbs_diffG1_sampling_bsurf == NULL)
					{
						qbs_diffG1_sampling_bsurf = bsurfBoundary->diffG1_sampling();
					}
					else
					{
						qbs_diffG1_sampling_bsurf = qbs_diffG1_sampling_bsurf->diffG1() < bsurfBoundary->diffG1_sampling()->diffG1() ? bsurfBoundary->diffG1_sampling() : qbs_diffG1_sampling_bsurf;
					}

					if (qbs_diff_duu_sampling_bsurf == NULL)
					{
						qbs_diff_duu_sampling_bsurf = bsurfBoundary->diff_duu_sampling();
					}
					else
					{
						qbs_diff_duu_sampling_bsurf = qbs_diff_duu_sampling_bsurf->diff_duu() < bsurfBoundary->diff_duu_sampling()->diff_duu() ? bsurfBoundary->diff_duu_sampling() : qbs_diff_duu_sampling_bsurf;
					}

					if (qbs_diff_duv_sampling_bsurf == NULL)
					{
						qbs_diff_duv_sampling_bsurf = bsurfBoundary->diff_duv_sampling();
					}
					else
					{
						qbs_diff_duv_sampling_bsurf = qbs_diff_duv_sampling_bsurf->diff_duv() < bsurfBoundary->diff_duv_sampling()->diff_duv() ? bsurfBoundary->diff_duv_sampling() : qbs_diff_duv_sampling_bsurf;
					}

					if (qbs_diff_dvv_sampling_bsurf == NULL)
					{
						qbs_diff_dvv_sampling_bsurf = bsurfBoundary->diff_dvv_sampling();
					}
					else
					{
						qbs_diff_dvv_sampling_bsurf = qbs_diff_dvv_sampling_bsurf->diff_dvv() < bsurfBoundary->diff_dvv_sampling()->diff_dvv() ? bsurfBoundary->diff_dvv_sampling() : qbs_diff_dvv_sampling_bsurf;
					}
				}
			}

			if (tempMark)
			{
				bsurf->diffG0_sampling() = qbs_diffG0_sampling_bsurf;
				bsurf->diffG1_sampling() = qbs_diffG1_sampling_bsurf;
				bsurf->diff_duu_sampling() = qbs_diff_duu_sampling_bsurf;
				bsurf->diff_duv_sampling() = qbs_diff_duv_sampling_bsurf;
				bsurf->diff_dvv_sampling() = qbs_diff_dvv_sampling_bsurf;

				// Update local thread-local maxima
				if (local_diffG0 == NULL || (qbs_diffG0_sampling_bsurf != NULL && local_diffG0->diffG0() < qbs_diffG0_sampling_bsurf->diffG0()))
					local_diffG0 = qbs_diffG0_sampling_bsurf;

				if (local_diffG1 == NULL || (qbs_diffG1_sampling_bsurf != NULL && local_diffG1->diffG1() < qbs_diffG1_sampling_bsurf->diffG1()))
					local_diffG1 = qbs_diffG1_sampling_bsurf;

				if (local_diff_duu == NULL || (qbs_diff_duu_sampling_bsurf != NULL && local_diff_duu->diff_duu() < qbs_diff_duu_sampling_bsurf->diff_duu()))
					local_diff_duu = qbs_diff_duu_sampling_bsurf;

				if (local_diff_duv == NULL || (qbs_diff_duv_sampling_bsurf != NULL && local_diff_duv->diff_duv() < qbs_diff_duv_sampling_bsurf->diff_duv()))
					local_diff_duv = qbs_diff_duv_sampling_bsurf;

				if (local_diff_dvv == NULL || (qbs_diff_dvv_sampling_bsurf != NULL && local_diff_dvv->diff_dvv() < qbs_diff_dvv_sampling_bsurf->diff_dvv()))
					local_diff_dvv = qbs_diff_dvv_sampling_bsurf;
			}
		}

		// Reduction for global maxima
#pragma omp critical
		{
			if (local_diffG0 != NULL && (qbs_diffG0_sampling_model == NULL || qbs_diffG0_sampling_model->diffG0() < local_diffG0->diffG0()))
				qbs_diffG0_sampling_model = local_diffG0;

			if (local_diffG1 != NULL && (qbs_diffG1_sampling_model == NULL || qbs_diffG1_sampling_model->diffG1() < local_diffG1->diffG1()))
				qbs_diffG1_sampling_model = local_diffG1;

			if (local_diff_duu != NULL && (qbs_diff_duu_sampling_model == NULL || qbs_diff_duu_sampling_model->diff_duu() < local_diff_duu->diff_duu()))
				qbs_diff_duu_sampling_model = local_diff_duu;

			if (local_diff_duv != NULL && (qbs_diff_duv_sampling_model == NULL || qbs_diff_duv_sampling_model->diff_duv() < local_diff_duv->diff_duv()))
				qbs_diff_duv_sampling_model = local_diff_duv;

			if (local_diff_dvv != NULL && (qbs_diff_dvv_sampling_model == NULL || qbs_diff_dvv_sampling_model->diff_dvv() < local_diff_dvv->diff_dvv()))
				qbs_diff_dvv_sampling_model = local_diff_dvv;
		}
	}

	model.diffG0_sampling() = qbs_diffG0_sampling_model;
	model.diffG1_sampling() = qbs_diffG1_sampling_model;
	model.diff_duu_sampling() = qbs_diff_duu_sampling_model;
	model.diff_duv_sampling() = qbs_diff_duv_sampling_model;
	model.diff_dvv_sampling() = qbs_diff_dvv_sampling_model;
}

void CCG_QMSLib::obtainSamplingPoints_bsurfBoundarySegments(NUBE_bsurfBoudnary* bBoundaryCurve, NUBE_bsurfSegment* bBoundaryCurveSeg, const std::vector<std::vector<CPoint>>& cp, const std::vector<double>& kontsu, const std::vector<double>& kontsv, int degree_u, int degree_v, double stepL, double u_max, double v_max)
{
	/*std::cout << "bBoundaryCurveSeg->bsurfId(): " << bBoundaryCurveSeg->bsurfId() << std::endl;
	std::cout << "bBoundaryCurveSeg->bsurfBoundaryLoacalId(): " << bBoundaryCurveSeg->bsurfBoundaryLoacalId() << std::endl;
	std::cout << "u_max_v_max: " << (bBoundaryCurve->halfedges().size()) * 1.0 << std::endl;
	std::cout << "u_max: " << u_max << std::endl;
	std::cout << "v_max: " << v_max << std::endl;*/
	double u, v;
	if (bBoundaryCurveSeg->bsurfBoundaryLoacalId() == 1) {
		u = 0.0;
		double v_max = (bBoundaryCurve->halfedges().size() - (bBoundaryCurveSeg->firstHalfedgeIndex() - 1)) * 1.0;
		double v_min = (bBoundaryCurve->halfedges().size() - (bBoundaryCurveSeg->endHalfedgeIndex())) * 1.0;
		for (v = v_max; v >= v_min; v -= stepL) {
			if (v > (bBoundaryCurve->halfedges().size()) * 1.0)
			{
				//std::cout << "boundary type 1 para v: " << v << std::endl;
				v = (bBoundaryCurve->halfedges().size()) * 1.0;
			}
			if (v < 0.0)
			{
				//std::cout << "boundary type 1 para v: " << v << std::endl;
				v = 0.0;
			}
			CPoint p = DeBoorAlgorithmSurface1(cp, kontsu, kontsv, degree_u, degree_v, u, v);
			QB_sampling* sample = new QB_sampling;
			sample->pos() = p;
			sample->uv() = { u, v };
			/*std::cout << "boundary type 1 (u,v) =  (" << u << "," << v << ")" << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[0] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[1] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[2] << std::endl;*/
			bBoundaryCurveSeg->samplings().push_back(sample);
		}
	}
	else if (bBoundaryCurveSeg->bsurfBoundaryLoacalId() == 2) {
		v = 0.0;
		double u_min = (bBoundaryCurveSeg->firstHalfedgeIndex() - 1) * 1.0;
		double u_max = (bBoundaryCurveSeg->endHalfedgeIndex()) * 1.0;
		for (u = u_min; u <= u_max; u += stepL) {
			if (u > (bBoundaryCurve->halfedges().size()) * 1.0)
			{
				//std::cout << "boundary type 2 para u: " << v << std::endl;
				u = (bBoundaryCurve->halfedges().size()) * 1.0;
			}
			if (u < 0.0)
			{
				//std::cout << "boundary type 2 para u: " << v << std::endl;
				u = 0.0;
			}
			CPoint p = DeBoorAlgorithmSurface1(cp, kontsu, kontsv, degree_u, degree_v, u, v);
			QB_sampling* sample = new QB_sampling;
			sample->pos() = p;
			sample->uv() = { u, v };
			/*std::cout << "boundary type 2 (u,v) =  (" << u << "," << v << ")" << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[0] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[1] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[2] << std::endl;*/
			bBoundaryCurveSeg->samplings().push_back(sample);
		}
	}
	else if (bBoundaryCurveSeg->bsurfBoundaryLoacalId() == 3) {
		u = u_max;
		double v_min = (bBoundaryCurveSeg->firstHalfedgeIndex() - 1) * 1.0;
		double v_max = (bBoundaryCurveSeg->endHalfedgeIndex()) * 1.0;
		for (v = v_min; v <= v_max; v += stepL) {
			if (v > (bBoundaryCurve->halfedges().size()) * 1.0)
			{
				//std::cout << "boundary type 3 para v: " << v << std::endl;
				v = (bBoundaryCurve->halfedges().size()) * 1.0;
			}
			if (v < 0.0)
			{
				//std::cout << "boundary type 3 para v: " << v << std::endl;
				v = 0.0;
			}
			CPoint p = DeBoorAlgorithmSurface1(cp, kontsu, kontsv, degree_u, degree_v, u, v);
			QB_sampling* sample = new QB_sampling;
			sample->pos() = p;
			sample->uv() = { u, v };
			/*std::cout << "boundary type 3 (u,v) =  (" << u << "," << v << ")" << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[0] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[1] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[2] << std::endl;*/
			bBoundaryCurveSeg->samplings().push_back(sample);
		}
	}
	else {
		v = v_max;
		double u_max = (bBoundaryCurve->halfedges().size() - (bBoundaryCurveSeg->firstHalfedgeIndex() - 1)) * 1.0;
		double u_min = (bBoundaryCurve->halfedges().size() - (bBoundaryCurveSeg->endHalfedgeIndex())) * 1.0;
		for (u = u_max; u >= u_min; u -= stepL) {
			if (u > (bBoundaryCurve->halfedges().size()) * 1.0)
			{
				//std::cout << "boundary type 4 para u: " << v << std::endl;
				u = (bBoundaryCurve->halfedges().size()) * 1.0;
			}
			if (u < 0.0)
			{
				//std::cout << "boundary type 4 para u: " << v << std::endl;
				u = 0.0;
			}
			CPoint p = DeBoorAlgorithmSurface1(cp, kontsu, kontsv, degree_u, degree_v, u, v);
			QB_sampling* sample = new QB_sampling;
			sample->pos() = p;
			sample->uv() = { u, v };
			/*std::cout << "boundary type 4 (u,v) =  (" << u << "," << v << ")" << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[0] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[1] << std::endl;
			std::cout << std::fixed << std::setprecision(19) << sample->pos()[2] << std::endl;*/
			bBoundaryCurveSeg->samplings().push_back(sample);
		}
	}
}

void CCG_QMSLib::computingSamplingInformation_bsurfBoundarySegment(NUBE_bsurfSegment* bBoundaryCurveSeg, std::shared_ptr<CCG_QMS_bsplineSurf> bsurf)
{
	for (auto sp : bBoundaryCurveSeg->samplings())
	{
		std::vector<std::vector<CPoint>> cpts = bsurf->cpts();   //control points	
		int p = bsurf->degree()[0];                //the degree of U
		int q = bsurf->degree()[1];                //the degree of V

		int n = p + 1;  //the number of control points of U
		int m = q + 1;  //the number of control points of V

		int d = 2;   //the highest partial derivative, setting to 2
		std::vector<std::vector<CPoint>> S(d + 1, std::vector<CPoint>(d + 1));  //saving partical derivative

		getMaincurvature_normal(sp->kappa1(), sp->kappa2(), sp->t1(), sp->t2(), sp->normal(), n, p, bsurf->knotsU(), m, q, bsurf->knotsV(), cpts, sp->uv()[0], sp->uv()[1], d, S);
		/*std::cout << "Sampling Point Face ID: " << sp1.fId() << std::endl;
		std::cout << "Sampling Point ID: " << sp1.vId() << std::endl;
		std::cout << "Principal Curvature k1: " << sp1.kappa1() << std::endl;
		std::cout << "Principal Curvature k2: " << sp1.kappa2() << std::endl;
		std::cout << "Principal Direction m1: (" << sp1.t1()[0] << ", " << sp1.t1()[1] << ", " << sp1.t1()[2] << ")" << std::endl;
		std::cout << "Principal Direction m2: (" << sp1.t2()[0] << ", " << sp1.t2()[1] << ", " << sp1.t2()[2] << ")" << std::endl;
		std::cout << "Normal Vector: (" << sp1.normal()[0] << ", " << sp1.normal()[1] << ", " << sp1.normal()[2] << ")" << std::endl;*/
		sp->du() = S[1][0];
		sp->dv() = S[0][1];
		sp->duu() = S[2][0];
		sp->duv() = S[1][1];
		sp->dvv() = S[0][2];
	}
}

void CCG_QMSLib::obtainingDualSamplingPoint_samplingPt_segment(NUBE_bsurfSegment* bBoundaryCurveSeg)
{
	for (int i = 0; i < bBoundaryCurveSeg->samplings().size(); i++)
	{
		if (bBoundaryCurveSeg->symBsurfSeg() != NULL)
		{
			bBoundaryCurveSeg->samplings()[i]->dualSampling() = bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i];

			/*std::cout << "faceId: " << bBoundaryCurveSeg->bsurfId() << " " << bBoundaryCurveSeg->symBsurfSeg()->bsurfId()
				<< "  localId: " << bBoundaryCurveSeg->bsurfBoundaryLoacalId() << " " << bBoundaryCurveSeg->symBsurfSeg()->bsurfBoundaryLoacalId()
				<< "  (" << bBoundaryCurveSeg->samplings()[i]->uv()[0] << "," << bBoundaryCurveSeg->samplings()[i]->uv()[1] << "),  (" << bBoundaryCurveSeg->samplings()[i]->pos()[0] << "," << bBoundaryCurveSeg->samplings()[i]->pos()[1] << "," << bBoundaryCurveSeg->samplings()[i]->pos()[2] << "),  (" << bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i]->uv()[0] << "," << bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i]->uv()[1] << "), ("
				<< bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i]->pos()[0] << "," << bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i]->pos()[1] << "," << bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i]->pos()[2] << ")" << "  sp->diffG0(): " << bBoundaryCurveSeg->samplings()[i]->diffG0() << " diff: " << (bBoundaryCurveSeg->samplings()[i]->pos() - bBoundaryCurveSeg->symBsurfSeg()->samplings()[bBoundaryCurveSeg->samplings().size() - 1 - i]->pos()).norm() << std::endl;*/

				//G0
				//double _diffG0 = (bBoundaryCurveSeg->samplings()[i]->pos() - bBoundaryCurveSeg->samplings()[i]->dualSampling()->pos()).norm();
				//bBoundaryCurveSeg->samplings()[i]->diffG0() = _diffG0;
				//std::cout << " _diffG0: " << _diffG0 << std::endl;
				//std::cout << " diffG0(): " << bBoundaryCurveSeg->samplings()[i]->diffG0() << std::endl;
			bBoundaryCurveSeg->samplings()[i]->diffG0() = (bBoundaryCurveSeg->samplings()[i]->pos() - bBoundaryCurveSeg->samplings()[i]->dualSampling()->pos()).norm();
			//G1
			double dotProduct = bBoundaryCurveSeg->samplings()[i]->normal() * bBoundaryCurveSeg->samplings()[i]->dualSampling()->normal();
			double magnitude1 = bBoundaryCurveSeg->samplings()[i]->normal().norm();
			double magnitude2 = bBoundaryCurveSeg->samplings()[i]->dualSampling()->normal().norm();
			// Check magnitude of normal is valid (is not 0)
			if (magnitude1 == 0.0 || magnitude2 == 0.0) {
				bBoundaryCurveSeg->samplings()[i]->diffG1() = 0.0;
				std::cout << "normal is 0!!!" << std::endl;
				continue;
			}
			// cos value of 2 normal
			double cosTheta = dotProduct / (magnitude1 * magnitude2);
			cosTheta = std::max(-1.0, std::min(1.0, cosTheta));  // verify cosθ ∈ [-1, 1]
			//radian value
			double angle = std::acos(cosTheta);
			bBoundaryCurveSeg->samplings()[i]->diffG1() = angle;
			//std::cout << "Angle between normals: " << angle << " radians" << std::endl;

			//duu
			bBoundaryCurveSeg->samplings()[i]->diff_duu() = (bBoundaryCurveSeg->samplings()[i]->duu() - bBoundaryCurveSeg->samplings()[i]->dualSampling()->duu()).norm();
			//duv
			bBoundaryCurveSeg->samplings()[i]->diff_duv() = (bBoundaryCurveSeg->samplings()[i]->duv() - bBoundaryCurveSeg->samplings()[i]->dualSampling()->duv()).norm();
			//dvv
			bBoundaryCurveSeg->samplings()[i]->diff_dvv() = (bBoundaryCurveSeg->samplings()[i]->dvv() - bBoundaryCurveSeg->samplings()[i]->dualSampling()->dvv()).norm();
		}

	}
}

int CCG_QMSLib::FindSpan(int degree, const std::vector<double>& knots, double u)
{
	int n = knots.size() - degree - 1;
	if (u >= knots[n]) return n - 1;
	if (u <= knots[degree]) return degree;
	int low = degree;
	int high = n;
	int mid = (low + high) / 2;
	while (u < knots[mid] || u >= knots[mid + 1]) {
		if (u < knots[mid]) high = mid;
		else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

CPoint CCG_QMSLib::DeBoorAlgorithmSurface1(const std::vector<std::vector<CPoint>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, double u, double v)
{
	// Search for the knot intervals of u and v
	int i = FindSpan(degreeU, knotVectorU, u);
	int j = FindSpan(degreeV, knotVectorV, v);

	// Recursion along the U direction
	std::vector<CPoint> temp(degreeV + 1);
	for (int col = 0; col <= degreeV; ++col) {
		std::vector<CPoint> column(degreeU + 1);
		for (int k = 0; k <= degreeU; ++k) {
			column[k] = controlPoints[i - degreeU + k][j - degreeV + col];
		}
		// Deboor algorithm along the U direction
		for (int r = 1; r <= degreeU; ++r) {
			for (int k = degreeU; k >= r; --k) {
				double alpha = (u - knotVectorU[i - degreeU + k]) /
					(knotVectorU[i + k - r + 1] - knotVectorU[i - degreeU + k]);
				for (int d = 0; d < 3; ++d) {
					column[k][d] = (1 - alpha) * column[k - 1][d] + alpha * column[k][d];
				}
			}
		}
		temp[col] = column[degreeU];
	}

	// Recursion along the V direction
	for (int r = 1; r <= degreeV; ++r) {
		for (int k = degreeV; k >= r; --k) {
			double alpha = (v - knotVectorV[j - degreeV + k]) /
				(knotVectorV[j + k - r + 1] - knotVectorV[j - degreeV + k]);
			for (int d = 0; d < 3; ++d) {
				temp[k][d] = (1 - alpha) * temp[k - 1][d] + alpha * temp[k][d];
			}
		}
	}

	return temp[degreeV];
}

void CCG_QMSLib::getMaincurvature_normal(double& k1, double& k2, CPoint& m1, CPoint& m2, CPoint& bpser_n, int n, int p, std::vector<double> U, int m, int q, std::vector<double> V, std::vector<std::vector<CPoint>> P, double u, double v, int d, std::vector<std::vector<CPoint>>& S)
{
	CPoint pu, puu, puv, pv, pvv;
	double pux, puy, puz, puux, puuy, puuz, puvx, puvy, puvz,
		pvx, pvy, pvz, pvvx, pvvy, pvvz, pnx, pny, pnz;
	double nlenth, lenth, ulenth, vlenth, k_delta, n_delta, tanpha1, tanpha2;
	double E, F, G, L, M, N;

	getsurfaceDerivat_allDegree(n, p, U, m, q, V, P, u, v, d, S);
	pu = S[1][0];	/*The first-order partial derivative in the u direction, no derivative in the v direction.*/
	pux = pu[0]; puy = pu[1]; puz = pu[2];
	ulenth = sqrt(pux * pux + puy * puy + puz * puz);
	puu = S[2][0];	/*u方向的二阶偏导，v方向不导*/
	puux = puu[0]; puuy = puu[1]; puuz = puu[2];
	puv = S[1][1];	/*u方向的一阶偏导，v方向的一阶偏导*/
	puvx = puv[0]; puvy = puv[1]; puvz = puv[2];
	pv = S[0][1];	/*u方向不导，v方向一阶偏导*/
	pvx = pv[0]; pvy = pv[1]; pvz = pv[2];
	vlenth = sqrt(pvx * pvx + pvy * pvy + pvz * pvz);
	pvv = S[0][2];	/*u方向不导，v方向二阶偏导*/
	pvvx = pvv[0]; pvvy = pvv[1]; pvvz = pvv[2];
	pnx = puy * pvz - puz * pvy;
	pny = -(pux * pvz - puz * pvx);
	pnz = pux * pvy - puy * pvx;
	nlenth = sqrt(pnx * pnx + pny * pny + pnz * pnz);
	pnx = pnx / nlenth; pny = pny / nlenth; pnz = pnz / nlenth;
	/* 单位法向量 */
	bpser_n[0] = pnx; bpser_n[1] = pny; bpser_n[2] = pnz;
	E = pux * pux + puy * puy + puz * puz;
	F = pux * pvx + puy * pvy + puz * pvz;
	G = pvx * pvx + pvy * pvy + pvz * pvz;
	L = pnx * puux + pny * puuy + pnz * puuz;
	M = pnx * puvx + pny * puvy + pnz * puvz;
	N = pnx * pvvx + pny * pvvy + pnz * pvvz;
	/* 判断方程是否有根 */
	n_delta = (G * L - E * N) * (G * L - E * N) - 4. * (E * M - F * L) * (F * N - G * M);
	k_delta = (E * N + G * L - 2 * F * M) * (E * N + G * L - 2 * F * M) - 4 * (E * G - F * F) * (L * N - M * M);
	if (k_delta != 0. && n_delta != 0.) {
		tanpha1 = ((G * L - E * N) + sqrt(n_delta)) / (2 * (F * N - G * M));
		tanpha2 = ((G * L - E * N) - sqrt(n_delta)) / (2 * (F * N - G * M));
		k1 = ((E * N + G * L - 2 * F * M) + sqrt(k_delta)) / (2 * (E * G - F * F));	/* (-b + sqrt(delta))/2a */
		k2 = ((E * N + G * L - 2 * F * M) - sqrt(k_delta)) / (2 * (E * G - F * F));	/* (-b - sqrt(delta))/2a */
		m1[0] = pux + tanpha1 * pvx;
		m1[1] = puy + tanpha1 * pvy;
		m1[2] = puz + tanpha1 * pvz;
		lenth = sqrt(m1[0] * m1[0] + m1[1] * m1[1] + m1[2] * m1[2]);
		m1[0] = m1[0] / lenth; m1[1] = m1[1] / lenth; m1[2] = m1[2] / lenth;
		m2[0] = pux + tanpha2 * pvx;
		m2[1] = puy + tanpha2 * pvy;
		m2[2] = puz + tanpha2 * pvz;
		lenth = sqrt(m2[0] * m2[0] + m2[1] * m2[1] + m2[2] * m2[2]);
		m2[0] = m2[0] / lenth; m2[1] = m2[1] / lenth; m2[2] = m2[2] / lenth;
	}
}

void CCG_QMSLib::getsurfaceDerivat_allDegree(int n, int p, std::vector<double> U, int m, int q, std::vector<double> V, std::vector<std::vector<CPoint>> P, double u, double v, int d, std::vector<std::vector<CPoint>>& S)
{
	std::vector<std::vector<CPoint>> controlPoints;
	controlPoints.resize(P[0].size(), std::vector<CPoint>(P.size()));
	std::vector<std::vector<double>> Nu;
	//Nu.resize(n + 1, std::vector<double>(p + 1));
	Nu.resize(n, std::vector<double>(p + 1));
	std::vector<std::vector<double>> Nv;
	//Nv.resize(m + 1, std::vector<double>(q + 1));
	Nv.resize(m, std::vector<double>(q + 1));
	std::vector<CPoint> temp;
	temp.resize(q + 1);
	int u_degree = std::min(d, p);
	int v_degree = std::min(d, q);
	n = U.size() - p - 2;
	m = V.size() - p - 2;
	/*对控制点矩阵进行转置*/
	/*for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[0].size(); j++)
			controlPoints[j][i] = P[i][j];
	}*/

	int uspan = FindSpan(n, p, u, U);
	DersBasisFuns(uspan, u, p, u_degree, U, Nu);
	int vspan = FindSpan(m, q, v, V);
	DersBasisFuns(vspan, v, q, v_degree, V, Nv);
	for (int k = 0; k <= u_degree; k++) {
		for (int s = 0; s <= q; s++) {
			temp[s][0] = 0.0;
			temp[s][1] = 0.0;
			temp[s][2] = 0.0;
			for (int r = 0; r <= p; r++) {
				temp[s] = temp[s] + (P[uspan - p + r][vspan - q + s] * Nu[k][r]);
				//temp[s] = temp[s] + (controlPoints[uspan - p + r][vspan - q + s] * Nu[k][r]);
			}
		}
		int dd = std::min(d - k, v_degree);
		for (int l = 0; l <= dd; l++) {
			S[k][l][0] = 0.0;
			S[k][l][1] = 0.0;
			S[k][l][2] = 0.0;
			for (int s = 0; s <= q; s++) {
				S[k][l] = S[k][l] + temp[s] * Nv[l][s];
			}
		}
	}
}

void CCG_QMSLib::DersBasisFuns(int i, double u, int p, int n, std::vector<double> U, std::vector<std::vector<double>>& Der)
{
	std::vector<std::vector<double>> ndu; 	/* 存储基函数和节点之差 */
	ndu.resize(p + 1, std::vector<double>(p + 1));
	std::vector<double> right;
	right.resize(p + 1);
	std::vector<double> left;
	left.resize(p + 1);
	std::vector<std::vector<double>> a;
	a.resize(2, std::vector<double>(p + 1));
	double saved;
	double temp;
	int j1;
	int j2;
	int j;
	int r;
	int k;

	ndu[0][0] = 1.0;
	for (j = 1; j <= p; j++) {
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++) {
			/* 下三角 */
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			/* 上三角 */
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
	for (j = 0; j <= p; j++) 	/* 载入基函数的值*/
		Der[0][j] = ndu[j][p];
	/* 下面计算导数 */
	for (r = 0; r <= p; r++) {
		/* 改变a数组的行*/
		int s1 = 0;
		int s2 = 1;
		a[0][0] = 1.0;
		/* 循环计算k阶导数, k = 1,2,...,n */
		for (k = 1; k <= n; k++) {
			double d = 0.0;
			int rk = r - k;
			int pk = p - k;
			if (r >= k) {
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)    j1 = 1;
			else		     j1 = -rk;
			if (r - 1 <= pk) j2 = k - 1;
			else             j2 = p - r;
			for (j = j1; j <= j2; j++) {
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk) {
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			Der[k][r] = d;
			/* 转换行 */
			j = s1;
			s1 = s2;
			s2 = j;
		}
	}
	/* 对结果乘以正确的因子 */
	r = p;
	for (k = 1; k <= n; k++) {
		for (j = 0; j <= p; j++)
			Der[k][j] *= r;
		r *= (p - k);
	}
}

int CCG_QMSLib::FindSpan(int n, int p, double u, std::vector<double> U)
{
	if (u == U[n + 1]) return n;	/* 特殊情况 */
	int low = p;					/* 进行二分搜索 */
	int high = n + 1;
	int mid = (low + high) / 2;
	while (u < U[mid] || u >= U[mid + 1]) {
		if (u < U[mid]) high = mid;
		else
		{
			low = mid;
		}
		mid = (low + high) / 2;
	}
	return mid;
}

void CCG_QMSLib::findMaxPositionDiffForEachEdge_T(CCG_QMS_model& model, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxPosDiff)
{
	for (auto bsurf : model.bsplineSurfs())
	{
		for (auto bsurfBoundary : bsurf->boundarys())
		{
			if (!bsurfBoundary->boundary())
			{
				maxPosDiff.insert(std::pair<NUBE_bsurfBoudnary*, double>(bsurfBoundary, bsurfBoundary->diffG0_sampling()->diffG0()));
			}
		}
	}
}

void CCG_QMSLib::findMaxAngleForEachEdge_T(CCG_QMS_model& model, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxAngles)
{
	for (auto bsurf : model.bsplineSurfs())
	{
		for (auto bsurfBoundary : bsurf->boundarys())
		{
			if (!bsurfBoundary->boundary())
			{
				maxAngles.insert(std::pair<NUBE_bsurfBoudnary*, double>(bsurfBoundary, bsurfBoundary->diffG1_sampling()->diffG1()));
			}
		}
	}
}

void CCG_QMSLib::findMaxSecondDerivativeDiffForEachEdge_T(CCG_QMS_model& model, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxDuuDiff, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxDuvDiff, std::unordered_map<NUBE_bsurfBoudnary*, double>& maxDvvDiff)
{
	for (auto bsurf : model.bsplineSurfs())
	{
		for (auto bsurfBoundary : bsurf->boundarys())
		{
			if (!bsurfBoundary->boundary())
			{
				maxDuuDiff.insert(std::pair<NUBE_bsurfBoudnary*, double>(bsurfBoundary, bsurfBoundary->diff_duu_sampling()->diff_duu()));
				maxDuvDiff.insert(std::pair<NUBE_bsurfBoudnary*, double>(bsurfBoundary, bsurfBoundary->diff_duv_sampling()->diff_duv()));
				maxDvvDiff.insert(std::pair<NUBE_bsurfBoudnary*, double>(bsurfBoundary, bsurfBoundary->diff_dvv_sampling()->diff_dvv()));
			}
		}
	}
}

void CCG_QMSLib::computeBoundaryEdgeRatio(const std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>>& allSurfs)
{
	for (auto surf : allSurfs) {

		// 先取出这片曲面上的所有边界（通常是 4 条，也可能更少/更多）
		const std::vector<NUBE_bsurfBoudnary*>& boundaries = surf->boundarys();
		// 如果没有边界，直接返回 0
		if (boundaries.empty()) {
			//std::cout << "曲面片没有边界，无法计算边界长度比。请检查输入数据。" << std::endl;
			surf->boundaryEdgeRatio() = 0.0;
		}

		// 初始化最小值和最大值
		double minLen = std::numeric_limits<double>::infinity();
		double maxLen = 0.0;

		// 遍历每一条边界
		for (NUBE_bsurfBoudnary* bnd : boundaries) {
			if (!bnd) continue;
			// 取出这条边界上的半边列表（按逆时针顺序）
			const std::vector<H*>& halfedges = bnd->halfedges();

			// 如果这条边界半边列表为空，就跳过
			if (halfedges.empty()) continue;

			// 遍历当前边界的所有半边
			for (H* h : halfedges) {
				if (!h) continue;

				// 取源顶点和目标顶点
				const CPoint& p_src = h->source()->point();
				const CPoint& p_tgt = h->target()->point();
				// 点到点的欧氏距离
				double dx = p_tgt[0] - p_src[0];
				double dy = p_tgt[1] - p_src[1];
				double dz = p_tgt[2] - p_src[2];
				double len = std::sqrt(dx * dx + dy * dy + dz * dz);

				// 更新最小值和最大值
				if (len < minLen) minLen = len;
				if (len > maxLen) maxLen = len;
			}
		}

		// 如果遍历后最小值还是无穷大，说明没遍历到任何边，直接返回 0
		if (minLen == std::numeric_limits<double>::infinity() || minLen <= 1e-16) {
			//std::cout << "曲面片没有有效边界，无法计算边界长度比。请检查输入数据。" << std::endl;
		}

		double ratio = maxLen / minLen;
		surf->boundaryEdgeRatio() = ratio;

		/*std::cout << "曲面片 ID = " << surf->id()
			<< " 的边界最长边/最短边长度比 = " << ratio
			<< std::endl;*/
	}
}

void CCG_QMSLib::computeControlPointMinDistance(const std::vector<std::shared_ptr<CCG_QMS_bsplineSurf>>& allSurfs)
{
	for (auto surf : allSurfs) {
		if (!surf) continue;

		// 取出控制点网格：二维数组 qbb_cpts
		std::vector<std::vector<CPoint>>& ctrlPtsGrid = surf->cpts();
		// 如果没有控制点，直接跳过
		if (ctrlPtsGrid.empty() || ctrlPtsGrid[0].empty()) {
			//std::cout << "曲面片 ID=" << surf->id() << " 没有控制点，无法计算最小距离。\n";
			surf->minCtrlPtDist() = 0.0;
			continue;
		}

		// 将 minDist 初始化为无穷大
		double minDist = std::numeric_limits<double>::infinity();

		// 获取行数和列数
		size_t rows = ctrlPtsGrid.size();
		size_t cols = ctrlPtsGrid[0].size();
		// 如果只有一个控制点，最小距离为 0
		if (rows * cols < 2) {
			minDist = 0.0;
			surf->minCtrlPtDist() = minDist;
			//std::cout << "曲面片 ID=" << surf->id() << " 只有一个控制点，最小距离设为 0。\n";
			continue;
		}

		// 遍历所有控制点对 i1,j1 与 i2,j2 （只枚举一次，i2,j2 从后面开始以免重复）
		for (size_t i1 = 0; i1 < rows; ++i1) {
			for (size_t j1 = 0; j1 < cols; ++j1) {
				const CPoint& p1 = ctrlPtsGrid[i1][j1];
				// 后续从 (i1,j1) 之后的索引开始，避免重复计算
				for (size_t i2 = i1; i2 < rows; ++i2) {
					// 如果 i2 == i1，则 j2 从 j1+1 开始；否则 j2 从 0 开始
					size_t startJ = (i2 == i1 ? j1 + 1 : 0);
					for (size_t j2 = startJ; j2 < cols; ++j2) {
						const CPoint& p2 = ctrlPtsGrid[i2][j2];

						// 计算 p1 与 p2 之间的欧氏距离
						double dx = p2[0] - p1[0];
						double dy = p2[1] - p1[1];
						double dz = p2[2] - p1[2];
						double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
						if (dist < minDist && dist > 1e-16) {
							// dist>1e-16 用于排除同一点或数值抖动导致的近零距离
							minDist = dist;
						}
					}
				}
			}
		}

		// 如果始终没有更新过 minDist（例如所有点几乎重合），就把它设为 0
		if (minDist == std::numeric_limits<double>::infinity()) {
			minDist = 0.0;
		}
		// 写入到对象里
		surf->minCtrlPtDist() = minDist;

		// 可选：打印出来检查
		/*std::cout << "曲面片 ID=" << surf->id()
			<< " 的控制点最小距离 = " << minDist
			<< std::endl;*/
	}
}

CPoint2 CCG_QMSLib::projectPointToTriangle(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c)
{
	// 计算重心坐标
	double u, v, w;
	computeBarycentricCoordinates(p, a, b, c, u, v, w);

	// 将重心坐标裁剪到合法范围
	u = std::max(0.0, std::min(1.0, u));
	v = std::max(0.0, std::min(1.0, v));
	w = std::max(0.0, std::min(1.0, w));

	double sum = u + v + w;
	if (sum == 0.0) sum = 1.0;

	u /= sum;
	v /= sum;
	w /= sum;

	return a * u + b * v + c * w;
}

CPoint2 CCG_QMSLib::projectPointToSegment(const CPoint2& p, const CPoint2& a, const CPoint2& b)
{
	CPoint2 ab = b - a;
	CPoint2 ap = p - a;

	double dot = ap[0] * ab[0] + ap[1] * ab[1];
	double abLenSq = ab[0] * ab[0] + ab[1] * ab[1];

	if (abLenSq == 0.0) {
		// a和b点重合，直接返回a点
		return a;
	}

	double proj = dot / abLenSq;
	proj = std::max(0.0, std::min(1.0, proj)); // 限制在[0,1]范围内

	CPoint2 result = a;
	result[0] += ab[0] * proj;
	result[1] += ab[1] * proj;

	return result;
}

double CCG_QMSLib::distanceToSegment(const CPoint2& p, const CPoint2& a, const CPoint2& b)
{
	CPoint2 ab = b - a;
	CPoint2 ap = p - a;

	double dot = ap[0] * ab[0] + ap[1] * ab[1];
	double abLenSq = ab[0] * ab[0] + ab[1] * ab[1];

	if (abLenSq == 0.0) {
		// a和b点重合，直接返回点到a的距离
		return sqrt(ap[0] * ap[0] + ap[1] * ap[1]);
	}

	double proj = dot / abLenSq;
	proj = std::max(0.0, std::min(1.0, proj)); // 限制在[0,1]范围内

	CPoint2 closest = a;
	closest[0] += ab[0] * proj;
	closest[1] += ab[1] * proj;

	double dx = p[0] - closest[0];
	double dy = p[1] - closest[1];

	return sqrt(dx * dx + dy * dy);
}

double CCG_QMSLib::distanceToTriangle(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c)
{
	// 计算点到三条边的最短距离
	double d1 = distanceToSegment(p, a, b);
	double d2 = distanceToSegment(p, b, c);
	double d3 = distanceToSegment(p, c, a);

	return std::min(std::min(d1, d2), d3);
}

void CCG_QMSLib::computeBarycentricCoordinates(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c, double& u, double& v, double& w)
{
	// 计算三角形的总面积
	double totalArea = std::abs(signedArea(a, b, c));
	if (totalArea < 1e-10) {
		// 处理退化三角形
		u = v = w = 1.0 / 3.0;
		return;
	}

	// 计算子三角形面积
	double area1 = std::abs(signedArea(p, b, c));
	double area2 = std::abs(signedArea(p, c, a));
	double area3 = std::abs(signedArea(p, a, b));

	// 计算重心坐标
	u = area1 / totalArea;
	v = area2 / totalArea;
	w = area3 / totalArea;

	// 处理数值误差，确保u+v+w=1
	double sum = u + v + w;
	if (std::abs(sum - 1.0) > 1e-10) {
		u /= sum;
		v /= sum;
		w /= sum;
	}
}

void CCG_QMSLib::computeQuadWeights(M* quadMesh, F* quadFace, const std::vector<H*>& fhs, const CPoint2& uv, QB_sampling& samplingPoint)
{
	const double diff_area = 1e-10; // 用于边界判断的容差

	// 四边形网格面第124个点组成的三角形的有向面积
	double area_tri_124 = signedArea(fhs[0]->uv(), fhs[1]->uv(), fhs[3]->uv());

	// 检查点是否在四边形的边上
	// 边12
	double area_tri_12s = signedArea(fhs[0]->uv(), fhs[1]->uv(), uv);
	if (std::abs(area_tri_12s) < diff_area) {
		double e_w = (uv - fhs[0]->uv()).norm() / ((fhs[1]->uv() - fhs[0]->uv()).norm());
		double f_w0 = 1 - e_w;
		double f_w1 = e_w;
		double f_w2 = 0.0;
		double f_w3 = 0.0;
		samplingPoint.ws().push_back(f_w0);
		samplingPoint.ws().push_back(f_w1);
		samplingPoint.ws().push_back(f_w2);
		samplingPoint.ws().push_back(f_w3);
		samplingPoint.boundary() = true;
		return;
	}

	// 计算w3权重
	double w3 = area_tri_12s / area_tri_124;

	// 计算w0权重 (三角形24s的面积)
	double area_tri_24s = signedArea(fhs[1]->uv(), fhs[3]->uv(), uv);
	double w0 = area_tri_24s / area_tri_124;

	// 边41
	double area_tri_41s = signedArea(fhs[3]->uv(), fhs[0]->uv(), uv);
	if (std::abs(area_tri_41s) < diff_area) {
		double e_w = (uv - fhs[3]->uv()).norm() / ((fhs[3]->uv() - fhs[0]->uv()).norm());
		double f_w0 = e_w;
		double f_w1 = 0.0;
		double f_w2 = 0.0;
		double f_w3 = 1 - e_w;
		samplingPoint.ws().push_back(f_w0);
		samplingPoint.ws().push_back(f_w1);
		samplingPoint.ws().push_back(f_w2);
		samplingPoint.ws().push_back(f_w3);
		samplingPoint.boundary() = true;
		return;
	}

	// 计算w1权重
	double w1 = area_tri_41s / area_tri_124;

	// 边23
	double area_tri_23s = signedArea(fhs[1]->uv(), fhs[2]->uv(), uv);
	if (std::abs(area_tri_23s) < diff_area) {
		double e_w = (uv - fhs[1]->uv()).norm() / ((fhs[2]->uv() - fhs[1]->uv()).norm());
		double f_w0 = 0.0;
		double f_w1 = 1 - e_w;
		double f_w2 = e_w;
		double f_w3 = 0.0;
		samplingPoint.ws().push_back(f_w0);
		samplingPoint.ws().push_back(f_w1);
		samplingPoint.ws().push_back(f_w2);
		samplingPoint.ws().push_back(f_w3);
		samplingPoint.boundary() = true;
		return;
	}

	// 边34
	double area_tri_34s = signedArea(fhs[2]->uv(), fhs[3]->uv(), uv);
	if (std::abs(area_tri_34s) < diff_area) {
		double e_w = (uv - fhs[2]->uv()).norm() / ((fhs[3]->uv() - fhs[2]->uv()).norm());
		double f_w0 = 0.0;
		double f_w1 = 0.0;
		double f_w2 = 1 - e_w;
		double f_w3 = e_w;
		samplingPoint.ws().push_back(f_w0);
		samplingPoint.ws().push_back(f_w1);
		samplingPoint.ws().push_back(f_w2);
		samplingPoint.ws().push_back(f_w3);
		samplingPoint.boundary() = true;
		return;
	}

	// 计算w2权重
	double w2 = 1.0 - w0 - w1 - w3;

	// 存储权重
	samplingPoint.ws().push_back(w0);
	samplingPoint.ws().push_back(w1);
	samplingPoint.ws().push_back(w2);
	samplingPoint.ws().push_back(w3);
	samplingPoint.boundary() = false;
}

double CCG_QMSLib::signedArea(CPoint2 p0, CPoint2 p1, CPoint2 p2)
{
	double area = ((p1[0] * p2[1] - p1[1] * p2[0]) - (p0[0] * p2[1] - p0[1] * p2[0]) + (p0[0] * p1[1] - p0[1] * p1[0])) * 0.5;
	return area;
}

bool CCG_QMSLib::isPointInTriangle(const CPoint2& p, const CPoint2& a, const CPoint2& b, const CPoint2& c)
{
	// 使用三个有向面积判断点是否在三角形内
	double area1 = signedArea(a, b, p);
	double area2 = signedArea(b, c, p);
	double area3 = signedArea(c, a, p);

	// 允许一定的数值误差
	const double EPSILON = 1e-10;

	// 如果三个面积都是同号的（都为正或都为负），则点在三角形内部
	return ((area1 >= -EPSILON && area2 >= -EPSILON && area3 >= -EPSILON) ||
		(area1 <= EPSILON && area2 <= EPSILON && area3 <= EPSILON));
}

bool CCG_QMSLib::locateTrianglePointByUV1(M* triMesh, M* quadMesh, F* quadFace, const CPoint2& uv, QB_sampling& samplingPoint)
{
	// --------------------------------------------------------------------------
	// Important: Ensure triMesh->classifyFacesByPatch() has been called BEFORE this function!
	// --------------------------------------------------------------------------

	// 1. Initialize sampling point info
	samplingPoint.uv() = uv;
	samplingPoint.uv_init() = uv;
	samplingPoint.fId() = quadFace->id();
	samplingPoint.patchId() = quadFace->patchIndex();

	// 2. Prepare Quad Face Halfedges (for weights calculation)
	H* fh = quadMesh->faceHalfedge(quadFace);
	std::vector<H*> fhs; // Or use std::array<H*, 4> if performance is critical
	fhs.reserve(4);
	H* currH = quadMesh->halfedgeNext(fh); // start from next
	fhs.push_back(currH);
	currH = quadMesh->halfedgeNext(currH);
	fhs.push_back(currH);
	currH = quadMesh->halfedgeNext(currH);
	fhs.push_back(currH);
	fhs.push_back(fh); // put original last as per your original logic

	// 3. Retrieve Candidate Triangles via Spatial Index (The Core Optimization)
	int pIdx = quadFace->patchIndex();
	int sIdx = quadFace->subPatchIndex();

	// Access the map from the mesh class
	const auto& fullMap = triMesh->patchMap_mesh();

	// Safe Lookup: Check if Patch exists
	auto itPatch = fullMap.find(pIdx);
	if (itPatch == fullMap.end()) {
		// Patch ID not found in triangle mesh
		return false;
	}

	// Safe Lookup: Check if SubPatch exists
	const auto& subMap = itPatch->second;
	auto itSub = subMap.find(sIdx);
	if (itSub == subMap.end()) {
		// SubPatch ID not found in triangle mesh
		return false;
	}

	// Get the reference to the specific list of faces
	const std::vector<F*>& candidates = itSub->second;

	if (candidates.empty()) return false;

	// 4. Iterate ONLY through candidate triangles
	bool found = false;
	F* closestTriangle = nullptr;
	double minDistance = std::numeric_limits<double>::max();

	// Cache for closest triangle data
	std::vector<V*> closestTriVertices;
	std::vector<CPoint2> closestTriUVs;

	// Reusable containers to avoid reallocation inside loop
	std::vector<V*> triVertices;
	std::vector<CPoint2> triUVs;
	triVertices.reserve(3);
	triUVs.reserve(3);

	for (F* tf : candidates) {
		// No need to check patchIndex/subPatchIndex again, the map guarantees it.

		// --- Extract Triangle Data ---
		triVertices.clear();
		triUVs.clear();

		H* triHE = triMesh->faceHalfedge(tf);
		if (!triHE) continue;

		H* currHE = triHE;
		do {
			V* v = triMesh->halfedgeTarget(currHE);
			triVertices.push_back(v);
			triUVs.push_back(currHE->uv());
			currHE = triMesh->halfedgeNext(currHE);
		} while (currHE != triHE);

		if (triUVs.size() != 3) continue;

		// --- Exact Match Check ---
		if (isPointInTriangle(uv, triUVs[0], triUVs[1], triUVs[2])) {
			// Found! Calculate Barycentric
			double u, v, w;
			computeBarycentricCoordinates(uv, triUVs[0], triUVs[1], triUVs[2], u, v, w);

			// Interpolate Position
			CPoint position = triVertices[0]->point() * u +
				triVertices[1]->point() * v +
				triVertices[2]->point() * w;

			samplingPoint.pos() = position;

			// Compute Weights
			computeQuadWeights(quadMesh, quadFace, fhs, uv, samplingPoint);

			return true; // Return immediately
		}

		// --- Distance Check (Fallback) ---
		// If not inside, keep track of the closest one
		double dist = distanceToTriangle(uv, triUVs[0], triUVs[1], triUVs[2]);
		if (dist < minDistance) {
			minDistance = dist;
			closestTriangle = tf;
			closestTriVertices = triVertices; // Copy data
			closestTriUVs = triUVs;           // Copy data
		}
	}

	// 5. Fallback: Use Closest Triangle (if exact match not found)
	if (closestTriangle != nullptr) {
		CPoint2 projectedUV = projectPointToTriangle(uv, closestTriUVs[0], closestTriUVs[1], closestTriUVs[2]);

		double u, v, w;
		computeBarycentricCoordinates(projectedUV, closestTriUVs[0], closestTriUVs[1], closestTriUVs[2], u, v, w);

		CPoint position = closestTriVertices[0]->point() * u +
			closestTriVertices[1]->point() * v +
			closestTriVertices[2]->point() * w;

		samplingPoint.pos() = position;

		computeQuadWeights(quadMesh, quadFace, fhs, uv, samplingPoint);

		return true;
	}

	return false;
}

void CCG_QMSLib::sampleQuadFace(M* triMesh, M* quadMesh, F* quadFace, const std::vector<CPoint2>& sampleUVs, std::vector<QB_sampling>& outSamplings)
{
	// 先清空输出容器
	outSamplings.clear();

	// 1. 先把 quadFace 对应的四条半边按顺序取出来，后面计算权重时会需要
	H* fh = quadMesh->faceHalfedge(quadFace);
	H* fhN1 = quadMesh->halfedgeNext(fh);
	H* fhN2 = quadMesh->halfedgeNext(fhN1);
	H* fhN3 = quadMesh->halfedgeNext(fhN2);
	std::vector<H*> fhs = { fhN1, fhN2, fhN3, fh };

	// 2. 对每个给定的 (u,v) 参数，调用 locateTrianglePointByUV1 去三角网格中定位
	for (const CPoint2& uv : sampleUVs)
	{
		QB_sampling st;
		// 先填充一些基本信息
		st.uv() = uv;
		st.uv_init() = uv;
		st.fId() = quadFace->id();          // 采样点归属的四边形面 ID
		st.patchId() = quadFace->patchIndex();  // 也可以把 patchId 设置为 quadFace 的 patchIndex

		// 2.1 先尝试“精确”在三角网格中找出包含此 UV 的那个三角形，
		//     并自动插值出三维坐标及重心权重。
		bool found = locateTrianglePointByUV1(triMesh, quadMesh, quadFace, uv, st);
		if (!found) {
			// 如果 locateTrianglePointByUV1 返回 false，说明该 UV 点不在对应 patch 的任何三角形内部。
			// 此时我们可以选择“跳过”或“继续用最近三角形投影”的方式，让 st.pos() 仍然有一个值。
			// 下面这里假设直接 skip，不把它加入输出：
			continue;
		}

		// 2.2 计算该采样点在四边形四个顶点（fhs）上的插值权重
		computeQuadWeights(quadMesh, quadFace, fhs, uv, st);

		// 2.3 把这一条采样信息加入输出容器
		outSamplings.push_back(st);
	}
}

void CCG_QMSLib::ensureQuadFaceSampling(M* triMesh, M* quadMesh, int minSamples, std::vector<QB_sampling>& sts)
{
	if (minSamples <= 0) return;

	// ==========================================================
	// 1. 优化统计：使用 Vector 代替 Map (从 O(logN) -> O(1))
	// ==========================================================

	// 假设 n_faces() 返回最大面索引范围。如果面索引不连续，需先进行垃圾回收(garbage collection)或整理
	int numFaces = quadMesh->numFaces();

	// 初始化计数器，利用连续内存提升缓存命中率
	std::vector<int> countPerFace(numFaces, 0);

	// 统计现有采样点
	for (auto& st : sts) {
		int _fid = st.fId();
		// 边界检查，防止脏数据导致越界
		if (_fid >= 0 && _fid < numFaces) {
			countPerFace[_fid]++;
		}
	}

	// ==========================================================
	// 2. OpenMP 并行计算补采样
	// ==========================================================

	// 获取最大线程数
	int max_threads = omp_get_max_threads();

	// 线程局部存储：每个线程拥有自己的结果容器，避免 push_back 竞争
	std::vector<std::vector<QB_sampling>> thread_buffers(max_threads);

	// 开启并行区域
#pragma omp parallel
	{
		int tid = omp_get_thread_num();

		// 【核心优化】将临时变量的内存分配提到循环外
		// 在该线程生命周期内复用内存，避免循环内反复 new/delete
		std::vector<CPoint2> cornerUVs(4);
		std::vector<CPoint2> sampleUVs;
		sampleUVs.reserve(minSamples); // 预留最大可能需要的空间
		std::vector<QB_sampling> newSamples;
		newSamples.reserve(minSamples);

		// 使用 dynamic 调度：因为不同面缺失的采样点数量差异可能很大，动态调度负载更均衡
#pragma omp for schedule(dynamic)
		for (int i = 0; i < numFaces; ++i) {

			// 快速检查：如果当前面采样足够，直接跳过
			if (countPerFace[i] >= minSamples) {
				continue;
			}

			// 获取面句柄 (假设 OpenMesh 风格接口，可以通过索引获取 handle)
			auto f = quadMesh->idFace(i+1);

			// 检查句柄有效性 (处理可能存在的已删除面)
			//if (!quadMesh->is_valid_handle(f)) continue;

			int needed = minSamples - countPerFace[i];

			// --- 获取 UV 坐标 ---
			H* fh = quadMesh->faceHalfedge(f);
			H* h1 = quadMesh->halfedgeNext(fh);
			H* h2 = quadMesh->halfedgeNext(h1);
			H* h3 = quadMesh->halfedgeNext(h2);

			// 直接赋值，无需重新分配 cornerUVs 大小
			cornerUVs[0] = h1->uv();
			cornerUVs[1] = h2->uv();
			cornerUVs[2] = h3->uv();
			cornerUVs[3] = fh->uv();

			// --- 计算采样点 UV (双线性插值) ---
			sampleUVs.clear(); // 仅重置 size，不释放 capacity

			// 计算网格密度
			int n = static_cast<int>(std::ceil(std::sqrt((double)needed)));
			double inv_denom = 1.0 / (n + 1); // 预计算除法

			// 生成 UV 采样网格
			for (int r = 1; r <= n && (int)sampleUVs.size() < needed; ++r) {
				double s = r * inv_denom;
				double one_minus_s = 1.0 - s;

				for (int c = 1; c <= n && (int)sampleUVs.size() < needed; ++c) {
					double t = c * inv_denom;
					double one_minus_t = 1.0 - t;

					// Bilinear Interpolation: 
					// P(s,t) = (1-s)(1-t)v0 + s(1-t)v1 + stv2 + (1-s)tv3
					CPoint2 uvSample =
						cornerUVs[0] * (one_minus_s * one_minus_t) +
						cornerUVs[1] * (s * one_minus_t) +
						cornerUVs[2] * (s * t) +
						cornerUVs[3] * (one_minus_s * t);

					sampleUVs.push_back(uvSample);
				}
			}

			// --- 调用原始采样逻辑 ---
			newSamples.clear();
			// 注意：此处假设 sampleQuadFace 是线程安全的（只读 Mesh，不写全局变量）
			sampleQuadFace(triMesh, quadMesh, f, sampleUVs, newSamples);

			// --- 将结果存入线程局部 Buffer ---
			int added = 0;
			for (auto& stNew : newSamples) {
				if (added >= needed) break;
				thread_buffers[tid].push_back(stNew);
				++added;
			}
		} // end for loop
	} // end parallel region

	// ==========================================================
	// 3. 合并结果
	// ==========================================================

	// 计算总共增加了多少点，以便一次性 reserve
	size_t total_new_samples = 0;
	for (const auto& buf : thread_buffers) {
		total_new_samples += buf.size();
	}

	// 预分配内存，避免 insert 过程中的多次 reallocate
	sts.reserve(sts.size() + total_new_samples);

	// 批量插入
	for (const auto& buf : thread_buffers) {
		sts.insert(sts.end(), buf.begin(), buf.end());
	}
}

/*-----------------------------------------------------------------------------------------------*/

