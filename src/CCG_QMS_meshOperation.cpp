#include"CCG_QMS_meshOperation.h"
#include<unordered_map>
#include<array>
#include <unordered_set>
#include <limits>
#include <omp.h>

//Attribute to quad mesh to manifold spline
/*-----------------------------------------------------------------------------------------------*/
void CCG_QMSLib::markFeatureEdgesBySharpEdges(M* pMesh)
{
	/*mark boundary edges to feature*/
	for (auto e : It::MEIterator(pMesh))
	{
		if (e->boundary())
		{
			e->feature() = true;
			e->sharp() = true;
		}
	}

	/*
	* Mark the feature points to facilitate the retention of
	* the original features when extracting the control mesh later.
	*/
	for (auto v : It::MVIterator(pMesh))
	{
		for (auto ve : It::VCcwEIterator(pMesh, v))
		{
			if (ve->feature())
			{
				v->feature() = true;
				break;
			}
		}
	}
}

void CCG_QMSLib::arrangeMeshVertexOrder(M* pMesh)
{
	int count = 0;
	for (auto v : It::MVIterator(pMesh))
	{
		count++;
		v->id() = count;
	}
	/*To prevent errors when using the idVertex function, the following updates need to be made*/
	pMesh->updateMapVertex();
}

void CCG_QMSLib::magnifyMeshPoint(M* pMesh, double scale)
{
	for (auto v : It::MVIterator(pMesh))
	{
		v->point() = v->point() * scale;
	}
}

void CCG_QMSLib::markTLayout_quadMesh(M* quadMesh)
{
	for (auto e : It::MEIterator(quadMesh))
	{
		if (e->boundary())
		{
			e->tLayout() = true;
		}
		else
		{
			F* ef1 = quadMesh->edgeFace1(e);
			F* ef2 = quadMesh->edgeFace2(e);
			/*
			* The subPatchIndex is different on faces linking
			* the TLayout edge
			*
			*/
			if (ef1->subPatchIndex() != ef2->subPatchIndex())
			{
				e->tLayout() = true;
			}
		}
	}
}

int CCG_QMSLib::statisticTLayoutPatchNum(M* quadMesh)
{
	int num_tLayout = 0;
	for (auto f : It::MFIterator(quadMesh))
	{
		if (f->subPatchIndex() > num_tLayout)
		{
			num_tLayout = f->subPatchIndex();
		}
	}
	num_tLayout++;
	quadMesh->tLayout_num() = num_tLayout;
	return num_tLayout;
}

void CCG_QMSLib::statisticTLayoutPatch_quadNum(M* quadMesh)
{
	/*
	* Define a vector to statistic the number of quad mesh faces on each TLayout
	*/
	std::vector<int> tLayout_quadNums(quadMesh->tLayout_num());
	/*
	* Initialize the element to 0
	*/
	for (int i = 0; i < tLayout_quadNums.size(); i++)
	{
		tLayout_quadNums[i] = 0;
	}

	for (auto f : It::MFIterator(quadMesh))
	{
		/*caution: default subpatchIndex starts from 0*/
		int f_tLayoutIndex = f->subPatchIndex();
		tLayout_quadNums[f_tLayoutIndex]++;
	}
	/*
	* update the result of statistics about the number of quad mesh faces on each TLayout
	*/
	for (int i = 0; i < tLayout_quadNums.size(); i++)
	{
		quadMesh->tLayoutPatches_quadNum().insert(std::pair<int, int>(i, tLayout_quadNums[i]));
	}

	/*check*/
	/*int totalQuadFaceNum = 0;
	for (auto tLPQM : quadMesh->tLayoutPatches_quadNum())
	{
		totalQuadFaceNum += tLPQM.second;
		std::cout <<" TLayout patch index: " << tLPQM.first << " # TLayout patch quad: " << tLPQM.second << std::endl;
	}
	std::cout << "totalQuadFaceNum: " << totalQuadFaceNum << std::endl;*/

}

void CCG_QMSLib::compute_mesh_normal(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		CPoint n(0, 0, 0);
		for (auto pF : It::VCcwFIterator(pMesh, v))
		{
			CPoint p[3];
			CHalfEdge* he = pF->halfedge();
			for (int k = 0; k < 3; k++)
			{
				p[k] = he->target()->point();
				he = he->he_next();
			}

			CPoint fn = (p[1] - p[0]) ^ (p[2] - p[0]);
			pF->normal() = fn / fn.norm();
			n += fn;
		}

		n = n / n.norm();
		v->normal() = n;
	}
}

void CCG_QMSLib::markExtraordinaryPts_cw(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		int n_vfs = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			n_vfs++;
		}
		if (v->boundary())
		{
			if (n_vfs == 2)
			{
				v->ifSingular() = false;
			}
			else
			{
				v->ifSingular() = false;
				//std::cout << " Boundary singular pt,degree " << n_vfs << std::endl;
			}
		}
		else
		{
			if (n_vfs == 4)
			{
				v->ifSingular() = false;
			}
			else
			{
				v->ifSingular() = true;
				//std::cout << " Inner singular pt,degree " << n_vfs << std::endl;

			}
		}
	}
}

void CCG_QMSLib::markExtraordinaryPts(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		int n_vfs = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			n_vfs++;
		}
		if (v->boundary())
		{
			if (n_vfs == 2)
			{
				v->ifSingular() = false;
			}
			else
			{
				v->ifSingular() = true;
			}
		}
		else
		{
			if (n_vfs == 4)
			{
				v->ifSingular() = false;
			}
			else
			{
				v->ifSingular() = true;
			}
		}
	}
}

void CCG_QMSLib::markHalfedgeLocalId(M* pMesh)
{
	for (auto f : It::MFIterator(pMesh))
	{
		H* f_h = pMesh->faceHalfedge(f);
		H* f_hNh = pMesh->halfedgeNext(f_h);
		H* f_hNhNh = pMesh->halfedgeNext(f_hNh);
		H* f_hPre = pMesh->halfedgePrev(f_h);
		f_h->localId() = 2;
		f_hNh->localId() = 3;
		f_hNhNh->localId() = 0;
		f_hPre->localId() = 1;
	}
}

void CCG_QMSLib::markCorner_quadMesh(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		v->ifCorner() = false;
		/*corner must be a boundary vertex and only linking 1 face
		* or all edges are feature
		*/
		//if (!v->boundary())
		//{
		//	bool ifFeature = true;
		//	for (auto ve : It::VClwEIterator(pMesh,v))
		//	{
		//		if (!ve->feature())
		//			ifFeature = false;
		//	}
		//	if (ifFeature)
		//	{
		//		v->ifCorner() = true;
		//		//std::cout << "------Find a corner on Quad Mesh inner!" << std::endl;
		//	}
		//}

		//int num_v_f = 0;
		//for (auto vf : It::VCcwFIterator(pMesh, v))
		//{
		//	num_v_f++;
		//}
		//if (num_v_f == 1)
		//{
		//	v->ifCorner() = true;
		//	//std::cout << "------Find a corner on Quad Mesh boundary!" << std::endl;
		//}
		bool ifFeature = true;
		for (auto ve : It::VClwEIterator(pMesh, v))
		{
			if (!(ve->feature()))
				ifFeature = false;
		}
		
		if (ifFeature)
		{
			v->ifCorner() = true;
			//std::cout << "------Find a corner on Quad Mesh inner!" << std::endl;
		}

		//adjust to singular boudary vertex
		int vNumFaces = 0;
		for (auto ve : It::VCcwFIterator(pMesh, v))
		{
			vNumFaces++;
		}

		if (vNumFaces != 2)
		{
			bool tempMark = true;
			for (auto ve : It::VClwEIterator(pMesh, v))
			{
				if (!(ve->feature() || ve->constrinedBoundary()))
					tempMark = false;
			}
			if (tempMark)
			{
				v->ifCorner() = true;
			}
		}
	}

}

void CCG_QMSLib::markCorner_feature_quadMesh(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		v->ifCorner() = false;
		/*corner must be a boundary vertex and only linking 1 face
		* or linking 2 feature edges both are on the one quad face
		*/
		//if (!v->boundary()) continue;
		int num_v_f = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			num_v_f++;
		}
		if (num_v_f == 1)
		{
			v->ifCorner() = true;
			//std::cout << "------Find a corner on Quad Mesh!" << std::endl;
		}
		else if (!v->boundary())
		{
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				if (pMesh->halfedgeEdge(vh)->feature() && pMesh->halfedgeEdge(pMesh->halfedgePrev(vh))->feature())
				{
					v->ifCorner() = true;
				}
			}
		}
	}
}

void CCG_QMSLib::remarkSingularPts(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->ifSingular() && v->feature())
		{
			v->ifSingular() = false;
		}
	}
}

void CCG_QMSLib::computeV_degree_mesh(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		int num_vf = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			num_vf++;
		}
		v->degree() = num_vf;
	}
}

void CCG_QMSLib::magnifyMeshUV(M* pMesh, double scale)
{
	for (auto f : It::MFIterator(pMesh))
	{
		for (auto fh : It::FHEIterator(pMesh, f))
		{
			fh->uv() = fh->uv() * scale;
		}
	}
}

void CCG_QMSLib::createTriQuadBoundaryVertexRealationByParameterization4(M* triMesh, M* quadMesh)
{
	/*obtaing all vertices from tri Mesh whose original id is valid*/
	std::vector<H*> triMeshBoundaryHEs;
	for (auto f : It::MFIterator(triMesh))
	{
		for (auto fh : It::FHEIterator(triMesh, f))
		{
			if (triMesh->halfedgeSource(fh)->originalId() != -1 || triMesh->halfedgeTarget(fh)->originalId() != -1)
			{
				triMeshBoundaryHEs.push_back(fh);
			}
		}

	}
	/*obtaining all boundary halfedges from quad mesh*/
	std::vector<H*> quadMeshBoundaryHEs;
	for (auto f : It::MFIterator(quadMesh))
	{
		for (auto fh : It::FHEIterator(quadMesh, f))
		{
			if (quadMesh->halfedgeEdge(fh)->boundary())
			{
				quadMeshBoundaryHEs.push_back(fh);
			}
		}
	}
	/*Find the quad mesh boundary face and halfedge for all triMeshBoundaryV*/
	for (auto tbhe : triMeshBoundaryHEs)
	{
		F* tbheF = triMesh->halfedgeFace(tbhe);
		for (auto qbhe : quadMeshBoundaryHEs)
		{
			F* qbheF = quadMesh->halfedgeFace(qbhe);
			if (tbheF->patchIndex() == qbheF->patchIndex() && tbheF->subPatchIndex() == qbheF->subPatchIndex() && tbhe->id_FeatureTrajectory() == qbhe->id_FeatureTrajectory())
			{
				tbhe->boundaryFaceId() = qbheF->id();
				tbhe->boundaryTargetVId() = qbhe->target()->id();
				break;
			}
		}
	}

	/*
	* check Make sure to find the half-edge corresponding to each point in the simplified triangular mesh
	* whose original Id is not -1 and whose id_FeatureTrajectory is not -1.
	*/
	//for (auto v : It::MVIterator(triMesh))
	//{
	//	v->ifSharp() = false;
	//	v->ifVisit() = false;
	//}

	//for (auto h : triMeshBoundaryHEs)
	//{
	//	if (triMesh->halfedgeSource(h)->originalId() != -1)
	//	{
	//		triMesh->halfedgeSource(h)->ifVisit() = true;
	//	}
	//	if (triMesh->halfedgeTarget(h)->originalId() != -1)
	//	{
	//		triMesh->halfedgeTarget(h)->ifVisit() = true;
	//	}
	//}

	//for (auto v : It::MVIterator(triMesh))
	//{
	//	if (v->originalId() != -1 && !v->ifVisit())
	//	{
	//		std::cout << "there is a  original Id != -1 vertex not on d_FeatureTrajectory ！= -1 halfedges " << std::endl;
	//		v->ifSharp() = true;
	//	}
	//	//v->ifSharp() = false;
	//	//v->ifVisit() = false;
	//}

	//for (auto v : It::MVIterator(triMesh))
	//{
	//	v->ifSharp() = false;
	//	v->ifVisit() = false;
	//}

	/*check*/
	/*for (auto tbhe : triMeshBoundaryHEs)
	{
		if (tbhe->boundaryFaceId() == -1)
		{
			std::cout << "there is an halfedge linking original vertex no nurbs boundary curve matching！" << std::endl;
		}
	}*/
}

void CCG_QMSLib::obtainTLayout_motorcycleGraph_featureFree(M* pMesh)
{
	/*1. Mark the boundary edges, ensuring that the boundary edges are definitely the edges on the QuadLayout, and also mark the boundary points*/
	for (auto e : It::MEIterator(pMesh))
	{
		e->sharp() = false;
		pMesh->edgeVertex1(e)->ifSharp() = false;
		pMesh->edgeVertex2(e)->ifSharp() = false;
		if (/*e->feature() ||*/ e->boundary())
		{
			e->sharp() = true;
			pMesh->edgeVertex1(e)->ifSharp() = true;
			pMesh->edgeVertex2(e)->ifSharp() = true;
		}

	}
	//2. Mark the singular points (including those on the boundaries and features)
	std::vector<V*> singularVs;
	for (auto v : It::MVIterator(pMesh))
	{
		//if (v->boundary())continue;
		int fNum = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			fNum++;
		}
		if (fNum != 4 && !v->boundary())
		{
			v->ifSingular() = true;
			v->ifSharp() = true;
			singularVs.push_back(v);
			//std::cout << "********" << std::endl;
		}
		else if (fNum != 2 && v->boundary())
		{
			v->ifSingular() = true;
			v->ifSharp() = true;
			singularVs.push_back(v);
		}

	}

	/*3. Starting from the singularity point, place the corresponding half-edges of all unmarked edges into the queue one by one*/
	std::queue<H*> hes;
	for (auto v : singularVs)
	{
		for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
		{
			E* vheE = pMesh->halfedgeEdge(vhe);
			if (vheE->sharp())continue;
			hes.push(vhe);
			//Mark the corresponding edge of the halfedge
			vheE->sharp() = true;
		}
	}
	while (!hes.empty())
	{
		/*4.Take halfedge of the items from the queue one by one*/
		H* he = hes.front();
		hes.pop();
		////(1).Mark the corresponding edge of the halfedge
		//E* heE = pMesh->halfedgeEdge(he);
		//heE->sharp() = true;
		/*(2).Determine whether one halfedge is pointing to the marked point.
		* If so, stop the operation; if not, place the corresponding halfedge of the next halfedge into the queue according to the direction of the halfedge,
		* and mark the point that the current half is pointing to.
		*/
		V* heTarget = pMesh->halfedgeTarget(he);
		if (heTarget->ifSharp()) continue;
		heTarget->ifSharp() = true;
		H* heN = pMesh->halfedgeNext(he);
		//(3).Determine whether the current heTarget point is a T node
		E* heNE = pMesh->halfedgeEdge(heN);
		H* heSym = pMesh->halfedgeSym(he);//heSymis not NULL
		H* heSymPre = pMesh->halfedgePrev(heSym);
		E* heSymPreE = pMesh->halfedgeEdge(heSymPre);
		if (heNE->sharp() && heSymPreE->sharp()) continue;
		H* heNSym = pMesh->halfedgeSym(heN);
		H* heNSymN = pMesh->halfedgeNext(heNSym);
		hes.push(heNSymN);
		E* heNsymNE = pMesh->halfedgeEdge(heNSymN);
		heNsymNE->sharp() = true;
	}

}

/*-----------------------------------------------------------------------------------------------*/

void CCG_QMSLib::obtainTLayout_motorcycleGraph(M* pMesh)
{
	/*1. Mark the boundary edges, ensuring that the boundary edges are definitely the edges on the QuadLayout, and also mark the boundary points*/
	for (auto e : It::MEIterator(pMesh))
	{
		e->sharp() = false;
		pMesh->edgeVertex1(e)->ifSharp() = false;
		pMesh->edgeVertex2(e)->ifSharp() = false;
		if (e->feature() || e->boundary())
		{
			e->sharp() = true;
			pMesh->edgeVertex1(e)->ifSharp() = true;
			pMesh->edgeVertex2(e)->ifSharp() = true;
		}

	}
	//2. Mark the singular points (including those on the boundaries and features)
	std::vector<V*> singularVs;
	for (auto v : It::MVIterator(pMesh))
	{
		//if (v->boundary())continue;
		int fNum = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			fNum++;
		}
		if (fNum != 4 && !v->boundary())
		{
			v->ifSingular() = true;
			v->ifSharp() = true;
			singularVs.push_back(v);
			//std::cout << "********" << std::endl;
		}
		else if (fNum != 2 && v->boundary())
		{
			v->ifSingular() = true;
			v->ifSharp() = true;
			singularVs.push_back(v);
		}
		else if (fNum == 4 && !v->boundary())
		{
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				if (pMesh->halfedgeEdge(vh)->feature() && pMesh->halfedgeEdge(pMesh->halfedgePrev(vh))->feature())
				{
					v->ifSharp() = true;
					singularVs.push_back(v);
				}
			}
		}

	}

	/*3. Starting from the singularity point, place the corresponding half-edges of all unmarked edges into the queue one by one*/
	std::queue<H*> hes;
	for (auto v : singularVs)
	{
		for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
		{
			E* vheE = pMesh->halfedgeEdge(vhe);
			if (vheE->sharp())continue;
			/**/

			hes.push(vhe);
			//Mark the corresponding edge of the halfedge
			vheE->sharp() = true;
		}
	}
	while (!hes.empty())
	{
		/*4.Take halfedge of the items from the queue one by one*/
		H* he = hes.front();
		hes.pop();
		////(1).Mark the corresponding edge of the halfedge
		//E* heE = pMesh->halfedgeEdge(he);
		//heE->sharp() = true;
		/*(2).Determine whether one halfedge is pointing to the marked point.
		* If so, stop the operation; if not, place the corresponding halfedge of the next halfedge into the queue according to the direction of the halfedge,
		* and mark the point that the current half is pointing to.
		*/
		V* heTarget = pMesh->halfedgeTarget(he);
		if (heTarget->ifSharp()) continue;
		heTarget->ifSharp() = true;
		H* heN = pMesh->halfedgeNext(he);
		//(3).Determine whether the current heTarget point is a T node
		E* heNE = pMesh->halfedgeEdge(heN);
		H* heSym = pMesh->halfedgeSym(he);//heSymis not NULL
		H* heSymPre = pMesh->halfedgePrev(heSym);
		E* heSymPreE = pMesh->halfedgeEdge(heSymPre);
		if (heNE->sharp() && heSymPreE->sharp()) continue;
		H* heNSym = pMesh->halfedgeSym(heN);
		H* heNSymN = pMesh->halfedgeNext(heNSym);
		hes.push(heNSymN);
		E* heNsymNE = pMesh->halfedgeEdge(heNSymN);
		heNsymNE->sharp() = true;
	}
}

void CCG_QMSLib::obtainTLayout_motorcycleGraph_featureAware(M* pMesh)
{
	/*1. Mark the boundary edges, ensuring that the boundary edges are definitely the edges on the QuadLayout, and also mark the boundary points*/
	for (auto e : It::MEIterator(pMesh))
	{
		e->sharp() = false;
		pMesh->edgeVertex1(e)->ifSharp() = false;
		pMesh->edgeVertex2(e)->ifSharp() = false;
	}
	for (auto e : It::MEIterator(pMesh))
	{
		if (e->feature() || e->boundary() || e->constrinedBoundary())
		{
			e->sharp() = true;
			pMesh->edgeVertex1(e)->ifSharp() = true;
			pMesh->edgeVertex2(e)->ifSharp() = true;
		}

	}

	/*2. collect outhalfedges to trace the quadlyout*/
	std::queue<H*> hes;
	int numSharpV = 0;
	int numBoundaryV = 0;
	/*
	* counting the number of middle edges from one edge of an initial mark
	* to another edge of an initial mark around the suspected corner point,
	* and if the number is not equal to 1, taking the middle edge as the boundary edge of the B-spline patch
	*/
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->boundary())
		{
			H* firstOutH = NULL;
			for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
			{
				E* vheE = pMesh->halfedgeEdge(vhe);
				if (vheE->boundary())
				{
					firstOutH = vhe;
					break;
				}
			}
			H* firstOutHPre = pMesh->halfedgePrev(firstOutH);
			H* firstOutHPreSym = pMesh->halfedgeSym(firstOutHPre);
			std::vector<H*> alternativeHs;
			bool firstMark = false;
			while (firstOutHPreSym != NULL && firstOutHPreSym->face()->id() != firstOutH->face()->id())
			{
				E* firstOutHPreSymE = pMesh->halfedgeEdge(firstOutHPreSym);
				if (firstOutHPreSymE->feature() || firstOutHPreSymE->constrinedBoundary())
				{
					if (!firstOutHPreSymE->feature() && firstOutHPreSymE->constrinedBoundary())
					{
						firstMark = true;
					}
					if (alternativeHs.size() != 1 && alternativeHs.size() > 0)
					{
						//std::cout << " boundary alternativeHs.size(): " << alternativeHs.size() << std::endl;
						for (auto alternativeH : alternativeHs)
						{
							E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
							if (alternativeHE->sharp())continue;
							hes.push(alternativeH);
							//Mark the corresponding edge of the halfedge
							alternativeHE->sharp() = true;
						}
					}
					else if (alternativeHs.size() == 1 && v->ifSingular())
					{
						if (firstMark)
						{
							for (auto alternativeH : alternativeHs)
							{
								E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
								if (alternativeHE->sharp())continue;
								hes.push(alternativeH);
								//Mark the corresponding edge of the halfedge
								alternativeHE->sharp() = true;
							}
						}
						else
						{
							firstMark = true;
						}
					}
					alternativeHs.clear();
				}
				else
				{
					alternativeHs.push_back(firstOutHPreSym);
				}
				firstOutHPre = pMesh->halfedgePrev(firstOutHPreSym);
				firstOutHPreSym = pMesh->halfedgeSym(firstOutHPre);
			}

			{
				if (alternativeHs.size() != 1 && alternativeHs.size() > 0)
				{
					//std::cout << " boudnary alternativeHs.size(): " << alternativeHs.size() << std::endl;
					for (auto alternativeH : alternativeHs)
					{
						E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
						if (alternativeHE->sharp())continue;
						hes.push(alternativeH);
						//Mark the corresponding edge of the halfedge
						alternativeHE->sharp() = true;
					}
				}
				else if (alternativeHs.size() == 1 && v->ifSingular())
				{
					if (firstMark)
					{
						for (auto alternativeH : alternativeHs)
						{
							E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
							if (alternativeHE->sharp())continue;
							hes.push(alternativeH);
							//Mark the corresponding edge of the halfedge
							alternativeHE->sharp() = true;
						}
					}
					else
					{
						firstMark = true;
					}
				}
				alternativeHs.clear();
			}
		}
		else if (v->ifSharp())
		{
			numSharpV++;
			H* firstOutH = NULL;
			for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
			{
				E* vheE = pMesh->halfedgeEdge(vhe);
				if (vheE->feature() || vheE->constrinedBoundary())
				{
					firstOutH = vhe;
					break;
				}
			}
			H* firstOutHPre = pMesh->halfedgePrev(firstOutH);
			H* firstOutHPreSym = pMesh->halfedgeSym(firstOutHPre);
			std::vector<H*> alternativeHs;
			bool firstMark = false;
			while (firstOutHPreSym->face()->id() != firstOutH->face()->id())
			{
				E* firstOutHPreSymE = pMesh->halfedgeEdge(firstOutHPreSym);
				if (!firstOutHPreSymE->feature() && firstOutHPreSymE->constrinedBoundary())
				{
					firstMark = true;
				}
				if (firstOutHPreSymE->feature() || firstOutHPreSymE->constrinedBoundary())
				{
					//std::cout << "-------------------" << std::endl;
					if (alternativeHs.size() != 1 && alternativeHs.size() > 0)
					{
						//std::cout << v->id()  << " feature 1 alternativeHs.size(): " << alternativeHs.size() << std::endl;
						for (auto alternativeH : alternativeHs)
						{
							E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
							if (alternativeHE->sharp())continue;
							hes.push(alternativeH);
							//Mark the corresponding edge of the halfedge
							alternativeHE->sharp() = true;
						}
					}
					else if (alternativeHs.size() == 1 && v->ifSingular())
					{
						if (firstMark)
						{
							for (auto alternativeH : alternativeHs)
							{
								E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
								if (alternativeHE->sharp())continue;
								hes.push(alternativeH);
								//Mark the corresponding edge of the halfedge
								alternativeHE->sharp() = true;
							}
						}
						else
						{
							firstMark = true;
						}
					}
					alternativeHs.clear();
				}
				else
				{
					alternativeHs.push_back(firstOutHPreSym);
					//std::cout << v->id() << " feature alternativeHs.size(): " << alternativeHs.size() << std::endl;
				}
				firstOutHPre = pMesh->halfedgePrev(firstOutHPreSym);
				firstOutHPreSym = pMesh->halfedgeSym(firstOutHPre);
			}

			{
				if (alternativeHs.size() != 1 && alternativeHs.size() > 0)
				{
					//std::cout << v->id() << " feature 2 alternativeHs.size(): " << alternativeHs.size() << std::endl;
					for (auto alternativeH : alternativeHs)
					{
						E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
						if (alternativeHE->sharp())continue;
						hes.push(alternativeH);
						//Mark the corresponding edge of the halfedge
						alternativeHE->sharp() = true;
					}
				}
				else if (alternativeHs.size() == 1 && v->ifSingular())
				{
					if (firstMark)
					{
						for (auto alternativeH : alternativeHs)
						{
							E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
							if (alternativeHE->sharp())continue;
							hes.push(alternativeH);
							//Mark the corresponding edge of the halfedge
							alternativeHE->sharp() = true;
						}
					}
					else
					{
						firstMark = true;
					}
				}
				alternativeHs.clear();
			}
		}
		else
		{
			int fNum = 0;
			for (auto vf : It::VCcwFIterator(pMesh, v))
			{
				fNum++;
			}
			if (fNum != 4)
			{
				v->ifSingular() = true;
				v->ifSharp() = true;
				for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
				{
					E* vheE = pMesh->halfedgeEdge(vhe);
					if (vheE->sharp())continue;
					hes.push(vhe);
					//Mark the corresponding edge of the halfedge
					vheE->sharp() = true;
				}
			}
		}
	}

	/*
	* For interior singular suspected corners,
	* count the number of middle edges from one feature edge to another feature edge around the suspected corners,
	* and if the number is not equal to 1, take the middle edge as the boundary edge of the B-spline patch
	*/
	markExtraordinaryPts(pMesh);
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->boundary())continue;
		else if (v->ifSharp() && v->ifSingular())
		{
			numSharpV++;
			H* firstOutH = NULL;
			for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
			{
				E* vheE = pMesh->halfedgeEdge(vhe);
				if (vheE->feature())
				{
					firstOutH = vhe;
					break;
				}
			}
			if (firstOutH == NULL)
			{
				v->ifSingular() = true;
				v->ifSharp() = true;
				for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
				{
					E* vheE = pMesh->halfedgeEdge(vhe);
					if (vheE->sharp())continue;
					hes.push(vhe);
					//Mark the corresponding edge of the halfedge
					vheE->sharp() = true;
				}
				continue;
			}
			H* firstOutHPre = pMesh->halfedgePrev(firstOutH);
			H* firstOutHPreSym = pMesh->halfedgeSym(firstOutHPre);
			std::vector<H*> alternativeHs;
			bool firstMark = false;
			while (firstOutHPreSym->face()->id() != firstOutH->face()->id())
			{
				E* firstOutHPreSymE = pMesh->halfedgeEdge(firstOutHPreSym);
				if (firstOutHPreSymE->feature())
				{
					if (alternativeHs.size() != 1 && alternativeHs.size() > 0)
					{
						for (auto alternativeH : alternativeHs)
						{
							E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
							if (alternativeHE->sharp())continue;
							hes.push(alternativeH);
							//Mark the corresponding edge of the halfedge
							alternativeHE->sharp() = true;
						}
					}
					alternativeHs.clear();
				}
				else
				{
					alternativeHs.push_back(firstOutHPreSym);
				}
				firstOutHPre = pMesh->halfedgePrev(firstOutHPreSym);
				firstOutHPreSym = pMesh->halfedgeSym(firstOutHPre);
			}

			{
				if (alternativeHs.size() != 1 && alternativeHs.size() > 0)
				{
					for (auto alternativeH : alternativeHs)
					{
						E* alternativeHE = pMesh->halfedgeEdge(alternativeH);
						if (alternativeHE->sharp())continue;
						hes.push(alternativeH);
						//Mark the corresponding edge of the halfedge
						alternativeHE->sharp() = true;
					}
				}
				alternativeHs.clear();

			}
		}
	}
	//std::cout << "numSharpV： " << numSharpV << std::endl;
	/*3. Starting from the singularity point, place the corresponding half-edges of all unmarked edges into the queue one by one*/
	while (!hes.empty())
	{
		/*4.Take halfedge of the items from the queue one by one*/
		H* he = hes.front();
		hes.pop();
		////(1).Mark the corresponding edge of the halfedge
		//E* heE = pMesh->halfedgeEdge(he);
		//heE->sharp() = true;
		/*(2).Determine whether one halfedge is pointing to the marked point.
		* If so, stop the operation; if not, place the corresponding halfedge of the next halfedge into the queue according to the direction of the halfedge,
		* and mark the point that the current half is pointing to.
		*/
		V* heTarget = pMesh->halfedgeTarget(he);
		if (heTarget->ifSharp()) continue;
		heTarget->ifSharp() = true;
		H* heN = pMesh->halfedgeNext(he);
		//(3).Determine whether the current heTarget point is a T node
		E* heNE = pMesh->halfedgeEdge(heN);
		H* heSym = pMesh->halfedgeSym(he);//heSymis not NULL
		H* heSymPre = pMesh->halfedgePrev(heSym);
		E* heSymPreE = pMesh->halfedgeEdge(heSymPre);
		if (heNE->sharp() && heSymPreE->sharp()) continue;
		H* heNSym = pMesh->halfedgeSym(heN);
		H* heNSymN = pMesh->halfedgeNext(heNSym);
		hes.push(heNSymN);
		E* heNsymNE = pMesh->halfedgeEdge(heNSymN);
		heNsymNE->sharp() = true;
	}
}

void CCG_QMSLib::obtainLayout_consistencyQuadMeshFaceNum(M* pMesh)
{
	/*1. Mark the boundary edges, ensuring that the boundary edges are definitely the edges on the QuadLayout, and also mark the boundary points*/
	for (auto e : It::MEIterator(pMesh))
	{
		e->sharp() = true;
		pMesh->edgeVertex1(e)->ifSharp() = true;
		pMesh->edgeVertex2(e)->ifSharp() = true;
	}
}

//Attribute to boundary curve fitting
/*-----------------------------------------------------------------------------------------------*/
void CCG_QMSLib::copyInitialMeshVPos(M* pMesh)
{
	computeCos(pMesh);
	computeGaussCur(pMesh);
	for (auto v : It::MVIterator(pMesh))
	{
		v->initialPos() = v->point();
	}
}

void CCG_QMSLib::computeCos(M* pMesh)
{
	for (auto f : It::MFIterator(pMesh))
	{
		for (auto fhe : It::FHEIterator(pMesh, f))
		{
			H* fhen = pMesh->halfedgeNext(fhe);
			CPoint fhenVec = fhen->target()->point() - fhen->source()->point();
			CPoint fheVec = fhe->source()->point() - fhe->target()->point();
			double fheTarAngCos = fhenVec * fheVec / (fhenVec.norm() * fheVec.norm());
			fhe->angleCos() = fheTarAngCos;
			//std::cout << "fhe->angleCos() " << fhe->angleCos() << std::endl;
		}
	}
}

void CCG_QMSLib::computeGaussCur(M* pMesh)
{
	double _M_PI = std::acos(-1.0);
	//std::cout << "pi: " << _M_PI << std::endl;
	for (auto v : It::MVIterator(pMesh))
	{
		double totalAngel = 0.0;
		for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
		{
			H* vhePre = pMesh->halfedgePrev(vhe);
			totalAngel += std::acos(vhePre->angleCos());
		}
		if (v->boundary())
		{
			v->gaussCur() = _M_PI - totalAngel;
			//std::cout << "Gauss Curvature: " << v->gaussCur() << std::endl;
		}
		else
		{
			v->gaussCur() = 2 * _M_PI - totalAngel;
		}
		//std::cout << "Gauss Curvature: " << v->gaussCur() << std::endl;
	}

	//std::cout << "Computing gauss curvature is finished! " << std::endl;
}
void CCG_QMSLib::recoverFlatBoundaryVertexPoint(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->boundary())
		{
			/*Use the boundary discrete Gaussian curvature of the current point*/
			if (abs(v->gaussCur()) < 1e-6)
			{
				/*std::cout << "--------------------------------------------------------" << std::endl;*/
				v->point() = v->initialPos();
			}
		}
	}
}

/*判断平面上三个点的有向面积是否不小于0*/
double diff_area = 1e-8;
bool signedAreaNotNegative(CPoint2 p0, CPoint2 p1, CPoint2 p2)
{
	double area = ((p1[0] * p2[1] - p1[1] * p2[0]) - (p0[0] * p2[1] - p0[1] * p2[0]) + (p0[0] * p1[1] - p0[1] * p1[0])) * 0.5;
	if (area > 0)
	{
		return true;
	}
	else if (abs(area) < diff_area)
	{
		// 检查点 p2 是否在 p0 和 p1 之间的线段上
		double minX = std::min(p0[0], p1[0]);
		double maxX = std::max(p0[0], p1[0]);
		double minY = std::min(p0[1], p1[1]);
		double maxY = std::max(p0[1], p1[1]);

		// 检查 p2 是否在 p0 和 p1 的边界框内
		if (p2[0] >= minX && p2[0] <= maxX && p2[1] >= minY && p2[1] <= maxY)
		{
			//std::cout << area << std::endl;
			///*输出参数点坐标*/
			//std::cout << "------" << std::endl;
			//std::cout << "Quad edge vertex 1 uv: " << p0[0] << " " << p0[1] << std::endl;
			//std::cout << "Quad edge vertex 2 uv: " << p1[0] << " " << p1[1] << std::endl;
			//std::cout << "Tri vertex         uv: " << p2[0] << " " << p2[1] << std::endl;
			//std::cout << "------" << std::endl;
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
double signedArea(CPoint2 p0, CPoint2 p1, CPoint2 p2)
{
	double area = ((p1[0] * p2[1] - p1[1] * p2[0]) - (p0[0] * p2[1] - p0[1] * p2[0]) + (p0[0] * p1[1] - p0[1] * p1[0])) * 0.5;
	return area;
}

void CCG_QMSLib::obtainSamplingFromTriMeshByParameterization(M* triMesh, M* quadMesh, std::vector<QB_sampling>& sts)
{
	// --------------------------------------------------------------------------
	// 1. 确保索引已构建
	// --------------------------------------------------------------------------
	// 这一步必须在并行区域之前串行执行
	// 重新构建索引，确保数据是最新的
	if (triMesh->patchMap_mesh().empty())
	{
		triMesh->classifyFacesByPatch();
	}

	// 获取只读引用，方便后续访问
	const auto& fullPatchMap = triMesh->patchMap_mesh();

	// --------------------------------------------------------------------------
	// 2. 准备并行容器
	// --------------------------------------------------------------------------
	int max_threads = omp_get_max_threads();
	// 每个线程独立的输出 buffer，避免 push_back 竞争
	std::vector<std::vector<QB_sampling>> thread_sts(max_threads);

	// --------------------------------------------------------------------------
	// 3. 并行遍历四边形网格
	// --------------------------------------------------------------------------
#pragma omp parallel
	{
		int tid = omp_get_thread_num();

		// 线程局部变量：预分配内存，避免循环内反复 new/delete
		std::vector<H*> fhs;
		fhs.reserve(4);

		// 局部去重容器：记录当前四边形已经处理过的三角形顶点ID
		// 替代原代码中修改全局 ifVisit 的逻辑
		std::unordered_set<int> visitedVertexIds;

#pragma omp for schedule(dynamic)
		for (auto f:It::MFIterator(quadMesh))
		{
			// --- A. 快速获取候选三角形 (利用内置索引) ---
			int pIdx = f->patchIndex();
			int sIdx = f->subPatchIndex();

			// 【关键安全操作】使用 find 而不是 []，避免多线程修改 map 导致崩溃
			auto itPatch = fullPatchMap.find(pIdx);
			if (itPatch == fullPatchMap.end()) continue; // 该 Patch 不存在

			const auto& subMap = itPatch->second;
			auto itSub = subMap.find(sIdx);
			if (itSub == subMap.end()) continue; // 该 SubPatch 不存在

			// 获取候选三角形列表
			const std::vector<F*>& candidateTris = itSub->second;
			if (candidateTris.empty()) continue;

			// --- B. 准备四边形几何信息 ---
			fhs.clear();
			H* fh = quadMesh->faceHalfedge(f);
			H* fhN1 = quadMesh->halfedgeNext(fh);
			H* fhN2 = quadMesh->halfedgeNext(fhN1);
			H* fhN3 = quadMesh->halfedgeNext(fhN2);
			fhs = { fhN1, fhN2, fhN3, fh };

			// --- C. 计算 UV 包围盒 (AABB) 用于预筛选 ---
			double min_u = std::numeric_limits<double>::max();
			double max_u = -std::numeric_limits<double>::max();
			double min_v = std::numeric_limits<double>::max();
			double max_v = -std::numeric_limits<double>::max();

			for (auto h : fhs) {
				const auto& uv = h->uv();
				if (uv[0] < min_u) min_u = uv[0];
				if (uv[0] > max_u) max_u = uv[0];
				if (uv[1] < min_v) min_v = uv[1];
				if (uv[1] > max_v) max_v = uv[1];
			}
			// 增加微小容差
			double eps = 1e-5;
			min_u -= eps; max_u += eps;
			min_v -= eps; max_v += eps;

			// 清空去重记录
			visitedVertexIds.clear();

			// --- D. 遍历候选三角形 ---
			for (auto tf : candidateTris)
			{
				// 遍历三角形的三个顶点 (通过半边迭代)
				// 假设 It::FHEIterator 是某种标准迭代器，如果太慢可以直接写 next->next
				H* triHE = triMesh->faceHalfedge(tf);
				H* currHE = triHE;

				// 手动展开循环，微小优化
				do {
					V* tfhTarget = triMesh->halfedgeTarget(currHE);
					int vId = tfhTarget->id();

					// 1. 局部去重：如果这个点已经在当前四边形计算过了，跳过
					if (visitedVertexIds.count(vId)) {
						currHE = triMesh->halfedgeNext(currHE);
						continue;
					}

					const auto& ptUV = currHE->uv();

					// 2. AABB 快速排斥 
					if (ptUV[0] < min_u || ptUV[0] > max_u || ptUV[1] < min_v || ptUV[1] > max_v) {
						visitedVertexIds.insert(vId); // 标记已访问（虽然不在范围内）
						currHE = triMesh->halfedgeNext(currHE);
						continue;
					}

					// 3. 精确几何检查 (Signed Area)
					bool markSample = true;
					for (int k = 0; k < 4; k++)
					{
						if (!signedAreaNotNegative(fhs[k]->uv(), fhs[(k + 1) % 4]->uv(), ptUV))
						{
							markSample = false;
							break;
						}
					}

					if (markSample)
					{
						visitedVertexIds.insert(vId); // 标记为已处理

						QB_sampling st;
						if (tfhTarget->boundary()) st.boundary() = true;

						st.fId() = f->id();
						st.patchId() = st.fId();
						st.pos() = tfhTarget->point();
						st.uv() = ptUV;
						st.uv_init() = ptUV;

						// --- 权重计算逻辑 (保持您原有的数学逻辑) ---
						// 计算各三角形面积用于重心坐标插值
						double area_tri_124 = signedArea(fhs[0]->uv(), fhs[1]->uv(), fhs[3]->uv());
						double area_tri_12s = signedArea(fhs[0]->uv(), fhs[1]->uv(), ptUV);

						// 边界处理 1-2
						if (std::abs(area_tri_12s) < diff_area)
						{
							double len_edge = (fhs[1]->uv() - fhs[0]->uv()).norm();
							double e_w = (len_edge > 1e-9) ? (ptUV - fhs[0]->uv()).norm() / len_edge : 0.0;
							st.ws() = { 1 - e_w, e_w, 0.0, 0.0 };
							st.vId() = fhs[0]->target()->id();
							thread_sts[tid].push_back(st);

							currHE = triMesh->halfedgeNext(currHE);
							continue;
						}

						double w3 = area_tri_12s / area_tri_124;
						double area_tri_24s = signedArea(fhs[1]->uv(), fhs[3]->uv(), ptUV);
						double w0 = area_tri_24s / area_tri_124;
						double area_tri_41s = signedArea(fhs[3]->uv(), fhs[0]->uv(), ptUV);

						// 边界处理 4-1
						if (std::abs(area_tri_41s) < diff_area)
						{
							double len_edge = (fhs[3]->uv() - fhs[0]->uv()).norm();
							double e_w = (len_edge > 1e-9) ? (ptUV - fhs[3]->uv()).norm() / len_edge : 0.0;
							st.ws() = { e_w, 0.0, 0.0, 1 - e_w };
							st.vId() = fhs[0]->target()->id();
							thread_sts[tid].push_back(st);

							currHE = triMesh->halfedgeNext(currHE);
							continue;
						}

						double w1 = area_tri_41s / area_tri_124;
						double area_tri_23s = signedArea(fhs[1]->uv(), fhs[2]->uv(), ptUV);

						// 边界处理 2-3
						if (std::abs(area_tri_23s) < diff_area)
						{
							double len_edge = (fhs[2]->uv() - fhs[1]->uv()).norm();
							double e_w = (len_edge > 1e-9) ? (ptUV - fhs[1]->uv()).norm() / len_edge : 0.0;
							st.ws() = { 0.0, 1 - e_w, e_w, 0.0 };
							st.vId() = fhs[0]->target()->id();
							thread_sts[tid].push_back(st);

							currHE = triMesh->halfedgeNext(currHE);
							continue;
						}

						double area_tri_34s = signedArea(fhs[2]->uv(), fhs[3]->uv(), ptUV);

						// 边界处理 3-4
						if (std::abs(area_tri_34s) < diff_area)
						{
							double len_edge = (fhs[3]->uv() - fhs[2]->uv()).norm();
							double e_w = (len_edge > 1e-9) ? (ptUV - fhs[2]->uv()).norm() / len_edge : 0.0;
							st.ws() = { 0.0, 0.0, 1 - e_w, e_w };
							st.vId() = fhs[0]->target()->id();
							thread_sts[tid].push_back(st);

							currHE = triMesh->halfedgeNext(currHE);
							continue;
						}

						// 内部点
						double w2 = 1.0 - w0 - w1 - w3;
						st.ws() = { w0, w1, w2, w3 };
						st.vId() = fhs[0]->target()->id();

						thread_sts[tid].push_back(st);
					}

					currHE = triMesh->halfedgeNext(currHE);

				} while (currHE != triHE);
			} // End Candidate Triangles Loop
		} // End Quad Loop
	} // End Parallel Region

	// --------------------------------------------------------------------------
	// 4. 合并所有线程的结果
	// --------------------------------------------------------------------------
	size_t total_count = 0;
	for (const auto& buf : thread_sts) total_count += buf.size();

	sts.reserve(sts.size() + total_count);

	for (const auto& buf : thread_sts) {
		sts.insert(sts.end(), buf.begin(), buf.end());
	}
}

void CCG_QMSLib::updateMeshPosByControlMesh(M* pMesh, CCG_QMS_model& model)
{
	for (auto v : It::MVIterator(pMesh))
	{
		v->point() = model.controlPoints()[v->id() - 1];
	}
}

std::vector<CCG_QMSLib::F*> CCG_QMSLib::FindQuadsInsideTris_Brute_WithBuckets(M* quadMesh, M* triMesh, double eps)
{
	// ---------------------------
	// 1) (patch, sub) 作为桶 key
	// ---------------------------
	struct PatchSubKey
	{
		int patch = -1;
		int sub = -1;

		bool operator==(const PatchSubKey& rhs) const noexcept
		{
			return patch == rhs.patch && sub == rhs.sub;
		}
	};

	struct PatchSubKeyHash
	{
		std::size_t operator()(const PatchSubKey& k) const noexcept
		{
			std::uint64_t a = static_cast<std::uint32_t>(k.patch);
			std::uint64_t b = static_cast<std::uint32_t>(k.sub);
			std::uint64_t x = (a << 32) ^ b;
			return std::hash<std::uint64_t>{}(x);
		}
	};

	using FaceBucket = std::unordered_map<PatchSubKey, std::vector<F*>, PatchSubKeyHash>;

	// ---------------------------
	// 2) 分桶函数
	// ---------------------------
	auto bucketFacesByPatchSub = [](M* mesh) -> FaceBucket
		{
			FaceBucket buckets;
			buckets.reserve(mesh->numFaces() / 4 + 8);

			for (auto f : It::MFIterator(mesh)) {
				if (!f) continue;
				PatchSubKey key;
				key.patch = f->patchIndex();
				key.sub = f->subPatchIndex();
				buckets[key].push_back(f);
			}
			return buckets;
		};

	auto printBucketStats = [](const char* name, const FaceBucket& buckets)
		{
			//std::cout << "===== " << name << " buckets =====\n";
			//std::cout << "bucket count = " << buckets.size() << "\n";
			int printed = 0;
			for (const auto& kv : buckets) {
				const auto& key = kv.first;
				const auto& vec = kv.second;
				/*std::cout << "  (patch=" << key.patch << ", sub=" << key.sub
					<< ") faces=" << vec.size() << "\n";*/
				if (++printed >= 10) break;
			}
			//std::cout << "===========================\n";
		};

	// ---------------------------
	// 3) 取面 UV（用 halfedge uv）
	// ---------------------------
	auto getFaceUVs_Quad = [](M* mesh, F* f, std::array<CPoint2, 4>& uvOut) -> bool
		{
			int k = 0;
			H* heStart = mesh->faceHalfedge(f);
			if (!heStart) return false;

			H* he = heStart;
			do {
				if (k >= 4) return false;     // 不是四边形
				uvOut[k++] = he->uv();        // halfedge 上存的 uv
				he = mesh->halfedgeNext(he);
			} while (he != heStart);

			return (k == 4);
		};

	auto getFaceUVs_Tri = [](M* mesh, F* f, std::array<CPoint2, 3>& uvOut) -> bool
		{
			int k = 0;
			H* heStart = mesh->faceHalfedge(f);
			if (!heStart) return false;

			H* he = heStart;
			do {
				if (k >= 3) return false;     // 不是三角形
				uvOut[k++] = he->uv();
				he = mesh->halfedgeNext(he);
			} while (he != heStart);

			return (k == 3);
		};

	// ---------------------------
	// 4) 点在三角形内（2D）
	// ---------------------------
	auto cross2 = [](const CPoint2& a, const CPoint2& b, const CPoint2& c) -> double
		{
			return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
		};

	auto pointInTri2D = [&](const CPoint2& p,
		const CPoint2& a,
		const CPoint2& b,
		const CPoint2& c) -> bool
		{
			double c1 = cross2(a, b, p);
			double c2 = cross2(b, c, p);
			double c3 = cross2(c, a, p);

			bool hasNeg = (c1 < -eps) || (c2 < -eps) || (c3 < -eps);
			bool hasPos = (c1 > eps) || (c2 > eps) || (c3 > eps);

			// 允许在边上：>= -eps
			return !(hasNeg && hasPos);
		};

	// ---------------------------
	// 5) brute：某个 quad 是否落在该桶的任意 tri 内
	// ---------------------------
	auto quadInsideAnyTri_Brute = [&](F* qf, const std::vector<F*>& triFacesInSameBucket) -> bool
		{
			std::array<CPoint2, 4> qUV;
			if (!getFaceUVs_Quad(quadMesh, qf, qUV)) return false;

			for (F* tf : triFacesInSameBucket) {
				std::array<CPoint2, 3> tUV;
				if (!getFaceUVs_Tri(triMesh, tf, tUV)) continue;

				bool ok =
					pointInTri2D(qUV[0], tUV[0], tUV[1], tUV[2]) &&
					pointInTri2D(qUV[1], tUV[0], tUV[1], tUV[2]) &&
					pointInTri2D(qUV[2], tUV[0], tUV[1], tUV[2]) &&
					pointInTri2D(qUV[3], tUV[0], tUV[1], tUV[2]);

				if (ok) return true;
			}
			return false;
		};

	// ---------------------------
	// 6) 分桶 + 执行 brute
	// ---------------------------
	FaceBucket quadBuckets = bucketFacesByPatchSub(quadMesh);
	FaceBucket triBuckets = bucketFacesByPatchSub(triMesh);

	printBucketStats("QuadMesh", quadBuckets);
	printBucketStats("TriMesh", triBuckets);

	std::vector<F*> hitQuads;
	hitQuads.reserve(quadMesh->numFaces() / 2);

	for (const auto& kv : quadBuckets) {
		const PatchSubKey& key = kv.first;
		const std::vector<F*>& qFaces = kv.second;

		auto itT = triBuckets.find(key);
		if (itT == triBuckets.end()) continue;

		const std::vector<F*>& tFaces = itT->second;

		for (F* qf : qFaces) {
			if (!qf) continue;
			if (quadInsideAnyTri_Brute(qf, tFaces)) {
				hitQuads.push_back(qf);
			}
		}
	}

	return hitQuads;
}
void CCG_QMSLib::calculateScaledJacobian_quadMesh(M* pMesh)
{
	for (auto f : It::MFIterator(pMesh))
	{
		// Store vertices in an array for easy indexing in the loop
		std::vector<CPoint> f_pts;
		for(auto fv:It::FVIterator(pMesh,f))
		{
			f_pts.push_back(fv->point());
		}
		double min_sj = 1.0;

		for (int i = 0; i < 4; ++i) {
			// Get the current vertex and its two adjacent vertices
			const CPoint& curr = f_pts[i];
			const CPoint& prev = f_pts[(i + 3) % 4]; // Previous point in the loop
			const CPoint& next = f_pts[(i + 1) % 4]; // Next point in the loop

			// Define two edge vectors originating from the current vertex
			CPoint L1 = next - curr;
			CPoint L2 = prev - curr;

			// Calculate the magnitude (length) of the edge vectors
			double l1_mag = L1.norm();
			double l2_mag = L2.norm();

			// Prevent division by zero for degenerate elements (edges with zero length)
			if (l1_mag < 1e-12 || l2_mag < 1e-12) {
				min_sj = 0.0;
			}

			// Calculate the cross product. 
			// In 3D space, the magnitude of the cross product represents the area of the parallelogram.
			CPoint cross_prod = L1 ^ L2;
			double area = cross_prod.norm();

			// Normalization (Scaled Jacobian at vertex i).
			// Geometric interpretation: This is equivalent to sin(theta), 
			// where theta is the angle between the two edges.
			double sj_i = area / (l1_mag * l2_mag);

			// Update the minimum value for the element
			if (sj_i < min_sj) {
				min_sj = sj_i;
			}
		}
		f->quality_scaledJocbian() = min_sj;
	}
}
void CCG_QMSLib::calculateEdgeRatio_quadMesh(M* pMesh)
{
	for (auto f : It::MFIterator(pMesh))
	{
		// Store edge lengths in an array for easy indexing in the loop
		std::vector<double> edge_lengths;
		for (auto fe : It::FEIterator(pMesh, f))
		{
			V* v1 = pMesh->edgeVertex1(fe);
			V* v2 = pMesh->edgeVertex2(fe);
			double length = (v1->point() - v2->point()).norm();
			edge_lengths.push_back(length);
		}
		// Calculate edge ratio
		double min_length = *std::min_element(edge_lengths.begin(), edge_lengths.end());
		double max_length = *std::max_element(edge_lengths.begin(), edge_lengths.end());
		// Prevent division by zero for degenerate elements (edges with zero length)
		if (max_length < 1e-12) {
			f->quality_edgeRatio() = 0.0;
		}
		else {
			double edge_ratio = max_length / min_length;
			f->quality_edgeRatio() = edge_ratio;
		}
	}
}
/*-----------------------------------------------------------------------------------------------*/
