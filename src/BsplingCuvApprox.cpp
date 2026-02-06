#include "BsplingCuvApprox.h"
#include <iostream>
#include <Eigen/Dense>
#include<Eigen/Sparse>
#include <Eigen/SVD>
#include <cmath> 
#include <Eigen/SparseLU>
#define _USE_MATH_DEFINES  //Regarding the MSVC compiler
#define M_PI 3.14159265358979323846
#include <omp.h>  //Enable OpenMP support


CCG_QMSLib::CFeatureLoop::CFeatureLoop(M* pMesh)
{
	m_pMesh = pMesh;
	int ehBoundrynum = 0;
	int ehSharpNum = 0;
	/* Mark boundary edges as sharp */
	for (auto e : It::MEIterator(pMesh))
	{
		if (e->boundary())
		{
			e->sharp() = true;
			ehBoundrynum++;
			//std::cout << "1111111111111111111" << std::endl;
		}
		else
		{
			e->sharp() = false;
			//std::cout << "2222222222222222222" << std::endl;
		}
	}
	//for (auto e : It::MEIterator(pMesh))
	//{
	//    if (e->sharp())
	//    {
	//        ehSharpNum++;
	//        //std::cout << "3333333333333333333333333" << std::endl;
	//    }
	//}
	/*std::cout << "Number of boundary halfedges (ehBoundrynum): " << ehBoundrynum << std::endl;
	std::cout << "Number of feature edges (ehSharpNum): " << ehSharpNum << std::endl;*/
	// Collect all feature edges
	std::set<E*> feature_hes;
	for (auto eiter : It::MEIterator(pMesh))
	{
		E* e = eiter;
		// Check if this is a feature edge
		if (!e->sharp()) continue;
		// Get the halfedge corresponding to the feature edge
		//H* he = m_pMesh->edgeHalfedge(e, 0);
		feature_hes.insert(e);
	}
	//std::cout << "Number of feature edges (feature_hes.size): " << feature_hes.size() << std::endl;

	// 1. First process feature arcs
	while (!feature_hes.empty()) {
		bool foundFeaturePoint = false;
		std::set<E*>::iterator eSiter = feature_hes.begin();

		// Find an edge whose source vertex is a feature point
		for (; eSiter != feature_hes.end(); ++eSiter) {
			E* e = *eSiter;
			H* he = m_pMesh->edgeHalfedge(e, 0);
			V* sourceV = m_pMesh->halfedgeSource(he); // Get source vertex
			V* targetV = m_pMesh->halfedgeTarget(he);
			// Determine if the source vertex is a feature point
			int sourceVEdgeSharpCount = 0;
			int targetVEdgeSharpCount = 0;

			for (auto ve : It::VCcwEIterator(m_pMesh, sourceV)) {
				if (ve->sharp()) sourceVEdgeSharpCount++;
			}
			for (auto ve : It::VCcwEIterator(m_pMesh, targetV)) {
				if (ve->sharp()) targetVEdgeSharpCount++;
			}
			/*for (auto vh : It::VCcwOutHEIterator(m_pMesh, sourceV)) {
				if (m_pMesh->halfedgeEdge(vh)->sharp()) sourceVEdgeSharpCount++;
			}*/

			/*for (auto vh : It::VCcwOutHEIterator(m_pMesh, targetV)) {
				if (m_pMesh->halfedgeEdge(vh)->sharp()) targetVEdgeSharpCount++;
			}*/
			if (sourceVEdgeSharpCount != 2 || targetVEdgeSharpCount != 2) { // Source vertex is a feature point
				foundFeaturePoint = true;

				// Ensure we use the halfedge direction where the source vertex is the feature point
				/*H* startHe = he;*/
				H* startHe = (sourceVEdgeSharpCount != 2) ? he : m_pMesh->halfedgeSym(he);
				if (startHe == NULL) continue;
				CLoopSegment* pArc = new CLoopSegment(m_pMesh, startHe);
				assert(pArc);
				m_arcs.push_back(pArc);
				// Remove processed edges from the container
				for (auto he : pArc->halfedges()) {
					E* e = m_pMesh->halfedgeEdge(he);
					eSiter = feature_hes.find(e);
					if (eSiter == feature_hes.end()) continue;
					feature_hes.erase(eSiter);
				}
				break;
			}
		}

		if (!foundFeaturePoint) break; // No feature points left; remaining edges form loops
	}

	// 2. Process feature loops
	while (!feature_hes.empty()) {
		// Get the halfedge of the first feature edge
		std::set<E*>::iterator eSiter = feature_hes.begin();
		H* he = m_pMesh->edgeHalfedge(*eSiter, 0);
		//std::cout << "11111" << std::endl;
		// Trace along this feature halfedge to form a loop
		CLoop* pL = new CLoop(m_pMesh, he);
		assert(pL);
		m_loops.push_back(pL);

		// Remove all edges belonging to this loop from feature_hes
		for (auto hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter++) {
			H* he = *hiter;
			E* e = m_pMesh->halfedgeEdge(he);
			eSiter = feature_hes.find(e);
			if (eSiter == feature_hes.end()) continue;
			feature_hes.erase(eSiter);
		}
	}
}

CCG_QMSLib::CFeatureLoop::~CFeatureLoop()
{
}

CCG_QMSLib::CLoopSegment::CLoopSegment(M* pMesh, H* pH)
{
	static int TriArcCounter = 0;
	static int QuadArcCounter = 0;
	// Count the number of vertices on the face
	F* f = *pMesh->faces().begin();
	int vertexCount = 0;
	H* he1 = pMesh->faceHalfedge(f);
	while (he1 != nullptr) {
		vertexCount++;
		he1 = pMesh->halfedgeNext(he1);
		if (he1 == pMesh->faceHalfedge(f)) break;
	}
	if (vertexCount == 3) triArc_id = TriArcCounter++;
	else if (vertexCount == 4) quadArc_id = QuadArcCounter++;
	else
	{
		//std::cout << "This is a polygon with " << vertexCount << " vertices." << std::endl;
		//std::cout << "--------------------  " << vertexCount << std::endl;
	}
	m_pMesh = pMesh;
	m_halfedges.clear();
	auto isFeaturePoint = [&](V* v) {
		int sharpCount = 0;
		for (auto ve : It::VCcwEIterator(pMesh, v)) {
			if (ve->sharp()) sharpCount++;
		}
		return sharpCount != 2; // If the number of sharp edges at the vertex is not equal to 2, it's a feature point
		};

	H* currentHe = pH; // Starting halfedge at a feature point
	m_halfedges.push_back(currentHe);
	// Vertex tracing logic for arcs
	H* he = pH;                      // The starting halfedge
	m_vertexLength = 0;              // Reset the arc length
	m_vertexs.clear();               // Clear the vertex list
	std::set<V*> visitedVertices;    // Record visited vertices

	// Add starting vertex (the arc must have a clearly defined starting point)
	V* startV = m_pMesh->halfedgeSource(he);
	m_vertexs.push_back(startV);
	visitedVertices.insert(startV);
	V* endV = m_pMesh->halfedgeTarget(he);
	m_vertexs.push_back(endV);
	visitedVertices.insert(endV);
	/*std::cout << "Tracing arc through vertex " << startV->id()
		<< " at (" << startV->point()[0] << "," << startV->point()[1] << "," << startV->point()[3] << ")" << std::endl;
	std::cout << "Tracing arc through vertex " << endV->id()
		<< " at (" << endV->point()[0] << "," << endV->point()[1] << "," << endV->point()[3] << ")" << std::endl;*/
	while (true) {
		// 1. Search for the next feature halfedge
		bool foundNext = false;
		H* nextHe = nullptr;
		E* currentEdge = pMesh->halfedgeEdge(currentHe);
		V* currentV = pMesh->halfedgeTarget(currentHe);
		// Termination condition 1: encounter duplicate vertices (prevent loops)
		for (auto vh : It::VCcwOutHEIterator(pMesh, currentV)) {
			E* edge = pMesh->halfedgeEdge(vh);
			// Skip the opposite halfedge
			if (edge == currentEdge) continue;
			// Stop at the first unvisited feature edge
			if (edge->sharp() && edge->visit() == false) {
				nextHe = vh;
				foundNext = true;
				edge->visit() = true;
				V* v = pMesh->halfedgeTarget(vh);
				if (visitedVertices.find(v) == visitedVertices.end()) {
					// Record the new vertex
					m_vertexs.push_back(v);
					visitedVertices.insert(v);
					// Accumulate arc length (using edge length)
					m_vertexLength += m_pMesh->edgeLength(m_pMesh->halfedgeEdge(he));
					// Print path during tracing
					/*std::cout << "Tracing arc through vertex " << v->id()
						<< " at (" << v->point()[0] << "," << v->point()[1] << "," << v->point()[3] << ")" << std::endl;*/
				}
				// Termination condition 2: safety limit (prevent infinite loops)
				if (visitedVertices.size() >= 10000) {
					//std::cerr << "Warning: Arc vertex tracing reached safety limit (10000 vertices)" << std::endl;
					break;
				}
				break;
			}
		}
		// 2. Handle the search result
		if (!foundNext) {
			//std::cout << "No next feature edge found, arc segment ends" << std::endl;
			break;
		}

		// 3. Add the confirmed feature halfedge to the segment
		m_halfedges.push_back(nextHe);

		// 4. Update the current halfedge
		currentHe = nextHe;

		// 5. Check if the new target vertex is a feature point
		V* nextV = pMesh->halfedgeTarget(currentHe);
		if (isFeaturePoint(nextV)) {
			//std::cout << "Reached feature point, arc segment complete" << std::endl;
			break;
		}
	}

	// Validate arc segment (must contain at least 1 halfedge)
	if (m_halfedges.empty()) {
		//std::cerr << "Error: No valid feature arc line!" << std::endl;
		return;
	}

	// Ensure the starting halfedge is included at the front if missing
	if (!m_halfedges.empty() && m_halfedges.front() != pH) {
		m_halfedges.insert(m_halfedges.begin(), pH);
	}
	/*std::cout << "Successfully created feature arc containing " << m_halfedges.size()
		<< " halfedges (start: " << pMesh->halfedgeSource(m_halfedges.front())->id()
		<< ", end: " << pMesh->halfedgeTarget(m_halfedges.back())->id() << ")" << std::endl;*/
}

inline CCG_QMSLib::CLoop::CLoop(M* pMesh, H* pH) {
	static int TriLoopCounter = 0;
	static int QuadLoopCounter = 0;

	// Count the number of vertices of the face
	F* f = *pMesh->faces().begin();
	int vertexCount = 0;
	H* he1 = pMesh->faceHalfedge(f);
	while (he1 != nullptr) {
		vertexCount++;
		he1 = pMesh->halfedgeNext(he1);
		if (he1 == pMesh->faceHalfedge(f)) break;
	}
	if (vertexCount == 3) triLoop_id = TriLoopCounter++;
	else if (vertexCount == 4) quadLoop_id = QuadLoopCounter++;
	else {
		//std::cout << "This is a polygon with " << vertexCount << " vertices." << std::endl;
	}

	m_pMesh = pMesh;
	H* pHalfedge = pH;
	m_length = 0;
	H* he = pH;

	// Fix feature edge tracing logic
	std::set<E*> visitedEdges; // Record visited edges
	do {
		V* v = m_pMesh->halfedgeTarget(he);
		bool foundSharpEdge = false;

		// Iterate over all outgoing half-edges of the current vertex in CCW order
		for (auto vh : It::VCcwOutHEIterator(pMesh, v)) {
			E* edge = m_pMesh->halfedgeEdge(vh);
			if (edge->sharp() && visitedEdges.find(edge) == visitedEdges.end()) {
				m_halfedges.push_back(vh);
				visitedEdges.insert(edge);
				m_length += m_pMesh->edgeLength(edge);
				he = vh; // Move to the next half-edge
				foundSharpEdge = true;
				break;
			}
		}
		// Check if we have returned to the starting half-edge
		/*if (he == pHalfedge) {
			std::cout << "The feature line forms a closed loop." << std::endl;
			break;
		}*/
		if (!foundSharpEdge) {
			//std::cout << "The feature line forms an open arc." << std::endl;
			break;
		}
		if (he->target()->id() == pHalfedge->target()->id()) {
			//std::cout << "The feature line forms a closed loop." << std::endl;
			break;
		}

		// Prevent infinite loop: force exit if exceeding expected number of edges
		if (m_halfedges.size() > pMesh->numEdges()) {
			//std::cerr << "Warning: Possible infinite loop in feature line tracing." << std::endl;
			break;
		}

	} while (true);

	// Vertex tracing logic
	he = pH;
	m_vertexLength = 0;
	m_vertexs.clear();
	std::set<V*> visitedVertices; // Record visited vertices

	do {
		V* v = m_pMesh->halfedgeTarget(he);
		if (visitedVertices.find(v) != visitedVertices.end()) break; // Avoid duplicate vertices
		m_vertexs.push_back(v);
		visitedVertices.insert(v);
		m_vertexLength += m_pMesh->edgeLength(m_pMesh->halfedgeEdge(he));
		he = m_pMesh->vertexMostClwOutHalfEdge(v);
	} while (he != pHalfedge && visitedVertices.size() < (pMesh->numVertices() + 1)); // Safety limit
}

inline CCG_QMSLib::CLoop::~CLoop()
{
	m_halfedges.clear();
	// Clear vertex list
	m_vertexs.clear();
	for (size_t i = 0; i < m_segments.size(); i++)
	{
		delete m_segments[i];
	}
	m_segments.clear();
}

void CCG_QMSLib::BE_curve_model::normalize_triQuadMesh(M* triMesh, M* quadMesh)
{
	CPoint s(0, 0, 0);
	for (auto v : It::MVIterator(triMesh))
	{
		s = s + v->point();
	}
	for (auto v : It::MVIterator(quadMesh))
	{
		s = s + v->point();
	}
	s = s / (triMesh->numVertices() + quadMesh->numVertices());

	for (auto v : It::MVIterator(triMesh))
	{
		CPoint p = v->point();
		p = p - s;
		v->point() = p;
	}

	for (auto v : It::MVIterator(quadMesh))
	{
		CPoint p = v->point();
		p = p - s;
		v->point() = p;
	}

	double d = 0;
	for (auto v : It::MVIterator(triMesh))
	{
		CPoint p = v->point();
		for (int k = 0; k < 3; k++)
		{
			d = (d > fabs(p[k])) ? d : fabs(p[k]);
		}
	}
	for (auto v : It::MVIterator(quadMesh))
	{
		CPoint p = v->point();
		for (int k = 0; k < 3; k++)
		{
			d = (d > fabs(p[k])) ? d : fabs(p[k]);
		}
	}

	for (auto v : It::MVIterator(triMesh))
	{
		CPoint p = v->point();
		p = p / d;
		v->point() = p;
	}
	for (auto v : It::MVIterator(quadMesh))
	{
		CPoint p = v->point();
		p = p / d;
		v->point() = p;
	}
}


void CCG_QMSLib::BE_curve_model::setControlPoints(M* pMesh)
{
	// Clear existing data
	for (auto curve : be_bcs) delete curve;
	be_bcs.clear();
	sharpFeaturePoints.clear();
	//originalSharpFeaturePoints.clear();

	// Iterate through all quadrilateral feature loops
	for (auto loop : be_quadFeatureLoops->loops()) {
		if (!loop) continue;

		// Iterate through each half-edge in the loop
		for (auto pH : loop->halfedges()) {
			// Get adjacent half-edges
			H* pHPrev = nullptr;
			H* pHNext = nullptr;

			// Get the previous half-edge
			V* vSource = pMesh->halfedgeSource(pH);
			for (auto vh : It::VClwInHEIterator(pMesh, vSource)) {
				if (std::find(loop->halfedges().begin(), loop->halfedges().end(), vh) != loop->halfedges().end()) {
					pHPrev = vh;
					break;
				}
			}

			// Get the next half-edge
			V* vTarget = pMesh->halfedgeTarget(pH);
			for (auto vh : It::VCcwOutHEIterator(pMesh, vTarget)) {
				if (std::find(loop->halfedges().begin(), loop->halfedges().end(), vh) != loop->halfedges().end()) {
					pHNext = vh;
				}
			}
			if (!pHPrev)
			{
				//std::cout << "Previous half-edge not found" << std::endl;
			}
			if (!pHNext)
			{
				//std::cout << "Next half-edge not found" << std::endl;
			}
			// Retrieve the four key points
			CPoint V0 = pHPrev ? pHPrev->source()->point() : pH->source()->point();
			CPoint V1 = pH->source()->point();
			CPoint V2 = pH->target()->point();
			CPoint V3 = pHNext ? pHNext->target()->point() : pH->target()->point();

			// Compute intermediate control points
			CPoint V4 = V0 * (1.0 / 3) + V1 * (2.0 / 3);
			CPoint V5 = V1 * (2.0 / 3) + V2 * (1.0 / 3);
			CPoint V6 = V1 * (1.0 / 3) + V2 * (2.0 / 3);
			CPoint V7 = V2 * (2.0 / 3) + V3 * (1.0 / 3);

			// Calculate quadrilateral mesh angles
			double quadAngle1 = calculateAngle(V0, V1, V2);
			double quadAngle2 = calculateAngle(V1, V2, V3);

			// Check corresponding triangle mesh points
			bool isTriCorner1 = false;
			bool isTriCorner2 = false;

			// Find the triangle mesh point corresponding to the current quad half-edge
			for (const auto& kv : triPointToQuadHalfedgeMap) {
				V* triPoint = kv.first;
				H* quadHalfedge = kv.second;
				if (quadHalfedge != pH) continue;

				isTriCorner1 = isTriCorner1 || checkTriCorner(V1, triPoint);
				isTriCorner2 = isTriCorner2 || checkTriCorner(V2, triPoint);
			}

			// Create Bezier curve
			BE_bezierCurve_approxi* be_curve = new BE_bezierCurve_approxi;

			// Set control points
			be_curve->cpts()[0] = ((quadAngle1 < 130.0 || quadAngle1 > 230.0) && isTriCorner1)
				? V1 : (V4 * 0.5 + V5 * 0.5);
			be_curve->cpts()[3] = ((quadAngle2 < 130.0 || quadAngle2 > 230.0) && isTriCorner2)
				? V2 : (V6 * 0.5 + V7 * 0.5);
			be_curve->cpts()[1] = V5;
			be_curve->cpts()[2] = V6;

			// Record feature points
			if ((quadAngle1 < 130.0 || quadAngle1 > 230.0) && isTriCorner1) {
				pMesh->halfedgeSource(pH)->feature() = true;
				sharpFeaturePoints[pH->source()->id()] = V1;
			}
			if ((quadAngle2 < 130.0 || quadAngle2 > 230.0) && isTriCorner2) {
				pMesh->halfedgeTarget(pH)->feature() = true;
				sharpFeaturePoints[pH->target()->id()] = V2;
				//std::cout << "V2:" << V2[0] << " " << V2[1] << " " << V2[2] << std::endl;
				//std::cout << "sharpFeaturePoints[pH->target()->id()]" << sharpFeaturePoints[pH->target()->id()][0]
				//          << " " << sharpFeaturePoints[pH->target()->id()][1]
				//          << " " << sharpFeaturePoints[pH->target()->id()][2] << std::endl;
			}

			// Store the curve
			be_bcs.push_back(be_curve);
			quadHalfedgeToBezierCurveMap[pH] = be_curve;
		}
	}
	// Iterate through all quadrilateral feature arcs
	for (auto arc : be_quadFeatureLoops->arcs()) {
		// Ensure feature arc is valid
		if (!arc || arc->halfedges().empty()) continue;

		// Iterate over all half-edges in the arc
		for (auto pH : arc->halfedges()) {
			// Get adjacent half-edges
			H* pHPrev = nullptr;
			H* pHNext = nullptr;

			// Get the previous half-edge
			V* vSource = pMesh->halfedgeSource(pH);
			for (auto vh : It::VClwInHEIterator(pMesh, vSource)) {
				if (std::find(arc->halfedges().begin(), arc->halfedges().end(), vh) != arc->halfedges().end()) {
					pHPrev = vh;
					break;
				}
			}

			// Get the next half-edge
			V* vTarget = pMesh->halfedgeTarget(pH);
			for (auto vh : It::VCcwOutHEIterator(pMesh, vTarget)) {
				if (std::find(arc->halfedges().begin(), arc->halfedges().end(), vh) != arc->halfedges().end()) {
					pHNext = vh;
				}
			}

			// Retrieve the four key points
			CPoint V0 = pHPrev ? pHPrev->source()->point() : pH->source()->point();
			CPoint V1 = pH->source()->point();
			CPoint V2 = pH->target()->point();
			CPoint V3 = pHNext ? pHNext->target()->point() : pH->target()->point();

			// Compute intermediate control points
			CPoint V4 = V0 * (1.0 / 3) + V1 * (2.0 / 3);
			CPoint V5 = V1 * (2.0 / 3) + V2 * (1.0 / 3);
			CPoint V6 = V1 * (1.0 / 3) + V2 * (2.0 / 3);
			CPoint V7 = V2 * (2.0 / 3) + V3 * (1.0 / 3);

			// Create Bezier curve
			BE_bezierCurve_approxi* be_curve = new BE_bezierCurve_approxi;

			// Set intermediate control points
			be_curve->cpts()[1] = V5;
			be_curve->cpts()[2] = V6;
			// Set control points
			be_curve->cpts()[0] = (V0 == V1) ? V1 : (V4 * 0.5 + V5 * 0.5);
			be_curve->cpts()[3] = (V2 == V3) ? V2 : (V6 * 0.5 + V7 * 0.5);
			// Store the curve
			be_bcs.push_back(be_curve);
			quadHalfedgeToBezierCurveMap[pH] = be_curve;
		}
	}
}

void CCG_QMSLib::BE_curve_model::markFeaturePts_input(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->feature() && v->boundary())
		{
			sharpFeaturePoints[v->id()] = v->point();
		}
	}
}

double CCG_QMSLib::BE_curve_model::calculateAngle(const CPoint& a, const CPoint& b, const CPoint& c)
{
	CPoint v1 = a - b;
	CPoint v2 = c - b;
	double cosAngle = (v1 * v2) / (v1.norm() * v2.norm());
	cosAngle = std::max(-1.0, std::min(1.0, cosAngle)); // Ensure value is within valid range
	return acos(cosAngle) * 180.0 / M_PI;               // Convert to degrees
}

bool CCG_QMSLib::BE_curve_model::checkTriCorner(const CPoint& quadPoint, V* triPoint) const
{
	if (!triPoint) return false;

	// Get adjacent half-edges of the triangle point (feature edges only)
	std::vector<H*> triHalfedges;
	M* triMesh = be_triFeatureLoops->getMesh();

	// Collect incoming boundary half-edges
	for (auto vh : It::VClwInHEIterator(triMesh, triPoint)) {
		if (vh->edge()->boundary()) {
			triHalfedges.push_back(vh);
		}
	}

	// Collect outgoing boundary half-edges
	for (auto vh : It::VCcwOutHEIterator(triMesh, triPoint)) {
		if (vh->edge()->boundary()) {
			triHalfedges.push_back(vh);
		}
	}

	// Compute angle (requires two feature edges)
	if (triHalfedges.size() == 2) {
		CPoint p0 = triHalfedges[0]->source()->point();
		CPoint p1 = triPoint->point();
		CPoint p2 = triHalfedges[1]->target()->point();

		double triAngle = calculateAngle(p0, p1, p2);
		return (triAngle < 130.0 || triAngle > 230.0);
	}
	return false;
}

void CCG_QMSLib::BE_curve_model::resetControlPoints(M* pMesh)
{
	// Clear existing bcs data
	for (auto curve : be_bcs) {
		delete curve;
	}
	be_bcs.clear();

	// Iterate over feature loops of the quadrilateral mesh
	for (auto loop : be_quadFeatureLoops->loops()) {
		if (!loop) continue;

		// Iterate over all half-edges in the loop
		for (auto pH : loop->halfedges())
		{
			// Retrieve start and end vertices of the half-edge
			CPoint V0, V1, V2, V3, V4, V5, V6, V7;
			H* pHNh = NULL;   // Next half-edge
			H* pHPrev = NULL; // Previous half-edge

			// Get source vertex of the current half-edge
			V* vSource = pMesh->halfedgeSource(pH);
			for (auto vh : It::VClwInHEIterator(pMesh, vSource))
			{
				if (std::find(loop->halfedges().begin(), loop->halfedges().end(), vh) != loop->halfedges().end()) {
					pHPrev = vh;
					break;
				}
			}

			// Get target vertex of the current half-edge
			V* vTarget = pMesh->halfedgeTarget(pH);
			for (auto vh : It::VCcwOutHEIterator(pMesh, vTarget))
			{
				if (std::find(loop->halfedges().begin(), loop->halfedges().end(), vh) != loop->halfedges().end()) {
					pHNh = vh;
				}
			}

			// Get coordinates of adjacent vertices
			V0 = pHPrev->source()->point();   // Start point of previous half-edge
			V1 = pH->source()->point();       // Start point of current half-edge
			V2 = pH->target()->point();       // End point of current half-edge
			V3 = pHNh->target()->point();     // End point of next half-edge

			// Create a new Bezier curve
			BE_bezierCurve_approxi* be_curve = new BE_bezierCurve_approxi;

			// Compute intermediate control points
			V4 = V0 * (1.0 / 3) + V1 * (2.0 / 3);
			V5 = V1 * (2.0 / 3) + V2 * (1.0 / 3);
			V6 = V1 * (1.0 / 3) + V2 * (2.0 / 3);
			V7 = V2 * (2.0 / 3) + V3 * (1.0 / 3);

			// Check if sharp feature point
			int sourceId = pH->source()->id();
			int targetId = pH->target()->id();

			// Set starting control point
			if (sharpFeaturePoints.find(sourceId) != sharpFeaturePoints.end()) {
				be_curve->cpts()[0] = sharpFeaturePoints[sourceId]; // Use original coordinates
				/*std::cout << "Restored sharp feature at vertex " << sourceId
						  << " |be_curve->cpts Position: ("
						  << be_curve->cpts()[0][0] << ", "
						  << be_curve->cpts()[0][1] << ", "
						  << be_curve->cpts()[0][2] << ")\n";*/
			}
			else {
				be_curve->cpts()[0] = V4 * 0.5 + V5 * 0.5;
			}

			// Set ending control point
			if (sharpFeaturePoints.find(targetId) != sharpFeaturePoints.end()) {
				be_curve->cpts()[3] = sharpFeaturePoints[targetId]; // Use original coordinates
				/*std::cout << "Restored sharp feature at vertex " << targetId
						  << " |be_curve->cpts Position: ("
						  << be_curve->cpts()[3][0] << ", "
						  << be_curve->cpts()[3][1] << ", "
						  << be_curve->cpts()[3][2] << ")\n";*/
			}
			else {
				be_curve->cpts()[3] = V6 * 0.5 + V7 * 0.5;
			}

			// Set intermediate control points
			be_curve->cpts()[1] = V5;
			be_curve->cpts()[2] = V6;

			// Add the curve to the bcs container
			this->bcs().push_back(be_curve);

			// Store mapping relationship
			quadHalfedgeToBezierCurveMap[pH] = be_curve;
		}
	}
	//std::cout << "reset curve ControlPoints successfully!" << std::endl;

	for (auto arc : be_quadFeatureLoops->arcs()) {
		if (!arc || arc->halfedges().empty()) continue;

		// Get fixed endpoints of the arc
		CPoint arcStartPoint = arc->halfedges().front()->source()->point();
		CPoint arcEndPoint = arc->halfedges().back()->target()->point();

		for (size_t i = 0; i < arc->halfedges().size(); ++i) {
			H* pH = arc->halfedges()[i];
			CPoint V1 = (i == 0) ? arcStartPoint : pH->source()->point();
			CPoint V2 = (i == arc->halfedges().size() - 1) ? arcEndPoint : pH->target()->point();

			// Handle adjacent points (ensure V0/V3 coincide with V1/V2 at boundary)
			CPoint V0 = (i == 0) ? V1 : arc->halfedges()[i - 1]->source()->point();
			CPoint V3 = (i == arc->halfedges().size() - 1) ? V2 : arc->halfedges()[i + 1]->target()->point();

			// Create curve and set control points
			BE_bezierCurve_approxi* be_curve = new BE_bezierCurve_approxi;

			// Fixed endpoints (force using original coordinates)
			be_curve->cpts()[0] = V1;
			be_curve->cpts()[3] = V2;

			// Compute internal control points (consistent with feature loop)
			CPoint V5 = V1 * (2.0 / 3) + V2 * (1.0 / 3);
			CPoint V6 = V1 * (1.0 / 3) + V2 * (2.0 / 3);
			be_curve->cpts()[1] = V5;
			be_curve->cpts()[2] = V6;

			// Store the curve
			this->bcs().push_back(be_curve);
			quadHalfedgeToBezierCurveMap[pH] = be_curve;

			// Debug output
			if (i == 0) {
				//std::cout << "Arc START point: (" << V1[0] << "," << V1[1] << "," << V1[2] << ")\n";
			}
			if (i == arc->halfedges().size() - 1) {
				//std::cout << "Arc END point: (" << V2[0] << "," << V2[1] << "," << V2[2] << ")\n";
			}
		}
	}
}

void CCG_QMSLib::BE_curve_model::mapTriLoopToQuadLoop(M* quadMesh, M* triMesh)
{
	//std::cout << "=========================================================" << std::endl;
	if (!be_triFeatureLoops || !be_quadFeatureLoops) {
		//std::cerr << "Feature loops are not initialized!" << std::endl;
		return;
	}

	//std::cout << "Mapping triFeatureLoops points to quadFeatureLoops..." << std::endl;

	// Iterate over feature loops of the triangular mesh
	for (auto triLoop : be_triFeatureLoops->loops()) {
		double minDistanceSum = std::numeric_limits<double>::max(); // Initialize to maximum value
		int closestQuadLoopId = -1;

		// Iterate over quadrilateral mesh feature loops, sum distances
		for (size_t quadLoopIdx = 0; quadLoopIdx < be_quadFeatureLoops->loops().size(); ++quadLoopIdx) {
			auto quadLoop = be_quadFeatureLoops->loops()[quadLoopIdx];
			double distanceSum = 0.0;

			// Iterate over all points in the triangular mesh loop
			for (auto triHalfEdge : triLoop->halfedges()) {
				V* triVertex = triMesh->halfedgeSource(triHalfEdge);
				CPoint triPoint = triVertex->point();

				// Iterate over all points in the quadrilateral mesh loop, compute distances
				double minPointDistance = std::numeric_limits<double>::max();
				for (auto quadHalfEdge : quadLoop->halfedges()) {
					V* quadVertex = quadMesh->halfedgeSource(quadHalfEdge);
					CPoint quadPoint = quadVertex->point();

					double distance = (triPoint - quadPoint).norm(); // Compute Euclidean distance
					minPointDistance = std::min(minPointDistance, distance);
				}

				// Accumulate the minimum distance from the current point to the quadrilateral loop
				distanceSum += minPointDistance;
			}

			// Compare total distance for the current loop, update minimum
			if (distanceSum < minDistanceSum) {
				minDistanceSum = distanceSum;
				closestQuadLoopId = quadLoop->getQuadLoopId(); // Use getQuadLoopId to get loop ID
			}
		}

		// Output matching results
		/*std::cout << "TriFeatureLoop with ID: " << triLoop->getTriLoopId()
				  << " matched with QuadFeatureLoop ID: " << closestQuadLoopId
				  << " with distance sum: " << minDistanceSum << std::endl;*/

				  // Store the mapping relationship
		triToQuadLoopMap[triLoop->getTriLoopId()] = closestQuadLoopId;
	}

	// 2. Matching new feature arcs
	if (be_triFeatureLoops->arcs().empty() || be_quadFeatureLoops->arcs().empty())
	{
		//std::cerr << "Feature loops are not initialized!" << std::endl;
		return;
	}
	//std::cout << "\nMapping triFeatureArcs to quadFeatureArcs..." << std::endl;

	for (auto triArc : be_triFeatureLoops->arcs()) {
		double minDistanceSum = std::numeric_limits<double>::max();
		int closestQuadArcId = -1;

		for (size_t quadArcIdx = 0; quadArcIdx < be_quadFeatureLoops->arcs().size(); ++quadArcIdx) {
			auto quadArc = be_quadFeatureLoops->arcs()[quadArcIdx];
			double distanceSum = 0.0;

			// Compute sum of minimum distances from triangle arc points to quadrilateral arc
			for (auto triHalfEdge : triArc->halfedges()) {
				V* triVertex = triMesh->halfedgeSource(triHalfEdge);
				CPoint triPoint = triVertex->point();

				double minPointDistance = std::numeric_limits<double>::max();
				for (auto quadHalfEdge : quadArc->halfedges()) {
					V* quadVertex = quadMesh->halfedgeSource(quadHalfEdge);
					CPoint quadPoint = quadVertex->point();
					double distance = (triPoint - quadPoint).norm();
					minPointDistance = std::min(minPointDistance, distance);
				}
				distanceSum += minPointDistance;
			}

			// Update best match
			if (distanceSum < minDistanceSum) {
				minDistanceSum = distanceSum;
				closestQuadArcId = quadArc->getQuadArcId();
			}
		}

		/*std::cout << "TriFeatureArc ID: " << triArc->getTriArcId()
				  << " -> QuadFeatureArc ID: " << closestQuadArcId
				  << " (distance: " << minDistanceSum << ")" << std::endl;*/
		triToQuadArcMap[triArc->getTriArcId()] = closestQuadArcId;
	}
	//std::cout << "All mappings complete!" << std::endl;
	//std::cout << "=========================================================" << std::endl;
}

void CCG_QMSLib::BE_curve_model::triPointToQuadHalfedge(M* quadMesh, M* triMesh) {
	if (!be_triFeatureLoops || !be_quadFeatureLoops) {
		//std::cerr << "Feature loops are not initialized!" << std::endl;
		return;
	}

	// Initialize quad half-edge ID map
	initializeQuadHalfedgeIdMap(quadMesh);
	int i = 0;
	// Iterate over triangular feature loops
	for (auto triLoop : be_triFeatureLoops->loops()) {
		// Find the corresponding quad feature loop for the triangular feature loop
		int quadLoopId = triToQuadLoopMap[triLoop->getTriLoopId()];
		auto quadLoop = be_quadFeatureLoops->loops()[quadLoopId];

		// Iterate over each point in the triangular feature loop
		for (auto triHalfEdge : triLoop->halfedges()) {
			V* triVertex = triMesh->halfedgeSource(triHalfEdge); // Triangular vertex
			CPoint triPoint = triVertex->point();

			// Initialize minimum distance and closest half-edge
			double minDistance = std::numeric_limits<double>::max();
			H* closestHalfEdge = nullptr;
			CPoint start, end;
			// Iterate over all half-edges in the quad feature loop
			for (auto quadHalfEdge : quadLoop->halfedges()) {
				V* quadSource = quadMesh->halfedgeSource(quadHalfEdge); // Half-edge source
				V* quadTarget = quadMesh->halfedgeTarget(quadHalfEdge); // Half-edge target

				CPoint p0 = quadSource->point();
				CPoint p1 = quadTarget->point();

				// Compute distance from point to segment
				double distance = pointToSegmentDistance(triPoint, p0, p1);

				// Update minimum distance and closest half-edge
				if (distance < minDistance) {
					minDistance = distance;
					closestHalfEdge = quadHalfEdge;
					start = p0;
					end = p1;
				}
			}

			// Store the mapping relationship
			triPointToQuadHalfedgeMap[triVertex] = closestHalfEdge;
			/*CPoint dataVertex = triVertex->point();
			dataPointToQuadHalfedgeMap[dataVertex] = closestHalfEdge;*/
			if ((triPoint - start) * (end - start) < 0.0 ||
				(end - start) * (end - start) < (triPoint - start) * (end - start))
			{
				//std::cout << "11111" << std::endl;
				i++;
				continue;
			}
			// Create sampling point and record information
			QB_curve_sampling* sampling = new QB_curve_sampling();
			sampling->pos() = triPoint;                     // Set the position of the point
			sampling->quadLoopId() = quadLoopId;            // Set the corresponding quad loop ID
			sampling->quadHalfedgesId() = quadHalfedgeToIdMap[closestHalfEdge]; // Set the corresponding half-edge ID
			qbm_samplings.push_back(sampling); // Add sampling point to the sample collection
			dataPointToQuadHalfedgeMap[sampling] = closestHalfEdge;
		}
	}
	//std::cout << "LOOP: number of ignored data points: " << i << std::endl;
	//std::cout << "LOOP: Triangle points mapped to quad half-edges and samples created!" << std::endl;

	// 2. Handle feature arc mapping
	int ignoredPoints = 0;
	for (auto triArc : be_triFeatureLoops->arcs()) {
		int quadArcId = triToQuadArcMap[triArc->getTriArcId()];
		auto quadArc = be_quadFeatureLoops->arcs()[quadArcId];

		for (auto triHalfEdge : triArc->halfedges()) {
			V* triVertex = triMesh->halfedgeSource(triHalfEdge);
			CPoint triPoint = triVertex->point();

			double minDistance = std::numeric_limits<double>::max();
			H* closestHalfEdge = nullptr;
			CPoint start, end;

			for (auto quadHalfEdge : quadArc->halfedges()) {
				V* quadSource = quadMesh->halfedgeSource(quadHalfEdge);
				V* quadTarget = quadMesh->halfedgeTarget(quadHalfEdge);
				CPoint p0 = quadSource->point();
				CPoint p1 = quadTarget->point();

				double distance = pointToSegmentDistance(triPoint, p0, p1);
				if (distance < minDistance) {
					minDistance = distance;
					closestHalfEdge = quadHalfEdge;
					start = p0;
					end = p1;
				}
			}

			triPointToQuadHalfedgeMap[triVertex] = closestHalfEdge;
			if ((triPoint - start) * (end - start) < 0.0 ||
				(end - start) * (end - start) < (triPoint - start) * (end - start)) {
				ignoredPoints++;
				continue;
			}

			QB_curve_sampling* sampling = new QB_curve_sampling();
			sampling->pos() = triPoint;
			sampling->quadArcId() = quadArcId;  // Use arcId instead of loopId
			sampling->quadHalfedgesId() = quadHalfedgeToIdMap[closestHalfEdge];
			qbm_samplings.push_back(sampling);
		}
	}

	//std::cout << "Arc: number of ignored data points: " << ignoredPoints << std::endl;
	//std::cout << "Arc: Triangle points mapped to quad half-edges (both loops and arcs) and samples created!" << std::endl;
}

void CCG_QMSLib::BE_curve_model::triPointToQuadHalfedge2(M* quadMesh, M* triMesh)
{
	// Iterate over triangular feature loops
	for (auto triLoop : be_triFeatureLoops->loops()) {
		// Find the corresponding quad feature loop for the triangular feature loop
		int quadLoopId = triToQuadLoopMap[triLoop->getTriLoopId()];
		auto quadLoop = be_quadFeatureLoops->loops()[quadLoopId];

		// Iterate over each point in the triangular feature loop
		for (auto triHalfEdge : triLoop->halfedges()) {
			V* triVertex = triMesh->halfedgeSource(triHalfEdge); // Triangular vertex
			CPoint triPoint = triVertex->point();

			// Initialize minimum distance and closest half-edge
			double minDistance = std::numeric_limits<double>::max();
			H* closestHalfEdge = nullptr;
			CPoint start, end;
			// Iterate over all half-edges in the quad feature loop
			for (auto quadHalfEdge : quadLoop->halfedges()) {
				V* quadSource = quadMesh->halfedgeSource(quadHalfEdge); // Half-edge source
				V* quadTarget = quadMesh->halfedgeTarget(quadHalfEdge); // Half-edge target

				CPoint p0 = quadSource->point();
				CPoint p1 = quadTarget->point();

				// Compute distance from point to segment
				double distance = pointToSegmentDistance(triPoint, p0, p1);

				// Update minimum distance and closest half-edge
				if (distance < minDistance) {
					minDistance = distance;
					closestHalfEdge = quadHalfEdge;
					start = p0;
					end = p1;
				}
			}

			// Store the mapping relationship
			triPointToQuadHalfedgeMap[triVertex] = closestHalfEdge;
		}
	}
}


void CCG_QMSLib::BE_curve_model::triBoundaryHalfedgesMapping2QuadBoundaryHalfedges()
{
	/*check sequence of tri boundary halfedges on loop*/
	/*for (auto triLoop : this->triFeatureLoops()->loops())
	{
		for (auto triBoundaryH : triLoop->halfedges())
		{
			std::cout << triBoundaryH->source()->id() << " ---> " << triBoundaryH->target()->id() << std::endl;
		}
	}*/

	///*check tri loop id*/
	//for (auto triLoop : this->triFeatureLoops()->loops())
	//{
	//	std::cout << "triLoop Id: " << triLoop->getTriLoopId() <<" quadLoop Id: " << getTriToQuadLoopMap().at(triLoop->getTriLoopId()) << std::endl;
	//}
	///*check quad loop id*/
	//for (auto quadLoop : this->quadFeatureLoops()->loops())
	//{
	//	std::cout << "quadLoop Id: " << quadLoop->getQuadLoopId() << std::endl;
	//}

	/*Given threshold for specify double type is 0*/
	double threshhold_triBoundaryHalfedgesMapping2QuadBoundaryHalfedges = 1e-10;
	initializeQuadHalfedgeIdMap();
	initializeTriHalfedgeIdMap();
	/*Traverse all  boundary loops on quad mesh*/
	for (auto quadLoop : this->quadFeatureLoops()->loops())
	{
		/*Traverse all boundary halfedges on quad mesh loops*/
		for (auto quadBoundaryH : quadLoop->halfedges())
		{
			H* locatingSourceTriH = NULL;
			H* locatingTargetTriH = NULL;
			double minDistance_source = std::numeric_limits<double>::max();
			double minDistance_target = std::numeric_limits<double>::max();
			std::vector<QB_curve_sampling*> quadBoundaryH_samplingPts;
			/*create mapping between quadBoundaryH and triBoundaryH*/
			/*1. find the triLoop according the relation triToQuadLoopMap*/
			int triLoopId_quadBoundaryH = -1;
			for (auto triLoop : this->triFeatureLoops()->loops())
			{
				if (getTriToQuadLoopMap().at(triLoop->getTriLoopId()) == quadLoop->getQuadLoopId())
				{
					triLoopId_quadBoundaryH = triLoop->getTriLoopId();
					/*2. Traverse all boundary halfedges on this tri loop*/
					for (auto triBoundaryH : triLoop->halfedges())
					{
						/*3. Specify quadBoundaryH's source is or not on this triBoundaryH*/
						{
							CPoint v = quadBoundaryH->source()->point();
							CPoint tri_v0 = triBoundaryH->source()->point();
							CPoint tri_v1 = triBoundaryH->target()->point();
							/*if (((v - tri_v0) ^ (tri_v1 - tri_v0)).norm() < threshhold_triBoundaryHalfedgesMapping2QuadBoundaryHalfedges && ((v - tri_v0)*(tri_v1 - tri_v0)) > 0.0 && (v - tri_v0).norm() <= (tri_v1 - tri_v0).norm())
							{
								locatingSourceTriH = triBoundaryH;
								locatingSuccess_source = true;
							}*/
							double dist = pointToSegmentDistance(v, tri_v0, tri_v1);
							if (dist < minDistance_source) {
								minDistance_source = dist;
								locatingSourceTriH = triBoundaryH;
							}
						}
						/*4. Specify quadBoundaryH's target is or not on this triBoundaryH*/
						{
							CPoint v = quadBoundaryH->target()->point();
							CPoint tri_v0 = triBoundaryH->source()->point();
							CPoint tri_v1 = triBoundaryH->target()->point();
							/*if (((v - tri_v0) ^ (tri_v1 - tri_v0)).norm() < threshhold_triBoundaryHalfedgesMapping2QuadBoundaryHalfedges && ((v - tri_v0) * (tri_v1 - tri_v0)) > 0.0 && (v - tri_v0).norm() <= (tri_v1 - tri_v0).norm())
							{
								locatingTargetTriH = triBoundaryH;
								locatingSuccess_target = true;
							}*/
							double dist = pointToSegmentDistance(v, tri_v0, tri_v1);
							if (dist < minDistance_target) {
								minDistance_target = dist;
								locatingTargetTriH = triBoundaryH;
							}
						}
					}
				}
				else
				{
					continue;
					//std::cout << " No finding tri loop at mapping between quadBoundaryH and triBoundaryH!" << std::endl;
					//std::cout << quadBoundaryH->source()->id() << " ---> " << quadBoundaryH->target()->id() << std::endl;
					//return;
				}
				if (locatingSourceTriH != NULL && locatingTargetTriH != NULL)
				{
					mapQuadHalfedgesAndTriHalfedges.insert(std::pair<H*, std::pair<H*, H*>>(quadBoundaryH, std::pair<H*, H*>(locatingSourceTriH, locatingTargetTriH)));
					/*std::cout << "Locating quadBoundaryH "<< quadBoundaryH->source()->id() << " is successful!" << std::endl;
					std::cout << "source: " << mapQuadHalfedgesAndTriHalfedges.at(quadBoundaryH).first->target()->id() << std::endl;
					std::cout << "target: " << mapQuadHalfedgesAndTriHalfedges.at(quadBoundaryH).second->target()->id() << std::endl;*/
					/*5. build mapping between sampling points from tri mesh boundary and quad boundary halfedges*/
					if (locatingSourceTriH->source()->id() == locatingTargetTriH->source()->id() && locatingSourceTriH->target()->id() == locatingTargetTriH->target()->id())
					{
						//collect all sampling points position about the quadBoundaryH
						std::vector<CPoint> samplingPts_quadBoundaryH;
						samplingPts_quadBoundaryH.push_back(quadBoundaryH->source()->point());
						samplingPts_quadBoundaryH.push_back(quadBoundaryH->target()->point());

						//create all sampling points about the quadBoundaryH
						for (auto samplingPt_quadBoundaryH : samplingPts_quadBoundaryH)
						{
							QB_curve_sampling* sampling = new QB_curve_sampling();
							sampling->pos() = samplingPt_quadBoundaryH;
							sampling->quadLoopId() = quadLoop->getQuadLoopId();
							sampling->quadHalfedgesId() = quadHalfedgeToIdMap[quadBoundaryH];
							qbm_samplings.push_back(sampling);
							dataPointToQuadHalfedgeMap[sampling] = quadBoundaryH;
							quadBoundaryH_samplingPts.push_back(sampling);
						}
					}
					else
					{
						if (triLoopId_quadBoundaryH > -1 && triLoopId_quadBoundaryH < (this->triFeatureLoops()->loops().size()))
						{
							auto startIt = std::find(this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().begin(), this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().end(), locatingSourceTriH);
							auto endIt = std::find(this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().begin(), this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().end(), locatingTargetTriH);

							if (startIt == this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().end() || endIt == this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().end()) {
								std::cout << "One or both elements not found!" << std::endl;
								return;
							}
							//collect all sampling points position about the quadBoundaryH
							std::vector<CPoint> samplingPts_quadBoundaryH;
							samplingPts_quadBoundaryH.push_back(quadBoundaryH->source()->point());

							if (triHalfedgeToIdMap.at(locatingSourceTriH) > triHalfedgeToIdMap.at(locatingTargetTriH))
							{
								std::list<H*> subList1(startIt, this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().end());
								std::list<H*> subList2(this->triFeatureLoops()->loops()[triLoopId_quadBoundaryH]->halfedges().begin(), endIt);
								/*std::list<H*> subList(startIt, endIt);
								std::cout << "---------" << quadBoundaryH->source()->id() << " ---> " << quadBoundaryH->target()->id() << std::endl;
								for (auto e : subList1)
								{
									std::cout << e->source()->id() << " ---> " << e->target()->id() << std::endl;
								}
								for (auto e : subList2)
								{
									std::cout << e->source()->id() << " ---> " << e->target()->id() << std::endl;
								}*/
								for (auto e : subList1)
								{
									samplingPts_quadBoundaryH.push_back(e->target()->point());
								}
								for (auto e : subList2)
								{
									samplingPts_quadBoundaryH.push_back(e->target()->point());
								}
							}
							else
							{
								std::list<H*> subList(startIt, endIt);
								//std::cout << "---------" << quadBoundaryH->source()->id() << " ---> " << quadBoundaryH->target()->id() << std::endl;
								////std::cout << " distance(startIt, endIt): " << distance(startIt, endIt) << std::endl;
								//for (auto e : subList)
								//{
								//	std::cout << e->source()->id() << " ---> " << e->target()->id() << std::endl;
								//}
								for (auto e : subList)
								{
									samplingPts_quadBoundaryH.push_back(e->target()->point());
								}
							}
							samplingPts_quadBoundaryH.push_back(quadBoundaryH->target()->point());

							//create all sampling points about the quadBoundaryH
							for (auto samplingPt_quadBoundaryH : samplingPts_quadBoundaryH)
							{
								QB_curve_sampling* sampling = new QB_curve_sampling();
								sampling->pos() = samplingPt_quadBoundaryH;
								sampling->quadLoopId() = quadLoop->getQuadLoopId();
								sampling->quadHalfedgesId() = quadHalfedgeToIdMap[quadBoundaryH];
								qbm_samplings.push_back(sampling);
								dataPointToQuadHalfedgeMap[sampling] = quadBoundaryH;
								quadBoundaryH_samplingPts.push_back(sampling);
							}
						}
						else
						{
							std::cout << " No finding tri loop at building mapping between sampling points from tri mesh boundary and quad boundary halfedges!" << std::endl;
							return;
						}
					}
					break;
				}
				else
				{
					std::cout << " No locating quadBoundaryH at mapping between quadBoundaryH and triBoundaryH!" << std::endl;
					std::cout << quadBoundaryH->source()->id() << " ---> " << quadBoundaryH->target()->id() << std::endl;
					//return;
				}
			}
			samplingPts2QuadHalfedges.insert(std::pair<H*, std::vector<QB_curve_sampling*>>(quadBoundaryH, quadBoundaryH_samplingPts));
		}
	}
	/*check sequence of tri boundary halfedges on loop*/
	/*for (auto triLoop : this->triFeatureLoops()->loops())
	{
		for (auto triBoundaryH : triLoop->halfedges())
		{
			std::cout << triBoundaryH->source()->id() << " ---> " << triBoundaryH->target()->id() << std::endl;
		}
	}*/
	/*check initial sampling points on quad boundary halfedge*/
	/*for (auto samplingPts2QuadHalfedge : samplingPts2QuadHalfedges)
	{
		std::cout << "******" << samplingPts2QuadHalfedge.first->source()->id() << "--->" << samplingPts2QuadHalfedge.first->target()->id() << std::endl;
		std::cout << " # sampling pts: " << samplingPts2QuadHalfedge.second.size() << std::endl;
	}*/
}

void CCG_QMSLib::BE_curve_model::locatePoorFittingBoundaryEdges(M* quadMesh, int minPointsThreshold, double errorThreshold)
{
	poorFittingEdges.clear();
	endpointCountMap.clear();
	m_boundaryFittingArcs.clear();  // 清空历史数据
	if (!quadFeatureLoops()) {
		std::cerr << "Error: Quad feature loops not initialized!" << std::endl;
		return;
	}

	// 确保参数化已完成
	//parametrize2(quadMesh);

	// 1. 遍历所有四边形网格特征环
	for (auto loop : quadFeatureLoops()->loops()) {
		for (auto he : loop->halfedges()) {
			auto it = samplingPts2QuadHalfedges.find(he);
			if (it == samplingPts2QuadHalfedges.end()) continue;

			const auto& pointsOnEdge = it->second;

			// 1.1 检查待拟合点数量
			if (pointsOnEdge.size() < minPointsThreshold) {
				//poorFittingEdges[he] = { 0.0, pointsOnEdge };
				continue;
			}

			// 1.2 计算最大拟合误差（基于参数化位置）//修改为：计算边界三角形网格点到它所对应的四边形网格边界边的距离，并且将最长距离作为它们所在的四边形边界边的贴体性度量
			double maxError = 0.0;
			V* startV = quadMesh->halfedgeSource(he);
			V* endV = quadMesh->halfedgeTarget(he);
			CPoint p0 = startV->point();
			CPoint p1 = endV->point();
			V* markTriV = NULL;
			for (auto point : pointsOnEdge) {
				//// 获取该点的参数值
				//double u = dataPointIdToParamMap[point->vId()];
				////std::cout << "u: " << u << std::endl;
				//// 计算参数化位置（线性插值）
				//CPoint interpolatedPos = p0 * (1.0 - u) + p1 * u;
				//// 计算拟合误差（数据点位置与参数化位置的距离）
				//double error = (point->pos() - interpolatedPos).norm();
				//std::cout << "error1: " << error << std::endl;
				double error = pointToSegmentDistance(point->pos(), p0, p1);
				//std::cout << "error2: " << error << std::endl;
				maxError = std::max(maxError, error);
			}

			// 1.3 检查误差阈值
			if (maxError > errorThreshold) {
				poorFittingEdges[he] = { maxError, pointsOnEdge };
				//std::cout << "startV->id:" << startV->id() << std::endl;
				//std::cout << "maxError: " << maxError << std::endl;
			}
		}
	}

	/*check poorFittingEdges */
	/*for (auto poorFittingEdge : poorFittingEdges)
	{
		std::cout << "*****" << poorFittingEdge.second.first << std::endl;
	}*/

	// 2. 统计端点出现次数
	for (const auto& pair : poorFittingEdges) {
		H* he = pair.first;
		V* startV = quadMesh->halfedgeSource(he);
		V* endV = quadMesh->halfedgeTarget(he);
		endpointCountMap[startV]++;
		endpointCountMap[endV]++;
	}

	//// 3. 初始化未被记录的顶点（出现次数为0）
	//for (auto loop : quadFeatureLoops()->loops()) {
	//	for (auto he : loop->halfedges()) {
	//		V* startV = quadMesh->halfedgeSource(he);
	//		V* endV = quadMesh->halfedgeTarget(he);
	//		if (endpointCountMap.find(startV) == endpointCountMap.end()) {
	//			endpointCountMap[startV] = 0;
	//		}
	//		if (endpointCountMap.find(endV) == endpointCountMap.end()) {
	//			endpointCountMap[endV] = 0;
	//		}
	//	}
	//}

	// 打印端点出现次数（调试用）
	/*std::cout << "Vertex occurrence count in poor fitting edges:" << std::endl;
	for (const auto& pair : endpointCountMap) {
		std::cout << "Vertex ID: " << pair.first->id() << ", Count: " << pair.second << std::endl;
	}*/
	// 4. 构建边界特征弧
	std::set<H*> processedEdges;  // 记录已处理的半边（避免重复）

	for (const auto& pair : poorFittingEdges) {
		H* he = pair.first;  // 当前待处理的半边
		if (processedEdges.count(he)) continue;  // 跳过已处理的边

		// 获取当前边的起点和终点顶点
		V* startV = quadMesh->halfedgeSource(he);
		V* endV = quadMesh->halfedgeTarget(he);

		// 1. 判断是否为弧的起点/终点（端点计数为1）
		bool isStartPoint = (endpointCountMap[startV] == 1);  // 起点：计数=1
		bool isEndPoint = (endpointCountMap[endV] == 1);      // 终点：计数=1

		// 跳过孤立边（两个端点均非起点/终点，计数≠1）
		if (!isStartPoint && !isEndPoint) {
			continue;
		}

		// 2. 确定弧的延伸方向（从起点出发）
		H* startHe = isStartPoint ? he : quadMesh->halfedgeSym(he);  // 正向或反向
		BoundaryFittingArc arc;  // 存储当前弧的所有半边
		H* currentHe = startHe;  // 当前延伸中的半边

		// 3. 沿着 poorFittingEdges 延伸弧
		while (currentHe && !processedEdges.count(currentHe)) {
			// 3.1 添加当前半边到弧
			arc.orderedHalfedges.push_back(currentHe);
			processedEdges.insert(currentHe);  // 标记为已处理

			// 3.2 寻找下一条有效半边（必须属于 poorFittingEdges）
			H* nextHe = nullptr;  // 下一条候选半边
			V* currentVertex = quadMesh->halfedgeTarget(currentHe);
			for (auto outH : It::VCcwOutHEIterator(quadMesh, currentVertex)) {
				// 条件1：outH 必须属于 poorFittingEdges（待拟合边）
				if (!poorFittingEdges.count(outH)) {
					continue;
				}

				// 条件2：outH 不能是 currentHe 的对偶半边（反向边）
				H* currentHeSym = quadMesh->halfedgeSym(currentHe);  // currentHe 的对偶边
				if (outH == currentHeSym) {
					continue;  // 跳过对偶边，避免反向延伸
				}
				// 条件3：outH 的方向应与边界一致（可选，确保沿着边界延伸）
				if (quadMesh->halfedgeEdge(outH)->boundary()) {
					nextHe = outH;
					break;
				}

			}
			// 3.3 检查终止条件（无法延伸或到达终点）
			if (!nextHe) {
				//std ::cout << "111111111111" << std::endl;// 无下一条有效边，终止延伸
				break;
			}

			// 检查下一条边的终点是否为弧的终点（计数=1）
			V* nextTargetV = quadMesh->halfedgeTarget(nextHe);
			if (endpointCountMap[nextTargetV] == 1) {
				// 添加终点边并终止
				arc.orderedHalfedges.push_back(nextHe);
				processedEdges.insert(nextHe);
				break;
			}

			// 3.4 继续延伸到下一条边
			currentHe = nextHe;
		}

		// 保存有效的弧
		if (!arc.orderedHalfedges.empty()) {
			m_boundaryFittingArcs.push_back(arc);
			/*std::cout << "Arc from vertex " << quadMesh->halfedgeSource(arc.orderedHalfedges.front())->id()
				<< " to " << quadMesh->halfedgeTarget(arc.orderedHalfedges.back())->id()
				<< ", size: " << arc.orderedHalfedges.size() << std::endl;*/
		}
	}
}

void CCG_QMSLib::BE_curve_model::addInterpolatedPointsForUnmappedHalfedges(M* pMesh)
{
	const int TARGET_POINTS_PER_EDGE = 4; // Target number of data points per halfedge

	// 1. Collect all unmapped quadrilateral halfedges and their associated loops
	std::vector<std::pair<H*, CLoop*>> unMappedPairs;
	for (auto quadLoop : be_quadFeatureLoops->loops()) {
		for (auto quadHe : quadLoop->halfedges()) {
			// Count the number of existing sampling points for the current halfedge
			int existingPoints = std::count_if(dataPointToQuadHalfedgeMap.begin(),
				dataPointToQuadHalfedgeMap.end(),
				[quadHe](const auto& entry) { return entry.second == quadHe; });

			if (existingPoints < TARGET_POINTS_PER_EDGE) {
				unMappedPairs.emplace_back(quadHe, quadLoop);
				/*std::cout << "Unmapped Quad Halfedge ID: " << quadHe->localId()
					<< " (needs " << TARGET_POINTS_PER_EDGE - existingPoints
					<< " more points)" << std::endl;*/
			}
		}
	}

	// 2. Process each unmapped halfedge
	for (auto& entry : unMappedPairs) {
		H* quadHe = entry.first;
		CLoop* quadLoop = entry.second;

		// Calculate how many additional points are needed
		int existingPoints = std::count_if(dataPointToQuadHalfedgeMap.begin(),
			dataPointToQuadHalfedgeMap.end(),
			[quadHe](const auto& entry) { return entry.second == quadHe; });
		int pointsNeeded = TARGET_POINTS_PER_EDGE - existingPoints;

		V* v0 = pMesh->halfedgeSource(quadHe);
		V* v1 = pMesh->halfedgeTarget(quadHe);
		CPoint p0 = v0->point();
		CPoint p1 = v1->point();

		// 3. Find the corresponding triangle loop
		int quadLoopId = quadLoop->getQuadLoopId();
		int triLoopId = -1;
		for (const auto& pair : triToQuadLoopMap) {
			if (pair.second == quadLoopId) {
				triLoopId = pair.first;
				break;
			}
		}
		if (triLoopId == -1) continue;

		CLoop* triLoop = be_triFeatureLoops->loops()[triLoopId];

		// 4. Find the closest point in the triangle loop
		std::pair<V*, double> result = findClosestTriPoint(triLoop, p0, p1);
		V* closestVert = result.first;
		if (!closestVert) continue;

		// 5. Recursive interpolation logic
		int recursionDepth = 0;
		const int maxRecursion = 30;
		std::vector<CPoint> newSamples; // Store newly found sample points in this iteration

		do {
			newSamples.clear();
			// Get adjacent halfedges
			std::pair<H*, H*> adjacentHEs = findAdjacentHalfedges(triLoop, closestVert);
			H* prevHe = adjacentHEs.first;
			H* nextHe = adjacentHEs.second;

			// Number of interpolation samples increases with recursion depth (3, 5, 7, ...)
			int sampleCount = 3 + recursionDepth * 2;
			std::vector<CPoint> samples = interpolateHalfedges(prevHe, nextHe, sampleCount);

			// Find the closest quadrilateral halfedge for each interpolation sample
			for (auto& sample : samples) {
				H* closestQuadHe = findClosestQuadHalfedge(quadLoop, sample);
				if (closestQuadHe == quadHe) {
					newSamples.push_back(sample);
					/*std::cout << "Found new sample for Halfedge ID: "
						<< quadHe->localId() << std::endl;*/

						// If enough new points have been found, break
					if (newSamples.size() >= pointsNeeded)
						break;
				}
			}

			recursionDepth++;

			// Exit condition: enough points found or reached maximum recursion depth
			if (newSamples.size() >= pointsNeeded || recursionDepth >= maxRecursion)
				break;

		} while (true);

		// 6. Add sampled points to the collection (add only the required number)
		int pointsToAdd = std::min(static_cast<int>(newSamples.size()), pointsNeeded);
		for (int i = 0; i < pointsToAdd; i++) {
			QB_curve_sampling* sampling = new QB_curve_sampling();
			sampling->pos() = newSamples[i];
			sampling->quadLoopId() = quadLoopId;
			sampling->quadHalfedgesId() = quadHalfedgeToIdMap[quadHe];
			qbm_samplings.push_back(sampling);
			dataPointToQuadHalfedgeMap[sampling] = quadHe;
		}

		/*std::cout << "Added " << pointsToAdd << " new points for Halfedge ID "
			<< quadHe->localId() << " (Total: "
			<< existingPoints + pointsToAdd << "/" << TARGET_POINTS_PER_EDGE
			<< ")" << std::endl;*/
	}
}

void CCG_QMSLib::BE_curve_model::addInterpolatedPointsForUnmappedHalfedges2()
{
	const int TARGET_POINTS_PER_EDGE = 8; // Target number of data points per halfedge
	for (auto& samplingPts2QuadHalfedge : samplingPts2QuadHalfedges)
	{
		if (TARGET_POINTS_PER_EDGE > samplingPts2QuadHalfedge.second.size())
		{
			std::vector<CPoint> newSamples;
			if (samplingPts2QuadHalfedge.second.size() > 1)
			{
				int sampleCount = (TARGET_POINTS_PER_EDGE - samplingPts2QuadHalfedge.second.size())
					/ (samplingPts2QuadHalfedge.second.size() - 1);
				sampleCount = sampleCount > 0 ? sampleCount : 1;
				for (int i = 0; i < samplingPts2QuadHalfedge.second.size() - 1; i++)
				{
					//obtaing new sampling points' position
					CPoint p0 = samplingPts2QuadHalfedge.second[i]->pos();
					CPoint p1 = samplingPts2QuadHalfedge.second[i + 1]->pos();
					for (int j = 1; j <= sampleCount; ++j) {
						double t = j / (double)(sampleCount + 1);
						newSamples.push_back(p0 * (1 - t) + p1 * t);
					}
				}
			}
			else
			{
				//std::cout << "error in initial quad halfedges sampling count!" << std::endl;
				return;
			}
			for (auto newSample : newSamples) {
				QB_curve_sampling* sampling = new QB_curve_sampling();
				sampling->pos() = newSample;
				sampling->quadLoopId() = samplingPts2QuadHalfedge.second[0]->quadLoopId();
				sampling->quadHalfedgesId() = samplingPts2QuadHalfedge.second[0]->quadHalfedgesId();
				qbm_samplings.push_back(sampling);
				dataPointToQuadHalfedgeMap[sampling] = samplingPts2QuadHalfedge.first;
				samplingPts2QuadHalfedge.second.push_back(sampling);
			}
		}
	}
	/*check initial sampling points on quad boundary halfedge*/
	/*for (auto samplingPts2QuadHalfedge : samplingPts2QuadHalfedges)
	{
		std::cout << "******" << samplingPts2QuadHalfedge.first->source()->id() << "--->"
				  << samplingPts2QuadHalfedge.first->target()->id() << std::endl;
		std::cout << " # sampling pts: " << samplingPts2QuadHalfedge.second.size() << std::endl;
	}*/
}

std::pair<CCG_QMSLib::V*, double> CCG_QMSLib::BE_curve_model::findClosestTriPoint(
	CLoop* triLoop,
	const MeshLib::CPoint& p0,
	const MeshLib::CPoint& p1)
{
	// Implementation of the method
	V* closestVertex = nullptr;
	double minDistance = std::numeric_limits<double>::max();

	for (auto triVert : triLoop->vertex()) {
		double dist = pointToSegmentDistance(triVert->point(), p0, p1);
		if (dist < minDistance) {
			minDistance = dist;
			closestVertex = triVert;
		}
	}

	return std::make_pair(closestVertex, minDistance);
}

// Find the previous and next halfedges of a triangle point
std::pair<CCG_QMSLib::H*, CCG_QMSLib::H*> CCG_QMSLib::BE_curve_model::findAdjacentHalfedges(CLoop* triLoop, V* vert) {
	if (!triLoop || !vert) {
		return { nullptr, nullptr };
	}

	H* prevHe = nullptr;
	H* nextHe = nullptr;

	auto& halfedges = triLoop->halfedges();
	for (auto it = halfedges.begin(); it != halfedges.end(); ++it) {
		if (m_pMesh->halfedgeTarget(*it) == vert) {
			prevHe = *it;
			nextHe = (std::next(it) == halfedges.end()) ? *halfedges.begin() : *std::next(it);
			break;
		}
	}
	return { prevHe, nextHe };
}

// Interpolate triangle halfedges
std::vector<CPoint> CCG_QMSLib::BE_curve_model::interpolateHalfedges(H* he1, H* he2, int count)
{
	std::vector<CPoint> samples;

	auto addSamples = [&](H* he) {
		if (!he) return;
		CPoint p0 = m_pMesh->halfedgeSource(he)->point();
		CPoint p1 = m_pMesh->halfedgeTarget(he)->point();
		for (int i = 1; i <= count; ++i) {
			double t = i / (double)(count + 1);
			samples.push_back(p0 * (1 - t) + p1 * t);
		}
		};

	addSamples(he1);
	addSamples(he2);
	return samples;
}

CCG_QMSLib::H* CCG_QMSLib::BE_curve_model::findClosestQuadHalfedge(CLoop* quadLoop, const CPoint& point)
{
	H* closestHe = nullptr;
	double minDist = std::numeric_limits<double>::max();

	for (auto quadHe : quadLoop->halfedges()) {
		CPoint p0 = m_pMesh->halfedgeSource(quadHe)->point();
		CPoint p1 = m_pMesh->halfedgeTarget(quadHe)->point();
		double dist = pointToSegmentDistance(point, p0, p1);
		if (dist < minDist) {
			minDist = dist;
			closestHe = quadHe;
		}
	}
	return closestHe;
}

void CCG_QMSLib::BE_curve_model::addSamplingPoint(const CPoint& point, int loopId, H* quadHe)
{
	QB_curve_sampling* sp = new QB_curve_sampling();
	sp->pos() = point;
	sp->quadLoopId() = loopId;
	sp->quadHalfedgesId() = quadHalfedgeToIdMap[quadHe];
	qbm_samplings.push_back(sp);
}

double CCG_QMSLib::BE_curve_model::pointToSegmentDistance(const CPoint& point, const CPoint& segStart, const CPoint& segEnd)
{
	CPoint v = segEnd - segStart;           // direction vector of the segment
	CPoint w = point - segStart;            // vector from segment start to point

	double c1 = w * v;                      // scalar projection of point onto segment
	if (c1 <= 0.0)
		return (point - segStart).norm();   // projection falls before the start point

	double c2 = v * v;                      // squared length of the segment direction vector
	if (c2 <= c1)
		return (point - segEnd).norm();     // projection falls after the end point

	CPoint pb = segStart + v * (c1 / c2);   // coordinates of the projected point
	return (point - pb).norm();             // return distance
}

void CCG_QMSLib::BE_curve_model::initializeQuadVertexIdMap(M* quadMesh)
{
	int id = 0; // 初始ID计数器

	// 1. 处理特征环顶点ID（原有代码）
	for (auto quadLoop : be_quadFeatureLoops->loops()) {
		//std::cout << "Loop vertices count: " << quadLoop->vertex().size() << std::endl;
		for (auto vertex : quadLoop->vertex()) {
			quadLoopVertexIdToMap[vertex] = id++;
			//std::cout << "Loop Vertex ID: " << quadLoopVertexIdToMap[vertex] << std::endl;
		}
	}
	//// 处理弧
	//for (auto arc : be_quadFeatureLoops->arcs()) {
	//	for (auto v : arc->vertex()) {
	//		if (quadLoopVertexIdToMap.find(v) == quadLoopVertexIdToMap.end()) {
	//			quadLoopVertexIdToMap[v] = id++;
	//		}
	//	}
	//}
	// 2. 处理特征弧顶点ID
	for (auto quadArc : be_quadFeatureLoops->arcs()) {
		//std::cout << "Arc vertices count: " << quadArc->vertex().size() << std::endl;

		// 创建临时集合避免重复顶点
		std::set<V*> uniqueVertices;
		for (auto vertex : quadArc->vertex()) {
			uniqueVertices.insert(vertex);
		}

		// 为弧中的顶点分配ID（跳过已分配ID的顶点）
		for (auto vertex : uniqueVertices) {
			if (quadLoopVertexIdToMap.find(vertex) == quadLoopVertexIdToMap.end()) {
				quadLoopVertexIdToMap[vertex] = id++;
				//std::cout << "Arc Vertex ID: " << quadLoopVertexIdToMap[vertex] << std::endl;
			}
		}
	}

	// --------------------------
 // 边界特征弧顶点ID映射（quadBoundaryArcVertexIdToMap）
 // --------------------------
	int boundary_id = 0; // 边界弧ID计数器（独立于id，避免冲突）
	quadBoundaryArcVertexIdToMap.clear(); // 清空边界弧映射表

	// 获取边界特征弧列表（假设通过BE_curve_model的成员函数获取）
	const std::vector<BoundaryFittingArc>& boundaryArcs = getBoundaryFittingArcs();

	// 遍历所有边界特征弧
	for (const auto& arc : boundaryArcs) {
		if (arc.orderedHalfedges.empty()) continue; // 跳过空弧

		std::set<V*> boundaryUniqueVertices; // 临时集合：边界弧顶点去重
		// 提取当前边界弧的所有顶点（遍历orderedHalfedges的起点和终点）
		for (H* he : arc.orderedHalfedges) {
			V* startV = quadMesh->halfedgeSource(he); // 半边起点顶点
			V* endV = quadMesh->halfedgeTarget(he);   // 半边终点顶点
			boundaryUniqueVertices.insert(startV);
			boundaryUniqueVertices.insert(endV);
		}

		// 为边界弧顶点分配独立ID（不与特征环/内部弧冲突）
		for (V* vertex : boundaryUniqueVertices) {
			// 1. 检查是否已在边界弧映射表中（避免重复分配）
			if (quadBoundaryArcVertexIdToMap.find(vertex) != quadBoundaryArcVertexIdToMap.end()) {
				continue;
			}
			// 2. 分配新ID并加入映射表
			quadBoundaryArcVertexIdToMap[vertex] = boundary_id++;
			// 调试输出：边界弧顶点ID分配情况
			/*std::cout << "Boundary Arc Vertex ID: " << quadBoundaryArcVertexIdToMap[vertex]
				<< " (Vertex ptr: " << vertex << ", ID: " << vertex->id() << ")" << std::endl;*/
		}
	}

	// --------------------------
	// 结果统计输出
	// --------------------------
	/*std::cout << "=== Vertex ID Map Summary ===" << std::endl;
	std::cout << "Loop + Internal Arc: " << quadLoopVertexIdToMap.size() << " vertices" << std::endl;
	std::cout << "Boundary Arc: " << quadBoundaryArcVertexIdToMap.size() << " vertices" << std::endl;
	std::cout << "Total unique vertices (all types): "
		<< quadLoopVertexIdToMap.size() + quadBoundaryArcVertexIdToMap.size() << std::endl;*/
}

void CCG_QMSLib::BE_curve_model::parametrize(M* quadMesh) {

	// Get mapping from triangle mesh points to quadrilateral halfedges
	auto& mapping = getTriPointToQuadHalfedgeMap();

	// Used to store groups of quadrilateral halfedges to triangle mesh points
	std::map<H*, std::vector<V*>> halfedgeToTriPoints;

	// Iterate over mapping and group triangle mesh points by corresponding quadrilateral halfedge
	for (auto& pair : mapping) {
		V* triPoint = pair.first;    // triangle mesh point
		H* quadHalfedge = pair.second;  // quadrilateral halfedge

		// Add triangle point to the corresponding quadrilateral halfedge group
		halfedgeToTriPoints[quadHalfedge].push_back(triPoint);
	}

	// Iterate over each quadrilateral halfedge to perform piecewise parameterization
	for (auto& group : halfedgeToTriPoints) {
		H* quadHalfedge = group.first;      // current quadrilateral halfedge
		std::vector<V*>& triPoints = group.second;  // triangle mesh points corresponding to the current halfedge

		// Get start and end vertices of the halfedge
		V* startVertex = quadMesh->halfedgeSource(quadHalfedge);
		V* endVertex = quadMesh->halfedgeTarget(quadHalfedge);

		// Create a new point list in the required order
		std::vector<V*> orderedPoints;
		orderedPoints.push_back(startVertex);  // store start vertex first
		orderedPoints.insert(orderedPoints.end(), triPoints.begin(), triPoints.end());  // then store the original triangle mesh points
		orderedPoints.push_back(endVertex);  // finally store end vertex

		// Sort the new point list by distance from the start vertex
		std::sort(orderedPoints.begin(), orderedPoints.end(), [startVertex](V* a, V* b) {
			return (a->point() - startVertex->point()).norm() <
				(b->point() - startVertex->point()).norm();
			});

		// Replace triPoints with the re-sorted point list (if original data needs modification)
		triPoints = std::move(orderedPoints);

		// Store triangle mesh points and chord lengths for this segment
		std::vector<double> lengths;  // chord lengths between points
		double segmentTotalLength = 0.0;

		// Calculate chord length between each point and the next
		for (size_t i = 1; i < triPoints.size(); ++i) {
			CPoint currPoint = triPoints[i]->point();
			CPoint prevPoint = triPoints[i - 1]->point();
			double distance = (currPoint - prevPoint).norm();  // distance calculation
			lengths.push_back(distance);
			segmentTotalLength += distance;
		}

		// Assign parameter values to each point
		double accumulatedLength = 0.0;
		for (size_t i = 0; i < triPoints.size(); ++i) {
			V* triPoint = triPoints[i];

			// Calculate parameter value
			double paramValue = (i == 0) ? 0.0 : (accumulatedLength / segmentTotalLength);

			// Store into parameterization mapping
			triPointIdToParamMap[triPoint->id()] = paramValue;
			//std::cout << "quadHalfedgeid" << quadHalfedgeToIdMap[quadHalfedge] << "   triPoint->id()" << triPoint->id() << "   paramValue" << paramValue << std::endl;

			// Update accumulated length
			if (i < lengths.size()) {
				accumulatedLength += lengths[i];
			}
		}
	}
}

void CCG_QMSLib::BE_curve_model::parametrize2(M* quadMesh) {

	// Get mapping from data points to quadrilateral halfedges
	auto& mapping = getDataPointToQuadHalfedgeMap();
	// Used to store groups of quadrilateral halfedges to data points
	std::map<H*, std::vector<QB_curve_sampling*>> halfedgeToTriPoints;

	// Iterate over mapping and group data points by corresponding quadrilateral halfedge
	for (auto& pair : mapping) {
		QB_curve_sampling* dataPoint = pair.first;  // data point on triangle mesh
		H* quadHalfedge = pair.second;              // quadrilateral halfedge

		// Add data point to the corresponding quadrilateral halfedge group
		halfedgeToTriPoints[quadHalfedge].push_back(dataPoint);
	}

	// Iterate over each quadrilateral halfedge to perform piecewise parameterization
	for (auto& pair : halfedgeToTriPoints) {
		H* quadHalfedge = pair.first;                     // current quadrilateral halfedge
		std::vector<QB_curve_sampling*>& dataPoints = pair.second;

		// Get start and end vertices of the halfedge
		V* startVertex = quadMesh->halfedgeSource(quadHalfedge);
		V* endVertex = quadMesh->halfedgeTarget(quadHalfedge);
		// Coordinates of start and end points
		CPoint startPoint = startVertex->point();
		CPoint endPoint = endVertex->point();

		// Calculate total length of the halfedge
		double TotalLength = (endPoint - startPoint).norm();

		// Sort data points by distance from start point
		std::sort(dataPoints.begin(), dataPoints.end(), [startPoint](QB_curve_sampling* a, QB_curve_sampling* b) {
			return (a->pos() - startPoint).norm() < (b->pos() - startPoint).norm();
			});

		// Assign parameter values to each data point
		for (QB_curve_sampling* dataPoint : dataPoints) {
			// Distance from the current point to the start point
			double distanceToStart = (dataPoint->pos() - startPoint).norm();
			// Parameter value = distance / total length
			double paramValue = distanceToStart / TotalLength;
			// Store parameterized value into mapping
			dataPointIdToParamMap[dataPoint->vId()] = paramValue;

			// Print debug information
			//std::cout << "quadHalfedgeid" << quadHalfedgeToIdMap[quadHalfedge] << "   dataPoint ID: " << dataPoint->vId()
			//          << ", Param Value: " << paramValue << std::endl;
		}
	}
	//std::cout << "# sampling on boundary: " << dataPointIdToParamMap.size() << std::endl;
}

void CCG_QMSLib::BE_curve_model::parametrize4(M* quadMesh)
{
	// Get mapping from triangle mesh points to quadrilateral halfedges
	auto& mapping = getDataPointToQuadHalfedgeMap();
	std::map<H*, std::vector<QB_curve_sampling*>> halfedgeToTriPoints;

	// Group data points by their associated halfedge
	for (auto& pair : mapping) {
		halfedgeToTriPoints[pair.second].push_back(pair.first);
	}

	// Process each halfedge
	for (auto& pair : halfedgeToTriPoints) {
		H* quadHalfedge = pair.first;
		std::vector<QB_curve_sampling*>& dataPoints = pair.second;

		// Get the start and end vertices of this halfedge
		V* startVertex = quadMesh->halfedgeSource(quadHalfedge);
		V* endVertex = quadMesh->halfedgeTarget(quadHalfedge);
		CPoint startPoint = startVertex->point();
		CPoint endPoint = endVertex->point();

		// 1. First, sort data points by distance from the start (preserving original order)
		std::sort(dataPoints.begin(), dataPoints.end(),
			[startPoint](QB_curve_sampling* a, QB_curve_sampling* b) {
				return (a->pos() - startPoint).norm() < (b->pos() - startPoint).norm();
			});

		// 2. Compute centripetal parameterization
		std::vector<double> params = computeCentripetalParams(startPoint, endPoint, dataPoints);

		// 3. Store the parameterization results
		for (size_t i = 0; i < dataPoints.size(); ++i) {
			dataPointIdToParamMap[dataPoints[i]->vId()] = params[i];

			// Debug output
			/*std::cout << "quadHalfedgeID: " << quadHalfedgeToIdMap[quadHalfedge]
					  << " dataPoint ID: " << dataPoints[i]->vId()
					  << ", Param Value: " << params[i] << std::endl;*/
		}
	}
	//std::cout << "# sampling on boundary: " << dataPointIdToParamMap.size() << std::endl;
}

// Newly added centripetal parameterization computation function
std::vector<double> CCG_QMSLib::BE_curve_model::computeCentripetalParams(
	const CPoint& startPoint,
	const CPoint& endPoint,
	const std::vector<QB_curve_sampling*>& dataPoints)
{
	std::vector<double> params(dataPoints.size() + 2); // including endpoints
	params[0] = 0.0;    // start parameter
	params.back() = 1.0; // end parameter

	// 1. Compute square-root distances between all adjacent points
	double total = 0.0;
	std::vector<double> sqrtDistances;

	// from start point to the first data point
	sqrtDistances.push_back(std::sqrt((dataPoints.front()->pos() - startPoint).norm()));
	total += sqrtDistances.back();

	// distances between intermediate data points
	for (size_t i = 1; i < dataPoints.size(); ++i) {
		double dist = std::sqrt((dataPoints[i]->pos() - dataPoints[i - 1]->pos()).norm());
		sqrtDistances.push_back(dist);
		total += dist;
	}

	// from the last data point to the end point
	sqrtDistances.push_back(std::sqrt((endPoint - dataPoints.back()->pos()).norm()));
	total += sqrtDistances.back();

	// 2. Compute cumulative parameter values
	double accum = 0.0;
	for (size_t i = 0; i < dataPoints.size(); ++i) {
		accum += sqrtDistances[i];
		params[i + 1] = accum / total; // i+1 because params[0] is the start
	}

	// 3. Return parameters corresponding to data points (excluding endpoints)
	return std::vector<double>(params.begin() + 1, params.end() - 1);
}

void CCG_QMSLib::BE_curve_model::parametrize3(M* quadMesh)
{
	const double EPSILON = 1e-10;  // Floating-point error threshold
	// Get mapping from triangle mesh points to quadrilateral halfedges
	auto& mapping = getTriPointToQuadHalfedgeMap();

	// Used to store groups of quadrilateral halfedges to triangle mesh points
	std::map<H*, std::vector<V*>> halfedgeToTriPoints;

	// Iterate over mapping and group triangle mesh points by corresponding quadrilateral halfedge
	for (auto& pair : mapping) {
		V* triPoint = pair.first;      // triangle mesh point
		H* quadHalfedge = pair.second; // quadrilateral halfedge

		// Add triangle point to the corresponding quadrilateral halfedge group
		halfedgeToTriPoints[quadHalfedge].push_back(triPoint);
	}

	// Process each halfedge group
	for (auto& entry : halfedgeToTriPoints) {
		H* quadHalfedge = entry.first;
		std::vector<V*>& triPoints = entry.second;
		V* startVertex = quadMesh->halfedgeSource(quadHalfedge);
		V* endVertex = quadMesh->halfedgeTarget(quadHalfedge);
		CPoint startPoint = startVertex->point();
		CPoint endPoint = endVertex->point();
		CPoint edgeVec = endPoint - startPoint;
		double totalLength = edgeVec.norm();

		// Sort by projection distance for precise ordering
		std::sort(triPoints.begin(), triPoints.end(),
			[&](V* a, V* b) {
				double proj_a = ((a->point() - startPoint) * edgeVec) / totalLength;
				double proj_b = ((b->point() - startPoint) * edgeVec) / totalLength;
				return proj_a < proj_b;
			});

		// Assign parameter values to each triangle point
		for (V* triPoint : triPoints) {
			CPoint vec = triPoint->point() - startPoint;
			double proj = (vec * edgeVec) / (totalLength * totalLength);

			// Handle boundary cases
			if (proj < EPSILON) {
				triPointIdToParamMap[triPoint->id()] = 0.0;
			}
			else if (proj > 1.0 - EPSILON) {
				triPointIdToParamMap[triPoint->id()] = 1.0;
			}
			else {
				triPointIdToParamMap[triPoint->id()] = proj;
			}
		}
	}
}



void CCG_QMSLib::BE_curve_model::compute_quad_vertex(M* quadMesh)
{

	using namespace Eigen;
	std::vector<double> v_x_vec, v_y_vec, v_z_vec;
	auto set_v_vaules = [](std::vector<double>& v_x_vec, std::vector<double>& v_y_vec, std::vector<double>& v_z_vec, double x, double y, double z)
		{
			v_x_vec.push_back(x);
			v_y_vec.push_back(y);
			v_z_vec.push_back(z);
		};
	std::vector<Eigen::Triplet<double>> tripletlists; // Define a container to store all coefficient pairs
	int total_equation_count = 0; // Indicates the equation index

	// Get mapping from quadrilateral halfedge to new ID
	auto& halfedgeToIdMap = getQuadHalfedgeToIdMap();
	int totalsize = 0;

	// Iterate over each feature loop
	for (auto loop : be_quadFeatureLoops->loops()) {

		// Ensure the current loop is valid
		if (!loop) continue;
		/*for (auto quadHalfEdge : loop->halfedges()) {
			std::cout << "Halfedge pointer: " << quadHalfEdge->localId() << std::endl;
		}*/
		// Get all halfedges in the loop
		auto& halfedges = loop->halfedges();
		std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end()); // Store the loop's halfedges into a vector
		//std::cout << "halfedgeList.size()" << halfedgeList.size() << std::endl;
		totalsize += halfedgeList.size();
		//std::cout << "totalsize" << totalsize << std::endl;
		H* startLHalfedge = halfedgeList[0]; // Initial halfedge
		int startLId = halfedgeToIdMap[startLHalfedge];
		// Iterate over each halfedge in the current loop
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i]; // Current halfedge

			// Get the ID of the current halfedge
			int currentId = halfedgeToIdMap[currentHalfedge];
			// Compute IDs of previous and next halfedges
			int prevId = (i == 0) ? (int)(totalsize - 1) : currentId - 1;
			int nextId = (currentId == (int)(totalsize - 1)) ? startLId : currentId + 1;

			// Get previous and next halfedges based on ID
			H* prevHalfedge = nullptr;
			H* nextHalfedge = nullptr;

			for (auto& pair : halfedgeToIdMap) {
				if (pair.second == prevId) prevHalfedge = pair.first;
				if (pair.second == nextId) nextHalfedge = pair.first;
			}

			// If both previous and next halfedges are found, proceed with logic
			if (prevHalfedge && nextHalfedge) {
				/*std::cout << "Current Halfedge ID: " << currentId << std::endl;
				std::cout << "Previous Halfedge ID: " << prevId << std::endl;
				std::cout << "Next Halfedge ID: " << nextId << std::endl;*/
				// Get new IDs for quadrilateral mesh vertices
				V* v0 = quadMesh->halfedgeSource(prevHalfedge);
				V* v1 = quadMesh->halfedgeTarget(prevHalfedge);
				V* v2 = quadMesh->halfedgeSource(nextHalfedge);
				V* v3 = quadMesh->halfedgeTarget(nextHalfedge);

				int v0Id = quadLoopVertexIdToMap[v0];
				int v1Id = quadLoopVertexIdToMap[v1];
				int v2Id = quadLoopVertexIdToMap[v2];
				int v3Id = quadLoopVertexIdToMap[v3];

				// Get triangle feature points associated with the current halfedge
				std::vector<QB_curve_sampling*> relatedTriPoints;
				for (auto& pair : dataPointToQuadHalfedgeMap) {
					if (pair.second == currentHalfedge) {
						relatedTriPoints.push_back(pair.first);
					}
				}

				// Iterate over triangle feature points, get parameter values and process
				for (auto dataPoint : relatedTriPoints) {
					int dataPointId = dataPoint->vId(); // Get the ID of the triangle feature point
					if (dataPointIdToParamMap.find(dataPointId) != dataPointIdToParamMap.end()) {
						double u = dataPointIdToParamMap[dataPointId];
						if (u < 0.0 || u > 1.0) {
							std::cerr << "Warning: u=" << u << " is out of range!" << std::endl;
							u = std::max(0.0, std::min(1.0, u)); // Clamp to [0, 1]
						}
						// Example: output the parameter value of the current triangle feature point
						//std::cout << "Halfedge ID: " << currentId << ", Triangle Point ID: " << triPointId << ", Param Value: " << u << std::endl;

						/*--------------------------Formula start------------------------------------------*/
						// Coefficient for V0
						int temp_index = v0Id;
						tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6);

						// Coefficient for V1
						temp_index = v1Id;
						tripletlists.emplace_back(total_equation_count, temp_index,
							(1 - u) * (1 - u) * (1 - u) * 2 / 3 +
							2 * u * (1 - u) * (1 - u) +
							u * u * (1 - u) +
							u * u * u / 6);

						// Coefficient for V2
						temp_index = v2Id;
						tripletlists.emplace_back(total_equation_count, temp_index,
							(1 - u) * (1 - u) * (1 - u) / 6 +
							u * (1 - u) * (1 - u) +
							u * u * (1 - u) * 2 +
							u * u * u * 2 / 3);

						// Coefficient for V3
						temp_index = v3Id;
						tripletlists.emplace_back(total_equation_count, temp_index, u * u * u / 6);
						/*--------------------------Formula end------------------------------------------*/

						double coeff_v0 = (1 - u) * (1 - u) * (1 - u) / 6;
						double coeff_v1 = (1 - u) * (1 - u) * (1 - u) * 2 / 3 +
							2 * u * (1 - u) * (1 - u) +
							u * u * (1 - u) +
							u * u * u / 6;
						double coeff_v2 = (1 - u) * (1 - u) * (1 - u) / 6 +
							u * (1 - u) * (1 - u) +
							u * u * (1 - u) * 2 +
							u * u * u * 2 / 3;
						double coeff_v3 = u * u * u / 6;
						double coeff_sum = coeff_v0 + coeff_v1 + coeff_v2 + coeff_v3;
						if (std::abs(coeff_sum - 1.0) > 1e-6) {
							//std::cerr << "Warning: B-spline coefficients do not sum to 1! Sum=" << coeff_sum << std::endl;
						}

						CPoint dataPointCoords = dataPoint->pos(); // Ensure triPoint->point() returns CPoint type
						set_v_vaules(v_x_vec, v_y_vec, v_z_vec,
							dataPointCoords[0],
							dataPointCoords[1],
							dataPointCoords[2]);
						total_equation_count++;
					}
				}
			}
		}
	}
	// Iterate over feature arcs (data point map not modified)
	for (auto arc : be_quadFeatureLoops->arcs()) {
		if (!arc) continue;

		// Get all halfedges in the arc
		auto& halfedges = arc->halfedges();
		std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end()); // Store the arc's halfedges into a vector
		totalsize += halfedgeList.size();
		//std::cout << "Arc halfedges count: " << halfedgeList.size() << std::endl;

		H* startAHalfedge = halfedgeList[0];
		int startAId = halfedgeToIdMap[startAHalfedge];
		//std::cout << "startAId" << startAId << std::endl;
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i];
			int currentId = halfedgeToIdMap[currentHalfedge];

			// Special handling for arcs: open paths only need to process valid adjacent edges
			int prevId = (i > 0) ? halfedgeToIdMap[halfedgeList[i - 1]] : -1;
			int nextId = (i < halfedgeList.size() - 1) ? halfedgeToIdMap[halfedgeList[i + 1]] : -1;

			// Get control vertices (same as loop processing)
			V* v0 = (prevId != -1) ? quadMesh->halfedgeSource(halfedgeList[i - 1]) : nullptr;
			V* v1 = quadMesh->halfedgeSource(currentHalfedge);
			V* v2 = quadMesh->halfedgeTarget(currentHalfedge);
			V* v3 = (nextId != -1)
				? quadMesh->halfedgeTarget(halfedgeList[i + 1])
				: nullptr;

			// Boundary protection
			//if (!v1) continue;

			// Get triangle feature points associated with the current halfedge
			std::vector<V*> relatedTriPoints;
			for (auto& pair : triPointToQuadHalfedgeMap) {
				if (pair.second == currentHalfedge) {
					relatedTriPoints.push_back(pair.first);
				}
			}

			// Process each feature point
			for (auto triPoint : relatedTriPoints) {
				int triPointId = triPoint->id();
				if (triPointIdToParamMap.find(triPointId) != triPointIdToParamMap.end()) {
					double u = triPointIdToParamMap[triPointId];
					CPoint triPointCoords = triPoint->point();

					// Special coefficient handling for arcs (using cubic B-splines)
					if (v0 && v3) { // Complete four-vertex case
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v0], (1 - u) * (1 - u) * (1 - u) / 6);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v1],
							(1 - u) * (1 - u) * (1 - u) * 2 / 3 +
							2 * u * (1 - u) * (1 - u) +
							u * u * (1 - u) +
							u * u * u / 6);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v2],
							(1 - u) * (1 - u) * (1 - u) / 6 +
							u * (1 - u) * (1 - u) +
							u * u * (1 - u) * 2 +
							u * u * u * 2 / 3);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v3], u * u * u / 6);
					}
					//else { // Boundary case (simplified handling)
					//    tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v1], 1.0);
					//}
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec,
						triPointCoords[0],
						triPointCoords[1],
						triPointCoords[2]);
					// Add endpoint constraints in the arc processing
					if (i == 0 || i == halfedgeList.size() - 1) { // First and last vertices
						int vid = quadLoopVertexIdToMap[v1];
						tripletlists.emplace_back(total_equation_count, vid, 1.0);
						set_v_vaules(v_x_vec, v_y_vec, v_z_vec,
							v1->point()[0],
							v1->point()[1],
							v1->point()[2]);
						total_equation_count++;
						//continue; // Skip B-spline computation
					}
					total_equation_count++;
				}
			}
		}
	}

	// Step 1: Build sparse matrix G and the right-hand vector g
	int G_row_num = total_equation_count;  // Total number of equations
	int G_col_num = quadLoopVertexIdToMap.size();     // Number of control points
	SparseMatrix<double> G(G_row_num, G_col_num);
	// Construct sparse matrix G from triplet list
	G.setFromTriplets(tripletlists.begin(), tripletlists.end());

	// Check equation count and number of control points
	//std::cout << "Total equations: " << total_equation_count << std::endl;
	//std::cout << "Control points: " << quadLoopVertexIdToMap.size() << std::endl;

	// After constructing all equations, convert std::vector to Eigen::VectorXd
	VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

	// Check data consistency
	//std::cout << "v_x_vec size: " << v_x_vec.size() << std::endl;
	//std::cout << "v_y_vec size: " << v_y_vec.size() << std::endl;
	//std::cout << "v_z_vec size: " << v_z_vec.size() << std::endl;
	//std::cout << "G rows: " << G.rows() << std::endl; // Rows

	// Add diagnostics before solving: output rank and column count
	//std::cout << "Rank of G: " << JacobiSVD<MatrixXd>(G).rank() << "/" << G.cols() << std::endl;
	// Step 2: Solve using least squares
	MatrixXd G_dense = G.toDense();

	// Compute the condition number of matrix G
	JacobiSVD<MatrixXd> svd(G_dense);
	double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//std::cout << "Condition number of G: " << cond << std::endl;

	// Solve using regularization method
	double lambda = 1e-6; // Regularization parameter
	MatrixXd G_reg = G_dense.transpose() * G_dense + lambda * MatrixXd::Identity(G_col_num, G_col_num);
	VectorXd v_x_reg = G_dense.transpose() * v_x;
	VectorXd v_y_reg = G_dense.transpose() * v_y;
	VectorXd v_z_reg = G_dense.transpose() * v_z;
	// Replace original solving approach

	VectorXd solution_x = G_reg.ldlt().solve(v_x_reg);
	VectorXd solution_y = G_reg.ldlt().solve(v_y_reg);
	VectorXd solution_z = G_reg.ldlt().solve(v_z_reg);

	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//solver.compute(G);  // G is sparse matrix
	//if (solver.info() != Eigen::Success) {
	//    std::cerr << "Decomposition failed!" << std::endl;
	//    return;
	//}
	//VectorXd solution_x = solver.solve(v_x);
	//VectorXd solution_y = solver.solve(v_y);
	//VectorXd solution_z = solver.solve(v_z);

	// Print non-zero elements after building the matrix
	//for (int k = 0; k < G.outerSize(); ++k) {
	//    for (SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
	//        if (it.row() >= total_equation_count - 10) { // Check last 10 equations
	//            std::cout << "Row " << it.row() << ", Col " << it.col()
	//                      << ", Value " << it.value() << std::endl;
	//        }
	//    }
	//}

	// Update quadrilateral mesh vertices
	for (auto& pair : quadLoopVertexIdToMap) {
		V* quadVertex = pair.first;
		int quadVertexId = pair.second;
		double new_x = solution_x[quadVertexId];
		double new_y = solution_y[quadVertexId];
		double new_z = solution_z[quadVertexId];
		CPoint& vertexPoint = quadVertex->point(); // Get a reference to the vertex
		vertexPoint[0] = new_x;
		vertexPoint[1] = new_y;
		vertexPoint[2] = new_z;
		/*std::cout << "111Vertex ID: " << quadVertexId
				  << ", Solution X: " << quadVertex->point()[0]
				  << ", Solution Y: " << quadVertex->point()[1]
				  << ", Solution Z: " << quadVertex->point()[2] << std::endl;*/
	}
	//std::cout << "Updated quad mesh vertices successfully." << std::endl;

	//for (auto arc : be_quadFeatureLoops->arcs()) {
	//    for (auto he : arc->halfedges()) {
	//        V* v = quadMesh->halfedgeTarget(he);
	//        if (quadLoopVertexIdToMap.find(v) == quadLoopVertexIdToMap.end()) {
	//            std::cerr << "Error: Arc vertex " << v << " not mapped!" << std::endl;
	//        }
	//        else
	//            std::cout << "Mapping correct" << std::endl;
	//    }
	//}

	//    // Step 1: Build sparse matrix G and the right-hand vector g
	//    int G_row_num = total_equation_count;  // Total number of equations
	//    int G_col_num = quadLoopVertexIdToMap.size();     // Number of control points
	//    SparseMatrix<double> G(G_row_num, G_col_num);
	//
	//    // Construct sparse matrix G from triplet list
	//    G.setFromTriplets(tripletlists.begin(), tripletlists.end());
	//    for (int i = 0; i < G.rows(); i++)
	//    {
	//        //std::cout << "G.row(i).sum(): " << G.row(i).sum() << std::endl;
	//        if (G.row(i).sum() < 0.9 || G.row(i).sum() > 1.1)
	//        {
	//            std::cout << "G.row(i).sum(): " << G.row(i).sum() << std::endl;
	//        }
	//    }
	//    // After constructing all equations, convert std::vector to Eigen::VectorXd
	//    VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	//    VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	//    VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());
	//
	//    std::cout << "Matrix G size: " << G.rows() << "x" << G.cols() << std::endl;
	//    std::cout << "v_x size: " << v_x.size() << std::endl;
	//
	//    // Step 2: Solve using least squares
	//    // Convert sparse matrix G to dense for SVD
	//    MatrixXd G_dense = G.toDense();
	//    for (int i = 0; i < G_dense.rows(); i++)
	//    {
	//        //std::cout << "G.row(i).sum(): " << G_dense.row(i).sum() << std::endl;
	//        if (G_dense.row(i).sum() < 0.9 || G_dense.row(i).sum() > 1.1)
	//        {
	//            std::cout << "G_dense.row(i).sum(): " << G_dense.row(i).sum() << std::endl;
	//        }
	//    }
	//    // Create an SVD solver
	//    JacobiSVD<MatrixXd> svd(G_dense, ComputeThinU | ComputeThinV);
	//
	//    // Compute pseudo-inverse
	//    MatrixXd G_pinv = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();
	//    // Solve for x, y, z respectively
	//    VectorXd solution_x = G_pinv * v_x;
	//    VectorXd solution_y = G_pinv * v_y;
	//    VectorXd solution_z = G_pinv * v_z;
	//    for (int i = 0; i < solution_x.size(); ++i) {
	//        std::cout << "Vertex ID: " << i
	//                  << ", Solution X: " << solution_x[i]
	//                  << ", Solution Y: " << solution_y[i]
	//                  << ", Solution Z: " << solution_z[i] << std::endl;
	//    }
	//
	//    // Update quadrilateral mesh vertices
	//    for (auto& pair : quadLoopVertexIdToMap) {
	//        V* quadVertex = pair.first;
	//        int quadVertexId = pair.second;
	//        double new_x = solution_x[quadVertexId];
	//        double new_y = solution_y[quadVertexId];
	//        double new_z = solution_z[quadVertexId];
	//        std::cout << "Vertex ID: " << quadVertexId
	//                  << ", Solution X: " << new_x
	//                  << ", Solution Y: " << new_y
	//                  << ", Solution Z: " << new_z << std::endl;
	//        CPoint& vertexPoint = quadVertex->point(); // Get a reference to the vertex
	//        vertexPoint[0] = new_x;
	//        vertexPoint[1] = new_y;
	//        vertexPoint[2] = new_z;
	//        std::cout << "111Vertex ID: " << quadVertexId
	//                  << ", Solution X: " << quadVertex->point()[0]
	//                  << ", Solution Y: " << quadVertex->point()[1]
	//                  << ", Solution Z: " << quadVertex->point()[2] << std::endl;
	//    }
	//    std::cout << "Updated quad mesh vertices successfully." << std::endl;
	//
} 
// Add smoothness constraints
void CCG_QMSLib::BE_curve_model::compute_quad_vertex1(M* quadMesh)
{
	using namespace Eigen;
	std::vector<double> v_x_vec, v_y_vec, v_z_vec;
	auto set_v_vaules = [](std::vector<double>& v_x_vec, std::vector<double>& v_y_vec, std::vector<double>& v_z_vec, double x, double y, double z)
		{
			v_x_vec.push_back(x);
			v_y_vec.push_back(y);
			v_z_vec.push_back(z);
		};
	std::vector<Eigen::Triplet<double>> tripletlists; // Container for storing all coefficient pairs
	int total_equation_count = 0; // Indicates which equation number

	// Get mapping from quadrilateral halfedges to new IDs
	auto& halfedgeToIdMap = getQuadHalfedgeToIdMap();
	int totalsize = 0;

	// Iterate over each feature loop
	for (auto loop : be_quadFeatureLoops->loops()) {

		// Ensure the current loop is valid
		if (!loop) continue;
		/*for (auto quadHalfEdge : loop->halfedges()) {
			std::cout << "Halfedge pointer: " << quadHalfEdge->localId() << std::endl;
		}*/
		// Get all halfedges in the loop
		auto& halfedges = loop->halfedges();
		std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end()); // Store loop halfedges in a vector
		//std::cout << "halfedgeList.size()" << halfedgeList.size() << std::endl;
		totalsize += halfedgeList.size();
		//std::cout << "totalsize" << totalsize << std::endl;
		H* startLHalfedge = halfedgeList[0]; // Initial halfedge
		int startLId = halfedgeToIdMap[startLHalfedge];
		// Iterate over each halfedge in the current loop
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i]; // Current halfedge

			// Get the ID of the current halfedge
			int currentId = halfedgeToIdMap[currentHalfedge];
			// Compute IDs for the previous and next halfedges
			int prevId = (i == 0) ? (int)(totalsize - 1) : currentId - 1;
			int nextId = (currentId == (int)(totalsize - 1)) ? startLId : currentId + 1;

			// Get the previous and next halfedges based on ID
			H* prevHalfedge = nullptr;
			H* nextHalfedge = nullptr;

			for (auto& pair : halfedgeToIdMap) {
				if (pair.second == prevId) prevHalfedge = pair.first;
				if (pair.second == nextId) nextHalfedge = pair.first;
			}

			// If both previous and next halfedges are found, proceed with further logic
			if (prevHalfedge && nextHalfedge) {
				/*std::cout << "Current Halfedge ID: " << currentId << std::endl;
				std::cout << "Previous Halfedge ID: " << prevId << std::endl;
				std::cout << "Next Halfedge ID: " << nextId << std::endl;*/
				// Get new IDs for quadrilateral mesh vertices
				V* v0 = quadMesh->halfedgeSource(prevHalfedge);
				V* v1 = quadMesh->halfedgeTarget(prevHalfedge);
				V* v2 = quadMesh->halfedgeSource(nextHalfedge);
				V* v3 = quadMesh->halfedgeTarget(nextHalfedge);

				int v0Id = quadLoopVertexIdToMap[v0];
				int v1Id = quadLoopVertexIdToMap[v1];
				int v2Id = quadLoopVertexIdToMap[v2];
				int v3Id = quadLoopVertexIdToMap[v3];

				// Get triangle feature points associated with the current halfedge
				std::vector<QB_curve_sampling*> relatedTriPoints;
				for (auto& pair : dataPointToQuadHalfedgeMap) {
					if (pair.second == currentHalfedge) {
						relatedTriPoints.push_back(pair.first);
					}
				}

				// Iterate over triangle feature points, retrieve parameter values, and process
				for (auto dataPoint : relatedTriPoints) {
					int dataPointId = dataPoint->vId(); // Get the ID of the triangle feature point
					if (dataPointIdToParamMap.find(dataPointId) != dataPointIdToParamMap.end()) {
						double u = dataPointIdToParamMap[dataPointId];
						if (u < 0.0 || u > 1.0) {
							std::cerr << "Warning: u=" << u << " is out of range!" << std::endl;
							u = std::max(0.0, std::min(1.0, u)); // Clamp to [0, 1]
						}
						// Example: print the parameter value of the current triangle feature point
						//std::cout << "Halfedge ID: " << currentId << ", Triangle Point ID: " << triPointId << ", Param Value: " << u << std::endl;
						/*--------------------------Start of formula------------------------------------------*/
						// Coefficient for V0
						int temp_index = v0Id;
						tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6);

						// Coefficient for V1
						temp_index = v1Id;
						tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6);

						// Coefficient for V2
						temp_index = v2Id;
						tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3);

						// Coefficient for V3
						temp_index = v3Id;
						tripletlists.emplace_back(total_equation_count, temp_index, u * u * u / 6);
						/*--------------------------End of formula------------------------------------------*/

						double coeff_v0 = (1 - u) * (1 - u) * (1 - u) / 6;
						double coeff_v1 = (1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6;
						double coeff_v2 = (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3;
						double coeff_v3 = u * u * u / 6;
						double coeff_sum = coeff_v0 + coeff_v1 + coeff_v2 + coeff_v3;
						if (std::abs(coeff_sum - 1.0) > 1e-6) {
							//std::cerr << "Warning: B-spline coefficients do not sum to 1! Sum=" << coeff_sum << std::endl;
						}
						if (u < 0.0 || u > 1.0) {
							//std::cerr << "Warning: u=" << u << " is out of range!" << std::endl;
							u = std::max(0.0, std::min(1.0, u)); // Clamp to [0, 1]
						}
						CPoint dataPointCoords = dataPoint->pos(); // Ensure triPoint->point() returns CPoint type
						set_v_vaules(v_x_vec, v_y_vec, v_z_vec, dataPointCoords[0], dataPointCoords[1], dataPointCoords[2]);
						total_equation_count++;

						double lambda1 = 1; // Weight for smoothing constraint
						/*--------------------------Start of smoothing constraint formula------------------------------------------*/
						// Coefficient for V0
						temp_index = v0Id;
						tripletlists.emplace_back(total_equation_count, temp_index, 4.0 * lambda1);

						// Coefficient for V1
						temp_index = v1Id;
						tripletlists.emplace_back(total_equation_count, temp_index, -6.0 * lambda1);

						//// Coefficient for V2
						//temp_index = v2Id;
						//tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3);

						// Coefficient for V3
						temp_index = v3Id;
						tripletlists.emplace_back(total_equation_count, temp_index, 2.0 * lambda1);
						// Target value is 0 (minimize curvature)
						v_x_vec.push_back(0.0);
						v_y_vec.push_back(0.0);
						v_z_vec.push_back(0.0);
						total_equation_count++;

						// Coefficient for V0
						temp_index = v0Id;
						tripletlists.emplace_back(total_equation_count, temp_index, -6.0 * lambda1);

						// Coefficient for V1
						temp_index = v1Id;
						tripletlists.emplace_back(total_equation_count, temp_index, 12.0 * lambda1);

						// Coefficient for V2
						temp_index = v2Id;
						tripletlists.emplace_back(total_equation_count, temp_index, -6.0 * lambda1);

						//// Coefficient for V3
						//temp_index = v3Id;
						//tripletlists.emplace_back(total_equation_count, temp_index, 2.0 * lambda1);
						// Target value is 0 (minimize curvature)
						v_x_vec.push_back(0.0);
						v_y_vec.push_back(0.0);
						v_z_vec.push_back(0.0);
						total_equation_count++;

						//// Coefficient for V0
						//int temp_index = v0Id;
						//tripletlists.emplace_back(total_equation_count, temp_index, -1.0 * lambda1);

						// Coefficient for V1
						temp_index = v1Id;
						tripletlists.emplace_back(total_equation_count, temp_index, -6.0 * lambda1);

						// Coefficient for V2
						temp_index = v2Id;
						tripletlists.emplace_back(total_equation_count, temp_index, 12.0 * lambda1);

						// Coefficient for V3
						temp_index = v3Id;
						tripletlists.emplace_back(total_equation_count, temp_index, -6.0 * lambda1);
						// Target value is 0 (minimize curvature)
						v_x_vec.push_back(0.0);
						v_y_vec.push_back(0.0);
						v_z_vec.push_back(0.0);
						total_equation_count++;

						// Coefficient for V0
						temp_index = v0Id;
						tripletlists.emplace_back(total_equation_count, temp_index, 2.0 * lambda1);

						//// Coefficient for V1
						//temp_index = v1Id;
						//tripletlists.emplace_back(total_equation_count, temp_index, 2.0 * lambda1);

						// Coefficient for V2
						temp_index = v2Id;
						tripletlists.emplace_back(total_equation_count, temp_index, -6.0 * lambda1);

						// Coefficient for V3
						temp_index = v3Id;
						tripletlists.emplace_back(total_equation_count, temp_index, 4.0 * lambda1);
						// Target value is 0 (minimize curvature)
						v_x_vec.push_back(0.0);
						v_y_vec.push_back(0.0);
						v_z_vec.push_back(0.0);
						total_equation_count++;
						/*--------------------------End of smoothing constraint formula------------------------------------------*/
					}
				}
			}
		}
	}

	// Iterate over feature arcs (data point map unchanged)
	for (auto arc : be_quadFeatureLoops->arcs()) {
		if (!arc) continue;

		// Get all halfedges in the arc
		auto& halfedges = arc->halfedges();
		std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end());
		totalsize += halfedgeList.size();
		//std::cout << "Arc halfedges count: " << halfedgeList.size() << std::endl;

		H* startAHalfedge = halfedgeList[0];
		int startAId = halfedgeToIdMap[startAHalfedge];
		//std::cout << "startAId: " << startAId << std::endl;
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i];
			int currentId = halfedgeToIdMap[currentHalfedge];

			// Special handling for arcs: open paths only need to process valid adjacent edges
			int prevId = (i > 0) ? halfedgeToIdMap[halfedgeList[i - 1]] : -1;
			int nextId = (i < halfedgeList.size() - 1) ? halfedgeToIdMap[halfedgeList[i + 1]] : -1;

			// Get control vertices (same as loop processing)
			V* v0 = (prevId != -1) ? quadMesh->halfedgeSource(halfedgeList[i - 1]) : nullptr;
			V* v1 = quadMesh->halfedgeSource(currentHalfedge);
			V* v2 = quadMesh->halfedgeTarget(currentHalfedge);
			V* v3 = (nextId != -1) ? quadMesh->halfedgeTarget(halfedgeList[i + 1]) : nullptr;

			// Boundary protection
			//if (!v1) continue;

			// Get associated triangle feature points
			std::vector<V*> relatedTriPoints;
			for (auto& pair : triPointToQuadHalfedgeMap) {
				if (pair.second == currentHalfedge) {
					relatedTriPoints.push_back(pair.first);
				}
			}

			// Process each feature point
			for (auto triPoint : relatedTriPoints) {
				int triPointId = triPoint->id();
				if (triPointIdToParamMap.find(triPointId) != triPointIdToParamMap.end()) {
					double u = triPointIdToParamMap[triPointId];
					CPoint triPointCoords = triPoint->point();

					// Special coefficient handling for arcs (using cubic B-spline)
					if (v0 && v3) { // Full four-vertex case
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v0], (1 - u) * (1 - u) * (1 - u) / 6);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v1],
							(1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v2],
							(1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v3], u * u * u / 6);
					}
					//else { // Boundary case (simplified)
					//    tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v1], 1.0);
					//}
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, triPointCoords[0], triPointCoords[1], triPointCoords[2]);
					// Add endpoint constraints in arc processing
					if (i == 0 || i == halfedgeList.size() - 1) { // First or last vertex
						int vid = quadLoopVertexIdToMap[v1];
						tripletlists.emplace_back(total_equation_count, vid, 1.0);
						set_v_vaules(v_x_vec, v_y_vec, v_z_vec, v1->point()[0], v1->point()[1], v1->point()[2]);
						total_equation_count++;
						//continue; // Skip B-spline calculation
					}
					total_equation_count++;
				}
			}
		}
	}

	// Step 1: Build sparse matrix G and right-hand side vectors g
	int G_row_num = total_equation_count;  // Total number of equations
	int G_col_num = quadLoopVertexIdToMap.size();     // Number of control points
	SparseMatrix<double> G(G_row_num, G_col_num);
	// Construct sparse matrix G from triplet list
	G.setFromTriplets(tripletlists.begin(), tripletlists.end());

	// Check equation count and control point count
	/*std::cout << "Total equations: " << total_equation_count << std::endl;
	std::cout << "Control points: " << quadLoopVertexIdToMap.size() << std::endl;*/

	// After constructing all equations, convert std::vector to Eigen::VectorXd
	VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

	// Check data consistency
	//std::cout << "v_x_vec size: " << v_x_vec.size() << std::endl;
	//std::cout << "v_y_vec size: " << v_y_vec.size() << std::endl;
	//std::cout << "v_z_vec size: " << v_z_vec.size() << std::endl;
	//std::cout << "G rows: " << G.rows() << std::endl; // Number of rows
	// Diagnostics before solving: output rank and column count
	//std::cout << "Rank of G: " << JacobiSVD<MatrixXd>(G).rank() << "/" << G.cols() << std::endl;
	// Step 2: Solve using least squares
	MatrixXd G_dense = G.toDense();

	// Compute condition number of matrix G
	JacobiSVD<MatrixXd> svd(G_dense);
	double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//std::cout << "Condition number of G: " << cond << std::endl;

	// Solve using regularization
	double lambda = 1e-4; // Regularization parameter

	//------------------------------------------------------------------------------------
	//const int num_integration_points = 7; // Increase integration points
	//double lambda1 = 1.0;
	//// 2. 7-point Gaussian integration
	//const std::vector<double> gauss_points = {
	//    0.025446, 0.129234, 0.297077, 0.5,
	//    0.702923, 0.870766, 0.974554
	//};
	//const std::vector<double> gauss_weights = {
	//    0.064742, 0.139853, 0.190915, 0.208979,
	//    0.190915, 0.139853, 0.064742
	//};

	//// Iterate over each feature loop
	//for (auto loop : be_quadFeatureLoops->loops()) {

	//    // Ensure the current loop is valid
	//    if (!loop) continue;
	//    /*for (auto quadHalfEdge : loop->halfedges()) {
	//        std::cout << "Halfedge pointer: " << quadHalfEdge->localId() << std::endl;
	//    }*/
	//    // Get all halfedges in the loop
	//    auto& halfedges = loop->halfedges();
	//    std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end()); // Store loop halfedges in a vector
	//    //std::cout << "halfedgeList.size(): " << halfedgeList.size() << std::endl;
	//    totalsize += halfedgeList.size();
	//    std::cout << "totalsize: " << totalsize << std::endl;
	//    H* startLHalfedge = halfedgeList[0]; // Initial halfedge
	//    int startLId = halfedgeToIdMap[startLHalfedge];
	//    // Iterate over each halfedge in the current loop
	//    for (size_t i = 0; i < halfedgeList.size(); ++i) {
	//        H* currentHalfedge = halfedgeList[i]; // Current halfedge

	//        // Get the ID of the current halfedge
	//        int currentId = halfedgeToIdMap[currentHalfedge];
	//        // Compute IDs for previous and next halfedges
	//        int prevId = (i == 0) ? (int)(totalsize - 1) : currentId - 1;
	//        int nextId = (currentId == (int)(totalsize - 1)) ? startLId : currentId + 1;

	//        // Get previous and next halfedges by ID
	//        H* prevHalfedge = nullptr;
	//        H* nextHalfedge = nullptr;

	//        for (auto& pair : halfedgeToIdMap) {
	//            if (pair.second == prevId) prevHalfedge = pair.first;
	//            if (pair.second == nextId) nextHalfedge = pair.first;
	//        }
	//        // Retrieve vertex IDs
	//        int vid[4];
	//        // If both prev and next halfedges found, proceed
	//        if (prevHalfedge && nextHalfedge) {
	//            /*std::cout << "Current Halfedge ID: " << currentId << std::endl;
	//            std::cout << "Previous Halfedge ID: " << prevId << std::endl;
	//            std::cout << "Next Halfedge ID: " << nextId << std::endl;*/
	//            // Get new IDs for quadrilateral mesh vertices
	//            // Retrieve 4 control points (needed for cubic B-spline)
	//            V* v[4];
	//            v[0] = quadMesh->halfedgeSource(prevHalfedge);
	//            v[1] = quadMesh->halfedgeTarget(prevHalfedge);
	//            v[2] = quadMesh->halfedgeSource(nextHalfedge);
	//            v[3] = quadMesh->halfedgeTarget(nextHalfedge);
	//            for (int i = 0; i < 4; ++i) vid[i] = quadLoopVertexIdToMap[v[i]];
	//            std::cout << "vid[0]: " << vid[0] << std::endl;
	//        }

	//        double len01 = std::sqrt(std::pow(v[0][0] - p1[0], 2) + std::pow(p0[1] - p1[1], 2) + std::pow(p0[2] - p1[2], 2));
	//        double len12 = std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1], 2) + std::pow(p1[2] - p2[2], 2));
	//        double len23 = std::sqrt(std::pow(p2[0] - p3[0], 2) + std::pow(p2[1] - p3[1], 2) + std::pow(p2[2] - p3[2], 2));

	//        // Compute normalization factor (based on the inverse cube of average length)
	//        double avg_length = (len01 + len12 + len23) / 3.0;
	//        double normalization = 1.0 / (avg_length * avg_length * avg_length);

	//        for (int k = 0; k < num_integration_points; ++k) {
	//            double u = gauss_points[k];
	//            double w = gauss_weights[k] * normalization;
	//            // 4. Numeric integration to compute curvature energy
	//            for (int k = 0; k < num_integration_points; ++k) {
	//                double u = gauss_points[k];
	//                double w = gauss_weights[k];

	//                // Compute second derivatives of cubic B-spline basis functions
	//                double N_pp[4];
	//                N_pp[0] = (1 - u);
	//                N_pp[1] = (3 * u - 2);
	//                N_pp[2] = (-3 * u + 1);
	//                N_pp[3] = u;

	//                // Construct curvature penalty term (λ * N''(u)^T * N''(u))
	//                for (int i = 0; i < 4; ++i) {
	//                    for (int j = 0; j < 4; ++j) {
	//                        double val = lambda1 * w * N_pp[i] * N_pp[j];
	//                        tripletlists.emplace_back(total_equation_count, vid[i], val);
	//                        tripletlists.emplace_back(total_equation_count + 1, vid[i], val);
	//                        tripletlists.emplace_back(total_equation_count + 2, vid[i], val);
	//                        std::cout << "----------------------------------------------" << std::endl;
	//                    }
	//                }

	//                // Target value is 0 (minimize curvature)
	//                v_x_vec.push_back(0);
	//                v_y_vec.push_back(0);
	//                v_z_vec.push_back(0);
	//            }
	//            total_equation_count += 3; // x, y, z components
	//        }
	//    }
	//}
	//------------------------------------------------------------------------------------

	MatrixXd G_reg = G_dense.transpose() * G_dense + lambda * MatrixXd::Identity(G_col_num, G_col_num);
	VectorXd v_x_reg = G_dense.transpose() * v_x;
	VectorXd v_y_reg = G_dense.transpose() * v_y;
	VectorXd v_z_reg = G_dense.transpose() * v_z;
	// Replace original solving method

	VectorXd solution_x = G_reg.ldlt().solve(v_x_reg);
	VectorXd solution_y = G_reg.ldlt().solve(v_y_reg);
	VectorXd solution_z = G_reg.ldlt().solve(v_z_reg);

	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//solver.compute(G);  // G is a sparse matrix
	//if (solver.info() != Eigen::Success) {
	//    std::cerr << "Decomposition failed!" << std::endl;
	//    return;
	//}
	//VectorXd solution_x = solver.solve(v_x);
	//VectorXd solution_y = solver.solve(v_y);
	//VectorXd solution_z = solver.solve(v_z);

	// Print non-zero elements after constructing the matrix
	//for (int k = 0; k < G.outerSize(); ++k) {
	//    for (SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
	//        if (it.row() >= total_equation_count - 10) { // Check last 10 equations
	//            std::cout << "Row " << it.row() << ", Col " << it.col()
	//                << ", Value " << it.value() << std::endl;
	//        }
	//    }
	//}

	// Print solution results
	/*for (int i = 0; i < solution_x.size(); ++i) {
		std::cout << "Vertex ID: " << i
			<< ", Solution X: " << solution_x[i]
			<< ", Solution Y: " << solution_y[i]
			<< ", Solution Z: " << solution_z[i] << std::endl;
	}*/
	// Update quadrilateral mesh vertices
	for (auto& pair : quadLoopVertexIdToMap) {
		V* quadVertex = pair.first;
		int quadVertexId = pair.second;
		double new_x = solution_x[quadVertexId];
		double new_y = solution_y[quadVertexId];
		double new_z = solution_z[quadVertexId];
		CPoint& vertexPoint = quadVertex->point(); // Get reference to vertex
		vertexPoint[0] = new_x;
		vertexPoint[1] = new_y;
		vertexPoint[2] = new_z;
		/*std::cout << "111Vertex ID: " << quadVertexId
			<< ", Solution X: " << quadVertex->point()[0]
			<< ", Solution Y: " << quadVertex->point()[1]
			<< ", Solution Z: " << quadVertex->point()[2] << std::endl;*/
	}
	//std::cout << "Updated quad mesh vertices successfully." << std::endl;
	// Add checks in arc processing
	/*for (auto arc : be_quadFeatureLoops->arcs()) {
		for (auto he : arc->halfedges()) {
			V* v = quadMesh->halfedgeTarget(he);
			if (quadLoopVertexIdToMap.find(v) == quadLoopVertexIdToMap.end()) {
				std::cerr << "Error: Arc vertex " << v << " not mapped!" << std::endl;
			}
			else
				std::cout << "Mapping correct" << std::endl;
		}
	}*/
	// Step 1: build sparse matrix G and right-hand side vector g
	//int G_row_num = total_equation_count;  // Total number of equations
	//int G_col_num = quadLoopVertexIdToMap.size();     // Number of control points
	//SparseMatrix<double> G(G_row_num, G_col_num);
	//
	// Construct sparse matrix G from tripletlist
	//G.setFromTriplets(tripletlists.begin(), tripletlists.end());
	//for (int i = 0; i < G.rows(); i++)
	//{
	//    //std::cout << "G.row(i).sum(): " << G.row(i).sum() << std::endl;
	//    if (G.row(i).sum() < 0.9 || G.row(i).sum() > 1.1)
	//    {
	//        std::cout << "G.row(i).sum(): " << G.row(i).sum() << std::endl;
	//    }
	//}
	// After constructing all equations, convert std::vector to Eigen::VectorXd
	//VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	//VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	//VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());
	//
	//std::cout << "Matrix G size: " << G.rows() << "x" << G.cols() << std::endl;
	//std::cout << "v_x size: " << v_x.size() << std::endl;
	//
	// Step 2: Solve using least squares
	// Convert sparse matrix G to dense to use SVD
	//MatrixXd G_dense = G.toDense();
	//for (int i = 0; i < G_dense.rows(); i++)
	//{
	//    //std::cout << "G_dense.row(i).sum(): " << G_dense.row(i).sum() << std::endl;
	//    if (G_dense.row(i).sum() < 0.9 || G_dense.row(i).sum() > 1.1)
	//    {
	//        std::cout << "G_dense.row(i).sum(): " << G_dense.row(i).sum() << std::endl;
	//    }
	//}
	// Create an SVD solver
	//JacobiSVD<MatrixXd> svd(G_dense, ComputeThinU | ComputeThinV);
	//
	// Compute pseudo-inverse
	//MatrixXd G_pinv = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();
	// Solve x, y, z separately
	//VectorXd solution_x = G_pinv * v_x;
	//VectorXd solution_y = G_pinv * v_y;
	//VectorXd solution_z = G_pinv * v_z;
	//for (int i = 0; i < solution_x.size(); ++i) {
	//    std::cout << "Vertex ID: " << i
	//        << ", Solution X: " << solution_x[i]
	//        << ", Solution Y: " << solution_y[i]
	//        << ", Solution Z: " << solution_z[i] << std::endl;
	//}
	//
	// Update quadrilateral mesh vertices
	//for (auto& pair : quadLoopVertexIdToMap) {
	//    V* quadVertex = pair.first;
	//    int quadVertexId = pair.second;
	//    double new_x = solution_x[quadVertexId];
	//    double new_y = solution_y[quadVertexId];
	//    double new_z = solution_z[quadVertexId];
	//    std::cout << "Vertex ID: " << quadVertexId
	//        << ", Solution X: " << new_x
	//        << ", Solution Y: " << new_y
	//        << ", Solution Z: " << new_z << std::endl;
	//    CPoint& vertexPoint = quadVertex->point(); // Get reference to vertex
	//    vertexPoint[0] = new_x;
	//    vertexPoint[1] = new_y;
	//    vertexPoint[2] = new_z;
	//    std::cout << "111Vertex ID: " << quadVertexId
	//        << ", Solution X: " << quadVertex->point()[0]
	//        << ", Solution Y: " << quadVertex->point()[1]
	//        << ", Solution Z: " << quadVertex->point()[2] << std::endl;
	//}
	//
	//std::cout << "Updated quad mesh vertices successfully." << std::endl;
	//
	//
}
//// Improved version: construct subsystems by loop/arc and solve in parallel (local solve + OpenMP)
// Each loop independently constructs its own G_local and v_local, then solves in parallel and merges results

void CCG_QMSLib::BE_curve_model::compute_quad_vertex2(M* quadMesh) {
	using namespace Eigen;
	using namespace std;

	// Total number of control points
	int total_control_points = quadLoopVertexIdToMap.size();

	// Final solution results
	VectorXd solution_x = VectorXd::Zero(total_control_points);
	VectorXd solution_y = VectorXd::Zero(total_control_points);
	VectorXd solution_z = VectorXd::Zero(total_control_points);

	// Number of times each control point is updated, used for averaging
	vector<int> control_count(total_control_points, 0);

	// Process all loops together (parallel)
	auto& loops = be_quadFeatureLoops->loops();
#pragma omp parallel for
	for (int loopIdx = 0; loopIdx < (int)loops.size(); ++loopIdx) {
		auto loop = loops[loopIdx];
		if (!loop) continue;

		auto& halfedges = loop->halfedges();
		vector<H*> halfedgeList(halfedges.begin(), halfedges.end());
		vector<Triplet<double>> triplets_local;
		vector<double> v_x_vec, v_y_vec, v_z_vec;
		std::map<V*, int> localIdMap;
		vector<V*> localVertices;
		int eq_id = 0;

		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i];
			H* prevHalfedge = halfedgeList[(i == 0) ? halfedgeList.size() - 1 : i - 1];
			H* nextHalfedge = halfedgeList[(i + 1) % halfedgeList.size()];

			V* v0 = quadMesh->halfedgeSource(prevHalfedge);
			V* v1 = quadMesh->halfedgeTarget(prevHalfedge);
			V* v2 = quadMesh->halfedgeSource(nextHalfedge);
			V* v3 = quadMesh->halfedgeTarget(nextHalfedge);
			vector<V*> verts = { v0, v1, v2, v3 };

			// Local indexing mapping
			for (V* v : verts) {
				if (localIdMap.find(v) == localIdMap.end()) {
#pragma omp critical
					{
						int lid = (int)localIdMap.size();
						localIdMap[v] = lid;
						localVertices.push_back(v);
					}
				}
			}

			// Retrieve feature points
			vector<QB_curve_sampling*> relatedTriPoints;
			for (auto& pair : dataPointToQuadHalfedgeMap) {
				if (pair.second == currentHalfedge) {
					relatedTriPoints.push_back(pair.first);
				}
			}

			for (auto dataPoint : relatedTriPoints) {
				int pid = dataPoint->vId();
				if (dataPointIdToParamMap.find(pid) == dataPointIdToParamMap.end()) continue;
				double u = dataPointIdToParamMap[pid];
				double coeffs[4] = {
					pow(1 - u, 3) / 6.0,
					((1 - u) * (1 - u) * (1 - u)) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6,
					pow(1 - u, 3) / 6.0 + u * (1 - u) * (1 - u) + 2 * u * u * (1 - u) + 2 * pow(u, 3) / 3.0,
					pow(u, 3) / 6.0
				};
				V* vertsArr[4] = { v0, v1, v2, v3 };

				for (int j = 0; j < 4; ++j) {
					int lid = localIdMap[vertsArr[j]];
					triplets_local.emplace_back(eq_id, lid, coeffs[j]);
				}
				auto pos = dataPoint->pos();
				v_x_vec.push_back(pos[0]);
				v_y_vec.push_back(pos[1]);
				v_z_vec.push_back(pos[2]);
				eq_id++;
			}
		}

		int n_eq = eq_id;
		int n_var = localVertices.size();
		SparseMatrix<double> G_local(n_eq, n_var);
		G_local.setFromTriplets(triplets_local.begin(), triplets_local.end());
		VectorXd vx = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
		VectorXd vy = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
		VectorXd vz = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

		// Regularized solution: GᵀG x = Gᵀ b + λ I
		double lambda = 1e-6;
		SparseMatrix<double> I(n_var, n_var);
		I.setIdentity();
		SparseMatrix<double> GtG = G_local.transpose() * G_local + lambda * I;
		VectorXd bx = G_local.transpose() * vx;
		VectorXd by = G_local.transpose() * vy;
		VectorXd bz = G_local.transpose() * vz;

		SimplicialLDLT<SparseMatrix<double>> solver;
		solver.compute(GtG);
		VectorXd sol_x = solver.solve(bx);
		VectorXd sol_y = solver.solve(by);
		VectorXd sol_z = solver.solve(bz);

		// Merge into global solution
#pragma omp critical
		{
			for (int i = 0; i < localVertices.size(); ++i) {
				V* v = localVertices[i];
				int gid = quadLoopVertexIdToMap[v];
				solution_x[gid] += sol_x[i];
				solution_y[gid] += sol_y[i];
				solution_z[gid] += sol_z[i];
				control_count[gid]++;
			}
		}
	}

	// Normalize by average and update quadMesh
	for (auto& pair : quadLoopVertexIdToMap) {
		V* v = pair.first;
		int id = pair.second;
		if (control_count[id] == 0) continue;
		v->point()[0] = solution_x[id] / control_count[id];
		v->point()[1] = solution_y[id] / control_count[id];
		v->point()[2] = solution_z[id] / control_count[id];
	}
	//std::cout << "Updated quad mesh vertices (loop local solve) successfully." << std::endl;
}


void CCG_QMSLib::BE_curve_model::compute_quad_vertex3(M* quadMesh)
{
	using namespace Eigen;
	std::vector<double> v_x_vec, v_y_vec, v_z_vec;
	auto set_v_vaules = [](std::vector<double>& v_x_vec, std::vector<double>& v_y_vec, std::vector<double>& v_z_vec, double x, double y, double z)
		{
			v_x_vec.push_back(x);
			v_y_vec.push_back(y);
			v_z_vec.push_back(z);
		};
	std::vector<Eigen::Triplet<double>> tripletlists; // Container for storing all coefficient pairs
	int total_equation_count = 0; // Indicates which equation number

	// Get mapping from quadrilateral halfedges to new IDs
	auto& halfedgeToIdMap = getQuadHalfedgeToIdMap();
	int totalsize = 0;

	// Iterate over each feature loop
	for (auto loop : be_quadFeatureLoops->loops()) {

		// Ensure the current loop is valid
		if (!loop) continue;
		/*for (auto quadHalfEdge : loop->halfedges()) {
			std::cout << "Halfedge pointer: " << quadHalfEdge->localId() << std::endl;
		}*/
		// Get all halfedges in the loop
		auto& halfedges = loop->halfedges();
		std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end()); // Store loop halfedges in a vector
		//std::cout << "halfedgeList.size(): " << halfedgeList.size() << std::endl;
		totalsize += halfedgeList.size();
		std::cout << "totalsize: " << totalsize << std::endl;
		H* startLHalfedge = halfedgeList[0]; // Initial halfedge
		int startLId = halfedgeToIdMap[startLHalfedge];
		// Iterate over each halfedge in the current loop
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i]; // Current halfedge

			// Get the ID of the current halfedge
			int currentId = halfedgeToIdMap[currentHalfedge];
			// Compute IDs for the previous and next halfedges
			int prevId = (i == 0) ? (int)(totalsize - 1) : currentId - 1;
			int nextId = (currentId == (int)(totalsize - 1)) ? startLId : currentId + 1;

			// Get the previous and next halfedges based on ID
			H* prevHalfedge = nullptr;
			H* nextHalfedge = nullptr;
			for (auto& pair : halfedgeToIdMap) {
				if (pair.second == prevId) prevHalfedge = pair.first;
				if (pair.second == nextId) nextHalfedge = pair.first;
			}

			// If both previous and next halfedges are found, proceed with further logic
			if (prevHalfedge && nextHalfedge) {
				/*std::cout << "Current Halfedge ID: " << currentId << std::endl;
				std::cout << "Previous Halfedge ID: " << prevId << std::endl;
				std::cout << "Next Halfedge ID: " << nextId << std::endl;*/
				// Get vertex pointers for quadrilateral mesh based on halfedges
				V* v0 = quadMesh->halfedgeSource(prevHalfedge);
				V* v1 = quadMesh->halfedgeTarget(prevHalfedge);
				V* v2 = quadMesh->halfedgeSource(nextHalfedge);
				V* v3 = quadMesh->halfedgeTarget(nextHalfedge);

				int v0Id = quadLoopVertexIdToMap[v0];
				int v1Id = quadLoopVertexIdToMap[v1];
				int v2Id = quadLoopVertexIdToMap[v2];
				int v3Id = quadLoopVertexIdToMap[v3];

				// Get triangle feature points associated with the current halfedge
				std::vector<QB_curve_sampling*> relatedTriPoints;
				for (auto& pair : dataPointToQuadHalfedgeMap) {
					if (pair.second == currentHalfedge) {
						relatedTriPoints.push_back(pair.first);
					}
				}

				// Iterate over triangle feature points, retrieve parameter values, and process
				for (auto dataPoint : relatedTriPoints) {
					int dataPointId = dataPoint->vId(); // Get the ID of the triangle feature point
					if (dataPointIdToParamMap.find(dataPointId) != dataPointIdToParamMap.end()) {
						double u = dataPointIdToParamMap[dataPointId];
						if (u < 0.0 || u > 1.0) {
							std::cerr << "Warning: u=" << u << " is out of range!" << std::endl;
							u = std::max(0.0, std::min(1.0, u)); // Clamp to [0, 1]
						}
						// Example: print the parameter value of the current triangle feature point
						//std::cout << "Halfedge ID: " << currentId << ", Triangle Point ID: " << triPointId << ", Param Value: " << u << std::endl;

						double coeff_v0 = 0.0;
						double coeff_v1 = 0.0;
						double coeff_v2 = 0.0;
						double coeff_v3 = 0.0;

						/*--------------------------Start of formula------------------------------------------*/
						int temp_index = 0;
						if (v1->feature()) {
							if (v2->feature()) {
								std::cout << "---v1->feature()---v2->feature()---" << std::endl;
								// Coefficient for V1
								temp_index = v1Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) + 2.0 * u * (1 - u) * (1 - u) + u * u * (1 - u));
								coeff_v1 = (1 - u) * (1 - u) * (1 - u) + 2.0 * u * (1 - u) * (1 - u) + u * u * (1 - u);
								// Coefficient for V2
								temp_index = v2Id;
								tripletlists.emplace_back(total_equation_count, temp_index, u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u);
								coeff_v2 = u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u;
							}
							else {
								std::cout << "---v1->feature()---" << std::endl;
								// Coefficient for V1
								temp_index = v1Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6);
								coeff_v1 = (1 - u) * (1 - u) * (1 - u) + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6;
								// Coefficient for V2
								temp_index = v2Id;
								tripletlists.emplace_back(total_equation_count, temp_index, u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2.0 / 3);
								coeff_v2 = u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2.0 / 3;
								// Coefficient for V3
								temp_index = v3Id;
								tripletlists.emplace_back(total_equation_count, temp_index, u * u * u / 6.0);
								coeff_v3 = u * u * u / 6.0;
							}
						}
						else {
							if (v2->feature()) {
								std::cout << "---v2->feature()---" << std::endl;
								// Coefficient for V0
								temp_index = v0Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6.0);
								coeff_v0 = (1 - u) * (1 - u) * (1 - u) / 6.0;
								// Coefficient for V1
								temp_index = v1Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u));
								coeff_v1 = (1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u);
								// Coefficient for V2
								temp_index = v2Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u);
								coeff_v2 = (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u;
							}
							else {
								std::cout << "------" << std::endl;
								// Coefficient for V0
								temp_index = v0Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6);
								coeff_v0 = (1 - u) * (1 - u) * (1 - u) / 6;
								// Coefficient for V1
								temp_index = v1Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6);
								coeff_v1 = (1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6;
								// Coefficient for V2
								temp_index = v2Id;
								tripletlists.emplace_back(total_equation_count, temp_index, (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3);
								coeff_v2 = (1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3;
								// Coefficient for V3
								temp_index = v3Id;
								tripletlists.emplace_back(total_equation_count, temp_index, u * u * u / 6);
								coeff_v3 = u * u * u / 6;
							}
						}
						/*--------------------------End of formula------------------------------------------*/

						double coeff_sum = coeff_v0 + coeff_v1 + coeff_v2 + coeff_v3;
						if (std::abs(coeff_sum - 1.0) > 1e-6) {
							std::cerr << "Warning: B-spline coefficients do not sum to 1! Sum=" << coeff_sum << std::endl;
						}

						CPoint dataPointCoords = dataPoint->pos(); // Ensure dataPoint->pos() returns CPoint type
						set_v_vaules(v_x_vec, v_y_vec, v_z_vec, dataPointCoords[0], dataPointCoords[1], dataPointCoords[2]);
						total_equation_count++;
					}
				}
			}
		}
	}

	// Iterate over feature arcs (data point map unchanged)
	for (auto arc : be_quadFeatureLoops->arcs()) {
		if (!arc) continue;

		// Get all halfedges in the arc
		auto& halfedges = arc->halfedges();
		std::vector<H*> halfedgeList(halfedges.begin(), halfedges.end());
		totalsize += halfedgeList.size();
		std::cout << "Arc halfedges count: " << halfedgeList.size() << std::endl;

		H* startAHalfedge = halfedgeList[0];
		int startAId = halfedgeToIdMap[startAHalfedge];
		std::cout << "startAId: " << startAId << std::endl;
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i];
			int currentId = halfedgeToIdMap[currentHalfedge];

			// Special handling for arcs: open paths only need to process valid adjacent edges
			int prevId = (i > 0) ? halfedgeToIdMap[halfedgeList[i - 1]] : -1;
			int nextId = (i < halfedgeList.size() - 1) ? halfedgeToIdMap[halfedgeList[i + 1]] : -1;

			// Get control vertices (same as loop processing)
			V* v0 = (prevId != -1) ? quadMesh->halfedgeSource(halfedgeList[i - 1]) : nullptr;
			V* v1 = quadMesh->halfedgeSource(currentHalfedge);
			V* v2 = quadMesh->halfedgeTarget(currentHalfedge);
			V* v3 = (nextId != -1) ? quadMesh->halfedgeTarget(halfedgeList[i + 1]) : nullptr;

			// Boundary protection
			//if (!v1) continue;

			// Get associated triangle feature points
			std::vector<V*> relatedTriPoints;
			for (auto& pair : triPointToQuadHalfedgeMap) {
				if (pair.second == currentHalfedge) {
					relatedTriPoints.push_back(pair.first);
				}
			}

			// Process each feature point
			for (auto triPoint : relatedTriPoints) {
				int triPointId = triPoint->id();
				if (triPointIdToParamMap.find(triPointId) != triPointIdToParamMap.end()) {
					double u = triPointIdToParamMap[triPointId];
					CPoint triPointCoords = triPoint->point();

					// Special coefficient handling for arcs (using cubic B-spline)
					if (v0 && v3) { // Full four-vertex case
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v0], (1 - u) * (1 - u) * (1 - u) / 6);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v1],
							(1 - u) * (1 - u) * (1 - u) * 2 / 3 + 2 * u * (1 - u) * (1 - u) + u * u * (1 - u) + u * u * u / 6);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v2],
							(1 - u) * (1 - u) * (1 - u) / 6 + u * (1 - u) * (1 - u) + u * u * (1 - u) * 2 + u * u * u * 2 / 3);
						tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v3], u * u * u / 6);
					}
					//else { // Boundary case (simplified)
					//    tripletlists.emplace_back(total_equation_count, quadLoopVertexIdToMap[v1], 1.0);
					//}
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, triPointCoords[0], triPointCoords[1], triPointCoords[2]);
					// Add endpoint constraints in arc processing
					if (i == 0 || i == halfedgeList.size() - 1) { // First or last vertex
						int vid = quadLoopVertexIdToMap[v1];
						tripletlists.emplace_back(total_equation_count, vid, 1.0);
						set_v_vaules(v_x_vec, v_y_vec, v_z_vec, v1->point()[0], v1->point()[1], v1->point()[2]);
						total_equation_count++;
						//continue; // Skip B-spline calculation
					}
					total_equation_count++;
				}
			}
		}
	}

	// Step 1: Build sparse matrix G and right-hand side vectors g
	int G_row_num = total_equation_count;  // Total number of equations
	int G_col_num = quadLoopVertexIdToMap.size();     // Number of control points
	SparseMatrix<double> G(G_row_num, G_col_num);
	// Construct sparse matrix G from triplet list
	G.setFromTriplets(tripletlists.begin(), tripletlists.end());

	// Check equation count and control point count
	std::cout << "Total equations: " << total_equation_count << std::endl;
	std::cout << "Control points: " << quadLoopVertexIdToMap.size() << std::endl;

	// After constructing all equations, convert std::vector to Eigen::VectorXd
	VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

	// Check data consistency
	std::cout << "v_x_vec size: " << v_x_vec.size() << std::endl;
	std::cout << "v_y_vec size: " << v_y_vec.size() << std::endl;
	std::cout << "v_z_vec size: " << v_z_vec.size() << std::endl;
	std::cout << "G rows: " << G.rows() << std::endl; // Number of rows
	// Diagnostics before solving: output rank and column count
	std::cout << "Rank of G: " << JacobiSVD<MatrixXd>(G).rank() << "/" << G.cols() << std::endl;

	// Step 2: Solve using least squares
	MatrixXd G_dense = G.toDense();

	// Compute condition number of matrix G
	JacobiSVD<MatrixXd> svd(G_dense);
	double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	std::cout << "Condition number of G: " << cond << std::endl;

	// Solve using regularization
	double lambda = 1e-6; // Regularization parameter
	MatrixXd G_reg = G_dense.transpose() * G_dense + lambda * MatrixXd::Identity(G_col_num, G_col_num);
	VectorXd v_x_reg = G_dense.transpose() * v_x;
	VectorXd v_y_reg = G_dense.transpose() * v_y;
	VectorXd v_z_reg = G_dense.transpose() * v_z;
	// Replace original solving method

	VectorXd solution_x = G_reg.ldlt().solve(v_x_reg);
	VectorXd solution_y = G_reg.ldlt().solve(v_y_reg);
	VectorXd solution_z = G_reg.ldlt().solve(v_z_reg);

	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//solver.compute(G);  // G is a sparse matrix
	//if (solver.info() != Eigen::Success) {
	//    std::cerr << "Decomposition failed!" << std::endl;
	//    return;
	//}
	//VectorXd solution_x = solver.solve(v_x);
	//VectorXd solution_y = solver.solve(v_y);
	//VectorXd solution_z = solver.solve(v_z);

	// Print non-zero elements for last 10 equations
	for (int k = 0; k < G.outerSize(); ++k) {
		for (SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
			if (it.row() >= total_equation_count - 10) { // Check last 10 equations
				std::cout << "Row " << it.row() << ", Col " << it.col()
					<< ", Value " << it.value() << std::endl;
			}
		}
	}

	// Print solution results
	for (int i = 0; i < solution_x.size(); ++i) {
		std::cout << "Vertex ID: " << i
			<< ", Solution X: " << solution_x[i]
			<< ", Solution Y: " << solution_y[i]
				<< ", Solution Z: " << solution_z[i] << std::endl;
	}

	// Update quadrilateral mesh vertices
	for (auto& pair : quadLoopVertexIdToMap) {
		V* quadVertex = pair.first;
		int quadVertexId = pair.second;
		double new_x = solution_x[quadVertexId];
		double new_y = solution_y[quadVertexId];
		double new_z = solution_z[quadVertexId];
		CPoint& vertexPoint = quadVertex->point(); // Get reference to vertex
		vertexPoint[0] = new_x;
		vertexPoint[1] = new_y;
		vertexPoint[2] = new_z;
		/*std::cout << "111Vertex ID: " << quadVertexId
				  << ", Solution X: " << quadVertex->point()[0]
				  << ", Solution Y: " << quadVertex->point()[1]
				  << ", Solution Z: " << quadVertex->point()[2] << std::endl;*/
	}
	std::cout << "Updated quad mesh vertices successfully." << std::endl;

	// Add checks in arc processing
	/*for (auto arc : be_quadFeatureLoops->arcs()) {
		for (auto he : arc->halfedges()) {
			V* v = quadMesh->halfedgeTarget(he);
			if (quadLoopVertexIdToMap.find(v) == quadLoopVertexIdToMap.end()) {
				std::cerr << "Error: Arc vertex " << v << " not mapped!" << std::endl;
			}
			else
				std::cout << "Mapping correct" << std::endl;
		}
	}*/
	// Improved version: construct subsystems by loop/arc and solve in parallel (local solve + OpenMP)
	// Each loop independently constructs its own G_local and v_local, then solves in parallel and merges results
}


void CCG_QMSLib::BE_curve_model::compute_quad_vertex4(M* quadMesh)
{
	using namespace Eigen;
	using namespace std;

	// Total number of control points
	int total_control_points = quadLoopVertexIdToMap.size();

	// Final solved results
	VectorXd solution_x = VectorXd::Zero(total_control_points);
	VectorXd solution_y = VectorXd::Zero(total_control_points);
	VectorXd solution_z = VectorXd::Zero(total_control_points);

	// Number of times each control point is assigned, for averaging
	vector<int> control_count(total_control_points, 0);

	// Process all loops together (parallel)
	auto& loops = be_quadFeatureLoops->loops();
#pragma omp parallel for
	for (int loopIdx = 0; loopIdx < (int)loops.size(); ++loopIdx) {
		auto loop = loops[loopIdx];
		if (!loop) continue;

		auto& halfedges = loop->halfedges();
		vector<H*> halfedgeList(halfedges.begin(), halfedges.end());
		vector<Triplet<double>> triplets_local;
		vector<double> v_x_vec, v_y_vec, v_z_vec;
		std::map<V*, int> localIdMap;
		vector<V*> localVertices;
		int eq_id = 0;

		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHalfedge = halfedgeList[i];
			H* prevHalfedge = halfedgeList[(i == 0) ? halfedgeList.size() - 1 : i - 1];
			H* nextHalfedge = halfedgeList[(i + 1) % halfedgeList.size()];

			V* v0 = quadMesh->halfedgeSource(prevHalfedge);
			V* v1 = quadMesh->halfedgeTarget(prevHalfedge);
			V* v2 = quadMesh->halfedgeSource(nextHalfedge);
			V* v3 = quadMesh->halfedgeTarget(nextHalfedge);
			vector<V*> verts = { v0, v1, v2, v3 };

			// Local index mapping
			for (V* v : verts) {
				if (localIdMap.find(v) == localIdMap.end()) {
#pragma omp critical
					{
						int lid = (int)localIdMap.size();
						localIdMap[v] = lid;
						localVertices.push_back(v);
					}
				}
			}

			// Retrieve feature sampling points
			vector<QB_curve_sampling*> relatedTriPoints;
			for (auto& pair : dataPointToQuadHalfedgeMap) {
				if (pair.second == currentHalfedge) {
					relatedTriPoints.push_back(pair.first);
				}
			}

			for (auto dataPoint : relatedTriPoints) {
				int pid = dataPoint->vId();
				if (dataPointIdToParamMap.find(pid) == dataPointIdToParamMap.end()) continue;
				double u = dataPointIdToParamMap[pid];

				if (u < 0.0 || u > 1.0) {
					std::cerr << "Warning: u=" << u << " is out of range!" << std::endl;
					u = std::max(0.0, std::min(1.0, u)); // Clamp to [0, 1]
				}
				double coeffs[4] = { 0.0, 0.0, 0.0, 0.0 };
				if (v1->feature())
				{
					if (v2->feature())
					{
						//std::cout << "--- v1->feature() and v2->feature() ---" << std::endl;
						// Coefficient for V1
						coeffs[1] = (1 - u) * (1 - u) * (1 - u)
							+ 2.0 * u * (1 - u) * (1 - u)
							+ u * u * (1 - u);
						// Coefficient for V2
						coeffs[2] = u * (1 - u) * (1 - u)
							+ u * u * (1 - u) * 2
							+ u * u * u;
					}
					else
					{
						//std::cout << "--- v1->feature() only ---" << std::endl;
						// Coefficient for V1
						coeffs[1] = (1 - u) * (1 - u) * (1 - u)
							+ 2 * u * (1 - u) * (1 - u)
							+ u * u * (1 - u)
							+ u * u * u / 6;

						// Coefficient for V2
						coeffs[2] = u * (1 - u) * (1 - u)
							+ u * u * (1 - u) * 2
							+ u * u * u * 2.0 / 3;

						// Coefficient for V3
						coeffs[3] = u * u * u / 6.0;
					}
				}
				else
				{
					if (v2->feature())
					{
						//std::cout << "--- v2->feature() only ---" << std::endl;
						// Coefficient for V0
						coeffs[0] = (1 - u) * (1 - u) * (1 - u) / 6.0;

						// Coefficient for V1
						coeffs[1] = (1 - u) * (1 - u) * (1 - u) * 2 / 3
							+ 2 * u * (1 - u) * (1 - u)
							+ u * u * (1 - u);

						// Coefficient for V2
						coeffs[2] = (1 - u) * (1 - u) * (1 - u) / 6
							+ u * (1 - u) * (1 - u)
							+ u * u * (1 - u) * 2
							+ u * u * u;
					}
					else
					{
						//std::cout << "--- no features on v1 or v2 ---" << std::endl;
						// Coefficient for V0
						coeffs[0] = (1 - u) * (1 - u) * (1 - u) / 6;

						// Coefficient for V1
						coeffs[1] = (1 - u) * (1 - u) * (1 - u) * 2 / 3
							+ 2 * u * (1 - u) * (1 - u)
							+ u * u * (1 - u)
							+ u * u * u / 6;

						// Coefficient for V2
						coeffs[2] = (1 - u) * (1 - u) * (1 - u) / 6
							+ u * (1 - u) * (1 - u)
							+ u * u * (1 - u) * 2
							+ u * u * u * 2 / 3;

						// Coefficient for V3
						coeffs[3] = u * u * u / 6;
					}
				}
				V* verts[4] = { v0, v1, v2, v3 };

				for (int j = 0; j < 4; ++j) {
					int lid = localIdMap[verts[j]];
					triplets_local.emplace_back(eq_id, lid, coeffs[j]);
				}
				auto pos = dataPoint->pos();
				v_x_vec.push_back(pos[0]);
				v_y_vec.push_back(pos[1]);
				v_z_vec.push_back(pos[2]);
				eq_id++;
			}
		}

		int n_eq = eq_id;
		int n_var = localVertices.size();
		SparseMatrix<double> G_local(n_eq, n_var);
		G_local.setFromTriplets(triplets_local.begin(), triplets_local.end());
		VectorXd vx = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
		VectorXd vy = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
		VectorXd vz = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

		// Regularized solve: GᵗG x = Gᵗb + λI
		double lambda = 1e-6;
		SparseMatrix<double> I(n_var, n_var);
		I.setIdentity();
		SparseMatrix<double> GtG = G_local.transpose() * G_local + lambda * I;
		VectorXd bx = G_local.transpose() * vx;
		VectorXd by = G_local.transpose() * vy;
		VectorXd bz = G_local.transpose() * vz;

		SimplicialLDLT<SparseMatrix<double>> solver;
		solver.compute(GtG);
		VectorXd sol_x = solver.solve(bx);
		VectorXd sol_y = solver.solve(by);
		VectorXd sol_z = solver.solve(bz);

		// Merge into global solution
#pragma omp critical
		{
			for (int i = 0; i < localVertices.size(); ++i) {
				V* v = localVertices[i];
				int gid = quadLoopVertexIdToMap[v];
				solution_x[gid] += sol_x[i];
				solution_y[gid] += sol_y[i];
				solution_z[gid] += sol_z[i];
				control_count[gid]++;
			}
		}
		// Normalize by average and update to quadMesh
		for (auto& pair : quadLoopVertexIdToMap) {
			V* v = pair.first;
			int id = pair.second;
			if (control_count[id] == 0) continue;
			v->point()[0] = solution_x[id] / control_count[id];
			v->point()[1] = solution_y[id] / control_count[id];
			v->point()[2] = solution_z[id] / control_count[id];
		}

	for (auto& pair : quadLoopVertexIdToMap) {
		V* v = pair.first;
		int id = pair.second;
		if (control_count[id] == 0) continue;
		v->point()[0] = solution_x[id] / control_count[id];
		v->point()[1] = solution_y[id] / control_count[id];
		v->point()[2] = solution_z[id] / control_count[id];
	}
	//std::cout << "Updated quad mesh vertices (loop local solve) successfully." << std::endl;
}
}

void CCG_QMSLib::BE_curve_model::compute_quad_vertex5(M* quadMesh)
{
	using namespace Eigen;
	std::vector<double> v_x_vec, v_y_vec, v_z_vec;
	auto set_v_vaules = [](std::vector<double>& v_x_vec, std::vector<double>& v_y_vec, std::vector<double>& v_z_vec, double x, double y, double z)
	{
		v_x_vec.push_back(x);
		v_y_vec.push_back(y);
		v_z_vec.push_back(z);
	};
	std::vector<Eigen::Triplet<double>> tripletlists; // 定义一个用于存储所有系数对的容器
	int total_equation_count = 0;//表示的是第几个方程

	// 获取四边形半边到新 ID 的映射
	auto& halfedgeToIdMap = getQuadBoundaryArcHalfedgeToIdMap();
	int totalsize = 0;

	// 遍历边界特征弧
	// 1. 获取所有弧
	const std::vector<BoundaryFittingArc>& arcs = getBoundaryFittingArcs();

	for (const auto& arc : arcs) {
		if (arc.orderedHalfedges.empty()) continue;
		const auto& halfedgeList = arc.orderedHalfedges;
		totalsize += halfedgeList.size();

		// 遍历弧中的每条半边
		for (size_t i = 0; i < halfedgeList.size(); ++i) {
			H* currentHe = halfedgeList[i];
			int currentId = halfedgeToIdMap[currentHe];

			// --------------------------
			// 获取控制顶点v0/v1/v2/v3
			// --------------------------
			int prevId = (i > 0) ? halfedgeToIdMap[halfedgeList[i - 1]] : -1;
			int nextId = (i < halfedgeList.size() - 1) ? halfedgeToIdMap[halfedgeList[i + 1]] : -1;
			V* v0 = (prevId != -1) ? quadMesh->halfedgeSource(halfedgeList[i - 1]) : nullptr;
			V* v1 = quadMesh->halfedgeSource(currentHe);
			V* v2 = quadMesh->halfedgeTarget(currentHe);
			V* v3 = (nextId != -1) ? quadMesh->halfedgeTarget(halfedgeList[i + 1]) : nullptr;
			if (!v1 || !v2) continue;  // 边界保护（v1/v2必存在）

			// --------------------------
			// 查找已知量：前/后边界半边顶点（v0或v3缺失时补充）
			// --------------------------
			V* knownV0 = nullptr;  // 已知量v0（v0不存在时从边界查找）
			V* knownV3 = nullptr;  // 已知量v3（v3不存在时从边界查找）

			// 1. 若v0不存在（!v0），从当前弧起点（v1）查找前一条边界半边顶点
			if (!v0) {
				V* arcStartV = v1;
				for (auto outH : It::VClwInHEIterator(quadMesh, arcStartV)) {  // 顺时针遍历起点入边
					if (quadMesh->halfedgeEdge(outH)->boundary() && outH != currentHe) {  // 边界边且非当前半边
						knownV0 = quadMesh->halfedgeSource(outH);  // 前一条边界半边的起点（已知量）
						break;
					}
				}
				//for (auto outH : It::VCcwOutHEIterator(quadMesh, arcStartV)) {  // 逆时针遍历起点出边
				//	if (quadMesh->halfedgeEdge(outH)->boundary() && outH != currentHe) {  // 边界边且非当前半边
				//		knownV0 = quadMesh->halfedgeTarget(outH);  // 前一条边界半边的起点（已知量）
				//		break;
				//	}
				//}
			}

			// 2. 若v3不存在（!v3），从当前弧终点（v2）查找后一条边界半边顶点
			if (!v3) {
				V* arcEndV = v2;
				for (auto outH : It::VCcwOutHEIterator(quadMesh, arcEndV)) {  // 逆时针遍历终点出边
					if (quadMesh->halfedgeEdge(outH)->boundary() && outH != currentHe) {  // 边界边且非当前半边
						knownV3 = quadMesh->halfedgeTarget(outH);  // 后一条边界半边的终点（已知量）
						break;
					}
				}
			}
			if (!v0)
			{
				//if (!knownV0)
					//std::cout << "-------------0----------------------" << std::endl;
			}
			if (!v3)
			{
				//if (!knownV3)
					//std::cout << "--------------3---------------------" << std::endl;
			}
			// --------------------------
			// 获取关联的三角形特征点（原有逻辑）
			// --------------------------
			// 获取特征点
			std::vector<QB_curve_sampling*> relatedTriPoints;
			for (auto& pair : dataPointToQuadHalfedgeMap) {
				if (pair.second == currentHe) {
					relatedTriPoints.push_back(pair.first);
					//std::cout << "Found related tri point: " << pair.first->vId() << std::endl;
				}
			}

			// --------------------------
			// 处理每个特征点：分4种情况列方程（含v0/v3均不存在）
			// --------------------------
			for (auto triPoint : relatedTriPoints) {
				/*std::cout << "111111111111111111111111" << std::endl;*/
				int triPointId = triPoint->vId();
				auto paramIt = dataPointIdToParamMap.find(triPointId);
				if (paramIt == dataPointIdToParamMap.end()) continue;
				double u = paramIt->second;
				CPoint triPointCoords = triPoint->pos();  // 方程右端项（基础值）
				//std::cout << "111111111111111111111111" << std::endl;
				// --------------------------
				// 情况1：v0和v3均存在（完整四顶点，中间点）
				// --------------------------
				if (v0 && v3) {
					// 四变量均为弧内顶点（变量），使用边界弧ID映射表
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v0], pow(1 - u, 3) / 6.0);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v1],
						pow(1 - u, 3) * 2 / 3 + 2 * u * pow(1 - u, 2) + u * u * (1 - u) + pow(u, 3) / 6);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v2],
						pow(1 - u, 3) / 6 + u * pow(1 - u, 2) + 2 * u * u * (1 - u) + pow(u, 3) * 2 / 3);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v3], pow(u, 3) / 6.0);
					// 右端项：直接使用特征点坐标
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, triPointCoords[0], triPointCoords[1], triPointCoords[2]);
				}
				// --------------------------
				// 情况2：!v0（v0不存在，已知量代入v0项）
				// --------------------------
				else if (!v0 && knownV0 && v3) {
					// 已知量：knownV0（前一条边界顶点），变量：v1, v2, v3
					CPoint v0Fixed = knownV0->point();
					double coeffV0 = pow(1 - u, 3) / 6.0;
					CPoint v0Term = v0Fixed * coeffV0;  // 已知项：v0 * 系数

					// 变量项（v1, v2, v3）
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v1],
						pow(1 - u, 3) * 2 / 3 + 2 * u * pow(1 - u, 2) + u * u * (1 - u) + pow(u, 3) / 6);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v2],
						pow(1 - u, 3) / 6 + u * pow(1 - u, 2) + 2 * u * u * (1 - u) + pow(u, 3) * 2 / 3);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v3], pow(u, 3) / 6.0);
					// 右端项修正：特征点坐标 - 已知项v0Term
					CPoint correctedCoords = triPointCoords - v0Term;
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, correctedCoords[0], correctedCoords[1], correctedCoords[2]);
				}
				// --------------------------
				// 情况3：!v3（v3不存在，已知量代入v3项）
				// --------------------------
				else if (v0 && !v3 && knownV3) {
					// 已知量：knownV3（后一条边界顶点），变量：v0, v1, v2
					CPoint v3Fixed = knownV3->point();
					double coeffV3 = pow(u, 3) / 6.0;
					CPoint v3Term = v3Fixed * coeffV3;  // 已知项：v3 * 系数

					// 变量项（v0, v1, v2）
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v0], pow(1 - u, 3) / 6.0);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v1],
						pow(1 - u, 3) * 2 / 3 + 2 * u * pow(1 - u, 2) + u * u * (1 - u) + pow(u, 3) / 6);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v2],
						pow(1 - u, 3) / 6 + u * pow(1 - u, 2) + 2 * u * u * (1 - u) + pow(u, 3) * 2 / 3);
					// 右端项修正：特征点坐标 - 已知项v3Term
					CPoint correctedCoords = triPointCoords - v3Term;
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, correctedCoords[0], correctedCoords[1], correctedCoords[2]);
				}
				// --------------------------
				// 情况4：v0和v3均不存在（!v0 && !v3），双已知量代入
				// --------------------------
				else if (!v0 && !v3 && knownV0 && knownV3) {
					// 已知量：knownV0（前边界顶点）+ knownV3（后边界顶点），变量：v1, v2
					CPoint v0Fixed = knownV0->point();
					CPoint v3Fixed = knownV3->point();
					double coeffV0 = pow(1 - u, 3) / 6.0;
					double coeffV3 = pow(u, 3) / 6.0;
					CPoint knownTerms = v0Fixed * coeffV0 + v3Fixed * coeffV3;  // 总已知项

					// 变量项（仅v1, v2）
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v1],
						pow(1 - u, 3) * 2 / 3 + 2 * u * pow(1 - u, 2) + u * u * (1 - u) + pow(u, 3) / 6);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v2],
						pow(1 - u, 3) / 6 + u * pow(1 - u, 2) + 2 * u * u * (1 - u) + pow(u, 3) * 2 / 3);
					// 右端项修正：特征点坐标 - 总已知项
					CPoint correctedCoords = triPointCoords - knownTerms;
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, correctedCoords[0], correctedCoords[1], correctedCoords[2]);
				}
				// --------------------------
				// 异常处理：已知量未找到（跳过或按简化情况处理）
				// --------------------------
				else {
					// 若已知量查找失败（如边界不连续），可简化为双顶点线性插值（v1, v2）
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v1], 1 - u);
					tripletlists.emplace_back(total_equation_count, quadBoundaryArcVertexIdToMap[v2], u);
					set_v_vaules(v_x_vec, v_y_vec, v_z_vec, triPointCoords[0], triPointCoords[1], triPointCoords[2]);
					//std::cout << "Warning: Missing known vertices, fallback to linear interpolation for v1-v2." << std::endl;
					//std::cout << "v1->id" << v1->id() << std::endl;
					//std::cout << "v2->id" << v2->id() << std::endl;
				}
				total_equation_count++;
			}
		}
	}
	// Step 1: 构建稀疏矩阵 G 和右边的向量 g
	int G_row_num = total_equation_count;  // 方程总数
	int G_col_num = quadBoundaryArcVertexIdToMap.size();     // 控制点的数量
	SparseMatrix<double> G(G_row_num, G_col_num);
	// 从 tripletlist 构造稀疏矩阵 G
	G.setFromTriplets(tripletlists.begin(), tripletlists.end());

	// 检查方程数量和控制点数量
	//std::cout << "Total equations: " << total_equation_count << std::endl;
	//std::cout << "Control points: " << quadBoundaryArcVertexIdToMap.size() << std::endl;

	// 构造完所有方程后，将 std::vector 转换为 Eigen::VectorXd
	VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

	// 检查数据一致性
	//std::cout << "v_x_vec size: " << v_x_vec.size() << std::endl;
	//std::cout << "v_y_vec size: " << v_y_vec.size() << std::endl;
	//std::cout << "v_z_vec size: " << v_z_vec.size() << std::endl;
	//std::cout << "G rows: " << G.rows() << std::endl;//行数
	//在求解前添加诊断，输出秩和列数
	//std::cout << "Rank of G: " << JacobiSVD<MatrixXd>(G).rank() << "/" << G.cols() << std::endl;
	// Step 2: 使用最小二乘法求解
	MatrixXd G_dense = G.toDense();

	// 计算矩阵 G 的条件数
	JacobiSVD<MatrixXd> svd(G_dense);
	double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//std::cout << "Condition number of G: " << cond << std::endl;

	//使用正则化方法求解
	double lambda = 1e-6; // 正则化参数
	MatrixXd G_reg = G_dense.transpose() * G_dense + lambda * MatrixXd::Identity(G_col_num, G_col_num);
	VectorXd v_x_reg = G_dense.transpose() * v_x;
	VectorXd v_y_reg = G_dense.transpose() * v_y;
	VectorXd v_z_reg = G_dense.transpose() * v_z;
	// 替换原有求解方式

	VectorXd solution_x = G_reg.ldlt().solve(v_x_reg);
	VectorXd solution_y = G_reg.ldlt().solve(v_y_reg);
	VectorXd solution_z = G_reg.ldlt().solve(v_z_reg);

	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//solver.compute(G);  // G是稀疏矩阵
	//if (solver.info() != Eigen::Success) {
	//	std::cerr << "Decomposition failed!" << std::endl;
	//	return;
	//}
	//VectorXd solution_x = solver.solve(v_x);
	//VectorXd solution_y = solver.solve(v_y);
	//VectorXd solution_z = solver.solve(v_z);

	// 在构建矩阵后打印非零元素
	for (int k = 0; k < G.outerSize(); ++k) {
		for (SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
			if (it.row() >= total_equation_count - 10) { // 检查最后10个方程
				/*std::cout << "Row " << it.row() << ", Col " << it.col()
					<< ", Value " << it.value() << std::endl;*/
			}
		}
	}

	// 输出求解结果
	for (int i = 0; i < solution_x.size(); ++i) {
		/*std::cout << "Vertex ID: " << i
			<< ", Solution X: " << solution_x[i]
			<< ", Solution Y: " << solution_y[i]
			<< ", Solution Z: " << solution_z[i] << std::endl;*/
	}
	// 更新四边形网格点
	for (auto& pair : quadBoundaryArcVertexIdToMap) {
		V* quadVertex = pair.first;
		int quadVertexId = pair.second;
		double new_x = solution_x[quadVertexId];
		double new_y = solution_y[quadVertexId];
		double new_z = solution_z[quadVertexId];
		CPoint& vertexPoint = quadVertex->point(); // 获取顶点的引用
		vertexPoint[0] = new_x;
		vertexPoint[1] = new_y;
		vertexPoint[2] = new_z;
		/*std::cout << "111Vertex ID: " << quadVertexId
			<< ", Solution X: " << quadVertex->point()[0]
			<< ", Solution Y: " << quadVertex->point()[1]
			<< ", Solution Z: " << quadVertex->point()[2] << std::endl;*/
	}
	//std::cout << "Updated quad mesh vertices successfully." << std::endl;
	// 在弧处理部分添加检查
	/*for (auto arc : be_quadFeatureLoops->arcs()) {
		for (auto he : arc->halfedges()) {
			V* v = quadMesh->halfedgeTarget(he);
			if (quadLoopVertexIdToMap.find(v) == quadLoopVertexIdToMap.end()) {
				std::cerr << "Error: Arc vertex " << v << " not mapped!" << std::endl;
			}
			else
				std::cout << "map正确" << std::endl;
		}
	}*/
	//// Step 1: 构建稀疏矩阵 G 和右边的向量 g
	//int G_row_num = total_equation_count;  // 方程总数
	//int G_col_num = quadBoundaryArcVertexIdToMap.size();     // 控制点的数量
	//SparseMatrix<double> G(G_row_num, G_col_num);

	//// 从 tripletlist 构造稀疏矩阵 G
	//G.setFromTriplets(tripletlists.begin(), tripletlists.end());
	//for (int i = 0; i < G.rows(); i++)
	//{
	//	//std::cout << "G.row(i).sum(): " << G.row(i).sum() <<std::endl;
	//	if (G.row(i).sum() < 0.9 || G.row(i).sum() > 1.1)
	//	{
	//		std::cout << "G.row(i).sum(): " << G.row(i).sum() << std::endl;
	//	}
	//}
	//// 构造完所有方程后，将 std::vector 转换为 Eigen::VectorXd
	//VectorXd v_x = VectorXd::Map(v_x_vec.data(), v_x_vec.size());
	//VectorXd v_y = VectorXd::Map(v_y_vec.data(), v_y_vec.size());
	//VectorXd v_z = VectorXd::Map(v_z_vec.data(), v_z_vec.size());

	//std::cout << "Matrix G size: " << G.rows() << "x" << G.cols() << std::endl;
	//std::cout << "v_x size: " << v_x.size() << std::endl;

	//// Step 2: 使用最小二乘法求解
	//// 将稀疏矩阵G转换为密集矩阵，以便使用SVD
	//MatrixXd G_dense = G.toDense();
	//for (int i = 0; i < G_dense.rows(); i++)
	//{
	//	//std::cout << "G.row(i).sum(): " << G.row(i).sum() <<std::endl;
	//	if (G_dense.row(i).sum() < 0.9 || G_dense.row(i).sum() > 1.1)
	//	{
	//		std::cout << "G_dense.row(i).sum(): " << G_dense.row(i).sum() << std::endl;
	//	}
	//}
	//// 创建一个SVD求解器
	//JacobiSVD<MatrixXd> svd(G_dense, ComputeThinU | ComputeThinV);

	//// 计算伪逆
	//MatrixXd G_pinv = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();
	//// 分别求解 x, y, z
	//VectorXd solution_x = G_pinv * v_x;
	//VectorXd solution_y = G_pinv * v_y;
	//VectorXd solution_z = G_pinv * v_z;
	//for (int i = 0; i < solution_x.size(); ++i) {
	//	std::cout << "Vertex ID: " << i
	//		<< ", Solution X: " << solution_x[i]
	//		<< ", Solution Y: " << solution_y[i]
	//		<< ", Solution Z: " << solution_z[i] << std::endl;
	//}

	//// 更新四边形网格点
	//for (auto& pair : quadBoundaryArcVertexIdToMap) {
	//	V* quadVertex = pair.first;
	//	int quadVertexId = pair.second;
	//	double new_x = solution_x[quadVertexId];
	//	double new_y = solution_y[quadVertexId];
	//	double new_z = solution_z[quadVertexId];
	//	std::cout << "Vertex ID: " << quadVertexId
	//		<< ", Solution X: " << new_x
	//		<< ", Solution Y: " << new_y
	//		<< ", Solution Z: " << new_z << std::endl;
	//	CPoint& vertexPoint = quadVertex->point(); // 获取顶点的引用
	//	vertexPoint[0] = new_x;
	//	vertexPoint[1] = new_y;
	//	vertexPoint[2] = new_z;
	//}
	//std::cout << "Updated quad mesh vertices successfully." << std::endl;

}

double CCG_QMSLib::BE_curve_model::computeTotalError() {
	double totalError = 0.0; // Initialize total error to 0
	int pointCount = 0;     // Initialize point count to 0
	auto& dataToQuadMap = dataPointToQuadHalfedgeMap; // Get mapping from triangle sampling points to quad halfedges
	auto& paramMap = dataPointIdToParamMap;            // Get parameter mapping for triangle sampling points

	// Iterate over each triangle feature point
	for (auto& pair : dataToQuadMap) {
		QB_curve_sampling* dataV = pair.first; // Triangle sampling point
		H* quadH = pair.second;                // Corresponding quad halfedge

		// Find the parameter value for the current sampling point
		auto paramIt = paramMap.find(dataV->vId());
		if (paramIt == paramMap.end()) continue; // Skip if no parameter value found
		double u = paramIt->second;               // Get the local parameter value

		// Obtain the Bézier curve directly via quadHalfedgeToBezierCurveMap
		auto bezierIt = quadHalfedgeToBezierCurveMap.find(quadH);
		if (bezierIt == quadHalfedgeToBezierCurveMap.end()) continue; // Skip if no Bézier curve found
		BE_bezierCurve_approxi* bezierCurve = bezierIt->second;

		// Retrieve control points of the Bézier curve
		std::vector<CPoint>& controlPoints = bezierCurve->cpts();
		if (controlPoints.size() != 4) continue; // Ensure it's a cubic Bézier curve (4 control points)

		// Compute approximation point using the Bézier formula
		double u1 = 1.0 - u;
		CPoint approxPt =
			controlPoints[0] * (u1 * u1 * u1)
			+ controlPoints[1] * (3 * u1 * u1 * u)
			+ controlPoints[2] * (3 * u1 * u * u)
			+ controlPoints[3] * (u * u * u);

		// Get the original point position
		CPoint originalPt = dataV->pos();

		// Compute error (Euclidean distance)
		double error = (originalPt - approxPt).norm();
		//std::cout << "pointError " << error << std::endl;
		totalError += error; // Accumulate error
		pointCount++;        // Increment point count
	}

	// Compute average error
	double averageError = (pointCount > 0) ? (totalError / pointCount) : 0.0;

	//std::cout << "totalError: " << totalError << std::endl;
	//std::cout << "averageError: " << averageError << std::endl;

	return totalError; // Return total error
}

double CCG_QMSLib::BE_curve_model::computeTriTotalError() {
	double totalError = 0.0; // Initialize total error to 0
	int pointCount = 0;     // Initialize point count to 0
	auto& dataToQuadMap = triPointToQuadHalfedgeMap; // Get mapping from triangle points to quad halfedges
	auto& paramMap = triPointIdToParamMap;            // Get parameter mapping for triangle points

	// Iterate over each triangular feature point
	for (auto& pair : dataToQuadMap) {
		V* dataV = pair.first; // Triangle point
		H* quadH = pair.second; // Corresponding quad halfedge

		// Find the parameter value for the current triangle point
		auto paramIt = paramMap.find(dataV->id());
		if (paramIt == paramMap.end()) continue; // Skip if no parameter value found
		double u = paramIt->second;               // Get the local parameter value

		// Obtain the Bézier curve directly via quadHalfedgeToBezierCurveMap
		auto bezierIt = quadHalfedgeToBezierCurveMap.find(quadH);
		if (bezierIt == quadHalfedgeToBezierCurveMap.end()) continue; // Skip if no Bézier curve found
		BE_bezierCurve_approxi* bezierCurve = bezierIt->second;

		// Retrieve control points of the Bézier curve
		std::vector<CPoint>& controlPoints = bezierCurve->cpts();
		if (controlPoints.size() != 4) continue; // Ensure it's a cubic Bézier curve (4 control points)

		// Compute approximation point using the Bézier formula
		double u1 = 1.0 - u;
		CPoint approxPt =
			controlPoints[0] * (u1 * u1 * u1)
			+ controlPoints[1] * (3 * u1 * u1 * u)
			+ controlPoints[2] * (3 * u1 * u * u)
			+ controlPoints[3] * (u * u * u);

		// Get the original point position
		CPoint originalPt = dataV->point();

		// Compute error (Euclidean distance)
		double error = (originalPt - approxPt).norm();
		//std::cout << "pointError " << error << std::endl;
		totalError += error; // Accumulate error
		pointCount++;        // Increment point count
	}

	// Compute average error (unused here)
	double averageError = (pointCount > 0) ? (totalError / pointCount) : 0.0;
	/*std::cout << "----------------------------------------- " <<  std::endl;
	std::cout << "triTotalError: " << totalError << std::endl;
	std::cout << "triAverageError: " << averageError << std::endl;
	std::cout << "----------------------------------------- " << std::endl;*/
	return totalError; // Return total error
}

void CCG_QMSLib::BE_curve_model::preserveQuadCornerPoints(M* pMesh) {
	for (auto& entry : sharpFeaturePoints) {
		auto vid = entry.first;
		auto originalPos = entry.second;
		V* v = pMesh->idVertex(vid);
		if (v) {
			v->point() = originalPos;
			/*std::cout << "Restored sharp feature at vertex " << vid
					  << " to (" << originalPos[0] << ", "
					  << originalPos[1] << ", " << originalPos[2] << ")\n";*/
		}
	}

	// 2. Verify output
	//std::cout << "\n======== Sharp Feature Preservation ========\n"
	//          << "Restored " << sharpFeaturePoints.size() << " sharp features\n"
	//          /*<< "Restored " << originalArcStartEndPoints.size() << " arc endpoints\n"*/
	//          << "==========================================\n";
}


void CCG_QMSLib::BE_curve_model::CheckDegeneratedBoudnaryVertexAndUpdate(M* pMesh)
{
	/*obtain all boundary vertex*/
	std::vector<V*> boundaryVs;
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->boundary())
		{
			boundaryVs.push_back(v);
		}
		v->ifVisit() = false;
	}
	/*check degenerated boundary vertex*/
	bool mark_update = false;
	do
	{
		mark_update = false;
		for (auto bv : boundaryVs)
		{
			if (bv->ifVisit()) continue;
			bool mark_degenerate = false;
			for (auto inH : It::VClwInHEIterator(pMesh, bv))
			{
				if (pMesh->halfedgeEdge(inH)->boundary())
				{
					V* v1 = pMesh->halfedgeSource(inH);
					V* v2 = pMesh->halfedgeTarget(inH);
					V* v3 = pMesh->halfedgeTarget(pMesh->halfedgeNext(inH));
					V* v4 = pMesh->halfedgeSource(pMesh->halfedgePrev(inH));
					CPoint n1 = (v3->point() - v2->point()) ^ (v1->point() - v2->point());
					CPoint n2 = (v1->point() - v4->point()) ^ (v3->point() - v4->point());
					if (n1 * n2 < 0.0)
					{
						bv->point() = bv->initialPos();
						bv->ifVisit() = true;
						//std::cout << " ---degenerated boundary vertex bv->id(): " << bv->id() << std::endl;
						mark_update = true;
						mark_degenerate = true;
						break;
					}
				}
			}
			if (!mark_degenerate)
			{
				for (auto outH : It::VCcwOutHEIterator(pMesh, bv))
				{
					if (pMesh->halfedgeEdge(outH)->boundary())
					{
						V* v1 = pMesh->halfedgeSource(pMesh->halfedgePrev(outH));
						V* v2 = pMesh->halfedgeSource(outH);
						V* v3 = pMesh->halfedgeTarget(outH);
						V* v4 = pMesh->halfedgeTarget(pMesh->halfedgeNext(outH));
						CPoint n1 = (v3->point() - v2->point()) ^ (v1->point() - v2->point());
						CPoint n2 = (v1->point() - v4->point()) ^ (v3->point() - v4->point());
						if (n1 * n2 < 0.0)
						{
							bv->point() = bv->initialPos();
							bv->ifVisit() = true;
							//std::cout << " ---degenerated boundary vertex bv->id(): " << bv->id() << std::endl;
							mark_update = true;
							break;
						}
					}
				}
			}
		}
	} while (mark_update);
	//initialize the attribute about ifvisit
	for (auto bv : boundaryVs)
	{
		bv->ifVisit() = false;
	}
	//std::cout << "check degenerated boundary vertex finished!" << std::endl;
}

void CCG_QMSLib::BE_curve_model::postProcessing_checkAndModifyOrthogonality_boundaryFitting(M* pMesh)
{
	/*1. obtaing vertex set that is modified on boudnary fitting processing*/
	std::vector<V*> tempVs;
	for (auto poorFittingEdge : poorFittingEdges)
	{
		pMesh->halfedgeSource(poorFittingEdge.first)->ifVisit() = true;
		pMesh->halfedgeTarget(poorFittingEdge.first)->ifVisit() = true;
	}
	for (auto v : It::MVIterator(pMesh))
	{
		if (!v->boundary()) continue;
		/*Ignoring the extraordinary boundary vertex*/
		int num_vfs = 0;
		for (auto vf : It::VCcwFIterator(pMesh, v))
		{
			num_vfs++;
		}
		if (num_vfs != 2) continue;
		/*if (v->ifVisit())
		{
			tempVs.push_back(v);
		}*/
		tempVs.push_back(v);
	}
	/*for (auto tempV : tempVs)
	{
		tempV->ifVisit() = false;
	}*/
	std::cout << "# poor vertex: " << tempVs.size() << std::endl;
	/*2. check and modify*/
	for (auto tempV : tempVs)
	{
		/*2.1 find the boundary halfedge targeting tempV*/
		H* tempV_he = NULL;
		for (auto vhe : It::VClwInHEIterator(pMesh, tempV))
		{
			if (pMesh->halfedgeEdge(vhe)->boundary())
			{
				tempV_he = vhe;
				break;
			}
		}
		if (tempV_he == NULL)
		{
			std::cout << "error! find the boundary halfedge targeting tempV！" << std::endl;
		}
		/*2.2 calculating the difference in angles*/
		V* v1 = tempV;
		V* v2 = pMesh->halfedgeSource(tempV_he);
		V* v3 = pMesh->halfedgeTarget(pMesh->halfedgeNext(tempV_he));
		V* v4 = pMesh->halfedgeTarget(pMesh->halfedgeNext(pMesh->halfedgeSym(pMesh->halfedgeNext(tempV_he))));
		CPoint v12 = v2->point() - v1->point();
		CPoint v13 = v3->point() - v1->point();
		CPoint v14 = v4->point() - v1->point();
		double angle_213 = std::acos(v12 * v13 / (v12.norm() * v13.norm()));
		double angle_413 = std::acos(v14 * v13 / (v14.norm() * v13.norm()));
		CPoint v12_init = v2->initialPos() - v1->initialPos();
		CPoint v13_init = v3->initialPos() - v1->initialPos();
		CPoint v14_init = v4->initialPos() - v1->initialPos();
		double angle_213_init = std::acos(v12_init * v13_init / (v12_init.norm() * v13_init.norm()));
		double angle_413_init = std::acos(v14_init * v13_init / (v14_init.norm() * v13_init.norm()));
		double angleDiff = std::abs(angle_213 - angle_213_init) + std::abs(angle_413 - angle_413_init);

		//std::cout << "M_PI: " << std::acos(-1.0) << std::endl;
		if (angleDiff > std::acos(-1.0) / 3.0)
		{
			std::cout << " poor orthogonality vertex id: " << tempV->id() << std::endl;
			/* 2.3 modifying the vertex's position of poor orthogonality*/
			CPoint v31_init = v1->initialPos() - v3->point();
			CPoint v31 = v1->point() - v3->point();
			tempV->point() = v3->point() + (v31_init / v31_init.norm()) * (v31 * v31_init) / (v31_init.norm());
		}
	}
}

CPoint CCG_QMSLib::BE_bezierCurve_approxi::evaluate(double u) const
{
	// Extract control points (clear naming)
	const CPoint& P0 = be_cpts[0];
	const CPoint& P1 = be_cpts[1];
	const CPoint& P2 = be_cpts[2];
	const CPoint& P3 = be_cpts[3];

	// Explicit formula for cubic Bézier curve (more efficient)
	double u1 = 1.0 - u;
	double u1_sq = u1 * u1;
	double u_sq = u * u;

	return P0 * (u1_sq * u1) +
		P1 * (3 * u1_sq * u) +
		P2 * (3 * u1 * u_sq) +
		P3 * (u_sq * u);
}

void CCG_QMSLib::BE_curve_model::computeTriangleFeaturePointFitting()
{
	// Clear previous records
	originalToApproxPoints.clear();

	// Check required mappings exist
	if (!be_triFeatureLoops || triPointIdToParamMap.empty()) {
		//std::cerr << "Error: Missing required data for triangle feature point fitting\n";
		return;
	}

	// Iterate over all triangle feature loops
	for (auto loop : be_triFeatureLoops->loops()) {
		if (!loop) continue;

		// Iterate over halfedges on the loop
		for (auto he : loop->halfedges()) {
			V* triVertex = be_triFeatureLoops->getMesh()->halfedgeSource(he);
			int vid = triVertex->id();

			// Get original coordinates
			CPoint originalPos = triVertex->point();

			// Get associated quad halfedge
			H* quadHE = triPointToQuadHalfedgeMap[triVertex];
			if (!quadHE) continue;

			// Get parameter value
			double u = triPointIdToParamMap[vid];

			// Get associated Bézier curve
			BE_bezierCurve_approxi* curve = quadHalfedgeToBezierCurveMap[quadHE];
			if (!curve) continue;

			// Calculate fitted coordinates
			CPoint fittedPos = curve->evaluate(u);

			// Store in the mapping
			originalToApproxPoints[vid] = std::make_pair(originalPos, fittedPos);
		}
	}

	// Output statistics
	/*std::cout << "\n[Triangle Feature Point Fitting]\n"
			  << "Processed " << originalToApproxPoints.size()
			  << " triangle feature points\n";*/
}

void CCG_QMSLib::BE_curve_model::initializeQuadHalfedgeIdMap(M* quadMesh)
{
	int id = 0; // Initial ID
	for (auto quadLoop : be_quadFeatureLoops->loops()) {
		for (auto quadHalfEdge : quadLoop->halfedges()) {
			//std::cout << "22222 " << std::endl;
			quadHalfedgeToIdMap[quadHalfEdge] = id++;
			//std::cout <<"loopQuadHalfEdge id: " << quadHalfEdge->localId() << std::endl;
		}
	}
	for (auto quadArc : be_quadFeatureLoops->arcs()) {
		for (auto quadHalfEdge : quadArc->halfedges()) {
			//std::cout << "22222 " << std::endl;
			quadHalfedgeToIdMap[quadHalfEdge] = id++;
			//std::cout << "arcQuadHalfEdge id: " << quadHalfEdge->localId() << std::endl;
		}
	}
	//std::cout << "Quad halfedge ID map initialized!" << std::endl;
}

void CCG_QMSLib::BE_curve_model::initializeQuadHalfedgeIdMap()
{
	int id = 0; // Initial ID
	for (auto quadLoop : be_quadFeatureLoops->loops()) {
		for (auto quadHalfEdge : quadLoop->halfedges()) {
			quadHalfedgeToIdMap[quadHalfEdge] = id++;
			//std::cout << "loopQuadHalfEdge id: " << quadHalfEdge->localId() << std::endl;
		}
	}
	for (auto quadArc : be_quadFeatureLoops->arcs()) {
		for (auto quadHalfEdge : quadArc->halfedges()) {
			quadHalfedgeToIdMap[quadHalfEdge] = id++;
			//std::cout << "arcQuadHalfEdge id: " << quadHalfEdge->localId() << std::endl;
		}
	}
	//std::cout << "Quad halfedge ID map initialized!" << std::endl;
}

void CCG_QMSLib::BE_curve_model::initializeTriHalfedgeIdMap()
{
	int id = 0; // Initial ID
	for (auto triLoop : be_triFeatureLoops->loops()) {
		for (auto triHalfEdge : triLoop->halfedges()) {
			triHalfedgeToIdMap[triHalfEdge] = id++;
			//std::cout << "loopQuadHalfEdge id: " << quadHalfEdge->localId() << std::endl;
		}
	}
}


void CCG_QMSLib::BE_surface_model::computeSamplingPosForTriangleFeatures(M* pMesh, CFeatureLoop* trifeatureLoop, std::map<int, double> ttp, std::map<V*, H*> ttq)
{
	double error_nonFixedBoudnary = 0;

	// Iterate over feature loops
	for (auto loop : trifeatureLoop->loops()) {
		// Ensure current loop is valid
		if (!loop) continue;

		// Iterate over all vertices in the loop
		for (auto he : loop->halfedges()) {
			// Get coordinates of the current vertex
			V* currentPoint = (V*)he->target();
			int vId = currentPoint->id();
			int fId = ttq[currentPoint]->face()->id();
			F* fidf = pMesh->idFace(fId);
			bool tempMarkFaceBoundary = false;
			for (auto fidfe : It::FEIterator(pMesh, fidf))
			{
				if (fidfe->boundary())
				{
					tempMarkFaceBoundary = true;
				}
			}
			if (tempMarkFaceBoundary)
			{
				//std::cout << "123456789------------------------------------------------------------------------" << std::endl;
			}
			// Get UV parameters via function
			CPoint2 uv = getUVFromSamplingPoint(currentPoint, ttp, ttq);

			// Store control point indices and weights
			std::vector<int> index;    // IDs of control points associated with this sampling point
			std::vector<double> weight; // Weights of the current sampling point on its associated control points
			computeRealSamplingPos(fId, vId, uv, index, weight);

			// Compute new sampling point coordinates
			CPoint approx_pos;
			for (int i = 0; i < index.size(); i++) {
				// Accumulate the product of control point positions and weights
				approx_pos += (this->controlPoints()[index[i] - 1] * weight[i]);
			}

			// Output debug information
			/*std::cout << "Triangle feature point " << vId << " approximate position: "
				<< approx_pos[0] << " " << approx_pos[1] << " " << approx_pos[2] << std::endl;
			std::cout << "Triangle feature point " << vId << " original point: "
				<< currentPoint->point()[0] << " " << currentPoint->point()[1] << " " << currentPoint->point()[2] << std::endl;
			std::cout << "Distance to original position: "
				<< (currentPoint->point() - approx_pos).norm() << std::endl;*/

				// Accumulate error
			error_nonFixedBoudnary += (currentPoint->point() - approx_pos).norm();

			// Store the original and fitted points in the map
			originalToApproxPoints[vId] = std::make_pair(currentPoint->point(), approx_pos);

			// Print debug information at key steps
			//std::cout << "UV: " << uv[0] << ", " << uv[1] << std::endl;
		}
	}

	// Output total fitting error
	//std::cout << "Total fitting error: " << error_nonFixedBoudnary << std::endl;

	// Print stored original and fitted points
	for (const auto& entry : originalToApproxPoints) {
		int vId = entry.first;
		const CPoint& originalPoint = entry.second.first;
		const CPoint& approxPoint = entry.second.second;

		/*std::cout << "Vertex ID: " << vId << std::endl;
		std::cout << "  Original Point: " << originalPoint[0] << " " << originalPoint[1] << " " << originalPoint[2] << std::endl;
		std::cout << "  Approx Point:   " << approxPoint[0] << " " << approxPoint[1] << " " << approxPoint[2] << std::endl;*/
	}
}

CPoint2 CCG_QMSLib::BE_surface_model::getUVFromSamplingPoint(V* v, std::map<int, double> ttp, std::map<V*, H*> ttq)
{
	// Get UV parameters via function
	double u = ttp[v->id()];
	H* he = ttq[v];
	int id = he->localId();
	if (id == 2) {
		return CPoint2(u, 0.0);
	}
	else if (id == 3) {
		return CPoint2(1.0, u);
	}
	else if (id == 4) {
		return CPoint2(1.0 - u, 1.0);
	}
	else if (id == 1) {
		return CPoint2(0.0, 1.0 - u);
	}
}


