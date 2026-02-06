#include"CCG_QMS_quad2ManifoldSpline.h"
#include"CCG_QMS_meshOperation.h"
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<chrono>
#include"CCG_QMS_bsplineTool.h"
#include <vector>
#include <map>
#include <omp.h>
#include <Eigen/IterativeLinearSolvers>

//PI
#define _M_PI 3.141592653

/*Some functions ony are used here*/
bool bfIdCompare(std::shared_ptr<CCG_QMSLib::CCG_QMS_bezierSurf> be_f1, std::shared_ptr<CCG_QMSLib::CCG_QMS_bezierSurf> be_f2)
{
	return be_f1->id() < be_f2->id();
}
//Calculate the factorial
double fac(int n)
{
	if (n < 0)
	{
		std::cout << "ERROR! The factorial parameter n should be bigger than 0!" << std::endl;
	}
	double r = 1.0;
	for (int i = n; i >= 1; i--)
	{
		r *= i;
	}
	return r;
}

//Calculate the parameters of the binomial distribution
double binomial(int n, int k)
{
	if (k > n)
	{
		std::cout << "ERROR! The binomial's parameter k should be smaller than n!" << std::endl;
		return 0.0;
	}
	double r = fac(n) / (fac(k) * fac(n - k));
	return r;
}

/*Calculate the sum of the elements in the vector container*/
double vectorSum(std::vector<double> controlWeights)
{
	double s = 0.0;
	for (auto c : controlWeights)
	{
		s += c;
	}
	return s;
}

/**
 * Calculate the values of all n-degree Bernstein polynomials.
 * \param n degree
 * \param u parameter
 * \param B Array, used to store the values of all basis functions
 */
void allBernstein(int n, double u, double maxVal, std::vector<double>& B)
{
	B.resize(n + 1);
	B[0] = 1.0;
	double u1 = (maxVal - u) / (maxVal);
	for (int j = 1; j <= n; j++)
	{
		double saved = 0.0;
		for (int k = 0; k < j; k++)
		{
			double temp = B[k];
			B[k] = saved + u1 * temp;
			saved = u / maxVal * temp;
		}
		B[j] = saved;
	}
}

//Calculate the value of the basis function corresponding to a point on the Bezier surface at the control points.
std::vector<double> bezierSurfcontrolPtsWeights(std::vector<int>& degree, CPoint2 uv, double max_u, double max_v)
{
	/*Here, the default parameter range on the Bezier surface is [0,1]*[0,1]*/
	double umax = max_u, umin = 0.0, vmax = max_v, vmin = 0.0;
	double diff_error = 1e-8;
	if (degree.size() != 2)
	{
		std::cout << "ERROR! The size of parameter(degree) is not 2!" << std::endl;
	}

	/*The values of the basis functions corresponding to the control points in the U direction*/
	//std::cout << "-----------U-----------" << std::endl;
	//std::cout << "u: " << uv[0] << std::endl;
	std::vector<double> uw;
	allBernstein(degree[0], uv[0], umax, uw);
	//std::cout << std::endl;
	//std::cout << "-----------------------" << std::endl;
	//std::cout << "-----------V-----------" << std::endl;
	//std::cout << "v: " << uv[1] << std::endl;
	/*The values of the basis functions corresponding to the control points in the v direction*/
	std::vector<double> vw;
	allBernstein(degree[1], uv[1], vmax, vw);
	//std::cout << std::endl;
	//std::cout << "----------------------" << std::endl;
	//std::cout << "---------U x V-------------" << std::endl;
	//std::cout << "uv: " << uv[0] << " " << uv[1] << std::endl;
	std::vector<double> w;
	for (int i = 0; i <= degree[1]; i++)
	{
		for (int j = 0; j <= degree[0]; j++)
		{
			w.push_back(vw[i] * uw[j]);
			//std::cout << std::setw(12) << vw[i] * uw[j];
		}
		//std::cout << std::endl;
	}
	//std::cout << "----------------------" << std::endl;
	/*check w = 1 ?*/
	if (abs(vectorSum(w) - 1.0) > diff_error)
	{
		std::cout << "basis controlWeights sum: " << vectorSum(w) << std::endl;
	}
	return w;
}


/*--------------------------------functions using to output B spline in .iges format start----------------------------------------------------*/
/*Append the length of the string to the beginning of the current string.*/
std::string HString(std::string str)
{
	std::string hs = std::to_string(str.length()) + "H" + str;
	return hs;
}

/*Merge the string array into a fixed-length string array*/
std::vector<std::string> make_section(std::vector<std::string> fileds, int linewidth)
{
	std::vector<std::string> sec;

	int index = 1;
	std::string line = "";

	int num = fileds.size();

	for (int i = 0; i < num; i++)
	{
		std::string newitem;
		if (i + 1 < num)
		{
			newitem = fileds[i] + ",";
		}
		else
		{
			newitem = fileds[i] + ";";
		}

		int len = line.length() + newitem.length();

		if (len > linewidth)
		{
			sec.push_back(line);
			index++;
			line = "";
		}
		line = line + newitem;
	}
	sec.push_back(line);

	return sec;
}

//DIRECTORY ENTRY SECTION
struct DirectoryEntrySectionEle
{
	int type = 128;//128,NURBS
	int id;
	int p_start;
	int p_count;
};

//以字符串数组的形式获取一个NURBS曲面片的信息
void make_section_array(std::shared_ptr <CCG_QMSLib::CCG_QMS_bsplineSurf> nurbs, double USpanMin, double USpanMax, double VSpanMin, double VSpanMax, std::vector<std::string>& P)
{
	int dim = nurbs->degree().size();

	//in IGES the control points are stored in the format [x, y, z, w]
	//instead of [w*x, w*y, w*z, w]
	//default w = 1.0，

	std::vector<std::vector<CPoint>>cp = nurbs->cpts();
	int degU = nurbs->degree()[0];
	int degV = nurbs->degree()[1];
	std::vector<double> knotsU = nurbs->knotsU();
	std::vector<double> knotsV = nurbs->knotsV();
	std::vector<double> uspan = { USpanMin ,USpanMax };
	std::vector<double> vspan = { VSpanMin ,VSpanMax };
	P.push_back("128"); //NURBS surface
	//P.push_back(std::to_string(nurbs->cpts()[0].size() - 1));//Number of control points in U
	//P.push_back(std::to_string(nurbs->cpts().size() - 1));//Number of control points in V

	/*Modify so that the surface direction is consistent with the direction of the quadrilateral mesh,2025/3/1*/
	P.push_back(std::to_string(nurbs->cpts().size() - 1));//Number of control points in U
	P.push_back(std::to_string(nurbs->cpts()[0].size() - 1));//Number of control points in V

	P.push_back(std::to_string(degU)); //Degree in U
	P.push_back(std::to_string(degV));//Degree in V
	/*U is closed? 0(No),1(Yes)*/
	if (nurbs->closedAlongU())
	{
		P.push_back("1");
	}
	else
	{
		P.push_back("0");
	}
	/*V is closed? 0(No),1(Yes)*/
	if (nurbs->closedAlongV())
	{
		P.push_back("1");
	}
	else
	{
		P.push_back("0");
	}
	/*0 is a rational number, and 1 is a polynomial.*/
	P.push_back("0");
	/*U is period? 0(No),1(Yes)*/
	if (nurbs->periodicAlongU())
	{
		P.push_back("1");
	}
	else
	{
		P.push_back("0");
	}

	/*V is period? 0(No),1(Yes)*/
	if (nurbs->periodicAlongV())
	{
		P.push_back("1");
	}
	else
	{
		P.push_back("0");
	}
	/*knot vector*/
	for (int i = 0; i < knotsU.size(); i++)
	{
		P.push_back(std::to_string(knotsU[i])+"D0");
	}
	for (int i = 0; i < knotsV.size(); i++)
	{
		P.push_back(std::to_string(knotsV[i]) + "D0");
	}
	/*weight w*/
	for (int i = 0; i < nurbs->cpts().size(); i++)
	{
		for (int j = 0; j < nurbs->cpts()[0].size(); j++)
		{
			P.push_back(std::to_string(1.0) + "D0");
		}
	}
	/*control points coordinate x,y,z*/
	/*for (int i = 0; i < nurbs->cpts().size(); i++)
	{
		for (int j = 0; j < nurbs->cpts()[0].size(); j++)
		{
			for (int k = 0; k < 3; k++)
			{
				P.push_back(std::to_string(nurbs->cpts()[i][j][k]));
			}
		}
	}*/
	/*control points coordinate x,y,z,output in the order of first U direction and then V direction.2025/3/1*/
	for (int j = 0; j < nurbs->cpts()[0].size(); j++)
	{
		for (int i = 0; i < nurbs->cpts().size(); i++)
		{
			for (int k = 0; k < 3; k++)
			{
				std::ostringstream oss;
				oss << std::fixed << std::setprecision(19) << nurbs->cpts()[i][j][k];
				//P.push_back(std::to_string(nurbs->cpts()[i][j][k]));
				P.push_back(oss.str() + "D0");
			}
		}
	}

	/*U parameterization domain*/
	P.push_back(std::to_string(uspan[0]) + "D0");
	P.push_back(std::to_string(uspan[1]) + "D0");
	/*V parameterization domain*/
	P.push_back(std::to_string(vspan[0]) + "D0");
	P.push_back(std::to_string(vspan[1]) + "D0");
	/*P.push_back("0");
	P.push_back("0");*/

}


/*--------------------------------functions using to output B spline in .iges format end----------------------------------------------------*/

/*-------------------function for CCG_QMS_model start------------------------------------*/
void CCG_QMSLib::CCG_QMS_model::obtainCtlPts(M* pMesh)
{
	//intialize the number of control points
	this->controlPoints().resize(pMesh->numVertices());
	for(int i = 0; i < pMesh->numVertices(); i++)
	{
		this->controlPoints()[i] = pMesh->idVertex(i+1)->point();
	}
}

void CCG_QMSLib::CCG_QMS_model::buildBezierSurfs_cw(M* pMesh)
{
	/*Marking singularities*/
	markExtraordinaryPts_cw(pMesh);

	for (auto f : It::MFIterator(pMesh))
	{
		/*
		* build bezier surface according to quad face type
		* First, build face points:(9)(10)(11)(12)
		*/
		auto be_f = std::make_shared<CCG_QMS_bezierSurf>();
		if (!be_f) {
			// Handle allocation failure
			return;
		}
		/*
		* specify bezier surface type
		*/
		for (auto fv : It::FVIterator(pMesh, f))
		{
			if (fv->ifSingular())
			{
				be_f->type() = SingularSurf;
			}
		}
		
		be_f->degree(0) = 3;
		be_f->degree(1) = 3;

		/*
		* Give the bezier surface id,
		* identical to the quad face id
		*/
		be_f->id() = f->id();

		/*
		* Compute the face control points on bezier surface
		*/
		/*
		* Face control points' local id is 6,7,10,11
		* ///////////////////////////////////////////
		*     |                     |
		*     |                     |
		*     |<--faceHalfedge(f)---|
		* ---PC--------------------PD---
		*     |  10(B23)  11B(33)   |
		*     |                     |
		*     |  6(B22)   7(B32)    |
		* ---PA--------------------PB---
		*     |                     |
		*     |                     |
		*     |                     |
		* ///////////////////////////////////////////
		*/
		/*Initialize the face control points*/
		CPoint B32, B22, B23, B33;
		/*Obtain the vertex on the quad face*/
		CPoint pA, pB, pC, pD;
		V* PA; V* PB; V* PC; V* PD;
		/*Obtain the halfedges on the quad face*/
		H* f_h = pMesh->faceHalfedge(f);
		pC = f_h->target()->point();
		PC = pMesh->halfedgeTarget(f_h);
		pD = f_h->source()->point();
		PD = pMesh->halfedgeSource(f_h);
		H* f_hNh = pMesh->halfedgeNext(f_h);
		H* f_hNhNh = pMesh->halfedgeNext(f_hNh);
		pB = f_hNhNh->target()->point();
		PB = pMesh->halfedgeTarget(f_hNhNh);
		pA = f_hNhNh->source()->point();
		PA = pMesh->halfedgeSource(f_hNhNh);
		/*Computing face control points according to equation 9-12*/
		B22 = pA * 4.0 / 9 + pB * 2.0 / 9 + pC * 2.0 / 9 + pD * 1.0 / 9;
		B32 = pA * 2.0 / 9 + pB * 4.0 / 9 + pC * 1.0 / 9 + pD * 2.0 / 9;
		B23 = pA * 2.0 / 9 + pB * 1.0 / 9 + pC * 4.0 / 9 + pD * 2.0 / 9;
		B33 = pA * 1.0 / 9 + pB * 2.0 / 9 + pC * 2.0 / 9 + pD * 4.0 / 9;

		BE_linkingCtlPoints linkCW_B22, linkCW_B32, linkCW_B23, linkCW_B33;
		/*-----------------------------------------------------------------------------------------------------------------------------*/
		/*B22*/
		/*PA*/
		{
			std::pair<int, double> templinkCW(PA->id(), 4.0 / 9);
			if (linkCW_B22.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B22.linkingCtlPId_weights.end())
			{
				linkCW_B22.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f0It = linkCW_B22.linkingCtlPId_weights.find(templinkCW.first);
				f0It->second += templinkCW.second;
			}
		}
		/*PC*/
		{
			std::pair<int, double> templinkCW(PC->id(), 2.0 / 9);
			if (linkCW_B22.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B22.linkingCtlPId_weights.end())
			{
				linkCW_B22.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f0It = linkCW_B22.linkingCtlPId_weights.find(templinkCW.first);
				f0It->second += templinkCW.second;
			}
		}
		/*PB*/
		{
			std::pair<int, double> templinkCW(PB->id(), 2.0 / 9);
			if (linkCW_B22.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B22.linkingCtlPId_weights.end())
			{
				linkCW_B22.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f0It = linkCW_B22.linkingCtlPId_weights.find(templinkCW.first);
				f0It->second += templinkCW.second;
			}
		}
		/*PD*/
		{
			std::pair<int, double> templinkCW(PD->id(), 1.0 / 9);
			if (linkCW_B22.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B22.linkingCtlPId_weights.end())
			{
				linkCW_B22.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f0It = linkCW_B22.linkingCtlPId_weights.find(templinkCW.first);
				f0It->second += templinkCW.second;
			}
		}
		/*check B22*/
		/*{
			CPoint tempP;
			for (std::map<int, double>::iterator it = linkCW_B22.linkingCtlPId_weights.begin(); it != linkCW_B22.linkingCtlPId_weights.end(); it++)
			{
				tempP += pMesh->idVertex(it->first)->point() * it->second;
			}
			std::cout << " F0: " << F0[0] << " " << F0[1] << " " << F0[2] << std::endl;
			std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
		}*/
		/*-----------------------------------------------------------------------------------------------------------------------------*/
		/*B32*/
		/*PA*/
		{
			std::pair<int, double> templinkCW(PA->id(), 2.0 / 9);
			if (linkCW_B32.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B32.linkingCtlPId_weights.end())
			{
				linkCW_B32.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f2It = linkCW_B32.linkingCtlPId_weights.find(templinkCW.first);
				f2It->second += templinkCW.second;
			}
		}
		/*PC*/
		{
			std::pair<int, double> templinkCW(PC->id(), 1.0 / 9);
			if (linkCW_B32.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B32.linkingCtlPId_weights.end())
			{
				linkCW_B32.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f2It = linkCW_B32.linkingCtlPId_weights.find(templinkCW.first);
				f2It->second += templinkCW.second;
			}
		}
		/*PB*/
		{
			std::pair<int, double> templinkCW(PB->id(), 4.0 / 9);
			if (linkCW_B32.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B32.linkingCtlPId_weights.end())
			{
				linkCW_B32.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f2It = linkCW_B32.linkingCtlPId_weights.find(templinkCW.first);
				f2It->second += templinkCW.second;
			}
		}
		/*PD*/
		{
			std::pair<int, double> templinkCW(PD->id(), 2.0 / 9);
			if (linkCW_B32.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B32.linkingCtlPId_weights.end())
			{
				linkCW_B32.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f2It = linkCW_B32.linkingCtlPId_weights.find(templinkCW.first);
				f2It->second += templinkCW.second;
			}
		}
		/*check B32*/
		/*{
			CPoint tempP;
			for (std::map<int, double>::iterator it = linkCW_B32.linkingCtlPId_weights.begin(); it != linkCW_B32.linkingCtlPId_weights.end(); it++)
			{
				tempP += pMesh->idVertex(it->first)->point() * it->second;
			}
			std::cout << " F2: " << F2[0] << " " << F2[1] << " " << F2[2] << std::endl;
			std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
		}*/
		/*-----------------------------------------------------------------------------------------------------------------------------*/
		/*B23*/
		/*PA*/
		{
			std::pair<int, double> templinkCW(PA->id(), 2.0 / 9);
			if (linkCW_B23.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B23.linkingCtlPId_weights.end())
			{
				linkCW_B23.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f1It = linkCW_B23.linkingCtlPId_weights.find(templinkCW.first);
				f1It->second += templinkCW.second;
			}
		}
		/*PC*/
		{
			std::pair<int, double> templinkCW(PC->id(), 4.0 / 9);
			if (linkCW_B23.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B23.linkingCtlPId_weights.end())
			{
				linkCW_B23.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f1It = linkCW_B23.linkingCtlPId_weights.find(templinkCW.first);
				f1It->second += templinkCW.second;
			}
		}
		/*PB*/
		{
			std::pair<int, double> templinkCW(PB->id(), 1.0 / 9);
			if (linkCW_B23.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B23.linkingCtlPId_weights.end())
			{
				linkCW_B23.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f1It = linkCW_B23.linkingCtlPId_weights.find(templinkCW.first);
				f1It->second += templinkCW.second;
			}
		}
		/*PD*/
		{
			std::pair<int, double> templinkCW(PD->id(), 2.0 / 9);
			if (linkCW_B23.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B23.linkingCtlPId_weights.end())
			{
				linkCW_B23.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f1It = linkCW_B23.linkingCtlPId_weights.find(templinkCW.first);
				f1It->second += templinkCW.second;
			}
		}
		/*check B23*/
		/*{
			CPoint tempP;
			for (std::map<int, double>::iterator it = linkCW_B23.linkingCtlPId_weights.begin(); it != linkCW_B23.linkingCtlPId_weights.end(); it++)
			{
				tempP += pMesh->idVertex(it->first)->point() * it->second;
			}
			std::cout << " F1: " << F1[0] << " " << F1[1] << " " << F1[2] << std::endl;
			std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
		}*/
		/*-----------------------------------------------------------------------------------------------------------------------------*/
		/*B33*/
		/*PA*/
		{
			std::pair<int, double> templinkCW(PA->id(), 1.0 / 9);
			if (linkCW_B33.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B33.linkingCtlPId_weights.end())
			{
				linkCW_B33.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f3It = linkCW_B33.linkingCtlPId_weights.find(templinkCW.first);
				f3It->second += templinkCW.second;
			}
		}
		/*PC*/
		{
			std::pair<int, double> templinkCW(PC->id(), 2.0 / 9);
			if (linkCW_B33.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B33.linkingCtlPId_weights.end())
			{
				linkCW_B33.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f3It = linkCW_B33.linkingCtlPId_weights.find(templinkCW.first);
				f3It->second += templinkCW.second;
			}
		}
		/*PB*/
		{
			std::pair<int, double> templinkCW(PB->id(), 2.0 / 9);
			if (linkCW_B33.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B33.linkingCtlPId_weights.end())
			{
				linkCW_B33.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f3It = linkCW_B33.linkingCtlPId_weights.find(templinkCW.first);
				f3It->second += templinkCW.second;
			}
		}
		/*PD*/
		{
			std::pair<int, double> templinkCW(PD->id(), 4.0 / 9);
			if (linkCW_B33.linkingCtlPId_weights.find(templinkCW.first) == linkCW_B33.linkingCtlPId_weights.end())
			{
				linkCW_B33.linkingCtlPId_weights.insert(templinkCW);
			}
			else
			{
				std::map<int, double>::iterator f3It = linkCW_B33.linkingCtlPId_weights.find(templinkCW.first);
				f3It->second += templinkCW.second;
			}
		}
		/*check B33*/
		/*{
			CPoint tempP;
			for (std::map<int, double>::iterator it = linkCW_B33.linkingCtlPId_weights.begin(); it != linkCW_B33.linkingCtlPId_weights.end(); it++)
			{
				tempP += pMesh->idVertex(it->first)->point() * it->second;
			}
			std::cout << " F3: " << F3[0] << " " << F3[1] << " " << F3[2] << std::endl;
			std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
		}*/
		/*-----------------------------------------------------------------------------------------------------------------------------*/


		/*Face control points*/
		/*Check the number of face control points*/
		//std::cout << "be_f->cpts().size(): " << be_f->cpts().size() << std::endl;
		be_f->cpts()[5] = B22;
		be_f->cpts()[6] = B32;
		be_f->cpts()[9] = B23;
		be_f->cpts()[10] = B33;

		/*obtain mapping between face contro points and quad mesh vertex*/
		be_f->linkCptPWs()[5] = linkCW_B22;//6
		be_f->linkCptPWs()[6] = linkCW_B32;//7
		be_f->linkCptPWs()[9] = linkCW_B23;//10
		be_f->linkCptPWs()[10] = linkCW_B33;//11

		this->bfs().push_back(be_f);
	}
	//std::cout << "Obtaining faces' points succesfully!" << std::endl;

	/*Sort the created bezier surfaces*/
	sortBezierSurfaces();
	markHalfedgeLocalId(pMesh);
	/*Edge points*/
	for (auto e : It::MEIterator(pMesh))
	{
		BE_linkingCtlPoints ePW0, ePW1;
		/* Building edge points(13)(14)(16)(17)*/
		if (e->boundary())
		{
			/*Equation(16) (17)*/
			CPoint PA, PB;
			CPoint B21, B31;
			PA = e->halfedge(0)->source()->point();
			V* pt1 = pMesh->halfedgeSource(pMesh->edgeHalfedge(e, 0));
			PB = e->halfedge(0)->target()->point();
			V* pt2 = pMesh->halfedgeTarget(pMesh->edgeHalfedge(e, 0));
			H* eh = pMesh->edgeHalfedge(e, 0);
			F* ehf = pMesh->halfedgeFace(eh);
			auto ehf_beSurf = this->bfs()[ehf->id() - 1];
			//Equation 16,17 computing
			B21 = PA * 2.0 / 3 + PB * 1.0 / 3;
			B31 = PA * 1.0 / 3 + PB * 2.0 / 3;
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
				/*E0*/
				/*the first edge control point*/
			{
				std::pair<int, double> tempPW(pt1->id(), (2.0 / 3));
				if (ePW0.linkingCtlPId_weights.find(pt1->id()) == ePW0.linkingCtlPId_weights.end())
				{
					ePW0.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*the second edge control point*/
			{
				std::pair<int, double> tempPW(pt2->id(), (1.0 / 3));
				if (ePW0.linkingCtlPId_weights.find(pt2->id()) == ePW0.linkingCtlPId_weights.end())
				{
					ePW0.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*check*/
			/*{
				CPoint tempP;
				for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
				{
					tempP += pMesh->idVertex(it->first)->point() * it->second;
				}
				std::cout << " bfs()[f->id() - 1]->cpts()[4]: " << bfs()[f->id() - 1]->cpts()[4][0] << " " << bfs()[f->id() - 1]->cpts()[4][1] << " " << bfs()[f->id() - 1]->cpts()[4][2] << std::endl;
				std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
			}*/
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
			/*E1*/
			/*the first edge control point*/
			{
				std::pair<int, double> tempPW(pt1->id(), (1.0 / 3));
				if (ePW1.linkingCtlPId_weights.find(pt1->id()) == ePW1.linkingCtlPId_weights.end())
				{
					ePW1.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*the second edge control point*/
			{
				std::pair<int, double> tempPW(pt2->id(), (2.0 / 3));
				if (ePW1.linkingCtlPId_weights.find(pt2->id()) == ePW1.linkingCtlPId_weights.end())
				{
					ePW1.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*check*/
			/*{
				CPoint tempP;
				for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
				{
					tempP += pMesh->idVertex(it->first)->point() * it->second;
				}
				std::cout << " bfs()[f->id() - 1]->cpts()[8]: " << bfs()[f->id() - 1]->cpts()[8][0] << " " << bfs()[f->id() - 1]->cpts()[8][1] << " " << bfs()[f->id() - 1]->cpts()[8][2] << std::endl;
				std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
			}*/
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/

			switch (eh->localId())
			{
			case 0:
			{
				ehf_beSurf->cpts()[1] = B21;
				ehf_beSurf->cpts()[2] = B31;
				ehf_beSurf->linkCptPWs()[1] = ePW0;
				ehf_beSurf->linkCptPWs()[2] = ePW1;
				break;
			}
			case 1:
			{
				ehf_beSurf->cpts()[7] = B21;
				ehf_beSurf->cpts()[11] = B31;
				ehf_beSurf->linkCptPWs()[7] = ePW0;
				ehf_beSurf->linkCptPWs()[11] = ePW1;
				break;
			}
			case 2:
			{
				ehf_beSurf->cpts()[14] = B21;
				ehf_beSurf->cpts()[13] = B31;
				ehf_beSurf->linkCptPWs()[14] = ePW0;
				ehf_beSurf->linkCptPWs()[13] = ePW1;
				break;
			}
			case 3:
			{
				ehf_beSurf->cpts()[8] = B21;
				ehf_beSurf->cpts()[4] = B31;
				ehf_beSurf->linkCptPWs()[8] = ePW0;
				ehf_beSurf->linkCptPWs()[4] = ePW1;
				break;
			}
			default:
				break;
			}

		}
		else if (e->feature())
		{
			/*equation(16) (17)*/
			CPoint PA, PB;
			CPoint B21, B31;
			PA = e->halfedge(0)->source()->point();
			V* pt1 = pMesh->halfedgeSource(pMesh->edgeHalfedge(e, 0));
			PB = e->halfedge(0)->target()->point();
			V* pt2 = pMesh->halfedgeTarget(pMesh->edgeHalfedge(e, 0));
			H* eh = pMesh->edgeHalfedge(e, 0);
			F* ehf = pMesh->halfedgeFace(eh);
			auto ehf_beSurf = this->bfs()[ehf->id() - 1];
			//equation 16,17 computing
			B21 = PA * 2.0 / 3 + PB * 1.0 / 3;
			B31 = PA * 1.0 / 3 + PB * 2.0 / 3;
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
				/*E0*/
				/*the first edge control point*/
			{
				std::pair<int, double> tempPW(pt1->id(), (2.0 / 3));
				if (ePW0.linkingCtlPId_weights.find(pt1->id()) == ePW0.linkingCtlPId_weights.end())
				{
					ePW0.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*the second edge control point*/
			{
				std::pair<int, double> tempPW(pt2->id(), (1.0 / 3));
				if (ePW0.linkingCtlPId_weights.find(pt2->id()) == ePW0.linkingCtlPId_weights.end())
				{
					ePW0.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*check*/
			/*{
				CPoint tempP;
				for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
				{
					tempP += pMesh->idVertex(it->first)->point() * it->second;
				}
				std::cout << " bfs()[f->id() - 1]->cpts()[4]: " << bfs()[f->id() - 1]->cpts()[4][0] << " " << bfs()[f->id() - 1]->cpts()[4][1] << " " << bfs()[f->id() - 1]->cpts()[4][2] << std::endl;
				std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
			}*/
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
			/*E1*/
			/*the first edge control point*/
			{
				std::pair<int, double> tempPW(pt1->id(), (1.0 / 3));
				if (ePW1.linkingCtlPId_weights.find(pt1->id()) == ePW1.linkingCtlPId_weights.end())
				{
					ePW1.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*the second edge control point*/
			{
				std::pair<int, double> tempPW(pt2->id(), (2.0 / 3));
				if (ePW1.linkingCtlPId_weights.find(pt2->id()) == ePW1.linkingCtlPId_weights.end())
				{
					ePW1.linkingCtlPId_weights.insert(tempPW);
				}
				else
				{
					std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(pt1->id());
					eIt->second += tempPW.second;
				}
			}
			/*check*/
			/*{
				CPoint tempP;
				for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
				{
					tempP += pMesh->idVertex(it->first)->point() * it->second;
				}
				std::cout << " bfs()[f->id() - 1]->cpts()[8]: " << bfs()[f->id() - 1]->cpts()[8][0] << " " << bfs()[f->id() - 1]->cpts()[8][1] << " " << bfs()[f->id() - 1]->cpts()[8][2] << std::endl;
				std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
			}*/
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/

			switch (eh->localId())
			{
			case 0:
			{
				ehf_beSurf->cpts()[1] = B21;
				ehf_beSurf->cpts()[2] = B31;
				ehf_beSurf->linkCptPWs()[1] = ePW0;
				ehf_beSurf->linkCptPWs()[2] = ePW1;
				break;
			}
			case 1:
			{
				ehf_beSurf->cpts()[7] = B21;
				ehf_beSurf->cpts()[11] = B31;
				ehf_beSurf->linkCptPWs()[7] = ePW0;
				ehf_beSurf->linkCptPWs()[11] = ePW1;
				break;
			}
			case 2:
			{
				ehf_beSurf->cpts()[14] = B21;
				ehf_beSurf->cpts()[13] = B31;
				ehf_beSurf->linkCptPWs()[14] = ePW0;
				ehf_beSurf->linkCptPWs()[13] = ePW1;
				break;
			}
			case 3:
			{
				ehf_beSurf->cpts()[8] = B21;
				ehf_beSurf->cpts()[4] = B31;
				ehf_beSurf->linkCptPWs()[8] = ePW0;
				ehf_beSurf->linkCptPWs()[4] = ePW1;
				break;
			}
			default:
				break;
			}

			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
			/*The symmetric halfedge's control points*/
			eh = pMesh->halfedgeSym(eh);
			ehf = pMesh->halfedgeFace(eh);
			ehf_beSurf = this->bfs()[ehf->id() - 1];
			switch (eh->localId())
			{
			case 0:
			{
				ehf_beSurf->cpts()[2] = B21;
				ehf_beSurf->cpts()[1] = B31;
				ehf_beSurf->linkCptPWs()[2] = ePW0;
				ehf_beSurf->linkCptPWs()[1] = ePW1;
				break;
			}
			case 1:
			{
				ehf_beSurf->cpts()[11] = B21;
				ehf_beSurf->cpts()[7] = B31;
				ehf_beSurf->linkCptPWs()[11] = ePW0;
				ehf_beSurf->linkCptPWs()[7] = ePW1;
				break;
			}
			case 2:
			{
				ehf_beSurf->cpts()[13] = B21;
				ehf_beSurf->cpts()[14] = B31;
				ehf_beSurf->linkCptPWs()[13] = ePW0;
				ehf_beSurf->linkCptPWs()[14] = ePW1;
				break;
			}
			case 3:
			{
				ehf_beSurf->cpts()[4] = B21;
				ehf_beSurf->cpts()[8] = B31;
				ehf_beSurf->linkCptPWs()[4] = ePW0;
				ehf_beSurf->linkCptPWs()[8] = ePW1;
				break;
			}
			default:
				break;
			}
			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
		}
		else//not boundary/feature edges
		{
			//std::cout << "---regular edges---" << std::endl;
			/*equation (13)(14)*/
			H* eh = pMesh->edgeHalfedge(e, 0);
			F* ehf = pMesh->halfedgeFace(eh);
			auto ehf_beSurf = this->bfs()[ehf->id() - 1];
			H* ehSym = pMesh->halfedgeSym(eh);
			F* ehSymf = pMesh->halfedgeFace(ehSym);
			auto ehSymf_beSurf = this->bfs()[ehSymf->id() - 1];

			/*declare edge control points and preparing face control points*/
			CPoint B12, B13, B42, B43;
			CPoint B22, B23, B33, B32;
			switch (eh->localId())
			{
			case 0://halfedge local id is 0
			{
				switch (ehSym->localId())
				{
				case 0://symmetric halfedge local id 0
				{
					//local id 0,0
					//std::cout << "-----00-----" << std::endl;
					B32 = ehf_beSurf->cpts()[5];
					B22 = ehSymf_beSurf->cpts()[6];
					B33 = ehf_beSurf->cpts()[6];
					B23 = ehSymf_beSurf->cpts()[5];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[1] = ePW0;
					ehf_beSurf->linkCptPWs()[2] = ePW1;
					ehSymf_beSurf->linkCptPWs()[2] = ePW0;
					ehSymf_beSurf->linkCptPWs()[1] = ePW1;

					ehf_beSurf->cpts()[1] = B12;
					ehf_beSurf->cpts()[2] = B13;
					ehSymf_beSurf->cpts()[2] = B12;
					ehSymf_beSurf->cpts()[1] = B13;
					break;
				}
				case 1:
				{
					//std::cout << "-----01-----" << std::endl;
					B32 = ehf_beSurf->cpts()[5];
					B22 = ehSymf_beSurf->cpts()[10];
					B33 = ehf_beSurf->cpts()[6];
					B23 = ehSymf_beSurf->cpts()[6];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[1] = ePW0;
					ehf_beSurf->linkCptPWs()[2] = ePW1;
					ehSymf_beSurf->linkCptPWs()[11] = ePW0;
					ehSymf_beSurf->linkCptPWs()[7] = ePW1;

					ehf_beSurf->cpts()[1] = B12;
					ehf_beSurf->cpts()[2] = B13;
					ehSymf_beSurf->cpts()[11] = B12;
					ehSymf_beSurf->cpts()[7] = B13;
					break;
				}
				case 2:
				{
					//std::cout << "-----02-----" << std::endl;
					B32 = ehf_beSurf->cpts()[5];
					B22 = ehSymf_beSurf->cpts()[9];
					B33 = ehf_beSurf->cpts()[6];
					B23 = ehSymf_beSurf->cpts()[10];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[1] = ePW0;
					ehf_beSurf->linkCptPWs()[2] = ePW1;
					ehSymf_beSurf->linkCptPWs()[13] = ePW0;
					ehSymf_beSurf->linkCptPWs()[14] = ePW1;

					ehf_beSurf->cpts()[1] = B12;
					ehf_beSurf->cpts()[2] = B13;
					ehSymf_beSurf->cpts()[13] = B12;
					ehSymf_beSurf->cpts()[14] = B13;
					break;
				}
				case 3:
				{
					//std::cout << "-----03-----" << std::endl;
					B32 = ehf_beSurf->cpts()[5];
					B22 = ehSymf_beSurf->cpts()[5];
					B33 = ehf_beSurf->cpts()[6];
					B23 = ehSymf_beSurf->cpts()[9];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[1] = ePW0;
					ehf_beSurf->linkCptPWs()[2] = ePW1;
					ehSymf_beSurf->linkCptPWs()[4] = ePW0;
					ehSymf_beSurf->linkCptPWs()[8] = ePW1;

					ehf_beSurf->cpts()[1] = B12;
					ehf_beSurf->cpts()[2] = B13;
					ehSymf_beSurf->cpts()[4] = B12;
					ehSymf_beSurf->cpts()[8] = B13;
					break;
				}
				default:
					break;
				}
				break;
			}

			case 1:
			{
				switch (ehSym->localId())
				{
				case 0:
				{
					//std::cout << "-----10-----" << std::endl;
					B32 = ehf_beSurf->cpts()[6];
					B22 = ehSymf_beSurf->cpts()[6];
					B33 = ehf_beSurf->cpts()[10];
					B23 = ehSymf_beSurf->cpts()[5];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[7] = ePW0;
					ehf_beSurf->linkCptPWs()[11] = ePW1;
					ehSymf_beSurf->linkCptPWs()[2] = ePW0;
					ehSymf_beSurf->linkCptPWs()[1] = ePW1;

					ehf_beSurf->cpts()[7] = B12;
					ehf_beSurf->cpts()[11] = B13;
					ehSymf_beSurf->cpts()[2] = B12;
					ehSymf_beSurf->cpts()[1] = B13;
					break;
				}
				case 1:
				{
					//std::cout << "-----11-----" << std::endl;
					B32 = ehf_beSurf->cpts()[6];
					B22 = ehSymf_beSurf->cpts()[10];
					B33 = ehf_beSurf->cpts()[10];
					B23 = ehSymf_beSurf->cpts()[6];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[7] = ePW0;
					ehf_beSurf->linkCptPWs()[11] = ePW1;
					ehSymf_beSurf->linkCptPWs()[11] = ePW0;
					ehSymf_beSurf->linkCptPWs()[7] = ePW1;

					ehf_beSurf->cpts()[7] = B12;
					ehf_beSurf->cpts()[11] = B13;
					ehSymf_beSurf->cpts()[11] = B12;
					ehSymf_beSurf->cpts()[7] = B13;
					break;
				}
				case 2:
				{
					//std::cout << "-----12-----" << std::endl;
					B32 = ehf_beSurf->cpts()[6];
					B22 = ehSymf_beSurf->cpts()[9];
					B33 = ehf_beSurf->cpts()[10];
					B23 = ehSymf_beSurf->cpts()[10];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the fisrt face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[7] = ePW0;
					ehf_beSurf->linkCptPWs()[11] = ePW1;
					ehSymf_beSurf->linkCptPWs()[13] = ePW0;
					ehSymf_beSurf->linkCptPWs()[14] = ePW1;

					ehf_beSurf->cpts()[7] = B12;
					ehf_beSurf->cpts()[11] = B13;
					ehSymf_beSurf->cpts()[13] = B12;
					ehSymf_beSurf->cpts()[14] = B13;
					break;
				}
				case 3:
				{
					//std::cout << "-----13-----" << std::endl;
					B32 = ehf_beSurf->cpts()[6];
					B22 = ehSymf_beSurf->cpts()[5];
					B33 = ehf_beSurf->cpts()[10];
					B23 = ehSymf_beSurf->cpts()[9];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[7] = ePW0;
					ehf_beSurf->linkCptPWs()[11] = ePW1;
					ehSymf_beSurf->linkCptPWs()[4] = ePW0;
					ehSymf_beSurf->linkCptPWs()[8] = ePW1;

					ehf_beSurf->cpts()[7] = B12;
					ehf_beSurf->cpts()[11] = B13;
					ehSymf_beSurf->cpts()[4] = B12;
					ehSymf_beSurf->cpts()[8] = B13;
					break;
				}
				default:
					break;
				}
				break;
			}

			case 2:
			{
				switch (ehSym->localId())
				{
				case 0:
				{
					//break;
					//std::cout << "-----20-----" << std::endl;
					B32 = ehf_beSurf->cpts()[10];
					B22 = ehSymf_beSurf->cpts()[6];
					B33 = ehf_beSurf->cpts()[9];
					B23 = ehSymf_beSurf->cpts()[5];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[14] = ePW0;
					ehf_beSurf->linkCptPWs()[13] = ePW1;
					ehSymf_beSurf->linkCptPWs()[2] = ePW0;
					ehSymf_beSurf->linkCptPWs()[1] = ePW1;

					ehf_beSurf->cpts()[14] = B12;
					ehf_beSurf->cpts()[13] = B13;
					ehSymf_beSurf->cpts()[2] = B12;
					ehSymf_beSurf->cpts()[1] = B13;
					break;
				}
				case 1:
				{
					//std::cout << "-----21-----" << std::endl;
					B32 = ehf_beSurf->cpts()[10];
					B22 = ehSymf_beSurf->cpts()[10];
					B33 = ehf_beSurf->cpts()[9];
					B23 = ehSymf_beSurf->cpts()[6];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[14] = ePW0;
					ehf_beSurf->linkCptPWs()[13] = ePW1;
					ehSymf_beSurf->linkCptPWs()[11] = ePW0;
					ehSymf_beSurf->linkCptPWs()[7] = ePW1;

					ehf_beSurf->cpts()[14] = B12;
					ehf_beSurf->cpts()[13] = B13;
					ehSymf_beSurf->cpts()[11] = B12;
					ehSymf_beSurf->cpts()[7] = B13;
					break;
				}
				case 2:
				{
					//std::cout << "-----22-----" << std::endl;
					B32 = ehf_beSurf->cpts()[10];
					B22 = ehSymf_beSurf->cpts()[9];
					B33 = ehf_beSurf->cpts()[9];
					B23 = ehSymf_beSurf->cpts()[10];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[14] = ePW0;
					ehf_beSurf->linkCptPWs()[13] = ePW1;
					ehSymf_beSurf->linkCptPWs()[13] = ePW0;
					ehSymf_beSurf->linkCptPWs()[14] = ePW1;

					ehf_beSurf->cpts()[14] = B12;
					ehf_beSurf->cpts()[13] = B13;
					ehSymf_beSurf->cpts()[13] = B12;
					ehSymf_beSurf->cpts()[14] = B13;
					break;
				}
				case 3:
				{
					//std::cout << "-----23-----" << std::endl;
					B32 = ehf_beSurf->cpts()[10];
					B22 = ehSymf_beSurf->cpts()[5];
					B33 = ehf_beSurf->cpts()[9];
					B23 = ehSymf_beSurf->cpts()[9];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[14] = ePW0;
					ehf_beSurf->linkCptPWs()[13] = ePW1;
					ehSymf_beSurf->linkCptPWs()[4] = ePW0;
					ehSymf_beSurf->linkCptPWs()[8] = ePW1;

					ehf_beSurf->cpts()[14] = B12;
					ehf_beSurf->cpts()[13] = B13;
					ehSymf_beSurf->cpts()[4] = B12;
					ehSymf_beSurf->cpts()[8] = B13;
					break;
				}
				default:
					break;
				}
				break;
			}

			case 3:
			{

				switch (ehSym->localId())
				{
				case 0:
				{
					//std::cout << "-----30-----" << std::endl;
					B32 = ehf_beSurf->cpts()[9];
					B22 = ehSymf_beSurf->cpts()[6];
					B33 = ehf_beSurf->cpts()[5];
					B23 = ehSymf_beSurf->cpts()[5];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[8] = ePW0;
					ehf_beSurf->linkCptPWs()[4] = ePW1;
					ehSymf_beSurf->linkCptPWs()[2] = ePW0;
					ehSymf_beSurf->linkCptPWs()[1] = ePW1;

					ehf_beSurf->cpts()[8] = B12;
					ehf_beSurf->cpts()[4] = B13;
					ehSymf_beSurf->cpts()[2] = B12;
					ehSymf_beSurf->cpts()[1] = B13;
					break;
				}
				case 1:
				{
					//std::cout << "-----31-----" << std::endl;
					B32 = ehf_beSurf->cpts()[9];
					B22 = ehSymf_beSurf->cpts()[10];
					B33 = ehf_beSurf->cpts()[5];
					B23 = ehSymf_beSurf->cpts()[6];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[6].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[8] = ePW0;
					ehf_beSurf->linkCptPWs()[4] = ePW1;
					ehSymf_beSurf->linkCptPWs()[11] = ePW0;
					ehSymf_beSurf->linkCptPWs()[7] = ePW1;

					ehf_beSurf->cpts()[8] = B12;
					ehf_beSurf->cpts()[4] = B13;
					ehSymf_beSurf->cpts()[11] = B12;
					ehSymf_beSurf->cpts()[7] = B13;
					break;
				}
				case 2:
				{
					//std::cout << "-----32-----" << std::endl;
					B32 = ehf_beSurf->cpts()[9];
					B22 = ehSymf_beSurf->cpts()[9];
					B33 = ehf_beSurf->cpts()[5];
					B23 = ehSymf_beSurf->cpts()[10];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[10].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[8] = ePW0;
					ehf_beSurf->linkCptPWs()[4] = ePW1;
					ehSymf_beSurf->linkCptPWs()[13] = ePW0;
					ehSymf_beSurf->linkCptPWs()[14] = ePW1;

					ehf_beSurf->cpts()[8] = B12;
					ehf_beSurf->cpts()[4] = B13;
					ehSymf_beSurf->cpts()[13] = B12;
					ehSymf_beSurf->cpts()[14] = B13;
					break;
				}
				case 3:
				{
					//std::cout << "-----33-----" << std::endl;
					B32 = ehf_beSurf->cpts()[9];
					B22 = ehSymf_beSurf->cpts()[5];
					B33 = ehf_beSurf->cpts()[5];
					B23 = ehSymf_beSurf->cpts()[9];
					B12 = B32 * 0.5 + B22 * 0.5;
					B13 = B33 * 0.5 + B23 * 0.5;

					/*--------------------------------------------------------------------------------------------------------------------------------------------*/
					/*E1*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW0.linkingCtlPId_weights.find(it->first) == ePW0.linkingCtlPId_weights.end())
						{
							ePW0.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW0.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW0.linkingCtlPId_weights.begin(); it != ePW0.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[4]: " << ehf_beSurf->cpts()[4][0] << " " << ehf_beSurf->cpts()[4][1] << " " << ehf_beSurf->cpts()[4][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					/*E2*/
					/*the first face control point*/
					for (std::map<int, double>::iterator it = ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.begin(); it != ehf_beSurf->linkCptPWs()[5].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*the second face control point*/
					for (std::map<int, double>::iterator it = ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.begin(); it != ehSymf_beSurf->linkCptPWs()[9].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempPW = *it;
						tempPW.second *= 0.5;
						if (ePW1.linkingCtlPId_weights.find(it->first) == ePW1.linkingCtlPId_weights.end())
						{
							ePW1.linkingCtlPId_weights.insert(tempPW);
						}
						else
						{
							std::map<int, double>::iterator eIt = ePW1.linkingCtlPId_weights.find(it->first);
							eIt->second += tempPW.second;
						}
					}
					/*check*/
					/*{
						CPoint tempP;
						for (std::map<int, double>::iterator it = ePW1.linkingCtlPId_weights.begin(); it != ePW1.linkingCtlPId_weights.end(); it++)
						{
							tempP += pMesh->idVertex(it->first)->point() * it->second;
						}
						std::cout << " ehf_beSurf->cpts()[8]: " << ehf_beSurf->cpts()[8][0] << " " << ehf_beSurf->cpts()[8][1] << " " << ehf_beSurf->cpts()[8][2] << std::endl;
						std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
					}*/
					/*--------------------------------------------------------------------------------------------------------------------------------------------*/

					ehf_beSurf->linkCptPWs()[8] = ePW0;
					ehf_beSurf->linkCptPWs()[4] = ePW1;
					ehSymf_beSurf->linkCptPWs()[4] = ePW0;
					ehSymf_beSurf->linkCptPWs()[8] = ePW1;

					ehf_beSurf->cpts()[8] = B12;
					ehf_beSurf->cpts()[4] = B13;
					ehSymf_beSurf->cpts()[4] = B12;
					ehSymf_beSurf->cpts()[8] = B13;
					break;
				}
				default:
					break;
				}
				break;
			}
			default:
				break;
			}
		}

	}
	//std::cout << "Obtaining edges' points succesfully!" << std::endl;

	markCorner_quadMesh(pMesh);
	//markCorner_feature_quadMesh(pMesh);
	for (auto v : It::MVIterator(pMesh))
	{
		/*create bezier surface corner control points*/
		if (v->boundary())//boundary corner
		{
			int vNumFeatrueEdges = 0;
			for (auto ve : It::VCcwEIterator(pMesh, v))
			{
				if (ve->feature())
				{
					vNumFeatrueEdges++;
				}
			}

			int vNumFaces = 0;
			for (auto ve : It::VCcwFIterator(pMesh, v))
			{
				vNumFaces++;
			}

			if (v->ifCorner())
			{
				//std::cout << 1 << std::endl;
				/*equation(19)*/
				CPoint B44;
				CPoint PA;
				H* eh = pMesh->vertexMostCcwOutHalfEdge(v);
				F* ehf = pMesh->halfedgeFace(eh);
				auto ehf_beSurf = this->bfs()[ehf->id() - 1];
				PA = eh->source()->point();
				V* pt1 = pMesh->halfedgeSource(eh);
				//std::cout << "-----0-----" << std::endl;
				B44 = PA;
				BE_linkingCtlPoints vPW0;
				{
					std::pair<int, double> tempPW(pt1->id(), 1.0);
					vPW0.linkingCtlPId_weights.insert(tempPW);
				}
				for (auto eh : It::VCcwOutHEIterator(pMesh, v))
				{
					ehf = pMesh->halfedgeFace(eh);
					ehf_beSurf = this->bfs()[ehf->id() - 1];
					switch (eh->localId())
					{
					case 0://halfedg local id 0
					{
						//std::cout << "-----e0-----" << std::endl;
						ehf_beSurf->cpts()[0] = B44;
						ehf_beSurf->linkCptPWs()[0] = vPW0;
						break;
					}
					case 1:
					{
						//std::cout << "-----e1-----" << std::endl;
						ehf_beSurf->cpts()[3] = B44;
						ehf_beSurf->linkCptPWs()[3] = vPW0;
						break;
					}
					case 2:
					{
						ehf_beSurf->cpts()[15] = B44;
						ehf_beSurf->linkCptPWs()[15] = vPW0;
						break;
					}
					case 3:
					{
						ehf_beSurf->cpts()[12] = B44;
						ehf_beSurf->linkCptPWs()[12] = vPW0;
						break;
					}
					default:
						break;
					}
				}

			}
			else if (vNumFeatrueEdges > 0 && vNumFeatrueEdges != 2)//special process at feature/singular vertex on boundary
			{
				/*find an out halfedge that is feature and its  pre halfedge is not feature*/
				H* featureH = NULL;
				for (auto vh : It::VCcwOutHEIterator(pMesh, v))
				{
					if (pMesh->halfedgeEdge(vh)->feature() && !(pMesh->halfedgeEdge(pMesh->halfedgePrev(vh))->feature()))
					{
						featureH = vh;
						break;
					}
				}
				/*find the most in halfedges that is feature*/
				H* featureHPre = pMesh->halfedgePrev(featureH);
				while (!(pMesh->halfedgeEdge(featureHPre)->feature()))
				{
					featureHPre = pMesh->halfedgePrev(pMesh->halfedgeSym(featureHPre));
				}
				
				CPoint sumPoints(0, 0, 0);

				std::vector<BE_linkingCtlPoints> v_cws;
				{
					H* vh = featureH;
					F* vhf = pMesh->halfedgeFace(vh);
					auto vhf_beSurf = this->bfs()[vhf->id() - 1];
					switch (vh->localId())
					{
					case 0:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[1];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[1]);
						break;
					}
					case 1:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[7];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[7]);
						break;
					}
					case 2:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[14];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[14]);
						break;
					}
					case 3:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[8];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[8]);
						break;
					}
					default:
						break;
					}
				}

				{
					H* vh = featureHPre;
					F* vhf = pMesh->halfedgeFace(vh);
					auto vhf_beSurf = this->bfs()[vhf->id() - 1];
					switch (vh->localId())
					{
					case 0:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[2];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[2]);
						break;
					}
					case 1:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[11];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[11]);
						break;
					}
					case 2:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[13];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[13]);
						break;
					}
					case 3:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[4];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[4]);
						break;
					}
					default:
						break;
					}
				}

				CPoint B11 = sumPoints / 2.0;

				BE_linkingCtlPoints v_tempCW;
				for (int i = 0; i < v_cws.size(); i++)
				{
					for (std::map<int, double>::iterator it = v_cws[i].linkingCtlPId_weights.begin(); it != v_cws[i].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempV = *it;
						tempV.second *= (1.0 / 2.0);
						if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
						{
							v_tempCW.linkingCtlPId_weights.insert(tempV);
						}
						else
						{
							std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
							vit->second += tempV.second;
						}
					}
				}

				for (auto vh : It::VCcwOutHEIterator(pMesh, v))
				{
					F* vhf = pMesh->halfedgeFace(vh);
					auto vhf_beSurf = this->bfs()[vhf->id() - 1];
					switch (vh->localId())
					{
					case 0:
					{
						vhf_beSurf->cpts()[0] = B11;
						vhf_beSurf->linkCptPWs()[0] = v_tempCW;
						break;
					}
					case 1:
					{
						vhf_beSurf->cpts()[3] = B11;
						vhf_beSurf->linkCptPWs()[3] = v_tempCW;
						break;
					}
					case 2:
					{
						vhf_beSurf->cpts()[15] = B11;
						vhf_beSurf->linkCptPWs()[15] = v_tempCW;
						break;
					}
					case 3:
					{
						vhf_beSurf->cpts()[12] = B11;
						vhf_beSurf->linkCptPWs()[12] = v_tempCW;
						break;
					}
					default:
						break;
					}
				}
			}
			else if (vNumFaces != 2)
			{
			/*find an out halfedge that is feature (constrained boundary edge) and its  pre halfedge is not feature (constrained boundary edge)*/
			H* featureH = NULL;
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				if ((pMesh->halfedgeEdge(vh)->feature()||(pMesh->halfedgeEdge(vh)->constrinedBoundary())) && !(pMesh->halfedgeEdge(pMesh->halfedgePrev(vh))->feature()||(pMesh->halfedgeEdge(pMesh->halfedgePrev(vh))->constrinedBoundary())))
				{
					featureH = vh;
					break;
				}
			}
			/*find the most in halfedges that is feature*/
			H* featureHPre = pMesh->halfedgePrev(featureH);
			while (!(pMesh->halfedgeEdge(featureHPre)->feature()|| pMesh->halfedgeEdge(featureHPre)->constrinedBoundary()))
			{
				featureHPre = pMesh->halfedgePrev(pMesh->halfedgeSym(featureHPre));
			}

			CPoint sumPoints(0, 0, 0);

			std::vector<BE_linkingCtlPoints> v_cws;
			{
				H* vh = featureH;
				F* vhf = pMesh->halfedgeFace(vh);
				auto vhf_beSurf = this->bfs()[vhf->id() - 1];
				switch (vh->localId())
				{
				case 0:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[1];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[1]);
					break;
				}
				case 1:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[7];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[7]);
					break;
				}
				case 2:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[14];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[14]);
					break;
				}
				case 3:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[8];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[8]);
					break;
				}
				default:
					break;
				}
			}

			{
				H* vh = featureHPre;
				F* vhf = pMesh->halfedgeFace(vh);
				auto vhf_beSurf = this->bfs()[vhf->id() - 1];
				switch (vh->localId())
				{
				case 0:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[2];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[2]);
					break;
				}
				case 1:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[11];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[11]);
					break;
				}
				case 2:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[13];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[13]);
					break;
				}
				case 3:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[4];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[4]);
					break;
				}
				default:
					break;
				}
			}

			CPoint B11 = sumPoints / 2.0;

			BE_linkingCtlPoints v_tempCW;
			for (int i = 0; i < v_cws.size(); i++)
			{
				for (std::map<int, double>::iterator it = v_cws[i].linkingCtlPId_weights.begin(); it != v_cws[i].linkingCtlPId_weights.end(); it++)
				{
					std::pair<int, double> tempV = *it;
					tempV.second *= (1.0 / 2.0);
					if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
					{
						v_tempCW.linkingCtlPId_weights.insert(tempV);
					}
					else
					{
						std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
						vit->second += tempV.second;
					}
				}
			}

			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				F* vhf = pMesh->halfedgeFace(vh);
				auto vhf_beSurf = this->bfs()[vhf->id() - 1];
				switch (vh->localId())
				{
				case 0:
				{
					vhf_beSurf->cpts()[0] = B11;
					vhf_beSurf->linkCptPWs()[0] = v_tempCW;
					break;
				}
				case 1:
				{
					vhf_beSurf->cpts()[3] = B11;
					vhf_beSurf->linkCptPWs()[3] = v_tempCW;
					break;
				}
				case 2:
				{
					vhf_beSurf->cpts()[15] = B11;
					vhf_beSurf->linkCptPWs()[15] = v_tempCW;
					break;
				}
				case 3:
				{
					vhf_beSurf->cpts()[12] = B11;
					vhf_beSurf->linkCptPWs()[12] = v_tempCW;
					break;
				}
				default:
					break;
				}
			}
			}
			else
			{
				/*equation(18)*/
				H* eh = pMesh->vertexMostClwInHalfEdge(v);
				H* eh1 = pMesh->vertexMostClwOutHalfEdge(v);
				for (auto vh : It::VClwInHEIterator(pMesh, v))
				{
					if (vh->edge()->boundary()) {
						eh = vh;
					}
					else {
						eh1 = pMesh->halfedgeNext(vh);
					}
				}

				F* ehf = pMesh->halfedgeFace(eh);
				F* eh1f = pMesh->halfedgeFace(eh1);
				auto ehf_beSurf = this->bfs()[ehf->id() - 1];
				auto eh1f_beSurf = this->bfs()[eh1f->id() - 1];
				CPoint B41, B11, B31, B21;
				switch (eh->localId())
				{
				case 0:
				{
					//break;
					//std::cout << "-----V0-----" << std::endl;
					switch (eh1->localId())
					{
					case 0:
					{
						B31 = ehf_beSurf->cpts()[2];
						B21 = eh1f_beSurf->cpts()[1];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[2].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[1].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[3] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[0] = v_tempCW;

						ehf_beSurf->cpts()[3] = B41;
						eh1f_beSurf->cpts()[0] = B41;
						break;
					}
					case 1:
					{
						//std::cout << "-----V01-----" << std::endl;
						B31 = ehf_beSurf->cpts()[2];
						B21 = eh1f_beSurf->cpts()[7];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[2].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[7].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[3] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[3] = v_tempCW;

						ehf_beSurf->cpts()[3] = B41;
						eh1f_beSurf->cpts()[3] = B41;
						break;
					}
					case 2:
					{
						B31 = ehf_beSurf->cpts()[2];
						B21 = eh1f_beSurf->cpts()[14];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[2].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[14].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[3] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[15] = v_tempCW;

						ehf_beSurf->cpts()[3] = B41;
						eh1f_beSurf->cpts()[15] = B41;
						break;
					}
					case 3:
					{
						B31 = ehf_beSurf->cpts()[2];
						B21 = eh1f_beSurf->cpts()[8];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[2].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[8].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[3] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[12] = v_tempCW;

						ehf_beSurf->cpts()[3] = B41;
						eh1f_beSurf->cpts()[12] = B41;
						break;
					}
					default:
						break;
					}
					break;
				}
				case 1:
				{
					//break;
					//std::cout << "-----V1-----" << std::endl;
					switch (eh1->localId())
					{
					case 0:
					{
						B31 = ehf_beSurf->cpts()[11];
						B21 = eh1f_beSurf->cpts()[1];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[11].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[1].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[15] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[0] = v_tempCW;

						ehf_beSurf->cpts()[15] = B41;
						eh1f_beSurf->cpts()[0] = B41;
						break;
					}
					case 1:
					{
						B31 = ehf_beSurf->cpts()[11];
						B21 = eh1f_beSurf->cpts()[7];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[11].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[7].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[15] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[3] = v_tempCW;

						ehf_beSurf->cpts()[15] = B41;
						eh1f_beSurf->cpts()[3] = B41;
						break;
					}
					case 2:
					{
						B31 = ehf_beSurf->cpts()[11];
						B21 = eh1f_beSurf->cpts()[14];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[11].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[14].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[15] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[15] = v_tempCW;

						ehf_beSurf->cpts()[15] = B41;
						eh1f_beSurf->cpts()[15] = B41;
						break;
					}
					case 3:
					{
						B31 = ehf_beSurf->cpts()[11];
						B21 = eh1f_beSurf->cpts()[8];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[11].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[8].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[15] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[12] = v_tempCW;

						ehf_beSurf->cpts()[15] = B41;
						eh1f_beSurf->cpts()[12] = B41;
						break;
					}
					default:
						break;
					}
					break;
				}
				case 2:
				{
					//break;
					//std::cout << "-----V2-----" << std::endl;
					switch (eh1->localId())
					{
					case 0:
					{
						B31 = ehf_beSurf->cpts()[13];
						B21 = eh1f_beSurf->cpts()[1];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[13].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[1].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[12] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[0] = v_tempCW;

						ehf_beSurf->cpts()[12] = B41;
						eh1f_beSurf->cpts()[0] = B41;
						break;
					}
					case 1:
					{
						B31 = ehf_beSurf->cpts()[13];
						B21 = eh1f_beSurf->cpts()[7];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[13].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[7].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[12] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[3] = v_tempCW;

						ehf_beSurf->cpts()[12] = B41;
						eh1f_beSurf->cpts()[3] = B41;
						break;
					}
					case 2:
					{
						B31 = ehf_beSurf->cpts()[13];
						B21 = eh1f_beSurf->cpts()[14];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[13].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[14].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[12] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[15] = v_tempCW;

						ehf_beSurf->cpts()[12] = B41;
						eh1f_beSurf->cpts()[15] = B41;
						break;
					}
					case 3:
					{
						B31 = ehf_beSurf->cpts()[13];
						B21 = eh1f_beSurf->cpts()[8];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[13].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[8].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[12] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[12] = v_tempCW;

						ehf_beSurf->cpts()[12] = B41;
						eh1f_beSurf->cpts()[12] = B41;
						break;
					}
					default:
						break;
					}
					break;
				}
				case 3:
				{
					//break;
					//std::cout << "-----V3-----" << std::endl;
					switch (eh1->localId())
					{
					case 0:
					{
						//std::cout << "-----v30-----" << std::endl;
						B31 = ehf_beSurf->cpts()[4];
						B21 = eh1f_beSurf->cpts()[1];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[4].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[1].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[0] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[0] = v_tempCW;

						ehf_beSurf->cpts()[0] = B41;
						eh1f_beSurf->cpts()[0] = B41;
						break;
					}
					case 1:
					{
						B31 = ehf_beSurf->cpts()[4];
						B21 = eh1f_beSurf->cpts()[7];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[4].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[7].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[0] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[3] = v_tempCW;

						ehf_beSurf->cpts()[0] = B41;
						eh1f_beSurf->cpts()[3] = B41;
						break;
					}
					case 2:
					{
						B31 = ehf_beSurf->cpts()[4];
						B21 = eh1f_beSurf->cpts()[14];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[4].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[14].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[0] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[15] = v_tempCW;

						ehf_beSurf->cpts()[0] = B41;
						eh1f_beSurf->cpts()[15] = B41;
						break;
					}
					case 3:
					{
						B31 = ehf_beSurf->cpts()[4];
						B21 = eh1f_beSurf->cpts()[8];
						B41 = B31 * 0.5 + B21 * 0.5;

						/*--------------------------------------------------------------------------------------------------------------*/
						BE_linkingCtlPoints v_tempCW;
						/*f1*/
						for (auto cwit : ehf_beSurf->linkCptPWs()[4].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*f2*/
						for (auto cwit : eh1f_beSurf->linkCptPWs()[8].linkingCtlPId_weights)
						{
							std::pair<int, double> tempV = cwit;
							tempV.second *= 0.5;
							if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
							{
								v_tempCW.linkingCtlPId_weights.insert(tempV);
							}
							else
							{
								std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
								vit->second += tempV.second;
							}
						}
						/*check*/
						/*{
							CPoint tempP;
							for (auto it : v_tempCW.linkingCtlPId_weights)
							{
								tempP += (pMesh->idVertex(it.first)->point() * it.second);
							}
							std::cout << " vRegP: " << vRegP[0] << " " << vRegP[1] << " " << vRegP[2] << std::endl;
							std::cout << " tempP: " << tempP[0] << " " << tempP[1] << " " << tempP[2] << std::endl;
						}*/
						/*--------------------------------------------------------------------------------------------------------------*/
						ehf_beSurf->linkCptPWs()[0] = v_tempCW;
						eh1f_beSurf->linkCptPWs()[12] = v_tempCW;

						ehf_beSurf->cpts()[0] = B41;
						eh1f_beSurf->cpts()[12] = B41;
						break;
					}
					default:
						break;
					}
					break;
				}
				default:
					break;
				}
			}
		}
		else if (v->feature())//feature corner contrl points
		{
			/*statistic the number of feature edges incident to inner vertex */
			int vNumFeatrueEdges = 0;
			int vNumEdges = 0;
			for (auto ve : It::VCcwEIterator(pMesh, v))
			{
				vNumEdges++;
				if (ve->feature())
				{
					vNumFeatrueEdges++;
				}
			}
			if (vNumFeatrueEdges == vNumEdges)//singular feature point as the corner control point
			{
				//std::cout << 1 << std::endl;
				/*equation(19)*/
				CPoint B44;
				CPoint PA;
				H* eh = pMesh->vertexMostCcwOutHalfEdge(v);
				F* ehf = pMesh->halfedgeFace(eh);
				auto ehf_beSurf = this->bfs()[ehf->id() - 1];
				PA = eh->source()->point();
				V* pt1 = pMesh->halfedgeSource(eh);
				//std::cout << "-----0-----" << std::endl;
				B44 = PA;
				BE_linkingCtlPoints vPW0;
				{
					std::pair<int, double> tempPW(pt1->id(), 1.0);
					vPW0.linkingCtlPId_weights.insert(tempPW);
				}
				for (auto eh : It::VCcwOutHEIterator(pMesh, v))
				{
					ehf = pMesh->halfedgeFace(eh);
					ehf_beSurf = this->bfs()[ehf->id() - 1];
					switch (eh->localId())
					{
					case 0://halfedg local id 0
					{
						//std::cout << "-----e0-----" << std::endl;
						ehf_beSurf->cpts()[0] = B44;
						ehf_beSurf->linkCptPWs()[0] = vPW0;
						break;
					}
					case 1:
					{
						//std::cout << "-----e1-----" << std::endl;
						ehf_beSurf->cpts()[3] = B44;
						ehf_beSurf->linkCptPWs()[3] = vPW0;
						break;
					}
					case 2:
					{
						ehf_beSurf->cpts()[15] = B44;
						ehf_beSurf->linkCptPWs()[15] = vPW0;
						break;
					}
					case 3:
					{
						ehf_beSurf->cpts()[12] = B44;
						ehf_beSurf->linkCptPWs()[12] = vPW0;
						break;
					}
					default:
						break;
					}
				}
			}
			else
			{
				/*equation (15) update leveraging to the feature vertex*/
				//H* vh = pMesh->vertexMostClwInHalfEdge(v);
				std::vector<H*> vhs;
				/*find an out halfedge that is feature and its  pre halfedge is not feature*/
				H* featureH = NULL;
				for (auto vh : It::VCcwOutHEIterator(pMesh, v))
				{
					if (pMesh->halfedgeEdge(vh)->feature() && !(pMesh->halfedgeEdge(pMesh->halfedgePrev(vh))->feature()))
					{
						featureH = vh;
						break;
					}
				}
				vhs.push_back(featureH);
				H* featureHPre = pMesh->halfedgePrev(featureH);
				while (!(pMesh->halfedgeEdge(featureHPre)->feature()))
				{
					featureH = pMesh->halfedgeSym(featureHPre);
					featureHPre = pMesh->halfedgePrev(featureH);
				}
				featureH = pMesh->halfedgeSym(featureHPre);
				vhs.push_back(featureH);
				
				CPoint sumPoints(0, 0, 0);
				
				std::vector<BE_linkingCtlPoints> v_cws;
				for (int i = 0; i < vhs.size(); i++)
				{
					H* vh = vhs[i];
					F* vhf = pMesh->halfedgeFace(vh);
					auto vhf_beSurf = this->bfs()[vhf->id() - 1];
					switch (vh->localId())
					{
					case 0:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[1];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[1]);
						break;
					}
					case 1:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[7];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[7]);
						break;
					}
					case 2:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[14];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[14]);
						break;
					}
					case 3:
					{
						sumPoints = sumPoints + vhf_beSurf->cpts()[8];
						v_cws.push_back(vhf_beSurf->linkCptPWs()[8]);
						break;
					}
					default:
						break;
					}
				}
				CPoint B11 = sumPoints / vhs.size();
				
				BE_linkingCtlPoints v_tempCW;
				for (int i = 0; i < v_cws.size(); i++)
				{
					for (std::map<int, double>::iterator it = v_cws[i].linkingCtlPId_weights.begin(); it != v_cws[i].linkingCtlPId_weights.end(); it++)
					{
						std::pair<int, double> tempV = *it;
						tempV.second *= (1.0 / vhs.size());
						if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
						{
							v_tempCW.linkingCtlPId_weights.insert(tempV);
						}
						else
						{
							std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
							vit->second += tempV.second;
						}
					}
				}

				for (auto vh : It::VCcwOutHEIterator(pMesh, v))
				{
					F* vhf = pMesh->halfedgeFace(vh);
					auto vhf_beSurf = this->bfs()[vhf->id() - 1];
					switch (vh->localId())
					{
					case 0:
					{
						vhf_beSurf->cpts()[0] = B11;
						vhf_beSurf->linkCptPWs()[0] = v_tempCW;
						break;
					}
					case 1:
					{
						vhf_beSurf->cpts()[3] = B11;
						vhf_beSurf->linkCptPWs()[3] = v_tempCW;
						break;
					}
					case 2:
					{
						vhf_beSurf->cpts()[15] = B11;
						vhf_beSurf->linkCptPWs()[15] = v_tempCW;
						break;
					}
					case 3:
					{
						vhf_beSurf->cpts()[12] = B11;
						vhf_beSurf->linkCptPWs()[12] = v_tempCW;
						break;
					}
					default:
						break;
					}
				}
			}
		}
		else
		{
			//std::cout << "-------------regular vertex------------------" << std::endl;
			/*equation(15)*/
			//H* vh = pMesh->vertexMostClwInHalfEdge(v);
			std::vector<H*> vhs;
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				vhs.push_back(vh);
			}
			
			CPoint sumPoints(0, 0, 0);
			std::vector<BE_linkingCtlPoints> v_cws;
			for (int i = 0; i < vhs.size(); i++)
			{
				H* vh = vhs[i];
				F* vhf = pMesh->halfedgeFace(vh);
				auto vhf_beSurf = this->bfs()[vhf->id() - 1];
				switch (vh->localId())
				{
				case 0:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[5];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[5]);
					break;
				}
				case 1:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[6];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[6]);
					break;
				}
				case 2:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[10];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[10]);
					break;
				}
				case 3:
				{
					sumPoints = sumPoints + vhf_beSurf->cpts()[9];
					v_cws.push_back(vhf_beSurf->linkCptPWs()[9]);
					break;
				}
				default:
					break;
				}
			}
			//std::cout << "vhs.size():" << vhs.size() << std::endl;
			CPoint B11 = sumPoints / vhs.size();
			
			BE_linkingCtlPoints v_tempCW;
			for (int i = 0; i < v_cws.size(); i++)
			{
				for (std::map<int, double>::iterator it = v_cws[i].linkingCtlPId_weights.begin(); it != v_cws[i].linkingCtlPId_weights.end(); it++)
				{
					std::pair<int, double> tempV = *it;
					tempV.second *= (1.0 / vhs.size());
					if (v_tempCW.linkingCtlPId_weights.find(tempV.first) == v_tempCW.linkingCtlPId_weights.end())
					{
						v_tempCW.linkingCtlPId_weights.insert(tempV);
					}
					else
					{
						std::map<int, double>::iterator vit = v_tempCW.linkingCtlPId_weights.find(tempV.first);
						vit->second += tempV.second;
					}
				}
			}

			for (int i = 0; i < vhs.size(); i++)
			{
				H* vh = vhs[i];
				F* vhf = pMesh->halfedgeFace(vh);
				auto vhf_beSurf = this->bfs()[vhf->id() - 1];
				switch (vh->localId())
				{
				case 0:
				{
					vhf_beSurf->cpts()[0] = B11;
					vhf_beSurf->linkCptPWs()[0] = v_tempCW;
					break;
				}
				case 1:
				{
					vhf_beSurf->cpts()[3] = B11;
					vhf_beSurf->linkCptPWs()[3] = v_tempCW;
					break;
				}
				case 2:
				{
					vhf_beSurf->cpts()[15] = B11;
					vhf_beSurf->linkCptPWs()[15] = v_tempCW;
					break;
				}
				case 3:
				{
					vhf_beSurf->cpts()[12] = B11;
					vhf_beSurf->linkCptPWs()[12] = v_tempCW;
					break;
				}
				default:
					break;
				}
			}
		}
	}

	//std::cout << "Obtaining corners' points succesfully!" << std::endl;
	//std::cout << "----------------------------------------------" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::sortBezierSurfaces()
{
	/*sort according to ids*/
	std::sort(bfs().begin(), bfs().end(), bfIdCompare);
	/*check sorting result*/
	/*std::cout << "The consequence of bezier surfaces: " << std::endl;
	for (int i = 0; i < bfs().size(); i++)
	{
		std::cout << bfs()[i]->id() << std::endl;
	}*/
}

void CCG_QMSLib::CCG_QMS_model::raiseBezierSurfDegree3To5_all(M* pMesh)
{
	for (int i = 0; i < bfs().size(); i++)
	{
		bfs()[i]->raiseDegree3to5(pMesh);
		//bfs()[i]->raiseDegree(pMesh);
		bfs()[i]->degree(0) = 5;
		bfs()[i]->degree(1) = 5;
		//std::cout << bfs()[i]->id() << " Raising degree 3 to 5 finished!" << std::endl;
	}
}

void CCG_QMSLib::CCG_QMS_model::sortBezierCtlPts_cw(M* pMesh)
{
	/*Re-number the halfedges' local Id on the mesh faces,
	so that the half-edge with ID = 1 can point to both the first point
	and the first Bezier control point on the surface, 
	making it convenient for later acquisition of the B-spline surface.*/
	for (auto f : It::MFIterator(pMesh))
	{
		bool fB = false;
		for (auto fv : It::FVIterator(pMesh, f))
		{
			if (fv->ifSingular())
			{
				fB = true;
				break;
			}
		}
		if (fB)
		{
			/*Note: halfedge on the quad mesh face points to the fourth vertex*/
			//std::cout << "f->id(): " << f->id() << std::endl;
			H* fh = pMesh->faceHalfedge(f);
			//std::cout << "change before he id: " << fh->localId() << std::endl;
			fh->localId() = 4;
			//std::cout << "change after he id: " << fh->localId() << std::endl;
			H* fhn = pMesh->halfedgeNext(fh);
			//std::cout << "change before he id: " << fhn->localId() << std::endl;
			fhn->localId() = 1;
			//std::cout << "change after he id: " << fhn->localId() << std::endl;
			H* fhnn = pMesh->halfedgeNext(fhn);
			//std::cout << "change before he id: " << fhnn->localId() << std::endl;
			fhnn->localId() = 2;
			//std::cout << "change after he id: " << fhnn->localId() << std::endl;
			H* fhnnn = pMesh->halfedgeNext(fhnn);
			//std::cout << "change before he id: " << fhnnn->localId() << std::endl;
			fhnnn->localId() = 3;
			//std::cout << "change after he id: " << fhnnn->localId() << std::endl;
		}
		else
		{
			/*Note: halfedge on the quad mesh face points to the fourth vertex*/
			//std::cout << "f->id(): " << f->id() << std::endl;
			H* fh = pMesh->faceHalfedge(f);
			//std::cout << "change before he id: " << fh->localId();
			fh->localId() = 4;
			//std::cout << "change after he id: " << fh->localId();
			H* fhn = pMesh->halfedgeNext(fh);
			//std::cout << "change before he id: " << fhn->localId();
			fhn->localId() = 1;
			//std::cout << "change after he id: " << fhn->localId();
			H* fhnn = pMesh->halfedgeNext(fhn);
			//std::cout << "change before he id: " << fhnn->localId();
			fhnn->localId() = 2;
			//std::cout << "change after he id: " << fhnn->localId();
			H* fhnnn = pMesh->halfedgeNext(fhnn);
			//std::cout << "change before he id: " << fhnnn->localId();
			fhnnn->localId() = 3;
			//std::cout << "change after he id: " << fhnnn->localId();
		}
	}
}

void CCG_QMSLib::CCG_QMS_model::G1Constraints_optimization_check(M* pMesh)
{
	using namespace Eigen;

	auto set_g_vaules = [](std::vector<double>& g_x_vec, std::vector<double>& g_y_vec, std::vector<double>& g_z_vec, double x, double y, double z)
	{
		g_x_vec.push_back(x);
		g_y_vec.push_back(y);
		g_z_vec.push_back(z);
	};
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->ifSingular())
		{
			/*Reorder the control points of the singular Bezier surface*/
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				F* vh_f = pMesh->halfedgeFace(vh);
				auto vh_f_bf = this->bfs()[vh_f->id() - 1];
				vh_f_bf->cpts_singular().resize((vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1));
				switch (vh->localId())
				{
				case 2:
				{
					/*控制点的顺序不变*/
					for (int i = 0; i < vh_f_bf->cpts().size(); i++)
					{
						vh_f_bf->cpts_singular()[i] = vh_f_bf->cpts()[i];
						vh_f_bf->singular_linkCptPWs()[i] = vh_f_bf->linkCptPWs()[i];
					}
					break;
				}
				case 3:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							vh_f_bf->cpts_singular()[i * ((vh_f_bf->degree(0) + 1)) + j] = vh_f_bf->cpts()[vh_f_bf->degree(0) - i + j * ((vh_f_bf->degree(1) + 1))];
							/*check*/
							//std::cout << "before resort cpt index: " << vh_f_bf->degree(0) - i + j * ((vh_f_bf->degree(1) + 1)) << std::endl;
							//std::cout << "After resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							vh_f_bf->singular_linkCptPWs()[i * ((vh_f_bf->degree(0) + 1)) + j] = vh_f_bf->linkCptPWs()[vh_f_bf->degree(0) - i + j * ((vh_f_bf->degree(1) + 1))];
						}
					}
					break;
				}
				case 4:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->cpts()[(vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i];
							/*check*/
							//std::cout << "before resort cpt index: " << (vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i << std::endl;
							//std::cout << "After resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->linkCptPWs()[(vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i];
						}
					}

					break;
				}
				case 1:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->cpts()[vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j];
							/*check*/
							//std::cout << "before resort cpt index: " << vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j << std::endl;
							//std::cout << "After resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->linkCptPWs()[vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j];
						}
					}
					break;
				}
				default:
					break;
				}
			}

			std::vector<double> g_x_vec, g_y_vec, g_z_vec;
			/*Find the Bezier surface corresponding to the singular face that is associated with the current singularity point*/
			/*Locate the adjacent half side around the point*/
			std::vector<H*> vhs;
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				vhs.push_back(vh);
			}

			/*
			* Formulate a system of equations, 
			* optimize the control points on the Bezier surface
			* corresponding to the singularities around the singular points
			* and make them achieve G1 continuity
			*/
			/*Obtain all control points*/
			/*Assign an order to the control points*/
			/*
			* The i-th halfedge j-th control point's index = i*30+j，i=0,...,v->degree(),j=0,...,29
			* (No code display is required for the implementation. This is the default rule)
			*/
			std::vector<CPoint> v_all_cpts;
			for (int i = 0; i < vhs.size(); i++)
			{
				H* vh = vhs[i];
				F* vh_f = pMesh->halfedgeFace(vh);
				auto vh_f_bf = this->bfs()[vh_f->id() - 1];
				for (int j = 0; j < vh_f_bf->cpts_singular().size(); j++)
				{
					if (j % (vh_f_bf->degree(0) + 1) == 0)
					{
						//std::cout << std::endl;
						continue;
					}
					//std::cout << " " << j;
					v_all_cpts.push_back(vh_f_bf->cpts_singular()[j]);

				}
				if (i == vhs.size() - 1)
				{
					/*
					* Place the control points corresponding to the singularity into the current control point array.
					* Remember that this can only be done once, 
					* and they should be placed at the end of the current control point array. */
					v_all_cpts.push_back(vh_f_bf->cpts_singular()[0]);
				}
				//std::cout << "----------------------" << std::endl;
			}

			/*check*/
			//std::cout << "v_all_cpts.size(): " << v_all_cpts.size() << std::endl;

			/*Formulation*/
			std::vector<Eigen::Triplet<double>> tripletlists; // Define a container for storing all the coefficient pairs
			int total_equation_count = 0;//It indicates which equation it is
			for (int i = 0; i < vhs.size(); i++)
			{
				std::vector<Eigen::Triplet<double>> tripletlist; // Define a container for storing all the coefficient pairs
				H* vh = vhs[i];
				F* vh_f = pMesh->halfedgeFace(vh);
				auto vh_f_bf = this->bfs()[vh_f->id() - 1];
				int bf_id = i;//Corresponding to the id of face e in the paper's figure
				int pre_link_bf_id = ((i - 1) + vhs.size()) % vhs.size();//Corresponding to the id of face (e - 1) in the paper's figure
				int next_link_bf_id = ((i + 1) + vhs.size()) % vhs.size();//Corresponding to the id of face (e + 1) in the paper's figure
				//Formulation w1 and w2
				double w1 = 0.0, w2 = 0.0;
				if (v->boundary())
				{
					w1 = cos(1.0 * _M_PI / v->degree());
				}
				else
				{
					w1 = cos(2.0 * _M_PI / v->degree());
				}
				V* vhTarget = pMesh->halfedgeTarget(vh);
				if (vhTarget->boundary())
				{
					w2 = cos(1.0 * _M_PI / vhTarget->degree());
				}
				else
				{
					w2 = cos(2.0 * _M_PI / vhTarget->degree());
				}
				//std::cout << "w1: " << w1 << " w2: " << w2 << std::endl;
				int temp_index = -1;
				/*--------------(equation 27) start------------------------------------------*/
				//c_21^e-1 coefficient
				temp_index = pre_link_bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0);

				//c_21^e coefficient
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0 * w1);

				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w1);

				//c_12^e coefficient
				temp_index = next_link_bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				/*--------------(eqution 27) end------------------------------------------*/


				/*--------------(equation 28) start------------------------------------------*/
				//c_22^e-1 coefficient
				temp_index = pre_link_bf_id * 30 + 5;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_21^e coefficient
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w1 - 10.0);

				//c_31^e coefficient
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -8.0 * w1);

				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -2.0 * w1);

				//c_22^e coefficient
				temp_index = bf_id * 30 + 5;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				/*--------------(equation 28) end------------------------------------------*/


				/*--------------(equation 29) start------------------------------------------*/
				//c_23^e-1 coefficient
				temp_index = pre_link_bf_id * 30 + 10;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_31^e coefficient
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0);

				//c_51^e coefficient
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0 * w1);

				//c_41^e coefficient
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, 4.0 * w1);

				//c_61^e coefficient
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, w1);

				//c_21^e coefficient
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, w2);

				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -w2);

				//c_32^e coefficient
				temp_index = bf_id * 30 + 6;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				/*--------------(equation 29)end------------------------------------------*/


				/*--------------(equation 30)start------------------------------------------*/
				//c_24^e-1 coefficient
				temp_index = pre_link_bf_id * 30 + 15;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_41^e coefficient
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0);

				//c_61^e coefficient
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, -w1);

				//c_51^e coefficient
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, w1);

				//c_31^e coefficient
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, 4.0 * w2);

				//c_21^e coefficient
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0 * w2);

				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, w2);

				//c_42^e coefficient
				temp_index = bf_id * 30 + 7;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				/*--------------(equation 30) end------------------------------------------*/


				/*--------------(equation 31) start------------------------------------------*/
				//c_25^e-1 coefficient
				temp_index = pre_link_bf_id * 30 + 20;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_51^e coefficient
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w2 - 10.0);

				//c_41^e coefficient
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, -8.0 * w2);

				//c_61^e coefficient
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, -2.0 * w2);

				//c_52^e coefficient
				temp_index = bf_id * 30 + 8;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				/*--------------(equation 31) end------------------------------------------*/


				/*--------------(equation 32) start------------------------------------------*/
				//c_26^e-1 coefficient
				temp_index = pre_link_bf_id * 30 + 25;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				;
				//c_61^e coefficient
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w2 - 10.0);

				//c_51^e coefficient
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0 * w2);

				//c_62^e coefficient
				temp_index = bf_id * 30 + 9;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				/*--------------(equation 32) end------------------------------------------*/


				/*--------------(coefficient 33) start------------------------------------------*/
				//c_11^e coefficient
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -1.0);

				//c_21^e coefficient
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_31^e coefficient
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0);

				//c_41^e coefficient
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0);

				//c_51^e coefficient
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0);

				//c_61^e coefficient
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				total_equation_count++;
				//std::cout << "total_equation_count 33:" << total_equation_count << std::endl;
				////std::cout << " g_x_vec.size():" << g_x_vec.size() << std::endl;
				///*--------------(equation 33) end------------------------------------------*/

				///*--------------(equation 34-35) start------------------------------------------*/
				//c_51^e coefficient
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_52^e coefficient
				temp_index = bf_id * 30 + 8;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_53^e coefficient
				temp_index = bf_id * 30 + 13;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_54^e coefficient
				temp_index = bf_id * 30 + 18;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_55^e coefficient
				temp_index = bf_id * 30 + 23;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_56^e coefficient
				temp_index = bf_id * 30 + 28;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_61^e coefficient
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_62^e coefficient
				temp_index = bf_id * 30 + 9;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_63^e coefficient
				temp_index = bf_id * 30 + 14;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_64^e coefficient
				temp_index = bf_id * 30 + 19;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_65^e coefficient
				temp_index = bf_id * 30 + 24;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_66^e coefficient
				temp_index = bf_id * 30 + 29;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_26^e coefficient
				temp_index = bf_id * 30 + 25;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_36^e coefficient
				temp_index = bf_id * 30 + 26;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_46^e coefficient
				temp_index = bf_id * 30 + 27;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_25^e coefficient
				temp_index = bf_id * 30 + 20;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_35^e coefficient
				temp_index = bf_id * 30 + 21;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;

				//c_45^e coefficient
				temp_index = bf_id * 30 + 22;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				total_equation_count++;
				//std::cout << "total_equation_count:" << total_equation_count << std::endl;
				/*--------------(equation 34-35) end------------------------------------------*/

				//std::cout << " g_x_vec.size():" << g_x_vec.size() << std::endl;
				/*--------------(equation 36) start------------------------------------------*/
				//for (int j = 0; j < (vh_f_bf->degree(0) + 1) * vh_f_bf->degree(1); j++)
				//{
				//	if ((j + 1) % vh_f_bf->degree(0) == 0)
				//	{
				//		//c_2j^e coefficient
				//		temp_index = bf_id * 30 + j - 4;
				//		tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//		//c_1,j^e(c_i,1^e+1) coefficient
				//		if ((j + 1) / vh_f_bf->degree(0) == 1)
				//		{
				//			temp_index = v_all_cpts.size() - 1;
				//			tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//		}
				//		else
				//		{
				//			temp_index = next_link_bf_id * 30 + (j + 1) / vh_f_bf->degree(0) - 1;
				//			tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//		}
				//		set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0] - v_all_cpts[bf_id * 30 + j - 4][0], v_all_cpts[temp_index][1] - v_all_cpts[bf_id * 30 + j - 4][1], v_all_cpts[temp_index][2] - v_all_cpts[bf_id * 30 + j - 4][2]);
				//		total_equation_count++;
				//		//std::cout << "total_equation_count:" << total_equation_count << std::endl;
				//	}
				//	else
				//	{
				//		//c_ij^e coefficient
				//		temp_index = bf_id * 30 + j;
				//		tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//		//c_i+1,j^e coefficient
				//		temp_index = bf_id * 30 + j + 1;
				//		tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//		set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[bf_id * 30 + j][0] - v_all_cpts[temp_index][0], v_all_cpts[bf_id * 30 + j][1] - v_all_cpts[temp_index][1], v_all_cpts[bf_id * 30 + j][2] - v_all_cpts[temp_index][2]);
				//		total_equation_count++;
				//	}
				//}
				////std::cout << "total_equation_count 36:" << total_equation_count << std::endl;
				////c_11^e coefficient
				////c_21^e coefficient
				///*temp_index = v_all_cpts.size() - 1;
				//tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//temp_index = bf_id * 30;
				//tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[v_all_cpts.size() - 1][0] - v_all_cpts[temp_index][0], v_all_cpts[v_all_cpts.size() - 1][1] - v_all_cpts[temp_index][1], v_all_cpts[v_all_cpts.size() - 1][2] - v_all_cpts[temp_index][2]);
				//total_equation_count++;*/
				/*--------------(equation 36) end------------------------------------------*/

				/*--------------(equation 37) start------------------------------------------*/
				//for (int j = 0; j < (vh_f_bf->degree(0) + 1) * vh_f_bf->degree(1) - 5; j++)
				//{
				//	//c_ij^e coefficient
				//	temp_index = bf_id * 30 + j;
				//	tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//	//c_i,j+1^e coefficient
				//	temp_index = bf_id * 30 + j + 5;
				//	tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//	set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[bf_id * 30 + j][0] - v_all_cpts[temp_index][0], v_all_cpts[bf_id * 30 + j][1] - v_all_cpts[temp_index][1], v_all_cpts[bf_id * 30 + j][2] - v_all_cpts[temp_index][2]);
				//	total_equation_count++;
				//}
				////std::cout << "total_equation_count 37:" << total_equation_count << std::endl;
				/*--------------(equation 37) end------------------------------------------*/
			}
			//std::cout << "total_equation_count:" << total_equation_count << std::endl;
			// Step 1: Construct a sparse matrix G and on the right vector g
			int G_row_num = total_equation_count;  // Total number of equations
			int G_col_num = v_all_cpts.size();     // The number of control points
			SparseMatrix<double> G(G_row_num, G_col_num);

			//Construct the sparse matrix G from tripletlist
			G.setFromTriplets(tripletlists.begin(), tripletlists.end());

			// After constructing all equations, change std::vector to Eigen::VectorXd
			VectorXd g_x = VectorXd::Map(g_x_vec.data(), g_x_vec.size());
			VectorXd g_y = VectorXd::Map(g_y_vec.data(), g_y_vec.size());
			VectorXd g_z = VectorXd::Map(g_z_vec.data(), g_z_vec.size());

			//std::cout << "Matrix G size: " << G.rows() << " #x： " << G.cols() << std::endl;
			//std::cout << "g_x size: " << g_x.size() << std::endl;

			/*------------------------------------------------------------------------------------------------------------------------*/
			//// Step 2: Solve using the least squares method
			//// Convert the sparse matrix G into a dense matrix in order to use SVD
			//MatrixXd G_dense = G.toDense();

			//// Create a SVD solver
			//JacobiSVD<MatrixXd> svd(G_dense, ComputeThinU | ComputeThinV);

			//// Calculate the pseudo-inverse
			//MatrixXd G_pinv = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();

			//// Solve separately x, y, z
			//VectorXd solution_x = G_pinv * g_x;
			//VectorXd solution_y = G_pinv * g_y;
			//VectorXd solution_z = G_pinv * g_z;
			/*------------------------------------------------------------------------------------------------------------------------*/


			/*------------------------------------------------------------------------------------------------------------------------*/
			////Try a different approach to solve it
			//G.makeCompressed();
			//LeastSquaresConjugateGradient<SparseMatrix<double>> Solver_sparse;
			//Solver_sparse.setTolerance(1e-16);
			//Solver_sparse.compute(G);
			//VectorXd solution_x = Solver_sparse.solve(g_x);
			//VectorXd solution_y = Solver_sparse.solve(g_y);
			//VectorXd solution_z = Solver_sparse.solve(g_z);
			/*------------------------------------------------------------------------------------------------------------------------*/

			/*------------------------------------------------------------------------------------------------------------------------*/
			//Solution method 3:
			//Given initial values
			VectorXd solution_x0 = VectorXd::Constant(v_all_cpts.size(), 0.0);
			VectorXd solution_y0 = VectorXd::Constant(v_all_cpts.size(), 0.0);
			VectorXd solution_z0 = VectorXd::Constant(v_all_cpts.size(), 0.0);
			for (int i = 0; i < v_all_cpts.size(); i++)
			{
				solution_x0(i) = v_all_cpts[i][0];
				solution_y0(i) = v_all_cpts[i][1];
				solution_z0(i) = v_all_cpts[i][2];
			}
			G.makeCompressed();
			LeastSquaresConjugateGradient<SparseMatrix<double>> Solver_sparse;
			Solver_sparse.setTolerance(1e-16);
			Solver_sparse.compute(G);
			VectorXd solution_x = Solver_sparse.solveWithGuess(g_x, solution_x0);
			VectorXd solution_y = Solver_sparse.solveWithGuess(g_y, solution_y0);
			VectorXd solution_z = Solver_sparse.solveWithGuess(g_z, solution_z0);
			/*------------------------------------------------------------------------------------------------------------------------*/

			// Obtain the corresponding control point index from v_all_cpts
			int idx = 0;
			//std::cout << "Index\tv_all_cpts\t\t\tSolution_x\t\t\tSolution_y\t\t\tSolution_z\n";

			// Update the control point and redistribute the solution back to the control point
			for (int i = 0; i < vhs.size(); i++) {
				H* vh = vhs[i];
				F* vh_f = pMesh->halfedgeFace(vh);
				auto vh_f_bf = this->bfs()[vh_f->id() - 1];
				// Traverse the control points of the current Bezier surface (an irregular surface)
				for (int j = 0; j < vh_f_bf->cpts_singular().size(); j++)
				{
					// If the current control point coincides with the control point halfedge along the next surface, then skip it
					//FuweiChen update 2024.09.21
					if (j % (vh_f_bf->degree(0) + 1) == 0)
					{
						int id_special = 0;
						if (j == 0)
						{
							//id_special corresponding to the singularity
							id_special = solution_x.size() - 1;
						}
						else
						{
							id_special = j / (vh_f_bf->degree(0) + 1) - 1 + ((i + 1) % v->degree()) * ((vh_f_bf->degree(0)) * (vh_f_bf->degree(0) + 1));

						}
						//std::cout << "id_special: " << id_special << std::endl;
						vh_f_bf->cpts_singular()[j] = CPoint(solution_x[id_special], solution_y[id_special], solution_z[id_special]);
						continue;
					}

					// output intial  v_all_cpts and the control point values after the solution has been obtained
					//std::cout << idx << "\t(" << v_all_cpts[idx][0] << ", " << v_all_cpts[idx][1] << ", " << v_all_cpts[idx][2] << ")\t\t";
					//std::cout << "(" << solution_x[idx] << ", " << solution_y[idx] << ", " << solution_z[idx] << ")\n";

					vh_f_bf->cpts_singular()[j] = CPoint(solution_x[idx], solution_y[idx], solution_z[idx]);
					idx++;
				}

				//std::cout << "----------------------" << std::endl;
			}

			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				F* vh_f = pMesh->halfedgeFace(vh);
				auto vh_f_bf = this->bfs()[vh_f->id() - 1];
				vh_f_bf->cpts().resize((vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1));
				switch (vh->localId())
				{
				case 2:
				{
					//The sequence of the control points remains unchanged.
					for (int i = 0; i < vh_f_bf->cpts_singular().size(); i++)
					{
						vh_f_bf->cpts()[i] = vh_f_bf->cpts_singular()[i];
						vh_f_bf->linkCptPWs()[i] = vh_f_bf->singular_linkCptPWs()[i];
					}
					break;
				}
				case 3:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							int originalIndex = vh_f_bf->degree(0) - i + j * (vh_f_bf->degree(1) + 1);
							vh_f_bf->cpts()[originalIndex] = vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j];
							vh_f_bf->linkCptPWs()[originalIndex] = vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j];
						}
					}
					break;
				}
				case 4:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							int originalIndex = (vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i;
							vh_f_bf->cpts()[originalIndex] = vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j];
							vh_f_bf->linkCptPWs()[originalIndex] = vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j];
						}
					}
					break;
				}
				case 1:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							int originalIndex = vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j;
							vh_f_bf->cpts()[originalIndex] = vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j];
							//std::cout << "before resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							//std::cout << "After resort cpt index: " << vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j << std::endl;
							vh_f_bf->linkCptPWs()[originalIndex] = vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j];
						}
					}

					break;
				}
				default:
				{
					break;
				}
				}
			}

		}
	}
}

void CCG_QMSLib::CCG_QMS_model::G1Constraints_optimization_cws(M* pMesh)
{
	using namespace Eigen;
	auto set_g_vaules = [](std::vector<double>& g_x_vec, std::vector<double>& g_y_vec, std::vector<double>& g_z_vec, double x, double y, double z)
		{
			g_x_vec.push_back(x);
			g_y_vec.push_back(y);
			g_z_vec.push_back(z);
		};
	auto set_g_ws_vaules = [](std::vector<BE_linkingCtlPoints>& g_vec_ws, BE_linkingCtlPoints xyz_w)
		{
			g_vec_ws.push_back(xyz_w);
		};
	for (auto v : It::MVIterator(pMesh))
	{ 
		if (v->ifSingular())
		{
			/*对奇异Bezier曲面的控制点进行重新排序*/
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				F* vh_f = pMesh->halfedgeFace(vh);
				std::shared_ptr<CCG_QMS_bezierSurf> vh_f_bf = this->bfs()[vh_f->id() - 1];
				vh_f_bf->cpts_singular().resize((vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1));
				switch (vh->localId())
				{
				case 2:
				{
					/*控制点的顺序不变*/
					for (int i = 0; i < vh_f_bf->cpts().size(); i++)
					{
						vh_f_bf->cpts_singular()[i] = vh_f_bf->cpts()[i];
						vh_f_bf->singular_linkCptPWs()[i] = vh_f_bf->linkCptPWs()[i];
					}
					break;
				}
				case 3:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							vh_f_bf->cpts_singular()[i * ((vh_f_bf->degree(0) + 1)) + j] = vh_f_bf->cpts()[vh_f_bf->degree(0) - i + j * ((vh_f_bf->degree(1) + 1))];
							/*check*/
							//std::cout << "before resort cpt index: " << vh_f_bf->degree(0) - i + j * ((vh_f_bf->degree(1) + 1)) << std::endl;
							//std::cout << "After resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							vh_f_bf->singular_linkCptPWs()[i * ((vh_f_bf->degree(0) + 1)) + j] = vh_f_bf->linkCptPWs()[vh_f_bf->degree(0) - i + j * ((vh_f_bf->degree(1) + 1))];
						}
					}
					break;
				}
				case 4:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->cpts()[(vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i];
							/*check*/
							//std::cout << "before resort cpt index: " << (vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i << std::endl;
							//std::cout << "After resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->linkCptPWs()[(vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i];
						}
					}

					break;
				}
				case 1:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->cpts()[vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j];
							/*check*/
							//std::cout << "before resort cpt index: " << vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j << std::endl;
							//std::cout << "After resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j] = vh_f_bf->linkCptPWs()[vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j];
						}
					}
					break;
				}
				default:
					break;
				}
			}

			std::vector<double> g_x_vec, g_y_vec, g_z_vec;
			std::vector<BE_linkingCtlPoints> g_vec_ws;
			/*找到与当前奇异点关联的奇异面对应的Bezier曲面*/
			/*找到点周围邻接的半边*/
			std::vector<H*> vhs;
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				vhs.push_back(vh);
			}
			/*列方程组，优化奇异点周围奇异面对应的Beizer曲面上的控制点，让其达到G1连续*/

			/*获取所有控制点*/
			/*给控制点一个顺序*/
			/*第i个半边对应的第j个控制点的index=i*30+j，i=0,...,v->degree(),j=0,...,29(不需要代码显示实现，默认规则就是这样)*/
			std::vector<CPoint> v_all_cpts;
			std::vector<BE_linkingCtlPoints> v_all_cpts_ws;
			for (int i = 0; i < vhs.size(); i++)
			{
				H* vh = vhs[i];
				F* vh_f = pMesh->halfedgeFace(vh);
				std::shared_ptr<CCG_QMS_bezierSurf> vh_f_bf = this->bfs()[vh_f->id() - 1];
				for (int j = 0; j < vh_f_bf->cpts_singular().size(); j++)
				{
					if (j % (vh_f_bf->degree(0) + 1) == 0)
					{
						//std::cout << std::endl;
						continue;
					}
					//std::cout << " " << j;
					v_all_cpts.push_back(vh_f_bf->cpts_singular()[j]);
					v_all_cpts_ws.push_back(vh_f_bf->singular_linkCptPWs()[j]);
				}
				if (i == vhs.size() - 1)
				{
					//将奇异点对应的控制点放到当前控制点数组中，注意只能放一次，并且要放到当前控制点数组的最后
					v_all_cpts.push_back(vh_f_bf->cpts_singular()[0]);
					v_all_cpts_ws.push_back(vh_f_bf->singular_linkCptPWs()[0]);
				}
				//std::cout << "----------------------" << std::endl;
			}

			/*check*/
			//std::cout << "v_all_cpts.size(): " << v_all_cpts.size() << std::endl;

			/*列方程*/
			std::vector<Eigen::Triplet<double>> tripletlists; // 定义一个用于存储所有系数对的容器
			int total_equation_count = 0;//表示的是第几个方程
			for (int i = 0; i < vhs.size(); i++)
			{
				std::vector<Eigen::Triplet<double>> tripletlist; // 定义一个用于存储所有系数对的容器
				H* vh = vhs[i];
				F* vh_f = pMesh->halfedgeFace(vh);
				std::shared_ptr<CCG_QMS_bezierSurf> vh_f_bf = this->bfs()[vh_f->id() - 1];
				int bf_id = i;//对应到论文图中面e的id
				int pre_link_bf_id = ((i - 1) + vhs.size()) % vhs.size();//对应到论文图中面e-1的id
				int next_link_bf_id = ((i + 1) + vhs.size()) % vhs.size();//对应到论文图中面e+1的id
				/*列出w1,w2的公式*/
				double w1 = 0.0, w2 = 0.0;
				if (v->boundary())
				{
					w1 = cos(1.0 * _M_PI / v->degree());
				}
				else
				{
					w1 = cos(2.0 * _M_PI / v->degree());
				}
				V* vhTarget = pMesh->halfedgeTarget(vh);
				if (vhTarget->boundary())
				{
					w2 = cos(1.0 * _M_PI / vhTarget->degree());
				}
				else
				{
					w2 = cos(2.0 * _M_PI / vhTarget->degree());
				}
				//std::cout << "w1: " << w1 << " w2: " << w2 << std::endl;
				int temp_index = -1;
				/*--------------(公式27)开始------------------------------------------*/
				//c_21^e-1对应的系数
				temp_index = pre_link_bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_11^e对应的系数
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0);

				//c_21^e对应的系数
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0 * w1);

				//c_11^e对于的系数
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w1);

				//c_12^e对于的系数
				temp_index = next_link_bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_11^e对于的系数
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				/*--------------(公式27)结束------------------------------------------*/


				/*--------------(公式28)开始------------------------------------------*/
				//c_22^e-1对应的系数
				temp_index = pre_link_bf_id * 30 + 5;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_21^e对应的系数
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w1 - 10.0);

				//c_31^e对应的系数
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -8.0 * w1);

				//c_11^e对于的系数
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -2.0 * w1);

				//c_22^e对于的系数
				temp_index = bf_id * 30 + 5;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				/*--------------(公式28)结束------------------------------------------*/


				/*--------------(公式29)开始------------------------------------------*/
				//c_23^e-1对应的系数
				temp_index = pre_link_bf_id * 30 + 10;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_31^e对应的系数
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0);

				//c_51^e对应的系数
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0 * w1);

				//c_41^e对于的系数
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, 4.0 * w1);

				//c_61^e对于的系数
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, w1);

				//c_21^e对于的系数 
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, w2);

				//c_11^e对于的系数 
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -w2);

				//c_32^e对于的系数 
				temp_index = bf_id * 30 + 6;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				/*--------------(公式29)结束------------------------------------------*/


				/*--------------(公式30)开始------------------------------------------*/
				//c_24^e-1对应的系数
				temp_index = pre_link_bf_id * 30 + 15;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_41^e对应的系数
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0);

				//c_61^e对于的系数
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, -w1);

				//c_51^e对应的系数
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, w1);

				//c_31^e对于的系数
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, 4.0 * w2);

				//c_21^e对于的系数
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0 * w2);

				//c_11^e对应的系数
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, w2);

				//c_42^e对于的系数 
				temp_index = bf_id * 30 + 7;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				/*--------------(公式30)结束------------------------------------------*/


				/*--------------(公式31)开始------------------------------------------*/
				//c_25^e-1对应的系数
				temp_index = pre_link_bf_id * 30 + 20;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_51^e对应的系数
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w2 - 10.0);

				//c_41^e对应的系数
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, -8.0 * w2);

				//c_61^e对于的系数
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, -2.0 * w2);

				//c_52^e对应的系数
				temp_index = bf_id * 30 + 8;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				/*--------------(公式31)结束------------------------------------------*/


				/*--------------(公式32)开始------------------------------------------*/
				//c_26^e-1对应的系数
				temp_index = pre_link_bf_id * 30 + 25;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				;
				//c_61^e对于的系数
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0 * w2 - 10.0);

				//c_51^e对应的系数
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0 * w2);

				//c_62^e对应的系数
				temp_index = bf_id * 30 + 9;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				/*--------------(公式32)结束------------------------------------------*/


				/*--------------(公式33)开始------------------------------------------*/
				//c_11^e对应的系数
				temp_index = v_all_cpts.size() - 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -1.0);

				//c_21^e对于的系数
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 5.0);

				//c_31^e对于的系数
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, -10.0);

				//c_41^e对应的系数
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, 10.0);

				//c_51^e对应的系数
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, -5.0);

				//c_61^e对于的系数
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, 0.0, 0.0, 0.0);
				{
					BE_linkingCtlPoints xyz_w;
					xyz_w.linkingCtlPId_weights.insert(std::pair<int, double>(1, 0.0));
					set_g_ws_vaules(g_vec_ws, xyz_w);
				}
				total_equation_count++;
				//std::cout << "total_equation_count 33:" << total_equation_count << std::endl;
				////std::cout << " g_x_vec.size():" << g_x_vec.size() << std::endl;
				///*--------------(公式33)结束------------------------------------------*/

				///*--------------(公式34-35)开始------------------------------------------*/
				//c_51^e对应的系数
				/**/
				temp_index = bf_id * 30 + 3;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_52^e对应的系数
				temp_index = bf_id * 30 + 8;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_53^e对应的系数
				temp_index = bf_id * 30 + 13;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_54^e对应的系数
				temp_index = bf_id * 30 + 18;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_55^e对应的系数
				temp_index = bf_id * 30 + 23;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_56^e对应的系数
				temp_index = bf_id * 30 + 28;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_61^e对于的系数
				temp_index = bf_id * 30 + 4;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_62^e对于的系数
				temp_index = bf_id * 30 + 9;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_63^e对于的系数
				temp_index = bf_id * 30 + 14;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_64^e对于的系数
				temp_index = bf_id * 30 + 19;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_65^e对于的系数
				temp_index = bf_id * 30 + 24;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_66^e对于的系数
				temp_index = bf_id * 30 + 29;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_26^e对于的系数
				temp_index = bf_id * 30 + 25;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_36^e对于的系数
				temp_index = bf_id * 30 + 26;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_46^e对于的系数
				temp_index = bf_id * 30 + 27;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_25^e对于的系数
				temp_index = bf_id * 30 + 20;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_35^e对于的系数
				temp_index = bf_id * 30 + 21;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_45^e对于的系数
				temp_index = bf_id * 30 + 22;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;
				//std::cout << "total_equation_count:" << total_equation_count << std::endl;//检查正确，25个
				/*--------------(公式34-35)结束------------------------------------------*/

				//std::cout << " g_x_vec.size():" << g_x_vec.size() << std::endl;
				/*--------------(公式36)开始------------------------------------------*/
				//for (int j = 0; j < (vh_f_bf->degree(0) + 1) * vh_f_bf->degree(1); j++)//计算由vh_f_bf->cpts_singular()存储方式的同一个面上表示的控制点
				//{
				//	if ((j + 1) % vh_f_bf->degree(0) == 0)
				//	{
				//		//c_2j^e对于的系数
				//		temp_index = bf_id * 30 + j - 4;
				//		tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//		//c_1,j^e(c_i,1^e+1)对于的系数
				//		if ((j + 1) / vh_f_bf->degree(0) == 1)
				//		{
				//			temp_index = v_all_cpts.size() - 1;
				//			tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//		}
				//		else
				//		{
				//			temp_index = next_link_bf_id * 30 + (j + 1) / vh_f_bf->degree(0) - 1;
				//			tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//		}
				//		set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0] - v_all_cpts[bf_id * 30 + j - 4][0], v_all_cpts[temp_index][1] - v_all_cpts[bf_id * 30 + j - 4][1], v_all_cpts[temp_index][2] - v_all_cpts[bf_id * 30 + j - 4][2]);
				//		{
				//			set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index] - v_all_cpts_ws[bf_id * 30 + j - 4]);
				//		}
				//		total_equation_count++;
				//		//std::cout << "total_equation_count:" << total_equation_count << std::endl;
				//	}
				//	else
				//	{
				//		//c_ij^e对于的系数
				//		temp_index = bf_id * 30 + j;
				//		tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//		//c_i+1,j^e对于的系数
				//		temp_index = bf_id * 30 + j + 1;
				//		tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//		set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[bf_id * 30 + j][0] - v_all_cpts[temp_index][0], v_all_cpts[bf_id * 30 + j][1] - v_all_cpts[temp_index][1], v_all_cpts[bf_id * 30 + j][2] - v_all_cpts[temp_index][2]);
				//		{
				//			set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[bf_id * 30 + j] - v_all_cpts_ws[temp_index]);
				//		}
				//		total_equation_count++;
				//	}
				//}
				////std::cout << "total_equation_count 36:" << total_equation_count << std::endl;
				////c_11^e对应的系数
				////c_21^e对于的系数
				///*temp_index = v_all_cpts.size() - 1;
				//tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//temp_index = bf_id * 30;
				//tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[v_all_cpts.size() - 1][0] - v_all_cpts[temp_index][0], v_all_cpts[v_all_cpts.size() - 1][1] - v_all_cpts[temp_index][1], v_all_cpts[v_all_cpts.size() - 1][2] - v_all_cpts[temp_index][2]);
				//total_equation_count++;*/
				/*--------------(公式36)结束------------------------------------------*/

				/*--------------(公式37)开始------------------------------------------*/
				//for (int j = 0; j < (vh_f_bf->degree(0) + 1) * vh_f_bf->degree(1) - 5; j++)
				//{
				//	//c_ij^e对于的系数
				//	temp_index = bf_id * 30 + j;
				//	tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				//	//c_i,j+1^e对于的系数
				//	temp_index = bf_id * 30 + j + 5;
				//	tripletlists.emplace_back(total_equation_count, temp_index, -1.0);
				//	set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[bf_id * 30 + j][0] - v_all_cpts[temp_index][0], v_all_cpts[bf_id * 30 + j][1] - v_all_cpts[temp_index][1], v_all_cpts[bf_id * 30 + j][2] - v_all_cpts[temp_index][2]);
				//	{
				//		set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[bf_id * 30 + j] - v_all_cpts_ws[temp_index]);
				//	}
				//	total_equation_count++;
				//}
				////std::cout << "total_equation_count 37:" << total_equation_count << std::endl;
				/*--------------(公式37)结束------------------------------------------*/

				/*--------------(公式36-37)的修改版，开始------------------------------------------*/
				//c_33^e
				temp_index = bf_id * 30 + 11;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_43^e
				temp_index = bf_id * 30 + 12;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_34^e
				temp_index = bf_id * 30 + 16;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_44^e
				temp_index = bf_id * 30 + 17;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_21^e
				temp_index = bf_id * 30;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_31^e
				temp_index = bf_id * 30 + 1;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_41^e
				temp_index = bf_id * 30 + 2;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_22^e
				temp_index = bf_id * 30 + 5;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_32^e
				temp_index = bf_id * 30 + 6;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_42^e
				temp_index = bf_id * 30 + 7;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_23^e
				temp_index = bf_id * 30 + 10;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;

				//c_24^e
				temp_index = bf_id * 30 + 15;
				tripletlists.emplace_back(total_equation_count, temp_index, 1.0);
				set_g_vaules(g_x_vec, g_y_vec, g_z_vec, v_all_cpts[temp_index][0], v_all_cpts[temp_index][1], v_all_cpts[temp_index][2]);
				{
					set_g_ws_vaules(g_vec_ws, v_all_cpts_ws[temp_index]);
				}
				total_equation_count++;
				/*--------------(公式36-37)的修改版，结束------------------------------------------*/
			}
			//std::cout << "total_equation_count:" << total_equation_count << std::endl;
			// Step 1: 构建稀疏矩阵 G 和右边的向量 g
			int G_row_num = total_equation_count;  // 方程总数
			int G_col_num = v_all_cpts.size();     // 控制点的数量
			SparseMatrix<double> G(G_row_num, G_col_num);

			// 从 tripletlist 构造稀疏矩阵 G
			G.setFromTriplets(tripletlists.begin(), tripletlists.end());

			// 构造完所有方程后，将 std::vector 转换为 Eigen::VectorXd
			VectorXd g_x = VectorXd::Map(g_x_vec.data(), g_x_vec.size());
			VectorXd g_y = VectorXd::Map(g_y_vec.data(), g_y_vec.size());
			VectorXd g_z = VectorXd::Map(g_z_vec.data(), g_z_vec.size());

			/*std::cout << "Matrix G size: " << G.rows() << "x" << G.cols() << std::endl;
			std::cout << "g_x size: " << g_x.size() << std::endl;

			std::cout << "v->degree(): " << v->degree() << std::endl;*/
			/*------------------------------------------------------------------------------------------------------------------------*/
			// Step 2: 使用最小二乘法求解
			// 将稀疏矩阵G转换为密集矩阵，以便使用SVD
			MatrixXd G_dense = G.toDense();

			// 创建一个SVD求解器
			JacobiSVD<MatrixXd> svd(G_dense, ComputeThinU | ComputeThinV);

			// 计算伪逆
			MatrixXd G_pinv = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();

			// 分别求解 x, y, z
			VectorXd solution_x = G_pinv * g_x;
			VectorXd solution_y = G_pinv * g_y;
			VectorXd solution_z = G_pinv * g_z;
			/*------------------------------------------------------------------------------------------------------------------------*/
			std::vector<BE_linkingCtlPoints> v_all_cpts_solution_ws;
			for (int i = 0; i < G_pinv.rows(); i++)
			{
				BE_linkingCtlPoints temp_w = g_vec_ws[0] * G_pinv(i, 0);
				for (int j = 1; j < G_pinv.cols(); j++)
				{
					temp_w = temp_w + (g_vec_ws[j] * G_pinv(i, j));
				}
				v_all_cpts_solution_ws.push_back(temp_w);
			}


			/*------------------------------------------------------------------------------------------------------------------------*/
			//换种求解方式
			/*G.makeCompressed();
			LeastSquaresConjugateGradient<SparseMatrix<double>> Solver_sparse;
			Solver_sparse.setTolerance(1e-7);
			Solver_sparse.compute(G);
			VectorXd solution_x = Solver_sparse.solve(g_x);
			VectorXd solution_y = Solver_sparse.solve(g_y);
			VectorXd solution_z = Solver_sparse.solve(g_z);*/
			/*------------------------------------------------------------------------------------------------------------------------*/

			// 获取 v_all_cpts 中对应的控制点索引
			int idx = 0;
			// 输出标题
			//std::cout << "Index\tv_all_cpts\t\t\tSolution_x\t\t\tSolution_y\t\t\tSolution_z\n";

			// 更新控制点，将解分配回控制点
			for (int i = 0; i < vhs.size(); i++) {
				H* vh = vhs[i];
				F* vh_f = pMesh->halfedgeFace(vh); // 获取当前半边所在的面
				std::shared_ptr<CCG_QMS_bezierSurf> vh_f_bf = this->bfs()[vh_f->id() - 1]; // 获取当前半边所在的bezier面
				// 遍历当前bezier面（不规则面）的控制点
				for (int j = 0; j < vh_f_bf->cpts_singular().size(); j++)
				{
					// 如果当前的控制点与下一个面出半边上的控制点重合，则跳过
					//FuweiChen update 2024.09.21
					if (j % (vh_f_bf->degree(0) + 1) == 0)
					{
						int id_special = 0;
						if (j == 0)
						{
							//id_special 对应着奇异点
							id_special = solution_x.size() - 1;
						}
						else
						{
							id_special = j / (vh_f_bf->degree(0) + 1) - 1 + ((i + 1) % v->degree()) * ((vh_f_bf->degree(0)) * (vh_f_bf->degree(0) + 1));

						}
						//std::cout << "id_special: " << id_special << std::endl;
						vh_f_bf->cpts_singular()[j] = CPoint(solution_x[id_special], solution_y[id_special], solution_z[id_special]);
						vh_f_bf->singular_linkCptPWs()[j] = v_all_cpts_solution_ws[id_special];
						continue;
					}

					// 输出原始 v_all_cpts 和求解后的控制点值
					//std::cout << idx << "\t(" << v_all_cpts[idx][0] << ", " << v_all_cpts[idx][1] << ", " << v_all_cpts[idx][2] << ")\t\t";
					//std::cout << "(" << solution_x[idx] << ", " << solution_y[idx] << ", " << solution_z[idx] << ")\n";

					vh_f_bf->cpts_singular()[j] = CPoint(solution_x[idx], solution_y[idx], solution_z[idx]);
					vh_f_bf->singular_linkCptPWs()[j] = v_all_cpts_solution_ws[idx];
					idx++;
				}
				//if (i == vhs.size() - 1)
				//{
				//	// 最后处理奇异点
				//	int last_idx = v_all_cpts.size() - 1;

				//	// 输出最后的奇异点
				//	std::cout << last_idx << "\t(" << v_all_cpts[last_idx][0] << ", " << v_all_cpts[last_idx][1] << ", " << v_all_cpts[last_idx][2] << ")\t\t";
				//	std::cout << "(" << solution_x[last_idx] << ", " << solution_y[last_idx] << ", " << solution_z[last_idx] << ")\n";

				//	// 将奇异点对应的控制点放到当前控制点数组中，注意只能放一次，并且要放到当前控制点数组的最后
				//	vh_f_bf->cpts_singular()[0] = CPoint(solution_x[last_idx], solution_y[last_idx], solution_z[last_idx]);
				//}

				//std::cout << "----------------------" << std::endl;
			}
			/*奇异面上的控制点调整回和规则面上控制点顺序一致的状态*/
			for (auto vh : It::VCcwOutHEIterator(pMesh, v))
			{
				F* vh_f = pMesh->halfedgeFace(vh);
				std::shared_ptr<CCG_QMS_bezierSurf> vh_f_bf = this->bfs()[vh_f->id() - 1];
				vh_f_bf->cpts().resize((vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1));
				switch (vh->localId())
				{
				case 2:
				{
					// 控制点顺序不变
					for (int i = 0; i < vh_f_bf->cpts_singular().size(); i++)
					{
						vh_f_bf->cpts()[i] = vh_f_bf->cpts_singular()[i];
						vh_f_bf->linkCptPWs()[i] = vh_f_bf->singular_linkCptPWs()[i];
					}
					break;
				}
				case 3:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							int originalIndex = vh_f_bf->degree(0) - i + j * (vh_f_bf->degree(1) + 1);
							vh_f_bf->cpts()[originalIndex] = vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j];
							vh_f_bf->linkCptPWs()[originalIndex] = vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j];
						}
					}
					break;
				}
				case 4:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							int originalIndex = (vh_f_bf->degree(0) + 1) * (vh_f_bf->degree(1) + 1) - j - 1 - (vh_f_bf->degree(0) + 1) * i;
							vh_f_bf->cpts()[originalIndex] = vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j];
							vh_f_bf->linkCptPWs()[originalIndex] = vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j];
						}
					}
					break;
				}
				case 1:
				{
					for (int i = 0; i < vh_f_bf->degree(0) + 1; i++)
					{
						for (int j = 0; j < vh_f_bf->degree(1) + 1; j++)
						{
							int originalIndex = vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j;
							vh_f_bf->cpts()[originalIndex] = vh_f_bf->cpts_singular()[i * (vh_f_bf->degree(0) + 1) + j];
							//std::cout << "before resort cpt index: " << i * (vh_f_bf->degree(0) + 1) + j << std::endl;
							//std::cout << "After resort cpt index: " << vh_f_bf->degree(0) * (vh_f_bf->degree(1) + 1) + i - (vh_f_bf->degree(0) + 1) * j << std::endl;
							vh_f_bf->linkCptPWs()[originalIndex] = vh_f_bf->singular_linkCptPWs()[i * (vh_f_bf->degree(0) + 1) + j];
						}
					}

					break;
				}
				default:
				{
					break;
				}
				}
			}

		}
	}
}

void CCG_QMSLib::CCG_QMS_model::computeBezierSurfacesKnots(M* pMesh)
{
	for (int i = 0; i < this->bfs().size(); i++)
	{
		auto p = this->bfs()[i];
		F* pf = pMesh->idFace(i + 1);
		H* pfH = pMesh->faceHalfedge(pf);
		H* pfHN = pMesh->halfedgeNext(pfH);
		double u_ki = pMesh->halfedgeEdge(pfH)->knotInterval();
		double v_ki = pMesh->halfedgeEdge(pfHN)->knotInterval();
		//Obtain knot vector from U
		int deg1 = p->degree(0) + 1;
		std::vector<double> u_k;
		u_k.resize(2 * deg1);
		for (int i = 0; i < deg1; i++)
		{
			u_k[i] = 0.0;
			u_k[2 * deg1 - i - 1] = 1.0 * u_ki;
		}
		p->knotsU() = u_k;
		p->maxU() = p->knotsU().back();
		//Obtain knot vector from V
		int deg2 = p->degree(1) + 1;
		std::vector<double> v_k;
		v_k.resize(2 * deg2);
		for (int i = 0; i < deg2; i++)
		{
			v_k[i] = 0.0;
			v_k[2 * deg2 - i - 1] = 1.0 * v_ki;
		}
		p->knotsV() = v_k;
		p->maxV() = p->knotsV().back();
		/*---check bezier knots---*/
		/*std::cout << "---Bezier surface id: " << i << std::endl;
		std::cout << "---knots U: ";
		for (auto tuk : p->knotsU())
		{
			std::cout << " " << tuk ;
		}
		std::cout << std::endl;
		std::cout << "---knots V: ";
			for (auto tvk : p->knotsV())
			{
				std::cout << " " << tvk;
			}
		std::cout << std::endl;*/
	}
}

void CCG_QMSLib::CCG_QMS_model::sortBezierCtlPts(M* pMesh)
{
	for (auto f : It::MFIterator(pMesh))
	{
		/*Note: The halfedge of the surface points to the fourth point ID within the quadrilateral file surface*/
		H* fh = pMesh->faceHalfedge(f);
		H* fhn = pMesh->halfedgeNext(fh);
		auto bf = this->bfs()[f->id() - 1];
		switch (fhn->localId())
		{
		case 1:
		{
			//The order does not need to be changed.
			//std::cout << "------1------" << std::endl;
			break;
		}
		case 2:
		{
			//std::cout << "------2------" << std::endl;
			std::vector<BE_linkingCtlPoints> tempCWs;
			int deg = bf->degree(0);
			for (int i = 0; i <= deg; i++)
			{
				for (int j = 0; j <= deg; j++)
				{
					int index = (deg - i) + (deg + 1) * j;
					tempCWs.push_back(bf->linkCptPWs()[index]);
				}
			}
			bf->linkCptPWs() = tempCWs;
			break;
		}
		case 3:
		{
			//std::cout << "------3------" << std::endl;
			std::vector<BE_linkingCtlPoints> tempCWs;
			int deg = bf->degree(0);
			for (int i = 0; i <= deg; i++)
			{
				for (int j = 0; j <= deg; j++)
				{
					int index = (deg + 1) * (deg + 1) - 1 - (deg + 1) * i - j;
					tempCWs.push_back(bf->linkCptPWs()[index]);
				}
			}
			bf->linkCptPWs() = tempCWs;
			break;
		}
		case 4:
		{
			//std::cout << "------4------" << std::endl;
			std::vector<BE_linkingCtlPoints> tempCWs;
			int deg = bf->degree(0);
			for (int i = 0; i <= deg; i++)
			{
				for (int j = 0; j <= deg; j++)
				{
					int index = deg * (deg + 1) + i - (deg + 1) * j;
					tempCWs.push_back(bf->linkCptPWs()[index]);
				}
			}
			bf->linkCptPWs() = tempCWs;
			break;
		}
		default:
			break;
		}

		/*chek, Is the sum of the weights of each Bezier control point equal to 1.0*/
		//{
		//	for (auto p : bf->linkCptPWs())
		//	{
		//		double ws = 0.0;
		//		for (auto pw : p.linkingCtlPId_weights)
		//		{
		//			ws += pw.second;
		//		}
		//		double error = 1e-5;
		//		if (abs(ws - 1.0) > error)
		//		{
		//			/*for (auto it : tempV.linkingCtlPId_weights)
		//			{
		//				std::cout << " index:" << std::setw(8) << it.first << " weight:" << std::setw(8) << it.second << std::endl;
		//			}*/
		//			std::cout << "faceId: " << f->id() << "  firstVid: " << fhn->target()->id() << std::endl;
		//			std::cout << "controlWeights sum: " << ws << std::endl;
		//		}
		//	}
		//}
	}

	/*check, Is the sum of the weights of each Bezier control point equal to 1.0*/
	/*for (auto bf : bfs())
	{
		std::cout << "---------------bezier surface: " << bf->id() << std::endl;
		for (int i=0;i<bf->linkCptPWs().size(); i++)
		{
			std::cout << "*******bezier control point: "<< i << std::endl;
			double tws = 0.0;
			for (auto w : bf->linkCptPWs()[i].linkingCtlPId_weights)
			{
				std::cout << " index: " << w.first << " weight: " << w.second << std::endl;
				tws += w.second;
			}
			if (abs(tws - 1) > 1e-8)
			{
				std::cout << " total weight: " << tws << std::endl;
			}
		}

	}*/

	/*
	* Re-number the halfedge's local IDs on the grid surface,
	* so that the halfedge with ID = 1 can point to both the first point
	* and the first Bezier control point on the surface,
	* making it convenient for later acquisition of the B-spline surface.
	*/
	for (auto f : It::MFIterator(pMesh))
	{
		bool fB = false;
		for (auto fv : It::FVIterator(pMesh, f))
		{
			if (fv->ifSingular())
			{
				fB = true;
				break;
			}
		}
		if (fB)
		{
			/*Note: The halfedge of the surface points to the fourth point ID within the quadrilateral file surface.*/
			//std::cout << "f->id(): " << f->id() << std::endl;
			H* fh = pMesh->faceHalfedge(f);
			//std::cout << "change before he id: " << fh->localId() << std::endl;
			fh->localId() = 4;
			//std::cout << "change after he id: " << fh->localId() << std::endl;
			H* fhn = pMesh->halfedgeNext(fh);
			//std::cout << "change before he id: " << fhn->localId() << std::endl;
			fhn->localId() = 1;
			//std::cout << "change after he id: " << fhn->localId() << std::endl;
			H* fhnn = pMesh->halfedgeNext(fhn);
			//std::cout << "change before he id: " << fhnn->localId() << std::endl;
			fhnn->localId() = 2;
			//std::cout << "change after he id: " << fhnn->localId() << std::endl;
			H* fhnnn = pMesh->halfedgeNext(fhnn);
			//std::cout << "change before he id: " << fhnnn->localId() << std::endl;
			fhnnn->localId() = 3;
			//std::cout << "change after he id: " << fhnnn->localId() << std::endl;
		}
		else
		{
			/*Note: The halfedge of the surface points to the fourth point ID within the quadrilateral file surface.*/
			//std::cout << "f->id(): " << f->id() << std::endl;
			H* fh = pMesh->faceHalfedge(f);
			//std::cout << "change before he id: " << fh->localId();
			fh->localId() = 4;
			//std::cout << "change after he id: " << fh->localId();
			H* fhn = pMesh->halfedgeNext(fh);
			//std::cout << "change before he id: " << fhn->localId();
			fhn->localId() = 1;
			//std::cout << "change after he id: " << fhn->localId();
			H* fhnn = pMesh->halfedgeNext(fhn);
			//std::cout << "change before he id: " << fhnn->localId();
			fhnn->localId() = 2;
			//std::cout << "change after he id: " << fhnn->localId();
			H* fhnnn = pMesh->halfedgeNext(fhnn);
			//std::cout << "change before he id: " << fhnnn->localId();
			fhnnn->localId() = 3;
			//std::cout << "change after he id: " << fhnnn->localId();
		}
	}
}

void CCG_QMSLib::CCG_QMS_model::obtainBezierCptsBy_cpts(double scale)
{
	/*Initialize all Bezier control points to {0, 0, 0}*/
	for (auto bf : this->bfs())
	{
		/*CPoint tempP;*/
		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			//bf->cpts()[i] = tempP;
			bf->cpts()[i] = bf->cpts()[i] / scale;
		}
	}
	/*Calculate the Bezier control points through the initial control points and weights.*/
	/*for (auto bf : this->bfs())
	{
		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			CPoint tempP;
			for (auto cwIt : bf->linkCptPWs()[i].linkingCtlPId_weights)
			{
				tempP += ((this->controlPoints()[cwIt.first - 1] * cwIt.second));
			}
			bf->cpts()[i] = tempP/scale;
		}
	}*/

}

void CCG_QMSLib::CCG_QMS_model::obtainBSplineSurfacePatches_general_UVConsistency_topology_boundary(M* pMesh)
{
	/*Making sure the quad mesh marked the information of quadLayout*/
	/*Marking the boundary edges to sharp edges*/
	for (auto e : It::MEIterator(pMesh))
	{
		if (e->boundary())
		{
			e->sharp() = true;
		}
	}

	/*
	* Marking quadLayout's corner whose number of sharp edges is not equal to 2
	* (In fact, the corner on the boundary has 2 sharp edges and its degree is 1 )
	*/
	for (auto v : It::MVIterator(pMesh))
	{
		int numVSharpE = 0;
		for (auto ve : It::VCcwEIterator(pMesh, v))
		{
			if (ve->sharp())
			{
				numVSharpE++;
			}
		}
		if (numVSharpE > 2)
		{
			v->ifCorner() = true;
		}
		else if (numVSharpE == 2)
		{
			int numF = 0;
			for (auto vf : It::VCcwFIterator(pMesh, v))
			{
				numF++;
			}
			if (numF == 1)
			{
				v->ifCorner() = true;
			}
		}
	}

	/*Numbering the quad in different patch, and identify the faces that contain a single boundary edge on the surfaces without singular points.*/
	std::vector<F*> patches_noSingularities;//Store the faces that contain a single boundary edge on a surface without singular points.
	int patch_num = 0;//Count the obtained patches
	for (auto f : It::MFIterator(pMesh))
	{
		//Find one unvisited quad
		if (f->f_spline_patchId() == 0)
		{
			bool if_patch_noSingularities = true;
			F* patch_noSingularities_boundaryFace = NULL;
			//Increment the current patch counter by 1
			patch_num++;
			//update f_spline_patchId
			f->f_spline_patchId() = patch_num;
			for (auto fh : It::FHEIterator(pMesh, f))
			{
				V* fhv = pMesh->halfedgeTarget(fh);
				if (fhv->ifCorner())
				{
					if_patch_noSingularities = false;
					/*--------------------------------------------------------------*/
					//cosider the case including T node 2025/02/27,cfw
					int numFSharpEdges = 0;
					for (auto temp_fh : It::FHEIterator(pMesh, f))
					{
						if (pMesh->halfedgeEdge(temp_fh)->sharp())
						{
							numFSharpEdges++;
						}
					}
					if (numFSharpEdges == 1)//T node
					{
						if_patch_noSingularities = true;
					}
					/*--------------------------------------------------------------*/
				}
				E* fhe = pMesh->halfedgeEdge(fh);
				if (fhe->sharp() && patch_noSingularities_boundaryFace == NULL)
				{
					patch_noSingularities_boundaryFace = f;
				}
			}
			//Use breadth-first traversal to obtain other faces that share the same patch with the current face.
			std::queue<F*> qfs;
			qfs.push(f);
			while (!qfs.empty())
			{
				F* qf = qfs.front();
				qfs.pop();
				for (auto qfhe : It::FHEIterator(pMesh, qf))
				{
					E* qfhee = pMesh->halfedgeEdge(qfhe);
					if (qfhee->sharp())continue;
					H* qfheSym = pMesh->halfedgeSym(qfhe);
					F* qfheSymf = pMesh->halfedgeFace(qfheSym);
					if (qfheSymf->f_spline_patchId() != 0)continue;
					qfheSymf->f_spline_patchId() = patch_num;
					for (auto fh : It::FHEIterator(pMesh, qfheSymf))
					{
						V* fhv = pMesh->halfedgeTarget(fh);
						if (fhv->ifCorner())
						{
							if_patch_noSingularities = false;
							/*--------------------------------------------------------------*/
							//cosider the case including T node 2025/02/27,cfw
							int numFSharpEdges = 0;
							for (auto temp_fh : It::FHEIterator(pMesh, qfheSymf))
							{
								if (pMesh->halfedgeEdge(temp_fh)->sharp())
								{
									numFSharpEdges++;
								}
							}
							if (numFSharpEdges == 1)//T node
							{
								if_patch_noSingularities = true;
							}
							/*--------------------------------------------------------------*/
						}
						E* fhe = pMesh->halfedgeEdge(fh);
						if (fhe->sharp() && patch_noSingularities_boundaryFace == NULL)
						{
							patch_noSingularities_boundaryFace = qfheSymf;
						}
					}
					qfs.push(qfheSymf);
				}
			}
			//Regarding the 1-genus model, for example: torus
			if (patch_noSingularities_boundaryFace == NULL)
			{
				patch_noSingularities_boundaryFace = f;
			}
			if (if_patch_noSingularities)
			{
				patches_noSingularities.push_back(patch_noSingularities_boundaryFace);
			}
		}
	}
	//std::cout << "# patches: " << patch_num << std::endl;
	//std::cout << "#patches_noSingularities: " << patches_noSingularities.size() << std::endl;]

	/*Initialize the model's patches*/
	for (int i = 0; i < patch_num; i++)
	{
		auto m_p = std::make_shared<CCG_QMS_bsplineSurf>();
		m_p->id() = i + 1;
		this->bsplineSurfs().push_back(m_p);
	}
	//std::cout << "# initial model bsplineSurfs: " << bsplineSurfs().size() << std::endl;

	/*
	* 5.0 First, process the surfaces without singularities.
	*/
	for (auto regularF : patches_noSingularities)
	{
		//Find the halfedge corresponding to a sharp edge in the surface.
		H* vhe = NULL;
		for (auto fh : It::FHEIterator(pMesh, regularF))
		{
			E* fhe = pMesh->halfedgeEdge(fh);
			if (fhe->sharp())
			{
				vhe = fh;
				break;
			}
		}
		if (vhe == NULL)
		{
			vhe = pMesh->faceHalfedge(regularF);
			//For the case where there is a T-node, for the edges that are not marked as "sharp", they should be skipped.
			E* vheE = pMesh->halfedgeEdge(vhe);
			E* vPreheE = pMesh->halfedgeEdge(pMesh->halfedgePrev(vhe));
			//By using the patchId of the face where halfedge is located to determine whether the current patch has been accessed before
			F* vheF = pMesh->halfedgeFace(vhe);
			int vhef_spline_patchId = vheF->f_spline_patchId();
			//The counting of patchId starts from 1, while vector access starts from 0.
			vhef_spline_patchId--;
			auto currentBS = this->bsplineSurfs()[vhef_spline_patchId];
			currentBS->id() = vheF->f_spline_patchId();
			/*At this point, the UV directions of the surface patches are all periodic and closed.*/
			currentBS->closedAlongU() = true;
			currentBS->periodicAlongU() = true;
			currentBS->closedAlongV() = true;
			currentBS->periodicAlongV() = true;
			if (currentBS->cpts().size() != 0)continue;
			//Find the number of control points in the V direction.
			std::vector<H*> vHes;//The halfedge in the V boundary direction
			vHes.push_back(vhe);
			/*
			* The half section in the V direction starts from the source of vhe and ends at the source of the next section,
			* which is equal to the end of the source of vhe (Is it possible that this constraint is not satisfied? cfw 2024/7/29
			* (It won't be unsatisfied. If it is not satisfied, the current patch is not a circular structure. 2025/3/1))
			*/
			H* vheN = pMesh->halfedgeNext(vhe);
			V* vheNSource = pMesh->halfedgeSource(vheN);
			V* vheSource = pMesh->halfedgeSource(vhe);
			while (vheNSource->id() != vheSource->id())
			{
				H* vheNSym = pMesh->halfedgeSym(vheN);
				vhe = pMesh->halfedgeNext(vheNSym);
				vHes.push_back(vhe);
				vheN = pMesh->halfedgeNext(vhe);
				vheNSource = pMesh->halfedgeSource(vheN);
			}

			/*Update the current information of QB_bsplineSurf*/
			int f_bfId = vheF->id();
			f_bfId--;//The id of the quadrilateral mesh face starts from 1.
			currentBS->degree()[0] = this->bfs()[f_bfId]->degree(0);//Update degree
			currentBS->degree()[1] = this->bfs()[f_bfId]->degree(1);
			currentBS->cpts().resize(vHes.size() * currentBS->degree()[0] + 1);//Update the number of control points in the u direction
			/*Assign the face ID corresponding to the first halfedge of vHes to the current B-spline surface patch.*/
			currentBS->quadMeshFaceId() = vHes[0]->face()->id();

			/*knot vector assignment (applicable to general cases)*/
			//Iniitialize patch's u_knots
			double vKnotMax = 0.0;
			std::vector<double> temp_uknots;
			for (auto th : vHes)
			{
				E* the = pMesh->halfedgeEdge(th);
				vKnotMax += the->knotInterval();
				temp_uknots.push_back(vKnotMax);
			}
			currentBS->knotsU().resize(currentBS->degree()[0] * (vHes.size() + 1) + 2, 0.0);
			currentBS->knotsU()[0] = 0.0;
			for (int i = 0; i < currentBS->degree()[0]; i++)
			{
				currentBS->knotsU()[1 + i] = 0.0;
			}
			currentBS->knotsU()[currentBS->knotsU().size() - 1] = vKnotMax;
			for (int i = 0; i < vHes.size(); i++)
			{
				for (int j = 0; j < currentBS->degree()[1]; j++)
				{
					currentBS->knotsU()[currentBS->degree()[0] + 1 + currentBS->degree()[0] * i + j] = temp_uknots[i];
				}
			}
			/*check u knots*/
			/*std::cout << "U knots: ";
			for (int i = 0; i < currentBS->knotsU().size(); i++)
			{
				std::cout << " " << currentBS->knotsU()[i];
			}
			std::cout << std::endl;*/

			//Obtain the control points in the patch, and traverse the mesh vertices within the patch in the order of U first and then V.
			//Starting from one halfedge of the U boundary and then another halfedge, for each halfedge, a column of control points in the V direction will be selected.
			/*Assign local IDs to the surfaces, so that the quadrilateral mesh surfaces can be matched with the IDs of the B-spline surfaces.*/
			int count_localId = 0;
			for (int i = 0; i < vHes.size(); i++)
			{
				H* v_he = vHes[i];
				std::vector<H*> u_hes;//U-directional halfedge set
				//Place the second boundary point in the U direction into the vertices of the patch.
				H* v_heP = pMesh->halfedgePrev(v_he);
				V* v_heS = pMesh->halfedgeSource(v_he);
				u_hes.push_back(v_heP);
				F* u_hePF = pMesh->halfedgeFace(v_heP);
				u_hePF->localId() = ++count_localId;
				H* vhePP = pMesh->halfedgePrev(v_heP);
				E* vhePPE = pMesh->halfedgeEdge(vhePP);
				while (vhePP->target()->id() != v_heS->id())
				{
					v_he = pMesh->halfedgeSym(vhePP);
					v_heP = pMesh->halfedgePrev(v_he);
					u_hes.push_back(v_heP);
					u_hePF = pMesh->halfedgeFace(v_heP);
					u_hePF->localId() = ++count_localId;
					vhePP = pMesh->halfedgePrev(v_heP);
					vhePPE = pMesh->halfedgeEdge(vhePP);
				}

				//Resize the number of V-direction control points in the current bspline Surf
				for (int j = 0; j < currentBS->degree()[0]; j++)
				{
					currentBS->cpts()[i * currentBS->degree()[0] + j].resize(u_hes.size() * currentBS->degree()[1] + 1);
				}

				if (i == 0)
				{
					/*knot vector assignment (applicable to general cases)*/
					//Iniitialize patch's v_knots
					double uKnotMax = 0.0;
					std::vector<double> temp_vknots;
					for (auto th : u_hes)
					{
						E* the = pMesh->halfedgeEdge(th);
						uKnotMax += the->knotInterval();
						temp_vknots.push_back(uKnotMax);
					}
					currentBS->knotsV().resize(currentBS->degree()[1] * (u_hes.size() + 1) + 2, 0.0);
					currentBS->knotsV()[0] = 0.0;
					for (int i = 0; i < currentBS->degree()[1]; i++)
					{
						currentBS->knotsV()[1 + i] = 0.0;
					}
					currentBS->knotsV()[currentBS->knotsU().size() - 1] = uKnotMax;
					for (int i = 0; i < u_hes.size(); i++)
					{
						for (int j = 0; j < currentBS->degree()[0]; j++)
						{
							currentBS->knotsV()[currentBS->degree()[0] + 1 + currentBS->degree()[0] * i + j] = temp_vknots[i];
						}
					}

					/*check v knots*/
					/*std::cout << "V knots: ";
					for (int i = 0; i < currentBS->knotsV().size(); i++)
					{
						std::cout << " " << currentBS->knotsV()[i];
					}
					std::cout << std::endl;*/
				}

				/*Assign values to the control points*/
				for (int j = 0; j < u_hes.size(); j++)
				{
					H* u_he = u_hes[j];
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					auto h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[jv + ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(jv + 1) * (currentBS->degree()[0] + 1) - 1 - ju];
							}
						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv - ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1]) - jv * (currentBS->degree()[0] + 1) + ju];
							}
						}
						break;
					}
					default:
						break;
					}
				}
				/*Assign values to the last control point in the U direction for each one.*/
				{
					H* u_he = u_hes.back();
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					auto h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0]) + jv];

						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1]) - jv];

						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					default:
						break;
					}
				}

				if (i == vHes.size() - 1)
				{
					currentBS->cpts()[vHes.size() * currentBS->degree()[1]].resize(u_hes.size() * currentBS->degree()[0] + 1);
					/*Assign values to the control points*/
					for (int j = 0; j < u_hes.size(); j++)
					{
						H* u_he = u_hes[j];
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						auto h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[1] + 1) - 1 + (currentBS->degree()[1] + 1) * ju];
							}

							break;
						}
						case 2:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - ju];
							}

							break;
						}
						case 3:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[1] + 1) - ju * (currentBS->degree()[1] + 1)];
							}

							break;
						}
						case 4:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[ju];
							}

							break;
						}
						default:
							break;
						}
					}
					/*Assign values to the last control point in the U direction for each one*/
					{
						H* u_he = u_hes.back();
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						auto h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0] + 1) - 1];
							break;
						}
						case 2:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[0] + 1)];
							break;
						}
						case 3:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[0];
							break;
						}
						case 4:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0])];
							break;
						}
						default:
							break;
						}
					}
				}
			}
		}
		else
		{
			//For the case where there is a T-node, for the edges that are not marked as "sharp", they should be skipped.
			E* vheE = pMesh->halfedgeEdge(vhe);
			E* vPreheE = pMesh->halfedgeEdge(pMesh->halfedgePrev(vhe));
			//By using the patchId of the face where halfedge is located to determine whether the current patch has been accessed before
			F* vheF = pMesh->halfedgeFace(vhe);
			int vhef_spline_patchId = vheF->f_spline_patchId();
			//The counting of patchId starts from 1, while vector access starts from 0.
			vhef_spline_patchId--;
			auto currentBS = this->bsplineSurfs()[vhef_spline_patchId];
			currentBS->id() = vheF->f_spline_patchId();
			/*At this point,the U directions of the surface patches is periodic and closed*/
			currentBS->closedAlongU() = true;
			currentBS->periodicAlongU() = true;
			if (currentBS->cpts().size() != 0)continue;
			//Find the number of control points in the V direction.
			std::vector<H*> vHes;//The halfedge in the U boundary direction
			vHes.push_back(vhe);
			//The halfedge in the U direction starts from the source of vhe and ends at the source of the next section, which is equal to the end of the source of vhe.
			H* vheN = pMesh->halfedgeNext(vhe);
			V* vheNSource = pMesh->halfedgeSource(vheN);
			V* vheSource = pMesh->halfedgeSource(vhe);
			while (vheNSource->id() != vheSource->id())
			{
				H* vheNSym = pMesh->halfedgeSym(vheN);
				vhe = pMesh->halfedgeNext(vheNSym);
				vHes.push_back(vhe);
				vheN = pMesh->halfedgeNext(vhe);
				vheNSource = pMesh->halfedgeSource(vheN);
			}

			/*
			* The source of the first half of vHes is the first point of the quadrilateral mesh corresponding to the B-spline surface.
			* The last half of the target of vhes is the second point of the quadrilateral mesh corresponding to the B-spline surface.
			*/
			currentBS->quadIds()[0] = pMesh->vertexId(pMesh->halfedgeSource(vHes.front()));
			pMesh->halfedgeSource(vHes.front())->ifCorner() = true;
			currentBS->quadIds()[1] = pMesh->vertexId(pMesh->halfedgeTarget(vHes.back()));
			pMesh->halfedgeTarget(vHes.back())->ifCorner() = true;

			/*Update the current information of QB_bsplineSurf*/
			int f_bfId = vheF->id();
			f_bfId--;//The id of the quadrilateral mesh face starts from 1.
			currentBS->degree()[0] = this->bfs()[f_bfId]->degree(0);//Update degree
			currentBS->degree()[1] = this->bfs()[f_bfId]->degree(1);
			currentBS->cpts().resize(vHes.size() * currentBS->degree()[0] + 1);//Update the number of control points in the u direction
			/*Assign the face ID corresponding to the first halfedge of vHes to the current B-spline surface patch*/
			currentBS->quadMeshFaceId() = vHes[0]->face()->id();

			/*knot vector assignment (applicable to general cases)*/
			//Iniitialize patch's u_knots
			double vKnotMax = 0.0;
			std::vector<double> temp_uknots;
			for (auto th : vHes)
			{
				E* the = pMesh->halfedgeEdge(th);
				vKnotMax += the->knotInterval();
				temp_uknots.push_back(vKnotMax);
			}
			currentBS->knotsU().resize(currentBS->degree()[0] * (vHes.size() + 1) + 2, 0.0);
			currentBS->knotsU()[0] = 0.0;
			for (int i = 0; i < currentBS->degree()[0]; i++)
			{
				currentBS->knotsU()[1 + i] = 0.0;
			}
			currentBS->knotsU()[currentBS->knotsU().size() - 1] = vKnotMax;
			for (int i = 0; i < vHes.size(); i++)
			{
				for (int j = 0; j < currentBS->degree()[1]; j++)
				{
					currentBS->knotsU()[currentBS->degree()[1] + 1 + currentBS->degree()[1] * i + j] = temp_uknots[i];
				}
			}
			/*check u knots*/
			/*std::cout << "U knots: ";
			for (int i = 0; i < currentBS->knotsU().size(); i++)
			{
				std::cout << " " << currentBS->knotsU()[i];
			}
			std::cout << std::endl;*/

			//Obtain the control points in the patch, and traverse the mesh vertices within the patch in the order of U first and then V.
			//Starting from one halfedge of the U boundary and then another halfedge, for each halfedge, a column of control points in the V direction will be selected.
			for (int i = 0; i < vHes.size(); i++)
			{
				H* v_he = vHes[i];
				std::vector<H*> u_hes;//U-directional halfedge set
				//Place the second boundary point in the U direction into the vertices of the patch.
				H* v_heP = pMesh->halfedgePrev(v_he);
				u_hes.push_back(v_heP);
				H* vhePP = pMesh->halfedgePrev(v_heP);
				E* vhePPE = pMesh->halfedgeEdge(vhePP);
				while (!vhePPE->sharp())
				{
					v_he = pMesh->halfedgeSym(vhePP);
					v_heP = pMesh->halfedgePrev(v_he);
					u_hes.push_back(v_heP);
					vhePP = pMesh->halfedgePrev(v_heP);
					vhePPE = pMesh->halfedgeEdge(vhePP);
				}

				//Resize the number of V-direction control points in the current bspline Surf
				for (int j = 0; j < currentBS->degree()[1]; j++)
				{
					currentBS->cpts()[i * currentBS->degree()[0] + j].resize(u_hes.size() * currentBS->degree()[1] + 1);
				}

				if (i == 0)
				{
					/*The last half of u_hes's source is the fourth point of the quadrilateral mesh corresponding to the B-spline surface.*/
					currentBS->quadIds()[3] = pMesh->vertexId(pMesh->halfedgeSource(u_hes.back()));
					pMesh->halfedgeSource(u_hes.back())->ifCorner() = true;
					/*knot vector assignment (applicable to general cases)*/
					//Iniitialize patch's v_knots
					double uKnotMax = 0.0;
					std::vector<double> temp_vknots;
					for (auto th : u_hes)
					{
						E* the = pMesh->halfedgeEdge(th);
						uKnotMax += the->knotInterval();
						temp_vknots.push_back(uKnotMax);
					}
					currentBS->knotsV().resize(currentBS->degree()[1] * (u_hes.size() + 1) + 2, 0.0);
					currentBS->knotsV()[0] = 0.0;
					for (int i = 0; i < currentBS->degree()[1]; i++)
					{
						currentBS->knotsV()[1 + i] = 0.0;
					}
					currentBS->knotsV()[currentBS->knotsV().size() - 1] = uKnotMax;
					for (int i = 0; i < u_hes.size(); i++)
					{
						for (int j = 0; j < currentBS->degree()[1]; j++)
						{
							currentBS->knotsV()[currentBS->degree()[1] + 1 + currentBS->degree()[1] * i + j] = temp_vknots[i];
						}
					}

					/*check v knots*/
					/*std::cout << "V knots: ";
					for (int i = 0; i < currentBS->knotsV().size(); i++)
					{
						std::cout << " " << currentBS->knotsV()[i];
					}
					std::cout << std::endl;*/
				}

				/*Assign values to the control points*/
				for (int j = 0; j < u_hes.size(); j++)
				{
					H* u_he = u_hes[j];
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					auto h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[jv + ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(jv + 1) * (currentBS->degree()[0] + 1) - 1 - ju];
							}
						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv - ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1]) - jv * (currentBS->degree()[0] + 1) + ju];
							}
						}
						break;
					}
					default:
						break;
					}
				}
				/*Assign values to the last control point in the U direction for each one.*/
				{
					H* u_he = u_hes.back();
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					auto h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0]) + jv];

						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1]) - jv];

						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					default:
						break;
					}
				}

				//Move the last row of boundary points in the V direction into the "vertices" of the patch.
				if (i == vHes.size() - 1)
				{
					/*The source of the "prev" attribute on the last half of u_hes is the third point of the quadrilateral mesh corresponding to the B-spline surface.*/
					currentBS->quadIds()[2] = pMesh->vertexId(pMesh->halfedgeSource(pMesh->halfedgePrev(u_hes.back())));
					pMesh->halfedgeSource(pMesh->halfedgePrev(u_hes.back()))->ifCorner() = true;
					currentBS->cpts()[vHes.size() * currentBS->degree()[1]].resize(u_hes.size() * currentBS->degree()[0] + 1);
					/*Assign values to the control points*/
					for (int j = 0; j < u_hes.size(); j++)
					{
						H* u_he = u_hes[j];
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						auto h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[1] + 1) - 1 + (currentBS->degree()[1] + 1) * ju];
							}

							break;
						}
						case 2:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - ju];
							}

							break;
						}
						case 3:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[1] + 1) - ju * (currentBS->degree()[1] + 1)];
							}

							break;
						}
						case 4:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[ju];
							}

							break;
						}
						default:
							break;
						}
					}
					/*Assign values to the last control point in the U direction for each one*/
					{
						H* u_he = u_hes.back();
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						auto h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0] + 1) - 1];
							break;
						}
						case 2:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[0] + 1)];
							break;
						}
						case 3:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[0];
							break;
						}
						case 4:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0])];
							break;
						}
						default:
							break;
						}
					}
				}
			}
		}
	}

	/*Building patches starts from corner*/
	 /*Assign value to model's patches*/
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->ifCorner())
		{
			//Search for the unassigned patch from the outHalfedge of the singularity point
			for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
			{
				//For the case where there is a T-node, for the edges that are not marked as "sharp", they should be skipped.
				E* vheE = pMesh->halfedgeEdge(vhe);
				E* vPreheE = pMesh->halfedgeEdge(pMesh->halfedgePrev(vhe));
				if ((vheE->sharp() && (!vPreheE->sharp())) || (!vheE->sharp()))continue;
				F* vheF = pMesh->halfedgeFace(vhe);
				int vhef_spline_patchId = vheF->f_spline_patchId();
				//The counting of patchId starts from 1, while vector access starts from 0.
				vhef_spline_patchId--;
				auto currentBS = this->bsplineSurfs()[vhef_spline_patchId];
				currentBS->id() = vheF->f_spline_patchId();
				if (currentBS->cpts().size() != 0)continue;
				//Find the number of control points in the U direction.
				std::vector<H*> vHes;//The halfedge in the U boundary direction
				vHes.clear();
				vHes.push_back(vhe);
				//The halfedge in the U direction starts from one corner and ends at the other corner.
				//That is, until the edge where the next halfedge of the halfedge starting from the U direction is located becomes a sharp edge.
				H* vheN = pMesh->halfedgeNext(vhe);
				E* vheNE = pMesh->halfedgeEdge(vheN);
				while (!vheNE->sharp())
				{
					H* vheNSym = pMesh->halfedgeSym(vheN);
					vhe = pMesh->halfedgeNext(vheNSym);
					vHes.push_back(vhe);
					vheN = pMesh->halfedgeNext(vhe);
					vheNE = pMesh->halfedgeEdge(vheN);
				}
				/*
				* The source of the first half of vHes is the first point of the quadrilateral mesh corresponding to the B-spline surface.
				* The last half of the target of vhes is the second point of the quadrilateral mesh corresponding to the B-spline surface.
				*/
				currentBS->quadIds()[0] = pMesh->vertexId(pMesh->halfedgeSource(vHes.front()));
				currentBS->quadIds()[1] = pMesh->vertexId(pMesh->halfedgeTarget(vHes.back()));

				/*Update the current information of QB_bsplineSurf*/
				int f_bfId = vheF->id();
				f_bfId--;//The id of the quadrilateral mesh face starts from 1.
				currentBS->degree()[0] = this->bfs()[f_bfId]->degree(0);//Update degree
				currentBS->degree()[1] = this->bfs()[f_bfId]->degree(1);
				currentBS->cpts().resize(vHes.size() * currentBS->degree()[0] + 1);//Update the number of control points in the u direction
				/*Assign the face ID corresponding to the first halfedge of vHes to the current B-spline surface patch*/
				currentBS->quadMeshFaceId() = vHes[0]->face()->id();

				/*knot vector assignment (applicable to general cases)*/
				//Iniitialize patch's u_knots
				double vKnotMax = 0.0;
				std::vector<double> temp_uknots;
				for (auto th : vHes)
				{
					E* the = pMesh->halfedgeEdge(th);
					vKnotMax += the->knotInterval();
					temp_uknots.push_back(vKnotMax);
				}
				currentBS->knotsU().resize(currentBS->degree()[0] * (vHes.size() + 1) + 2, 0.0);
				currentBS->knotsU()[0] = 0.0;
				for (int i = 0; i < currentBS->degree()[0]; i++)
				{
					currentBS->knotsU()[1 + i] = 0.0;
				}
				currentBS->knotsU()[currentBS->knotsU().size() - 1] = vKnotMax;
				for (int i = 0; i < vHes.size(); i++)
				{
					for (int j = 0; j < currentBS->degree()[1]; j++)
					{
						currentBS->knotsU()[currentBS->degree()[0] + 1 + currentBS->degree()[0] * i + j] = temp_uknots[i];
					}
				}
				/*check u knots*/
				/*std::cout << "U knots: ";
				for (int i = 0; i < currentBS->knotsU().size(); i++)
				{
					std::cout << " " << currentBS->knotsU()[i];
				}
				std::cout << std::endl;*/

				//Obtain the control points in the patch, and traverse the mesh vertices within the patch in the order of U first and then V.
				//Starting from one halfedge of the U boundary and then another halfedge, for each halfedge, a column of control points in the V direction will be selected.
				for (int i = 0; i < vHes.size(); i++)
				{
					H* v_he = vHes[i];
					std::vector<H*> u_hes;//V directional halfedge set
					//Place the second boundary point in the V direction into the vertices of the patch.
					H* v_heP = pMesh->halfedgePrev(v_he);
					u_hes.push_back(v_heP);
					H* vhePP = pMesh->halfedgePrev(v_heP);
					E* vhePPE = pMesh->halfedgeEdge(vhePP);
					while (!vhePPE->sharp())
					{
						v_he = pMesh->halfedgeSym(vhePP);
						v_heP = pMesh->halfedgePrev(v_he);
						u_hes.push_back(v_heP);
						vhePP = pMesh->halfedgePrev(v_heP);
						vhePPE = pMesh->halfedgeEdge(vhePP);
					}

					//Resize the number of V-direction control points in the current bspline Surf
					for (int j = 0; j < currentBS->degree()[1]; j++)
					{
						currentBS->cpts()[i * currentBS->degree()[0] + j].resize(u_hes.size() * currentBS->degree()[1] + 1);
					}

					if (i == 0)
					{
						/*The last half of u_hes's source is the fourth point of the quadrilateral mesh corresponding to the B-spline surface.*/
						currentBS->quadIds()[3] = pMesh->vertexId(pMesh->halfedgeSource(u_hes.back()));
						/*knot vector assignment (applicable to general cases)*/
						//Iniitialize patch's v_knots
						double uKnotMax = 0.0;
						std::vector<double> temp_vknots;
						for (auto th : u_hes)
						{
							E* the = pMesh->halfedgeEdge(th);
							uKnotMax += the->knotInterval();
							temp_vknots.push_back(uKnotMax);
						}
						currentBS->knotsV().resize(currentBS->degree()[1] * (u_hes.size() + 1) + 2, 0.0);
						currentBS->knotsV()[0] = 0.0;
						for (int i = 0; i < currentBS->degree()[0]; i++)
						{
							currentBS->knotsV()[1 + i] = 0.0;
						}
						currentBS->knotsV()[currentBS->knotsV().size() - 1] = uKnotMax;
						for (int i = 0; i < u_hes.size(); i++)
						{
							for (int j = 0; j < currentBS->degree()[1]; j++)
							{
								currentBS->knotsV()[currentBS->degree()[1] + 1 + currentBS->degree()[1] * i + j] = temp_vknots[i];
							}
						}

						/*check v knots*/
						/*std::cout << "V knots: ";
						for (int i = 0; i < currentBS->knotsV().size(); i++)
						{
							std::cout << " " << currentBS->knotsV()[i];
						}
						std::cout << std::endl;*/
					}

					/*Assign values to the control points*/
					for (int j = 0; j < u_hes.size(); j++)
					{
						H* u_he = u_hes[j];
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						auto h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[jv + ju * (currentBS->degree()[1] + 1)];
								}
							}
							break;
						}
						case 2:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(jv + 1) * (currentBS->degree()[0] + 1) - 1 - ju];
								}
							}
							break;
						}
						case 3:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv - ju * (currentBS->degree()[1] + 1)];
								}
							}
							break;
						}
						case 4:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1]) - jv * (currentBS->degree()[0] + 1) + ju];
								}
							}
							break;
						}
						default:
							break;
						}
					}
					/*Assign values to the last control point in the U direction for each one*/
					{
						H* u_he = u_hes.back();
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						auto h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0]) + jv];

							}
							break;
						}
						case 2:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[jv * (currentBS->degree()[0] + 1)];

							}
							break;
						}
						case 3:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1]) - jv];

							}
							break;
						}
						case 4:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv * (currentBS->degree()[0] + 1)];

							}
							break;
						}
						default:
							break;
						}
					}

					//Move the last row of boundary points in the V direction into the "vertices" of the patch.
					if (i == vHes.size() - 1)
					{
						/*The source of the "prev" attribute on the last half of u_hes is the third point of the quadrilateral mesh corresponding to the B-spline surface.*/
						currentBS->quadIds()[2] = pMesh->vertexId(pMesh->halfedgeSource(pMesh->halfedgePrev(u_hes.back())));

						currentBS->cpts()[vHes.size() * currentBS->degree()[1]].resize(u_hes.size() * currentBS->degree()[0] + 1);
						/*Assign values to the control points*/
						for (int j = 0; j < u_hes.size(); j++)
						{
							H* u_he = u_hes[j];
							F* h_heF = pMesh->halfedgeFace(u_he);
							int h_heFId = h_heF->id();
							h_heFId--;
							auto h_heBf = this->bfs()[h_heFId];
							switch (u_he->localId())
							{
							case 1:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[1] + 1) - 1 + (currentBS->degree()[1] + 1) * ju];
								}

								break;
							}
							case 2:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - ju];
								}

								break;
							}
							case 3:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[1] + 1) - ju * (currentBS->degree()[1] + 1)];
								}

								break;
							}
							case 4:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[ju];
								}

								break;
							}
							default:
								break;
							}
						}
						/*Assign values to the last control point in the U direction for each one*/
						{
							H* u_he = u_hes.back();
							F* h_heF = pMesh->halfedgeFace(u_he);
							int h_heFId = h_heF->id();
							h_heFId--;
							auto h_heBf = this->bfs()[h_heFId];
							switch (u_he->localId())
							{
							case 1:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0] + 1) - 1];
								break;
							}
							case 2:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[0] + 1)];
								break;
							}
							case 3:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[0];
								break;
							}
							case 4:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0])];
								break;
							}
							default:
								break;
							}
						}
					}
				}
			}
		}
	}
	//std::cout << "Obtaining B spline surfaces finished!" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::obtainBSplineSurfacePatches_general_UVConsistency_recordQuadBoundaryNURBSBoundaryCurveRelation(M* pMesh)
{
	/*Making sure the quad mesh marked the information of quadLayout*/
	/*Marking the boundary edges to sharp edges*/
	for (auto e : It::MEIterator(pMesh))
	{
		if (e->boundary())
		{
			e->sharp() = true;
		}
	}

	/*
	* Marking quadLayout's corner whose number of sharp edges is not equal to 2
	* (In fact, the corner on the boundary has 2 sharp edges and its degree is 1 )
	*/
	for (auto v : It::MVIterator(pMesh))
	{
		int numVSharpE = 0;
		for (auto ve : It::VCcwEIterator(pMesh, v))
		{
			if (ve->sharp())
			{
				numVSharpE++;
			}
		}
		if (numVSharpE > 2)
		{
			v->ifCorner() = true;
		}
		else if (numVSharpE == 2)
		{
			int numF = 0;
			for (auto vf : It::VCcwFIterator(pMesh, v))
			{
				numF++;
			}
			if (numF == 1)
			{
				v->ifCorner() = true;
			}
		}
	}

	/*Numbering the quad in different patch, and identify the faces that contain a single boundary edge on the surfaces without singular points.*/
	std::vector<F*> patches_noSingularities;//Store the faces that contain a single boundary edge on a surface without singular points.
	int patch_num = 0;//Count the obtained patches
	for (auto f : It::MFIterator(pMesh))
	{
		//Find one unvisited quad
		if (f->f_spline_patchId() == 0)
		{
			bool if_patch_noSingularities = true;
			F* patch_noSingularities_boundaryFace = NULL;
			//Increment the current patch counter by 1
			patch_num++;
			//update f_spline_patchId
			f->f_spline_patchId() = patch_num;
			for (auto fh : It::FHEIterator(pMesh, f))
			{
				V* fhv = pMesh->halfedgeTarget(fh);
				if (fhv->ifCorner())
				{
					if_patch_noSingularities = false;
					/*--------------------------------------------------------------*/
					//cosider the case including T node 2025/02/27,cfw
					int numFSharpEdges = 0;
					for (auto temp_fh : It::FHEIterator(pMesh, f))
					{
						if (pMesh->halfedgeEdge(temp_fh)->sharp())
						{
							numFSharpEdges++;
						}
					}
					if (numFSharpEdges == 1)//T node
					{
						if_patch_noSingularities = true;
					}
					/*--------------------------------------------------------------*/
				}
				E* fhe = pMesh->halfedgeEdge(fh);
				if (fhe->sharp() && patch_noSingularities_boundaryFace == NULL)
				{
					patch_noSingularities_boundaryFace = f;
				}
			}
			//Use breadth-first traversal to obtain other faces that share the same patch with the current face.
			std::queue<F*> qfs;
			qfs.push(f);
			while (!qfs.empty())
			{
				F* qf = qfs.front();
				qfs.pop();
				for (auto qfhe : It::FHEIterator(pMesh, qf))
				{
					E* qfhee = pMesh->halfedgeEdge(qfhe);
					if (qfhee->sharp())continue;
					H* qfheSym = pMesh->halfedgeSym(qfhe);
					F* qfheSymf = pMesh->halfedgeFace(qfheSym);
					if (qfheSymf->f_spline_patchId() != 0)continue;
					qfheSymf->f_spline_patchId() = patch_num;
					for (auto fh : It::FHEIterator(pMesh, qfheSymf))
					{
						V* fhv = pMesh->halfedgeTarget(fh);
						if (fhv->ifCorner())
						{
							if_patch_noSingularities = false;
							/*--------------------------------------------------------------*/
							//cosider the case including T node 2025/02/27,cfw
							int numFSharpEdges = 0;
							for (auto temp_fh : It::FHEIterator(pMesh, qfheSymf))
							{
								if (pMesh->halfedgeEdge(temp_fh)->sharp())
								{
									numFSharpEdges++;
								}
							}
							if (numFSharpEdges == 1)//T node
							{
								if_patch_noSingularities = true;
							}
							/*--------------------------------------------------------------*/
						}
						E* fhe = pMesh->halfedgeEdge(fh);
						if (fhe->sharp() && patch_noSingularities_boundaryFace == NULL)
						{
							patch_noSingularities_boundaryFace = qfheSymf;
						}
					}
					qfs.push(qfheSymf);
				}
			}
			//Regarding the 1-genus model, for example: torus
			if (patch_noSingularities_boundaryFace == NULL)
			{
				patch_noSingularities_boundaryFace = f;
			}
			if (if_patch_noSingularities)
			{
				patches_noSingularities.push_back(patch_noSingularities_boundaryFace);
			}
		}
	}
	//std::cout << "# patches: " << patch_num << std::endl;
	//std::cout << "#patches_noSingularities: " << patches_noSingularities.size() << std::endl;

	/*Initialize the model's patches*/
	for (int i = 0; i < patch_num; i++)
	{
		auto m_p = std::make_shared<CCG_QMS_bsplineSurf>();
		this->bsplineSurfs().push_back(m_p);
	}
	//std::cout << "# initial model bsplineSurfs: " << bsplineSurfs().size() << std::endl;

	/*
	* 5.0 First, process the surfaces without singularities.
	*/
	for (auto regularF : patches_noSingularities)
	{
		//Find the halfedge corresponding to a sharp edge in the surface.
		H* vhe = NULL;
		for (auto fh : It::FHEIterator(pMesh, regularF))
		{
			E* fhe = pMesh->halfedgeEdge(fh);
			if (fhe->sharp())
			{
				vhe = fh;
				break;
			}
		}
		if (vhe == NULL)
		{
			vhe = pMesh->faceHalfedge(regularF);
			//For the case where there is a T-node, for the edges that are not marked as "sharp", they should be skipped.
			E* vheE = pMesh->halfedgeEdge(vhe);
			E* vPreheE = pMesh->halfedgeEdge(pMesh->halfedgePrev(vhe));
			//By using the patchId of the face where halfedge is located to determine whether the current patch has been accessed before
			F* vheF = pMesh->halfedgeFace(vhe);
			int vhef_spline_patchId = vheF->f_spline_patchId();
			//The counting of patchId starts from 1, while vector access starts from 0.
			vhef_spline_patchId--;
			std::shared_ptr<CCG_QMS_bsplineSurf> currentBS = this->bsplineSurfs()[vhef_spline_patchId];
			/*At this point, the UV directions of the surface patches are all periodic and closed.*/
			currentBS->closedAlongU() = true;
			currentBS->periodicAlongU() = true;
			currentBS->closedAlongV() = true;
			currentBS->periodicAlongV() = true;
			if (currentBS->cpts().size() != 0)continue;
			currentBS->id() = vhef_spline_patchId + 1;
			//Find the number of control points in the V direction.
			std::vector<H*> vHes;//The halfedge in the V boundary direction
			vHes.push_back(vhe);
			/*
			* The half section in the V direction starts from the source of vhe and ends at the source of the next section,
			* which is equal to the end of the source of vhe (Is it possible that this constraint is not satisfied? cfw 2024/7/29
			* (It won't be unsatisfied. If it is not satisfied, the current patch is not a circular structure. 2025/3/1))
			*/
			H* vheN = pMesh->halfedgeNext(vhe);
			V* vheNSource = pMesh->halfedgeSource(vheN);
			V* vheSource = pMesh->halfedgeSource(vhe);
			while (vheNSource->id() != vheSource->id())
			{
				H* vheNSym = pMesh->halfedgeSym(vheN);
				vhe = pMesh->halfedgeNext(vheNSym);
				vHes.push_back(vhe);
				vheN = pMesh->halfedgeNext(vhe);
				vheNSource = pMesh->halfedgeSource(vheN);
			}

			/*Update the current information of QB_bsplineSurf*/
			int f_bfId = vheF->id();
			f_bfId--;//The id of the quadrilateral mesh face starts from 1.
			currentBS->degree()[0] = this->bfs()[f_bfId]->degree(0);//Update degree
			currentBS->degree()[1] = this->bfs()[f_bfId]->degree(1);
			currentBS->cpts().resize(vHes.size() * currentBS->degree()[0] + 1);//Update the number of control points in the u direction
			/*Assign the face ID corresponding to the first halfedge of vHes to the current B-spline surface patch.*/
			currentBS->quadMeshFaceId() = vHes[0]->face()->id();

			/*knot vector assignment (applicable to general cases)*/
			//Iniitialize patch's u_knots
			double vKnotMax = 0.0;
			std::vector<double> temp_uknots;
			for (auto th : vHes)
			{
				E* the = pMesh->halfedgeEdge(th);
				vKnotMax += the->knotInterval();
				temp_uknots.push_back(vKnotMax);
			}
			currentBS->knotsU().resize(currentBS->degree()[0] * (vHes.size() + 1) + 2, 0.0);
			currentBS->knotsU()[0] = 0.0;
			for (int i = 0; i < currentBS->degree()[0]; i++)
			{
				currentBS->knotsU()[1 + i] = 0.0;
			}
			currentBS->knotsU()[currentBS->knotsU().size() - 1] = vKnotMax;
			for (int i = 0; i < vHes.size(); i++)
			{
				for (int j = 0; j < currentBS->degree()[1]; j++)
				{
					currentBS->knotsU()[currentBS->degree()[0] + 1 + currentBS->degree()[0] * i + j] = temp_uknots[i];
				}
			}
			/*check u knots*/
			/*std::cout << "U knots: ";
			for (int i = 0; i < currentBS->knotsU().size(); i++)
			{
				std::cout << " " << currentBS->knotsU()[i];
			}
			std::cout << std::endl;*/

			//Obtain the control points in the patch, and traverse the mesh vertices within the patch in the order of U first and then V.
			//Starting from one halfedge of the U boundary and then another halfedge, for each halfedge, a column of control points in the V direction will be selected.
			/*Assign local IDs to the surfaces, so that the quadrilateral mesh surfaces can be matched with the IDs of the B-spline surfaces.*/
			int count_localId = 0;
			for (int i = 0; i < vHes.size(); i++)
			{
				H* v_he = vHes[i];
				v_he->NURBS_patchId() = currentBS->id();
				v_he->isoParameterLineType() = 2;
				std::vector<H*> u_hes;//V-directional halfedge set
				//Place the second boundary point in the V direction into the vertices of the patch.
				H* v_heP = pMesh->halfedgePrev(v_he);
				V* v_heS = pMesh->halfedgeSource(v_he);
				u_hes.push_back(v_heP);
				F* u_hePF = pMesh->halfedgeFace(v_heP);
				u_hePF->localId() = ++count_localId;
				H* vhePP = pMesh->halfedgePrev(v_heP);
				E* vhePPE = pMesh->halfedgeEdge(vhePP);
				while (vhePP->target()->id() != v_heS->id())
				{
					v_he = pMesh->halfedgeSym(vhePP);
					v_heP = pMesh->halfedgePrev(v_he);
					u_hes.push_back(v_heP);
					u_hePF = pMesh->halfedgeFace(v_heP);
					u_hePF->localId() = ++count_localId;
					vhePP = pMesh->halfedgePrev(v_heP);
					vhePPE = pMesh->halfedgeEdge(vhePP);
				}

				//Resize the number of V-direction control points in the current bspline Surf
				for (int j = 0; j < currentBS->degree()[0]; j++)
				{
					currentBS->cpts()[i * currentBS->degree()[0] + j].resize(u_hes.size() * currentBS->degree()[1] + 1);
				}

				if (i == 0)
				{
					/*knot vector assignment (applicable to general cases)*/
					//Iniitialize patch's v_knots
					double uKnotMax = 0.0;
					std::vector<double> temp_vknots;
					for (auto th : u_hes)
					{
						E* the = pMesh->halfedgeEdge(th);
						uKnotMax += the->knotInterval();
						temp_vknots.push_back(uKnotMax);
						th->NURBS_patchId() = currentBS->id();
						th->isoParameterLineType() = 1;
					}
					currentBS->knotsV().resize(currentBS->degree()[1] * (u_hes.size() + 1) + 2, 0.0);
					currentBS->knotsV()[0] = 0.0;
					for (int i = 0; i < currentBS->degree()[1]; i++)
					{
						currentBS->knotsV()[1 + i] = 0.0;
					}
					currentBS->knotsV()[currentBS->knotsU().size() - 1] = uKnotMax;
					for (int i = 0; i < u_hes.size(); i++)
					{
						for (int j = 0; j < currentBS->degree()[0]; j++)
						{
							currentBS->knotsV()[currentBS->degree()[0] + 1 + currentBS->degree()[0] * i + j] = temp_vknots[i];
						}
					}

					/*check v knots*/
					/*std::cout << "V knots: ";
					for (int i = 0; i < currentBS->knotsV().size(); i++)
					{
						std::cout << " " << currentBS->knotsV()[i];
					}
					std::cout << std::endl;*/
				}

				/*Assign values to the control points*/
				for (int j = 0; j < u_hes.size(); j++)
				{
					H* u_he = u_hes[j];
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[jv + ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(jv + 1) * (currentBS->degree()[0] + 1) - 1 - ju];
							}
						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv - ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1]) - jv * (currentBS->degree()[0] + 1) + ju];
							}
						}
						break;
					}
					default:
						break;
					}
				}
				/*Assign values to the last control point in the V direction for each one.*/
				{
					H* u_he = u_hes.back();
					pMesh->halfedgePrev(u_he)->NURBS_patchId() = currentBS->id();
					pMesh->halfedgePrev(u_he)->isoParameterLineType() = 4;
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0]) + jv];

						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1]) - jv];

						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					default:
						break;
					}
				}

				//Move the last column of boundary points in the U direction into the control points of the patch.
				if (i == vHes.size() - 1)
				{
					currentBS->cpts()[vHes.size() * currentBS->degree()[1]].resize(u_hes.size() * currentBS->degree()[0] + 1);
					/*Assign values to the control points*/
					for (int j = 0; j < u_hes.size(); j++)
					{
						H* u_he = u_hes[j];
						pMesh->halfedgeNext((pMesh->halfedgeNext(u_he)))->NURBS_patchId() = currentBS->id();
						pMesh->halfedgeNext((pMesh->halfedgeNext(u_he)))->isoParameterLineType() = 3;
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[1] + 1) - 1 + (currentBS->degree()[1] + 1) * ju];
							}

							break;
						}
						case 2:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - ju];
							}

							break;
						}
						case 3:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[1] + 1) - ju * (currentBS->degree()[1] + 1)];
							}

							break;
						}
						case 4:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[ju];
							}

							break;
						}
						default:
							break;
						}
					}
					/*Assign values to the last control point in the U direction for each one*/
					{
						H* u_he = u_hes.back();
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0] + 1) - 1];
							break;
						}
						case 2:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[0] + 1)];
							break;
						}
						case 3:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[0];
							break;
						}
						case 4:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0])];
							break;
						}
						default:
							break;
						}
					}
				}
			}
		}
		else
		{
			//For the case where there is a T-node, for the edges that are not marked as "sharp", they should be skipped.
			E* vheE = pMesh->halfedgeEdge(vhe);
			E* vPreheE = pMesh->halfedgeEdge(pMesh->halfedgePrev(vhe));
			//By using the patchId of the face where halfedge is located to determine whether the current patch has been accessed before
			F* vheF = pMesh->halfedgeFace(vhe);
			int vhef_spline_patchId = vheF->f_spline_patchId();
			//The counting of patchId starts from 1, while vector access starts from 0.
			vhef_spline_patchId--;
			std::shared_ptr<CCG_QMS_bsplineSurf> currentBS = this->bsplineSurfs()[vhef_spline_patchId];
			/*At this point,the U directions of the surface patches is periodic and closed*/
			currentBS->closedAlongU() = true;
			currentBS->periodicAlongU() = true;
			if (currentBS->cpts().size() != 0)continue;
			currentBS->id() = vhef_spline_patchId + 1;
			//Find the number of control points in the V direction.
			std::vector<H*> vHes;//The halfedge in the U boundary direction
			vHes.push_back(vhe);
			//The halfedge in the U direction starts from the source of vhe and ends at the source of the next halfedge, which is equal to the end of the source of vhe.
			H* vheN = pMesh->halfedgeNext(vhe);
			V* vheNSource = pMesh->halfedgeSource(vheN);
			V* vheSource = pMesh->halfedgeSource(vhe);
			while (vheNSource->id() != vheSource->id())
			{
				H* vheNSym = pMesh->halfedgeSym(vheN);
				vHes.push_back(vhe);
				vheN = pMesh->halfedgeNext(vhe);
				vheNSource = pMesh->halfedgeSource(vheN);
			}

			/*Update the current information of QB_bsplineSurf*/
			int f_bfId = vheF->id();
			f_bfId--;//The id of the quadrilateral mesh face starts from 1.
			currentBS->degree()[0] = this->bfs()[f_bfId]->degree(0);//Update degree
			currentBS->degree()[1] = this->bfs()[f_bfId]->degree(1);
			currentBS->cpts().resize(vHes.size() * currentBS->degree()[0] + 1);//Update the number of control points in the u direction
			/*Assign the face ID corresponding to the first halfedge of vHes to the current B-spline surface patch*/
			currentBS->quadMeshFaceId() = vHes[0]->face()->id();

			/*knot vector assignment (applicable to general cases)*/
			//Iniitialize patch's u_knots
			double vKnotMax = 0.0;
			std::vector<double> temp_uknots;
			for (auto th : vHes)
			{
				E* the = pMesh->halfedgeEdge(th);
				vKnotMax += the->knotInterval();
				temp_uknots.push_back(vKnotMax);
			}
			currentBS->knotsU().resize(currentBS->degree()[0] * (vHes.size() + 1) + 2, 0.0);
			currentBS->knotsU()[0] = 0.0;
			for (int i = 0; i < currentBS->degree()[0]; i++)
			{
				currentBS->knotsU()[1 + i] = 0.0;
			}
			currentBS->knotsU()[currentBS->knotsU().size() - 1] = vKnotMax;
			for (int i = 0; i < vHes.size(); i++)
			{
				for (int j = 0; j < currentBS->degree()[1]; j++)
				{
					currentBS->knotsU()[currentBS->degree()[1] + 1 + currentBS->degree()[1] * i + j] = temp_uknots[i];
				}
			}
			/*check u knots*/
			/*std::cout << "U knots: ";
			for (int i = 0; i < currentBS->knotsU().size(); i++)
			{
				std::cout << " " << currentBS->knotsU()[i];
			}
			std::cout << std::endl;*/

			//Obtain the control points in the patch, and traverse the mesh vertices within the patch in the order of U first and then V.
			//Starting from one halfedge of the U boundary and then another halfedge, for each halfedge, a column of control points in the V direction will be selected.
			for (int i = 0; i < vHes.size(); i++)
			{
				H* v_he = vHes[i];
				v_he->NURBS_patchId() = currentBS->id();
				v_he->isoParameterLineType() = 2;
				std::vector<H*> u_hes;//U-directional halfedge set
				//Place the second boundary point in the U direction into the vertices of the patch.
				H* v_heP = pMesh->halfedgePrev(v_he);
				u_hes.push_back(v_heP);
				H* vhePP = pMesh->halfedgePrev(v_heP);
				E* vhePPE = pMesh->halfedgeEdge(vhePP);
				while (!vhePPE->sharp())
				{
					v_he = pMesh->halfedgeSym(vhePP);
					v_heP = pMesh->halfedgePrev(v_he);
					u_hes.push_back(v_heP);
					vhePP = pMesh->halfedgePrev(v_heP);
					vhePPE = pMesh->halfedgeEdge(vhePP);
				}

				//Resize the number of V-direction control points in the current bspline Surf
				for (int j = 0; j < currentBS->degree()[1]; j++)
				{
					currentBS->cpts()[i * currentBS->degree()[0] + j].resize(u_hes.size() * currentBS->degree()[1] + 1);
				}

				if (i == 0)
				{
					/*knot vector assignment (applicable to general cases)*/
					//Iniitialize patch's v_knots
					double uKnotMax = 0.0;
					std::vector<double> temp_vknots;
					for (auto th : u_hes)
					{
						E* the = pMesh->halfedgeEdge(th);
						uKnotMax += the->knotInterval();
						temp_vknots.push_back(uKnotMax);
						th->NURBS_patchId() = currentBS->id();
						th->isoParameterLineType() = 1;
					}
					currentBS->knotsV().resize(currentBS->degree()[1] * (u_hes.size() + 1) + 2, 0.0);
					currentBS->knotsV()[0] = 0.0;
					for (int i = 0; i < currentBS->degree()[1]; i++)
					{
						currentBS->knotsV()[1 + i] = 0.0;
					}
					currentBS->knotsV()[currentBS->knotsV().size() - 1] = uKnotMax;
					for (int i = 0; i < u_hes.size(); i++)
					{
						for (int j = 0; j < currentBS->degree()[1]; j++)
						{
							currentBS->knotsV()[currentBS->degree()[1] + 1 + currentBS->degree()[1] * i + j] = temp_vknots[i];
						}
					}

					/*check v knots*/
					/*std::cout << "V knots: ";
					for (int i = 0; i < currentBS->knotsV().size(); i++)
					{
						std::cout << " " << currentBS->knotsV()[i];
					}
					std::cout << std::endl;*/
				}

				/*Assign values to the control points*/
				for (int j = 0; j < u_hes.size(); j++)
				{
					H* u_he = u_hes[j];
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[jv + ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(jv + 1) * (currentBS->degree()[0] + 1) - 1 - ju];
							}
						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv - ju * (currentBS->degree()[1] + 1)];
							}
						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{
							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1]) - jv * (currentBS->degree()[0] + 1) + ju];
							}
						}
						break;
					}
					default:
						break;
					}
				}
				/*Assign values to the last control point in the V direction for each one.*/
				{
					H* u_he = u_hes.back();
					pMesh->halfedgePrev(u_he)->NURBS_patchId() = currentBS->id();
					pMesh->halfedgePrev(u_he)->isoParameterLineType() = 4;
					F* h_heF = pMesh->halfedgeFace(u_he);
					int h_heFId = h_heF->id();
					h_heFId--;
					std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
					switch (u_he->localId())
					{
					case 1:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0]) + jv];

						}
						break;
					}
					case 2:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					case 3:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1]) - jv];

						}
						break;
					}
					case 4:
					{
						for (int jv = 0; jv < currentBS->degree()[1]; jv++)
						{

							currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv * (currentBS->degree()[0] + 1)];

						}
						break;
					}
					default:
						break;
					}
				}

				//Move the last row of boundary points in the V direction into the "vertices" of the patch.
				if (i == vHes.size() - 1)
				{
					currentBS->cpts()[vHes.size() * currentBS->degree()[1]].resize(u_hes.size() * currentBS->degree()[0] + 1);
					/*Assign values to the control points*/
					for (int j = 0; j < u_hes.size(); j++)
					{
						H* u_he = u_hes[j];
						pMesh->halfedgeNext(pMesh->halfedgeNext(u_he))->NURBS_patchId() = currentBS->id();
						pMesh->halfedgeNext(pMesh->halfedgeNext(u_he))->isoParameterLineType() = 3;
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[1] + 1) - 1 + (currentBS->degree()[1] + 1) * ju];
							}

							break;
						}
						case 2:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - ju];
							}

							break;
						}
						case 3:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[1] + 1) - ju * (currentBS->degree()[1] + 1)];
							}

							break;
						}
						case 4:
						{

							for (int ju = 0; ju < currentBS->degree()[0]; ju++)
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[ju];
							}

							break;
						}
						default:
							break;
						}
					}
					/*Assign values to the last control point in the U direction for each one*/
					{
						H* u_he = u_hes.back();
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0] + 1) - 1];
							break;
						}
						case 2:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[0] + 1)];
							break;
						}
						case 3:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[0];
							break;
						}
						case 4:
						{
							currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0])];
							break;
						}
						default:
							break;
						}
					}
				}
			}
		}
	}

	/*Building patches starts from corner*/
	 /*Assign value to model's patches*/
	for (auto v : It::MVIterator(pMesh))
	{
		if (v->ifCorner())
		{
			//Search for the unassigned patch from the outHalfedge of the singularity point
			for (auto vhe : It::VCcwOutHEIterator(pMesh, v))
			{
				//For the case where there is a T-node, for the edges that are not marked as "sharp", they should be skipped.
				E* vheE = pMesh->halfedgeEdge(vhe);
				E* vPreheE = pMesh->halfedgeEdge(pMesh->halfedgePrev(vhe));
				if ((vheE->sharp() && (!vPreheE->sharp())) || (!vheE->sharp()))continue;
				F* vheF = pMesh->halfedgeFace(vhe);
				int vhef_spline_patchId = vheF->f_spline_patchId();
				//The counting of patchId starts from 1, while vector access starts from 0.
				vhef_spline_patchId--;
				std::shared_ptr<CCG_QMS_bsplineSurf> currentBS = this->bsplineSurfs()[vhef_spline_patchId];
				if (currentBS->cpts().size() != 0)continue;
				currentBS->id() = vhef_spline_patchId + 1;
				//Find the number of control points in the U direction.
				std::vector<H*> vHes;//The halfedge in the U boundary direction
				vHes.push_back(vhe);
				//The halfedge in the U direction starts from one corner and ends at the other corner.
				//That is, until the edge where the next halfedge of the halfedge starting from the U direction is located becomes a sharp edge.
				H* vheN = pMesh->halfedgeNext(vhe);
				E* vheNE = pMesh->halfedgeEdge(vheN);
				while (!vheNE->sharp())
				{
					H* vheNSym = pMesh->halfedgeSym(vheN);
					vhe = pMesh->halfedgeNext(vheNSym);
					vHes.push_back(vhe);
					vheN = pMesh->halfedgeNext(vhe);
					vheNE = pMesh->halfedgeEdge(vheN);
				}

				/*Update the current information of QB_bsplineSurf*/
				int f_bfId = vheF->id();
				f_bfId--;//The id of the quadrilateral mesh face starts from 1.
				currentBS->degree()[0] = this->bfs()[f_bfId]->degree(0);//Update degree
				currentBS->degree()[1] = this->bfs()[f_bfId]->degree(1);
				currentBS->cpts().resize(vHes.size() * currentBS->degree()[0] + 1);//Update the number of control points in the u direction
				/*Assign the face ID corresponding to the first halfedge of vHes to the current B-spline surface patch*/
				currentBS->quadMeshFaceId() = vHes[0]->face()->id();

				/*knot vector assignment (applicable to general cases)*/
				//Iniitialize patch's u_knots
				double vKnotMax = 0.0;
				std::vector<double> temp_uknots;
				for (auto th : vHes)
				{
					E* the = pMesh->halfedgeEdge(th);
					vKnotMax += the->knotInterval();
					temp_uknots.push_back(vKnotMax);
				}
				currentBS->knotsU().resize(currentBS->degree()[0] * (vHes.size() + 1) + 2, 0.0);
				currentBS->knotsU()[0] = 0.0;
				for (int i = 0; i < currentBS->degree()[0]; i++)
				{
					currentBS->knotsU()[1 + i] = 0.0;
				}
				currentBS->knotsU()[currentBS->knotsU().size() - 1] = vKnotMax;
				for (int i = 0; i < vHes.size(); i++)
				{
					for (int j = 0; j < currentBS->degree()[1]; j++)
					{
						currentBS->knotsU()[currentBS->degree()[0] + 1 + currentBS->degree()[0] * i + j] = temp_uknots[i];
					}
				}
				/*check u knots*/
				/*std::cout << "U knots: ";
				for (int i = 0; i < currentBS->knotsU().size(); i++)
				{
					std::cout << " " << currentBS->knotsU()[i];
				}
				std::cout << std::endl;*/

				//Obtain the control points in the patch, and traverse the mesh vertices within the patch in the order of U first and then V.
				//Starting from one halfedge of the U boundary and then another halfedge, for each halfedge, a column of control points in the V direction will be selected.
				for (int i = 0; i < vHes.size(); i++)
				{
					H* v_he = vHes[i];
					v_he->NURBS_patchId() = currentBS->id();
					v_he->isoParameterLineType() = 2;
					std::vector<H*> u_hes;//V directional halfedge set
					//Place the second boundary point in the V direction into the vertices of the patch.
					H* v_heP = pMesh->halfedgePrev(v_he);
					u_hes.push_back(v_heP);
					H* vhePP = pMesh->halfedgePrev(v_heP);
					E* vhePPE = pMesh->halfedgeEdge(vhePP);
					while (!vhePPE->sharp())
					{
						v_he = pMesh->halfedgeSym(vhePP);
						v_heP = pMesh->halfedgePrev(v_he);
						u_hes.push_back(v_heP);
						vhePP = pMesh->halfedgePrev(v_heP);
						vhePPE = pMesh->halfedgeEdge(vhePP);
					}

					//Resize the number of V-direction control points in the current bspline Surf
					for (int j = 0; j < currentBS->degree()[1]; j++)
					{
						currentBS->cpts()[i * currentBS->degree()[0] + j].resize(u_hes.size() * currentBS->degree()[1] + 1);
					}

					if (i == 0)
					{
						/*knot vector assignment (applicable to general cases)*/
						//Iniitialize patch's v_knots
						double uKnotMax = 0.0;
						std::vector<double> temp_vknots;
						for (auto th : u_hes)
						{
							E* the = pMesh->halfedgeEdge(th);
							uKnotMax += the->knotInterval();
							temp_vknots.push_back(uKnotMax);
							th->NURBS_patchId() = currentBS->id();
							th->isoParameterLineType() = 1;
						}
						currentBS->knotsV().resize(currentBS->degree()[1] * (u_hes.size() + 1) + 2, 0.0);
						currentBS->knotsV()[0] = 0.0;
						for (int i = 0; i < currentBS->degree()[0]; i++)
						{
							currentBS->knotsV()[1 + i] = 0.0;
						}
						currentBS->knotsV()[currentBS->knotsV().size() - 1] = uKnotMax;
						for (int i = 0; i < u_hes.size(); i++)
						{
							for (int j = 0; j < currentBS->degree()[1]; j++)
							{
								currentBS->knotsV()[currentBS->degree()[1] + 1 + currentBS->degree()[1] * i + j] = temp_vknots[i];
							}
						}

						/*check v knots*/
						/*std::cout << "V knots: ";
						for (int i = 0; i < currentBS->knotsV().size(); i++)
						{
							std::cout << " " << currentBS->knotsV()[i];
						}
						std::cout << std::endl;*/
					}

					/*Assign values to the control points*/
					for (int j = 0; j < u_hes.size(); j++)
					{
						H* u_he = u_hes[j];
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[jv + ju * (currentBS->degree()[1] + 1)];
								}
							}
							break;
						}
						case 2:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(jv + 1) * (currentBS->degree()[0] + 1) - 1 - ju];
								}
							}
							break;
						}
						case 3:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv - ju * (currentBS->degree()[1] + 1)];
								}
							}
							break;
						}
						case 4:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{
								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[i * currentBS->degree()[1] + jv][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1]) - jv * (currentBS->degree()[0] + 1) + ju];
								}
							}
							break;
						}
						default:
							break;
						}
					}
					/*Assign values to the last control point in the U direction for each one*/
					{
						H* u_he = u_hes.back();
						pMesh->halfedgePrev(u_he)->NURBS_patchId() = currentBS->id();
						pMesh->halfedgePrev(u_he)->isoParameterLineType() = 4;
						F* h_heF = pMesh->halfedgeFace(u_he);
						int h_heFId = h_heF->id();
						h_heFId--;
						std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
						switch (u_he->localId())
						{
						case 1:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0]) + jv];

							}
							break;
						}
						case 2:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[jv * (currentBS->degree()[0] + 1)];

							}
							break;
						}
						case 3:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[1]) - jv];

							}
							break;
						}
						case 4:
						{
							for (int jv = 0; jv < currentBS->degree()[1]; jv++)
							{

								currentBS->cpts()[i * currentBS->degree()[1] + jv].back() = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - jv * (currentBS->degree()[0] + 1)];

							}
							break;
						}
						default:
							break;
						}
					}

					//Move the last row of boundary points in the V direction into the "vertices" of the patch.
					if (i == vHes.size() - 1)
					{
						currentBS->cpts()[vHes.size() * currentBS->degree()[1]].resize(u_hes.size() * currentBS->degree()[0] + 1);
						/*Assign values to the control points*/
						for (int j = 0; j < u_hes.size(); j++)
						{
							H* u_he = u_hes[j];
							pMesh->halfedgeNext(pMesh->halfedgeNext(u_he))->NURBS_patchId() = currentBS->id();
							pMesh->halfedgeNext(pMesh->halfedgeNext(u_he))->isoParameterLineType() = 3;
							F* h_heF = pMesh->halfedgeFace(u_he);
							int h_heFId = h_heF->id();
							h_heFId--;
							std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
							switch (u_he->localId())
							{
							case 1:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[1] + 1) - 1 + (currentBS->degree()[1] + 1) * ju];
								}

								break;
							}
							case 2:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0] + 1) * (currentBS->degree()[1] + 1) - 1 - ju];
								}

								break;
							}
							case 3:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[1] + 1) - ju * (currentBS->degree()[1] + 1)];
								}

								break;
							}
							case 4:
							{

								for (int ju = 0; ju < currentBS->degree()[0]; ju++)
								{
									currentBS->cpts()[vHes.size() * currentBS->degree()[1]][j * currentBS->degree()[0] + ju] = h_heBf->cpts()[ju];
								}

								break;
							}
							default:
								break;
							}
						}
						/*Assign values to the last control point in the U direction for each one*/
						{
							H* u_he = u_hes.back();
							F* h_heF = pMesh->halfedgeFace(u_he);
							int h_heFId = h_heF->id();
							h_heFId--;
							std::shared_ptr<CCG_QMS_bezierSurf> h_heBf = this->bfs()[h_heFId];
							switch (u_he->localId())
							{
							case 1:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[1] + 1) * (currentBS->degree()[0] + 1) - 1];
								break;
							}
							case 2:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0]) * (currentBS->degree()[0] + 1)];
								break;
							}
							case 3:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[0];
								break;
							}
							case 4:
							{
								currentBS->cpts()[vHes.size() * currentBS->degree()[1]].back() = h_heBf->cpts()[(currentBS->degree()[0])];
								break;
							}
							default:
								break;
							}
						}
					}
				}
			}
		}
	}
	//std::cout << "Obtaining B spline surfaces finished!" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::removeMultiKnots_BSplines()
{
	for (auto p : bsplineSurfs())
	{
		//p->removeSurfaceMultiKnots();
		//p->RemoveCurveKnot();
		//p->RemoveCurveKnot_UVconsisi(1e-6);
		p->RemoveCurveKnot_UVConsistency2(1e-14);
	}
	//std::cout << "Removing knots successfully!" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::computeBSplineSurfaceUniqueKnotMultiNum()
{
	//int count = 0;
	for (auto p : bsplineSurfs())
	{
		p->computeUniqueKnotNum();
		//if (count == 20)
		//{ //p->removeSurfaceMultiKnots();
		//	std::cout << "-----patch num: " << count << std::endl;
		//	p->computeUniqueKnotNum();
		//	std::cout << "---U cpt num: " << p->cpts().size() << std::endl;
		//	for (auto pp : p->cpts())
		//	{
		//		std::cout << "---V cpt num: " << pp.size() << std::endl;
		//	}
		//}
		//count++;
	}
	//std::cout << "Computing Unique Knot Multi Num successfully!" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::outputBSpline_iges(const char* output)
{
	/*Time string, which requires including the header file <chrono>*/
	auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::stringstream ss;
	ss << std::put_time(std::localtime(&t), "%Y-%m-%d-%H-%M-%S");
	std::string dt = ss.str();

	/*Define the dimension: a curve is of dimension 1, and a surface is of dimension 2.*/
	int dim = 2;
	std::string filename = output;

	/*Retrieve the file name from a relative or absolute path*/
	std::size_t found = filename.find_last_of("/\\");
	std::string igesFilename = filename.substr(found + 1);
	/*To prevent the name from exceeding the character limit of one line, the file name should be uniformly written as "QUAD_IGES". */
	igesFilename = "QUAD_IGES";
	//std::cout << " path: " << filename.substr(0, found) << '\n';
	//std::cout << " file: " << filename.substr(found + 1) << '\n';


	// START SECTION
	std::string S[4];
	S[0] = "";
	S[1] = "IGES obtained from Nurbs toolbox.";
	S[2] = "See <http://octave.sourceforge.net/nurbs/>.";
	S[3] = "";

	// GLOBAL SECTION
	std::string G[24];
	G[0] = "1H,";                           // Parameter Deliminator Character
	G[1] = "1H;";                           //Record Delimiter Character
	G[2] = HString("Nurbs toolbox");       // Product ID from Sender
	G[3] = HString(igesFilename);              // File Name
	G[4] = HString("Quad Nurbs");        // System ID
	G[5] = HString("nrb2iges");            // Pre-processor Version
	G[6] = "32";                            // Number of Bits for Integers (No. of bits present in the integer representation of the sending system)
	G[7] = "75";                            // Single Precision Magnitude (Maximum power of 10 which may be represented as a single precision floating point number from the sending system)
	G[8] = "6";                             // Single Precision Significance (No. of significant digits of a single precision floating point number on the sending system)
	G[9] = "75";                            // Double Precision Magnitude (Maximum power of 10 which may be represented as a double precision floating point number from the sending system)
	G[10] = "15";                            // Double Precision Significance (No. of significant digits of a double precision floating point number on the sending system)
	G[11] = HString("Nurbs from QUAD");    // Product ID for Receiver
	G[12] = "1.0";                           // Model Space Scale
	G[13] = "6";                             // Unit Flag (1 = inches)(6=Meters)
	G[14] = HString("M");                    // Units  (inches = "IN/INCH")(Meters = "M")(Milimeter = "MM")
	G[15] = "1000";                          // Maximum Number of Line Weights
	G[16] = "1.0";                           // Size of Maximum Line Width
	G[17] = HString(dt);                     // Date and Time of file generation
	G[18] = "0.000001";                      // Minimum User-intended Resolution
	G[19] = "10000.0";                       // Approximate Maximum Coordinate
	G[20] = HString("cc ffww");         // Name of Author
	G[21] = HString("DLUT");  // Author's Organization
	G[22] = "3";                             // IGES Version Number (3 = IGES version 2.0)
	G[23] = "0";                             // Drafting Standard Code (0 = no standard)

	//Convert section array to lines(maximum lenght 72)
	std::vector<std::string> GG;
	for (int i = 0; i < sizeof(G) / sizeof(G[0]); i++)
	{
		GG.push_back(G[i]);
	}
	std::vector<std::string> SectionG = make_section(GG, 72);

	// DIRECTORY ENTRY SECTION
	// Each directory entry consists of two, 80 character, fixed formatted lines
	std::vector<DirectoryEntrySectionEle> D;
	for (int i = 0; i < this->bsplineSurfs().size(); i++)
	{
		DirectoryEntrySectionEle nD;
		nD.type = 128;
		nD.id = 2 * (i + 1) - 1;
		nD.p_start = 0;
		nD.p_count = 0;
		D.push_back(nD);
	}

	//PARAMETER DATA SECTION
	/*
	*   The structure is a free formatted data entry from columns 1 to 64.
	*	Each line of free formatted data consists of the entity type number
	*	followed by the parameter data describing the entity.
	*	Columns 65 to 72 are reserved for a parameter data index which is an
	*	odd number counter, right justified in the field, which begins at the
	*	number 1 and progresses in odd increments for each entity entered.
	*	Column 73 is reserved for the letter to indicate the data element
	*	belongs to the parameter data section.
	*	Columns 74 to 80 are reserved for the sequence number. Each line of
	*	data corresponds to the entity type as specified in the global section.
	*/
	std::vector<std::vector<std::string>> SectionP;
	for (int i = 0; i < this->bsplineSurfs().size(); i++)
	{
		auto nf = this->bsplineSurfs()[i];
		double USpanMin = nf->knotsU().front();
		double USpanMax = nf->knotsU().back();
		double VSpanMin = nf->knotsV().front();
		double VSpanMax = nf->knotsV().back();
		std::vector<std::string> P;
		make_section_array(nf, USpanMin, USpanMax, VSpanMin, VSpanMax, P);
		//Convert section array to lines
		std::vector<std::string> SP = make_section(P, 64);
		D[i].p_count = SP.size();
		if (i == 0)
			D[i].p_start = 1;
		else
		{
			D[i].p_start = D[i - 1].p_start + D[i - 1].p_count;
		}
		SectionP.push_back(SP);
	}

	//save
	std::fstream _os(output, std::fstream::out);
	//Save Start Section
	for (int i = 0; i < (sizeof(S) / sizeof(S[0])); i++)
	{
		_os << std::left << std::setw(72) << S[i] << std::right << "S" << std::right << std::setw(7) << i + 1 << std::endl;
	}
	//Save Global Section
	for (int i = 0; i < SectionG.size(); i++)
	{
		_os << std::left << std::setw(72) << SectionG[i] << std::right << "G" << std::right << std::setw(7) << i + 1 << std::endl;
	}
	//Save Directory Entry Section
	for (int i = 0; i < D.size(); i++)
	{
		_os << std::setw(8) << D[i].type << std::setw(8) << D[i].p_start << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << "D" << std::setw(7) << (i + 1) * 2 - 1 << std::endl;
		_os << std::setw(8) << D[i].type << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << D[i].p_count << std::setw(8) << 0 << std::setw(8) << " " << std::setw(8) << " " << std::setw(8) << " " << std::setw(8) << 0 << "D" << std::setw(7) << (i + 1) * 2 << std::endl;
	}

	//Save Parameter Data Section
	int lines_p = 0;
	for (int i = 0; i < SectionP.size(); i++)
	{
		std::vector<std::string> sec = SectionP[i];
		for (int j = 0; j < sec.size(); j++)
		{
			lines_p++;
			_os << std::left << std::setw(64) << sec[j] << " " << std::right << std::setw(7) << D[i].id << "P" << std::setw(7) << lines_p << std::endl;
		}
	}

	//Save Terminate Section
	std::string sec_t;
	std::string sec_t1 = std::to_string(sizeof(S) / sizeof(S[0])) + "S";
	for (int i = sec_t1.length(); i <= 7; i++)
	{
		sec_t1.insert(0, " ");
	}
	sec_t += sec_t1;
	sec_t1 = std::to_string(SectionG.size()) + "G";
	for (int i = sec_t1.length(); i <= 7; i++)
	{
		sec_t1.insert(0, " ");
	}
	sec_t += sec_t1;
	sec_t1 = std::to_string(2 * D.size()) + "D";
	for (int i = sec_t1.length(); i <= 7; i++)
	{
		sec_t1.insert(0, " ");
	}
	sec_t += sec_t1;
	sec_t1 = std::to_string(lines_p) + "P";
	for (int i = sec_t1.length(); i <= 7; i++)
	{
		sec_t1.insert(0, " ");
	}
	sec_t += sec_t1;
	for (int i = 0; i < 7; i++)
	{
		sec_t += " ";
	}
	_os << std::left << std::setw(72) << sec_t << "T" << std::right << std::setw(7) << 1 << std::endl;
	_os.close();
	//std::cout << "Ouput file " << output << " successfully!" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::outputBSpline_iges_5_3(const char* output)
{
	/*Time string, which requires including the header file <chrono>*/
	auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::stringstream ss;
	//ss << std::put_time(std::localtime(&t), "%Y-%m-%d-%H-%M-%S");
	ss << std::put_time(std::localtime(&t), "%Y%m%d.%H%M%S");
	std::string dt = ss.str();

	/*Define the dimension: a curve is of dimension 1, and a surface is of dimension 2.*/
	int dim = 2;
	std::string filename = output;

	/*Retrieve the file name from a relative or absolute path*/
	std::size_t found = filename.find_last_of("/\\");
	std::string igesFilename = filename.substr(found + 1);
	/*To prevent the name from exceeding the character limit of one line, the file name should be uniformly written as "QUAD_IGES". */
	igesFilename = "QUAD_IGES";
	//std::cout << " path: " << filename.substr(0, found) << '\n';
	//std::cout << " file: " << filename.substr(found + 1) << '\n';


	// START SECTION
	/*std::string S[4];
	S[0] = "";
	S[1] = "IGES obtained from Nurbs toolbox.";
	S[2] = "See <http://octave.sourceforge.net/nurbs/>.";
	S[3] = "";*/
	std::string S[3];
	S[0] = "IGES obtained from Nurbs toolbox.";
	S[1] = "See <http://octave.sourceforge.net/nurbs/>.";
	S[2] = "";

	// GLOBAL SECTION
	std::string G[24];
	G[0] = "1H,";                           // Parameter Deliminator Character
	G[1] = "1H;";                           //Record Delimiter Character
	G[2] = HString("Nurbs toolbox");       // Product ID from Sender
	G[3] = HString(igesFilename);              // File Name
	G[4] = HString("Quad Nurbs");        // System ID
	G[5] = HString("nrb2iges");            // Pre-processor Version
	G[6] = "32";                            // Number of Bits for Integers (No. of bits present in the integer representation of the sending system)
	G[7] = "38";                            // 
	G[8] = "6";                             // Single Precision Significance (No. of significant digits of a single precision floating point number on the sending system)
	G[9] = "308";                            //
	G[10] = "15";                            // Double Precision Significance (No. of significant digits of a double precision floating point number on the sending system)
	G[11] = HString("Nurbs from QUAD");    // Product ID for Receiver
	G[12] = "1.0";                           // Model Space Scale
	G[13] = "6";                             // Unit Flag (1 = inches)(6=Meters)
	G[14] = HString("M");                    // Units  (inches = "IN/INCH")(Meters = "M")(Milimeter = "MM")
	G[15] = "1000";                          // Maximum Number of Line Weights
	G[16] = "1.0";                           // Size of Maximum Line Width
	G[17] = HString(dt);                     // Date and Time of file generation
	G[18] = "0.000001";                      // Minimum User-intended Resolution
	G[19] = "10000.0";                       // Approximate Maximum Coordinate
	G[20] = HString("cc ffww");         // Name of Author
	G[21] = HString("DLUT");  // Author's Organization
	G[22] = "11";                             //
	G[23] = "0";                             // Drafting Standard Code (0 = no standard)

	//Convert section array to lines(maximum lenght 72)
	std::vector<std::string> GG;
	for (int i = 0; i < sizeof(G) / sizeof(G[0]); i++)
	{
		GG.push_back(G[i]);
	}
	std::vector<std::string> SectionG = make_section(GG, 72);

	// DIRECTORY ENTRY SECTION
	// Each directory entry consists of two, 80 character, fixed formatted lines
	std::vector<DirectoryEntrySectionEle> D;
	for (int i = 0; i < this->bsplineSurfs().size(); i++)
	{
		DirectoryEntrySectionEle nD;
		nD.type = 128;
		nD.id = 2 * (i + 1) - 1;
		nD.p_start = 0;
		nD.p_count = 0;
		D.push_back(nD);
	}

	//PARAMETER DATA SECTION
	/*
	*   The structure is a free formatted data entry from columns 1 to 64.
	*	Each line of free formatted data consists of the entity type number
	*	followed by the parameter data describing the entity.
	*	Columns 65 to 72 are reserved for a parameter data index which is an
	*	odd number counter, right justified in the field, which begins at the
	*	number 1 and progresses in odd increments for each entity entered.
	*	Column 73 is reserved for the letter to indicate the data element
	*	belongs to the parameter data section.
	*	Columns 74 to 80 are reserved for the sequence number. Each line of
	*	data corresponds to the entity type as specified in the global section.
	*/
	std::vector<std::vector<std::string>> SectionP;
	for (int i = 0; i < this->bsplineSurfs().size(); i++)
	{
		auto nf = this->bsplineSurfs()[i];
		double USpanMin = nf->knotsU().front();
		double USpanMax = nf->knotsU().back();
		double VSpanMin = nf->knotsV().front();
		double VSpanMax = nf->knotsV().back();
		std::vector<std::string> P;
		make_section_array(nf, USpanMin, USpanMax, VSpanMin, VSpanMax, P);
		//Convert section array to lines
		std::vector<std::string> SP = make_section(P, 64);
		D[i].p_count = SP.size();
		if (i == 0)
			D[i].p_start = 1;
		else
		{
			D[i].p_start = D[i - 1].p_start + D[i - 1].p_count;
		}
		SectionP.push_back(SP);
	}

	//save
	std::fstream _os(output, std::fstream::out);
	//Save Start Section
	for (int i = 0; i < (sizeof(S) / sizeof(S[0])); i++)
	{
		_os << std::left << std::setw(72) << S[i] << std::right << "S" << std::right << std::setw(7) << i + 1 << std::endl;
	}
	//Save Global Section
	for (int i = 0; i < SectionG.size(); i++)
	{
		_os << std::left << std::setw(72) << SectionG[i] << std::right << "G" << std::right << std::setw(7) << i + 1 << std::endl;
	}
	//Save Directory Entry Section
	for (int i = 0; i < D.size(); i++)
	{
		_os << std::setw(8) << D[i].type << std::setw(8) << D[i].p_start << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << 0 << "D" << std::setw(7) << (i + 1) * 2 - 1 << std::endl;
		_os << std::setw(8) << D[i].type << std::setw(8) << 0 << std::setw(8) << 0 << std::setw(8) << D[i].p_count << std::setw(8) << 0 << std::setw(8) << " " << std::setw(8) << " " << std::setw(8) << "TrimSrf" << std::setw(8) << 0 << "D" << std::setw(7) << (i + 1) * 2 << std::endl;
	}

	//Save Parameter Data Section
	int lines_p = 0;
	for (int i = 0; i < SectionP.size(); i++)
	{
		std::vector<std::string> sec = SectionP[i];
		for (int j = 0; j < sec.size(); j++)
		{
			lines_p++;
			_os << std::left << std::setw(64) << sec[j] << " " << std::right << std::setw(7) << D[i].id << "P" << std::setw(7) << lines_p << std::endl;
		}
	}

	//Save Terminate Section
	std::string sec_t;
	std::string sec_t1 = "S";
	std::string sec_t1_1 = std::to_string(sizeof(S) / sizeof(S[0]));
	for (int i = sec_t1.length(); i <= (7- sec_t1_1.length()); i++)
	{
		sec_t1.insert(1, "0");
	}
	sec_t += sec_t1;
	sec_t += sec_t1_1;
	sec_t1 = "G";
	sec_t1_1 = std::to_string(SectionG.size());
	for (int i = sec_t1.length(); i <= (7 - sec_t1_1.length()); i++)
	{
		sec_t1.insert(1, "0");
	}
	sec_t += sec_t1;
	sec_t += sec_t1_1;
	sec_t1 = "D";
	sec_t1_1 = std::to_string(2 * D.size());
	for (int i = sec_t1.length(); i <= (7 - sec_t1_1.length()); i++)
	{
		sec_t1.insert(1, "0");
	}
	sec_t += sec_t1;
	sec_t += sec_t1_1;
	sec_t1 = "P";
	sec_t1_1 = std::to_string(lines_p);
	for (int i = sec_t1.length(); i <= (7 - sec_t1_1.length()); i++)
	{
		sec_t1.insert(1, "0");
	}
	sec_t += sec_t1;
	sec_t += sec_t1_1;
	for (int i = 0; i < 7; i++)
	{
		sec_t += " ";
	}
	_os << std::left << std::setw(72) << sec_t << "T" << std::right << std::setw(7) << 1 << std::endl;
	_os.close();
	//std::cout << "Ouput file " << output << " successfully!" << std::endl;
}

void CCG_QMSLib::CCG_QMS_model::outputBSplineSurf_fixedFormat3(const char* output)
{
	std::fstream _os(output, std::fstream::out);
	if (_os.fail())
	{
		fprintf(stderr, "Error is opening file %s\n", output);
		return;
	}

	/*Traverse all the B-spline surface patches in the model*/
	_os << "TotalPatch " << this->bsplineSurfs().size() << std::endl;
	for (int i = 0; i < this->bsplineSurfs().size(); i++)
	{
		_os << "Patch " << i + 1 << std::endl;
		auto bf = this->bsplineSurfs()[i];
		/*Number_of_unique_knots_u (int)*/
		int Number_of_unique_knots_u = bf->knotsU_unique_num().size();
		_os << Number_of_unique_knots_u << ";" << std::endl;
		/*Number_of_unique_knots_v (int)*/
		int Number_of_unique_knots_v = bf->knotsV_unique_num().size();
		_os << Number_of_unique_knots_v << ";" << std::endl;
		/*Number_of_ctrl_pts_u (int)*/
		int Number_of_ctrl_pts_u = bf->knotsU().size() - 1 - bf->degree()[0];
		_os << Number_of_ctrl_pts_u << ";" << std::endl;
		/*Number_of_ctrl_pts_v (int)*/
		int Number_of_ctrl_pts_v = bf->knotsV().size() - 1 - bf->degree()[1];
		_os << Number_of_ctrl_pts_v << ";" << std::endl;
		/*Unique_knots_u vector (doubles)*/
		int tempCount = 0;
		for (auto kuu : bf->knotsU_unique_num())
		{
			tempCount++;
			_os << std::fixed;
			_os << std::setprecision(3);
			if (tempCount == bf->knotsU_unique_num().size())
			{
				_os << kuu.first << ";" << std::endl;
			}
			else
			{
				_os << kuu.first << ",";
			}
		}
		/*Unique_knots_v vector (doubles)*/
		tempCount = 0;
		for (auto kvu : bf->knotsV_unique_num())
		{
			tempCount++;
			_os << std::fixed;
			_os << std::setprecision(3);
			if (tempCount == bf->knotsV_unique_num().size())
			{
				_os << kvu.first << ";" << std::endl;
			}
			else
			{
				_os << kvu.first << ",";
			}
		}
		/*Multiplicity_knots_u vector (ints)*/
		tempCount = 0;
		for (auto kun : bf->knotsU_unique_num())
		{
			tempCount++;
			_os << std::fixed;
			_os << std::setprecision(3);
			if (tempCount == bf->knotsU_unique_num().size())
			{
				_os << kun.second << ";" << std::endl;
			}
			else
			{
				_os << kun.second << ",";
			}
		}
		/*Multiplicity_knots_v vector (ints)*/
		tempCount = 0;
		for (auto kvn : bf->knotsV_unique_num())
		{
			tempCount++;
			_os << std::fixed;
			_os << std::setprecision(3);
			if (tempCount == bf->knotsV_unique_num().size())
			{
				_os << kvn.second << ";" << std::endl;
			}
			else
			{
				_os << kvn.second << ",";
			}
		}
		/*Periodic along U? (1/0 means True/False)*/
		if (bf->periodicAlongU())
		{
			_os << 1 << std::endl;
		}
		else
		{
			_os << 0 << std::endl;
		}
		/*Periodic along V? (1/0 means True/False)*/
		if (bf->periodicAlongV())
		{
			_os << 1 << std::endl;
		}
		else
		{
			_os << 0 << std::endl;
		}
		/* Closed along U? (1/0 means True/False)*/
		if (bf->closedAlongU())
		{
			_os << 1 << std::endl;
		}
		else
		{
			_os << 0 << std::endl;
		}
		/* Closed along V? (1/0 means True/False)*/
		if (bf->closedAlongV())
		{
			_os << 1 << std::endl;
		}
		else
		{
			_os << 0 << std::endl;
		}

		/*control pts*/
		/*for (auto cpts : bf->cpts())
		{
			for (auto cpt : cpts)
			{
				_os << std::fixed;
				_os << std::setprecision(15);
				_os << cpt[0] << "," << cpt[1] << "," << cpt[2] << ";" << std::endl;
			}
		}*/
		for (int i = 0; i < bf->cpts()[0].size(); i++)
		{
			for (int j = 0; j < bf->cpts().size(); j++)
			{
				_os << std::fixed;
				_os << std::setprecision(15);
				_os << bf->cpts()[j][i][0] << "," << bf->cpts()[j][i][1] << "," << bf->cpts()[j][i][2] << ";" << std::endl;
			}
		}
	}
	_os.close();
}

void CCG_QMSLib::CCG_QMS_model::obtainNURBS_boundaryCurve_output3(M* triMesh, M* quadMesh, const char* input, const char* output)
{
	//struct NURBS_boundaryCurve
	//{
	//	int patchId;//nurbs patch id
	//	int curveType;//( 1: u=u_min; 2: v = v_min; 3: u = u_max; 4:  v=v_max )
	//	std::vector<int> originalTriMeshVIds;//record the original tri mesh vetex id
	//	std::vector<int> triMeshVIds;//record the tri mesh vetex id
	//};

	/*store all nurbs boundary curve*/
	//std::vector<NURBS_boundaryCurve> nbcs;

	/*Mark boundary halfedges isoParameterLineType and NURBS_patchId*/
	for (auto bsurf : this->bsplineSurfs())
	{
		for (auto bsurfBoundary : bsurf->boundarys())
		{
			if (bsurfBoundary->boundary())
			{
				for (auto boundaryH : bsurfBoundary->halfedges())
				{
					boundaryH->isoParameterLineType() = bsurfBoundary->localId();
					boundaryH->NURBS_patchId() = bsurfBoundary->bsurfId();
				}
			}
		}
	}

	/*obtaining all nurbs boundary curve by quad mesh boudnary halfedges*/
	for (auto e : It::MEIterator(quadMesh))
	{
		if (e->boundary())
		{
			H* eh = quadMesh->edgeHalfedge(e, 0);
			bool mark = false;
			if (this->nbcs().size() == 0)
			{
				mark = false;
			}
			else
			{
				for (int i = 0; i < this->nbcs().size(); i++)
				{
					if (this->nbcs()[i].curveType == eh->isoParameterLineType() && this->nbcs()[i].patchId == eh->NURBS_patchId())
					{
						mark = true;
						break;
					}
				}
			}
			if (!mark)
			{
				NURBS_boundaryCurve temp_nbc;
				temp_nbc.curveType = eh->isoParameterLineType();
				temp_nbc.patchId = eh->NURBS_patchId();
				//std::cout << "temp_nbc.curveType: " << temp_nbc.curveType << " temp_nbc.patchId: " << temp_nbc.patchId << std::endl;
				this->nbcs().push_back(temp_nbc);
			}
		}
	}
	//std::cout << "# nurbs boundary curve: " << this->nbcs().size() << std::endl;

	/*obtaing all vertices from tri Mesh whose original id is valid*/
	std::vector<V*> triMeshBoundaryVs;
	for (auto v : It::MVIterator(triMesh))
	{
		if (v->originalId() != -1)
		{
			triMeshBoundaryVs.push_back(v);
			v->ifVisit() = false;
		}
	}
	//std::cout << "triMeshBoundaryVs.size(): " << triMeshBoundaryVs.size() << std::endl;
	/*obtaing all halfedges whose original id is valid*/
	std::vector<H*> triMeshBoundaryHEs;
	for (auto f : It::MFIterator(triMesh))
	{
		for (auto fh : It::FHEIterator(triMesh, f))
		{
			if (triMesh->halfedgeSource(fh)->originalId() != -1 || triMesh->halfedgeTarget(fh)->originalId() != -1)
			{
				E* fhE = triMesh->halfedgeEdge(fh);
				//fhE->mark() = true;
				triMeshBoundaryHEs.push_back(fh);
			}
		}

	}
	//std::cout << "triMeshBoundaryHEs.size(): " << triMeshBoundaryHEs.size() << std::endl;
	/*obtaining the originalTriMeshVId and triMeshVId*/
	for (int i = 0; i < this->nbcs().size(); i++)
	{
		/*iterator all halfedges from tri Mesh whose original id is valid*/
		for (auto triMeshBoundaryHE : triMeshBoundaryHEs)
		{
			if (triMeshBoundaryHE->boundaryFaceId() == -1)continue;
			F* bvQuadF = quadMesh->idFace(triMeshBoundaryHE->boundaryFaceId());
			for (auto fh : It::FHEIterator(quadMesh, bvQuadF))
			{
				/*if (fh->target()->id() == bv->boundaryTargetVId())
				{
					std::cout << "fh->NURBS_patchId(): " << fh->NURBS_patchId()<< " fh->isoParameterLineType(): "<< fh->isoParameterLineType() << std::endl;
					std::cout << "this->nbcs()[i].patchId: " << this->nbcs()[i].patchId << " this->nbcs()[i].curveType: " << this->nbcs()[i].curveType << std::endl;
				}*/
				if (fh->target()->id() == triMeshBoundaryHE->boundaryTargetVId() && fh->NURBS_patchId() == this->nbcs()[i].patchId && fh->isoParameterLineType() == this->nbcs()[i].curveType)
				{
					if (triMesh->halfedgeSource(triMeshBoundaryHE)->originalId() != -1)
					{
						fh->originalIds().push_back(triMesh->halfedgeSource(triMeshBoundaryHE)->originalId());
						triMesh->halfedgeSource(triMeshBoundaryHE)->ifVisit() = true;
						triMesh->halfedgeSource(triMeshBoundaryHE)->boundaryFaceId() = triMeshBoundaryHE->boundaryFaceId();
						triMesh->halfedgeSource(triMeshBoundaryHE)->boundaryTargetVId() = triMeshBoundaryHE->boundaryTargetVId();
					}
					if (triMesh->halfedgeTarget(triMeshBoundaryHE)->originalId() != -1)
					{
						fh->originalIds().push_back(triMesh->halfedgeTarget(triMeshBoundaryHE)->originalId());
						triMesh->halfedgeTarget(triMeshBoundaryHE)->ifVisit() = true;
						triMesh->halfedgeTarget(triMeshBoundaryHE)->boundaryFaceId() = triMeshBoundaryHE->boundaryFaceId();
						triMesh->halfedgeTarget(triMeshBoundaryHE)->boundaryTargetVId() = triMeshBoundaryHE->boundaryTargetVId();
					}
				}
				//std::cout << "fh->NURBS_patchId(): " << fh->NURBS_patchId() << " fh->isoParameterLineType(): " << fh->isoParameterLineType()<<" fh->originalIds().size(): "<< fh->originalIds().size() << std::endl;
			}
		}

		/*Iterator all vertices from tri Mesh whose original id is valid*/
		for (auto bv : triMeshBoundaryVs)
		{
			if (bv->ifVisit())
			{
				//std::cout << " ID: " << bv->id() << " count: " << bv->count() << std::endl;
				this->nbcs()[i].triMeshVIds.push_back(bv->id());
				this->nbcs()[i].originalTriMeshVIds.push_back(bv->originalId());
			}
		}
		/*initialize visit state to false*/
		for (auto bv : triMeshBoundaryVs)
		{
			bv->ifVisit() = false;
		}
	}

	//obtaining map between simplified IDs and initial IDs
	//std::map<int, std::vector<int>> vertexMap;
	{
		std::ifstream file(input);
		std::string line;
		int currentVertex = -1;

		if (!file.is_open()) {
			std::cerr << "Error opening file: " << input << std::endl;
			return;
		}

		while (std::getline(file, line)) {
			// Skip empty lines
			if (line.empty()) continue;

			// Check if line starts with "Vertex"
			if (line.find("Vertex") == 0) {
				// Extract vertex number
				std::istringstream iss(line);
				std::string token;
				iss >> token >> currentVertex;
			}
			else {
				// Process numbers under current vertex
				if (currentVertex != -1) {
					std::istringstream iss(line);
					int number;
					while (iss >> number) {
						vertexMap()[currentVertex].push_back(number);
					}
				}
			}
		}

		file.close();
	}

	/*check result of obtaining map between simplified IDs and initial IDs*/
	//for (auto v : It::MVIterator(objTriMesh))
	//{
	//	v->ifSharp() = false;
	//}

	//for (auto vp : vertexMap())
	//{
	//	std::cout << "vp.first: " << vp.first << std::endl;
	//	objTriMesh->idVertex(vp.first)->ifSharp() = true;
	//	for (auto tempId : vp.second)
	//	{
	//		if (!objTriMesh->idVertex(tempId)->boundary())
	//		{
	//			std::cout << "VId: " << tempId << " is not on the obj Mesh boundary! " << std::endl;
	//		}
	//		//std::cout << tempId << std::endl;
	//	}
	//	for (auto tempId : vp.second)
	//	{
	//		//std::cout << tempId << std::endl;
	//		objTriMesh->idVertex(tempId)->ifVisit() = true;
	//	}
	//}

	/*obtaing all boudnary halfedges from quadMesh*/
	std::vector<H*> quadMeshBoundaryHEs;
	for (auto e : It::MEIterator(quadMesh))
	{
		if (e->boundary())
		{
			H* eh0 = quadMesh->edgeHalfedge(e, 0);
			quadMeshBoundaryHEs.push_back(eh0);
		}
	}


	/*output all nurbs boundary curve*/
	std::fstream _os(output, std::fstream::out);
	if (_os.fail())
	{
		fprintf(stderr, "Error is opening file %s\n", output);
		return;
	}

	/*statistic valid NURBS curve number*/
	int validNURBSCurveSize = 0;
	for (int k = 0; k < this->nbcs().size(); k++)
	{
		NURBS_boundaryCurve nbc = this->nbcs()[k];
		bool tempMark = false;
		if (nbc.originalTriMeshVIds.size() == 0)
		{
			tempMark = true;
			/* find a boundary halfedges from quadMesh locating this nurbs curve*/
			H* ncHe = NULL;
			for (auto quadMeshBoundaryHE : quadMeshBoundaryHEs)
			{
				if (quadMeshBoundaryHE->NURBS_patchId() == nbc.patchId && quadMeshBoundaryHE->isoParameterLineType() == nbc.curveType)
				{
					ncHe = quadMeshBoundaryHE;
				}
			}
			//std::cout << "-----------------------------------------------------------------------------" << std::endl;
			//std::cout << " nbc.patchId(): " << nbc.patchId << " nbc.curveType(): " << nbc.curveType << std::endl;
			/*if (ncHe != NULL)
			{
				std::cout << "ncHe->NURBS_patchId(): " << ncHe->NURBS_patchId() << " ncHe->isoParameterLineType(): " << ncHe->isoParameterLineType() << std::endl;
			}
			std::cout << "--------------------------------------" << std::endl;*/

			/*find the nearest boundary halfedges from previous direction of ncHe, satisfy whose orinigalIds.size()!=0*/
			H* ncHePre = NULL;
			H* tempNCHe = ncHe;
			while (ncHePre == NULL)
			{
				V* ncHeSourceV = quadMesh->halfedgeSource(tempNCHe);
				for (auto ncHeSourceVHe : It::VClwInHEIterator(quadMesh, ncHeSourceV))
				{
					if (quadMesh->halfedgeEdge(ncHeSourceVHe)->boundary())
					{
						tempNCHe = ncHeSourceVHe;
						if (tempNCHe->NURBS_patchId() != ncHe->NURBS_patchId() && tempNCHe->originalIds().size() > 1)
						{
							ncHePre = tempNCHe;
							break;
						}
					}
				}
			}
			/*if (ncHePre != NULL)
			{
				std::cout << " ncHePre->NURBS_patchId(): " << ncHePre->NURBS_patchId() << " ncHePre->isoParameterLineType(): " << ncHePre->isoParameterLineType()<<" ncHePre->originalIds().size(): "<<ncHePre->originalIds().size() << std::endl;

			}
			std::cout << "--------------------------------------" << std::endl;*/



			/*find the nearest boundary halfedges from next direction of ncHe, satisfy whose orinigalIds.size()!=0*/
			H* ncHeNext = NULL;
			tempNCHe = ncHe;
			while (ncHeNext == NULL)
			{
				V* ncHeTargetV = quadMesh->halfedgeTarget(tempNCHe);
				for (auto ncHeTargetVHe : It::VCcwOutHEIterator(quadMesh, ncHeTargetV))
				{
					if (quadMesh->halfedgeEdge(ncHeTargetVHe)->boundary())
					{
						tempNCHe = ncHeTargetVHe;
						if (tempNCHe->NURBS_patchId() != ncHe->NURBS_patchId() && tempNCHe->originalIds().size() > 1)
						{
							ncHeNext = tempNCHe;
							break;
						}
					}
				}
			}
			/*if (ncHeNext != NULL)
			{
				std::cout << " ncHeNext->NURBS_patchId(): " << ncHeNext->NURBS_patchId() << " ncHeNext->isoParameterLineType(): " << ncHeNext->isoParameterLineType() << " ncHeNext->originalIds().size(): " << ncHeNext->originalIds().size() << std::endl;

			}
			std::cout << "--------------------------------------" << std::endl;*/

			//std::cout << "Before: nbc.originalTriMeshVIds.size(): " << nbc.originalTriMeshVIds.size() << std::endl;
			/*choose smaller originalIds*/
			if (ncHeNext->originalIds().size() > ncHePre->originalIds().size())
			{
				for (auto originalId : ncHePre->originalIds())
				{
					nbc.originalTriMeshVIds.push_back(originalId);
				}
			}
			else
			{
				for (auto originalId : ncHeNext->originalIds())
				{
					nbc.originalTriMeshVIds.push_back(originalId);
				}
			}
			this->nbcs()[k] = nbc;
			//std::cout << "After: nbc.originalTriMeshVIds.size(): " << nbc.originalTriMeshVIds.size() << std::endl;
			//std::cout << "-----------------------------------------------------------------------------" << std::endl;
			//continue;
		}
		for (int i = 0; i < nbc.originalTriMeshVIds.size(); i++)
		{
			if (vertexMap().find(nbc.originalTriMeshVIds[i]) == vertexMap().end())
			{
				//std::cout << "------------------------------------" << std::endl;
				tempMark = true;
				continue;
			}
			tempMark = false;
			break;
		}

		if (!tempMark)
		{
			validNURBSCurveSize++;
		}
		/*else
		{
			std::cout << "Test check: \n"
				<< "NURBS_BoundaryCurve[" << k + 1 << "]: Patch " << nbc.patchId << " , Side " << nbc.curveType << "\n  originalTriMeshVIds.size"
				<< nbc.originalTriMeshVIds.size()
				<< std::endl;
		}*/
	}

	_os << "Total_NURBS_BoundaryCurve_number " << validNURBSCurveSize << std::endl;

	int validNURBSCurveCount = 0;
	for (int k = 0; k < this->nbcs().size(); k++)
	{
		NURBS_boundaryCurve nbc = this->nbcs()[k];

		bool tempMark = false;
		if (nbc.originalTriMeshVIds.size() == 0)
		{
			tempMark = true;
		}
		for (int i = 0; i < nbc.originalTriMeshVIds.size(); i++)
		{
			if (vertexMap().find(nbc.originalTriMeshVIds[i]) == vertexMap().end())
			{
				//std::cout << "------------------------------------" << std::endl;
				tempMark = true;
				continue;
			}
			tempMark = false;
			break;
		}

		if (tempMark) continue;

		if (!tempMark)
		{
			validNURBSCurveCount++;
		}

		/*statistic all valid originalTriMeshVIds.size */
		int valid_originalTriMeshVIds_size = 0;
		for (int i = 0; i < nbc.originalTriMeshVIds.size(); i++)
		{
			if (vertexMap().find(nbc.originalTriMeshVIds[i]) == vertexMap().end())
			{
				//std::cout << "------------------------------------" << std::endl;
				continue;
			}
			valid_originalTriMeshVIds_size++;
		}

		_os << "NURBS_BoundaryCurve[" << validNURBSCurveCount << "]: Patch " << nbc.patchId << ", Side " << nbc.curveType;

		/*if (valid_originalTriMeshVIds_size < nbc.originalTriMeshVIds.size())
		{
			std::cout << "Test check: \n"
				<< "NURBS_BoundaryCurve[" << k + 1 << "]: Patch " << nbc.patchId << " , Side " << nbc.curveType <<"\n  originalTriMeshVIds.size"
				<< nbc.originalTriMeshVIds.size() << ". (BCV# = " << valid_originalTriMeshVIds_size << ")"
				<< std::endl;

		}*/
		_os << ". (BCV# = " << valid_originalTriMeshVIds_size << ")" << std::endl;

		int countvalid_originalTriMeshVIds = 0;
		for (int i = 0; i < nbc.originalTriMeshVIds.size(); i++)
		{
			/*_os << std::fixed;
			_os << std::setprecision(15);
			_os << "Vertex " << nbc.triMeshVIds[i] << " " << triMesh->idVertex(nbc.triMeshVIds[i])->point()[0]
				<< " " << triMesh->idVertex(nbc.triMeshVIds[i])->point()[1]
				<< " " << triMesh->idVertex(nbc.triMeshVIds[i])->point()[2]
				<< " " << nbc.originalTriMeshVIds[i] << std::endl;*/
				/*std::cout<< nbc.originalTriMeshVIds[i] << ": ";
				for (int j = 0; j < vertexMap.at(nbc.originalTriMeshVIds[i]).size() - 1; j++)
				{
					std::cout << vertexMap.at(nbc.originalTriMeshVIds[i])[j] << ", ";
				}
				std::cout << vertexMap.at(nbc.originalTriMeshVIds[i]).back() << ";"<<std::endl;*/

			if (vertexMap().find(nbc.originalTriMeshVIds[i]) == vertexMap().end())
			{
				//std::cout << "------------------------------------" << std::endl;
				continue;
			}
			countvalid_originalTriMeshVIds++;
			_os << "  BCV[" << countvalid_originalTriMeshVIds << "] IMBV " << nbc.originalTriMeshVIds[i] << " (Collapsed# = " << vertexMap().at(nbc.originalTriMeshVIds[i]).size() << "):" << std::endl;
			_os << "        ";
			for (int j = 0; j < vertexMap().at(nbc.originalTriMeshVIds[i]).size() - 1; j++)
			{
				_os << vertexMap().at(nbc.originalTriMeshVIds[i])[j] << ", ";
				//objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i])[j])->ifVisit() = true;
				//objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i])[j])->count()++;
				/*if (objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i])[j])->count() > 1)
				{
					std::cout << "boundary common vertex id: " << objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i])[j])->id()<< " count: " <<
						objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i])[j])->count() <<std::endl;
				}*/

			}
			_os << vertexMap().at(nbc.originalTriMeshVIds[i]).back() << ";" << std::endl;
			//objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i]).back())->ifVisit() = true;
			//objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i]).back())->count()++;
			/*if(objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i]).back())->count()>1)
			{
				std::cout << "boundary common vertex id: " << objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i]).back())->id()<<" count: "<<
					objTriMesh->idVertex(vertexMap().at(nbc.originalTriMeshVIds[i]).back())->count() << std::endl;
			}*/
		}
	}
	_os.close();

}

void CCG_QMSLib::CCG_QMS_model::obtainNURBS_boundaryCurve(M* quadMesh)
{
	/*obtaining all nurbs boundary curve by quad mesh boudnary halfedges*/
	for (auto e : It::MEIterator(quadMesh))
	{
		if (e->boundary())
		{
			H* eh = quadMesh->edgeHalfedge(e, 0);
			bool mark = false;
			if (this->nbcs().size() == 0)
			{
				mark = false;
			}
			else
			{
				for (int i = 0; i < this->nbcs().size(); i++)
				{
					if (this->nbcs()[i].curveType == eh->isoParameterLineType() && this->nbcs()[i].patchId == eh->NURBS_patchId())
					{
						mark = true;
						break;
					}
				}
			}
			if (!mark)
			{
				NURBS_boundaryCurve temp_nbc;
				temp_nbc.curveType = eh->isoParameterLineType();
				temp_nbc.patchId = eh->NURBS_patchId();
				//std::cout << "temp_nbc.curveType: " << temp_nbc.curveType << " temp_nbc.patchId: " << temp_nbc.patchId << std::endl;
				this->nbcs().push_back(temp_nbc);
			}
		}
	}
	//std::cout << "# nurbs boundary curve: " << this->nbcs().size() << std::endl;
}

double CCG_QMSLib::CCG_QMS_model::computeHausdorffDistance_allNURBSBoundaryCurve(M* objTriMesh, double sampleStep, std::vector<std::pair<int, double>>& sortedDistances)
{
	std::cout << "=== start computing Hausdorff distance ===" << std::endl;
	std::cout << " # NURBS boundary curve: " << nbcs().size() << std::endl;
	std::cout << "sampling step length: " << sampleStep << std::endl;

	//clear and initialize the sampling points vector
	m_boundarySamples.clear();
	m_boundarySamples.resize(nbcs().size());

	int totalSamples = 0; // recoard the total number of sampling points

	// Step 1 & 2: obtaining sampling points per nurbs boundary curve
	std::cout << "start sampling on every nurbs boundary curve..." << std::endl;
	for (size_t curveIdx = 0; curveIdx < nbcs().size(); ++curveIdx) {
		std::cout << "process curve " << curveIdx + 1 << "/" << nbcs().size() << std::endl;

		const NURBS_boundaryCurve& curve = nbcs()[curveIdx];
		int patchId = curve.patchId;
		int curveType = curve.curveType;

		std::cout << "  PatchID: " << patchId << ", curve type: " << curveType << std::endl;

		// 获取对应的B样条曲面
		std::shared_ptr<CCG_QMS_bsplineSurf> bsurf = nullptr;
		for (auto& surf : bsplineSurfs()) {
			if (surf->id() == patchId) {
				bsurf = surf;
				break;
			}
		}

		if (!bsurf) {
			std::cerr << "  无法找到曲面ID: " << patchId << "，跳过该曲线" << std::endl;
			continue;
		}

		// 根据边界曲线类型确定采样参数范围
		double fixedParam = 0.0; // 固定的参数值
		double startParam = 0.0; // 可变参数的起始值
		double endParam = 1.0;   // 可变参数的结束值
		bool isUFixed = false;   // 标记是否u参数固定

		switch (curveType) {
		case 1: // u = u_min
			isUFixed = true;
			fixedParam = bsurf->knotsU()[bsurf->degree()[0]]; // u_min
			startParam = bsurf->knotsV()[bsurf->degree()[1]]; // v_min
			endParam = bsurf->knotsV()[bsurf->knotsV().size() - bsurf->degree()[1] - 1]; // v_max
			break;
		case 2: // v = v_min
			isUFixed = false;
			fixedParam = bsurf->knotsV()[bsurf->degree()[1]]; // v_min
			startParam = bsurf->knotsU()[bsurf->degree()[0]]; // u_min
			endParam = bsurf->knotsU()[bsurf->knotsU().size() - bsurf->degree()[0] - 1]; // u_max
			break;
		case 3: // u = u_max
			isUFixed = true;
			fixedParam = bsurf->knotsU()[bsurf->knotsU().size() - bsurf->degree()[0] - 1]; // u_max
			startParam = bsurf->knotsV()[bsurf->degree()[1]]; // v_min
			endParam = bsurf->knotsV()[bsurf->knotsV().size() - bsurf->degree()[1] - 1]; // v_max
			break;
		case 4: // v = v_max
			isUFixed = false;
			fixedParam = bsurf->knotsV()[bsurf->knotsV().size() - bsurf->degree()[1] - 1]; // v_max
			startParam = bsurf->knotsU()[bsurf->degree()[0]]; // u_min
			endParam = bsurf->knotsU()[bsurf->knotsU().size() - bsurf->degree()[0] - 1]; // u_max
			break;
		default:
			std::cerr << "  未知的曲线类型: " << curveType << "，跳过该曲线" << std::endl;
			continue;
		}

		std::cout << "  参数范围: " << (isUFixed ? "u" : "v") << "=" << fixedParam
			<< ", " << (isUFixed ? "v" : "u") << " ∈ [" << startParam << ", " << endParam << "]" << std::endl;

		// 检查参数范围是否合理
		if (startParam >= endParam) {
			std::cerr << "  错误: 参数范围无效，跳过该曲线" << std::endl;
			continue;
		}

		// 计算采样点数量，避免过多采样
		int numSamples = static_cast<int>((endParam - startParam) / sampleStep) + 1;

		std::cout << "  开始采样，预计采样点数量: " << numSamples << std::endl;
		int actualSamples = 0;

		// 在曲线上采样点
		for (double t = startParam; t <= endParam; t += sampleStep) {
			BoundaryCurveSample sample;
			if (isUFixed) {
				sample.uvParam = CPoint2(fixedParam, t);
			}
			else {
				sample.uvParam = CPoint2(t, fixedParam);
			}

			try {
				// 使用DeBoorAlgorithmSurface计算采样点坐标
				sample.point = bsurf->DeBoorAlgorithmSurface(
					bsurf->cpts(),
					bsurf->knotsU(),
					bsurf->knotsV(),
					bsurf->degree()[0],
					bsurf->degree()[1],
					sample.uvParam[0],
					sample.uvParam[1]
				);
				sample.curveIndex = curveIdx;
				m_boundarySamples[curveIdx].push_back(sample);
				actualSamples++;
			}
			catch (const std::exception& e) {
				std::cerr << "  计算采样点时出错: " << e.what() << std::endl;
				continue;
			}
		}

		std::cout << "  实际采样点数量: " << actualSamples << std::endl;
		totalSamples += actualSamples;

		// 确保端点被采样到
		if (actualSamples > 0 && std::abs(m_boundarySamples[curveIdx].back().uvParam[isUFixed ? 1 : 0] - endParam) > 1e-6) {
			try {
				BoundaryCurveSample sample;
				if (isUFixed) {
					sample.uvParam = CPoint2(fixedParam, endParam);
				}
				else {
					sample.uvParam = CPoint2(endParam, fixedParam);
				}
				// 使用DeBoorAlgorithmSurface计算端点坐标
				sample.point = bsurf->DeBoorAlgorithmSurface(
					bsurf->cpts(),
					bsurf->knotsU(),
					bsurf->knotsV(),
					bsurf->degree()[0],
					bsurf->degree()[1],
					sample.uvParam[0],
					sample.uvParam[1]
				);
				sample.curveIndex = curveIdx;
				m_boundarySamples[curveIdx].push_back(sample);
				actualSamples++;
				totalSamples++;
			}
			catch (const std::exception& e) {
				std::cerr << "  计算端点采样点时出错: " << e.what() << std::endl;
			}
		}
	}

	std::cout << "采样完成，总采样点数量: " << totalSamples << std::endl;

	// Step 3: 计算每个三角网格边界点到所有NURBS边界曲线的最小距离
	std::cout << "开始计算三角网格边界点到NURBS边界曲线的最小距离..." << std::endl;

	// 计算三角网格边界点数量
	int boundaryPointCount = 0;
	for (auto v :It::MVIterator(objTriMesh)) {
		if (v->boundary()) boundaryPointCount++;
	}
	std::cout << "三角网格边界点数量: " << boundaryPointCount << std::endl;

	double maxMinDistance = 0.0; // 最大的最小距离（Hausdorff距离）
	std::vector<std::pair<int, double>> pointDistances;
	int processedPoints = 0;
	int matchedPoints = 0;

	for (auto v :It::MVIterator(objTriMesh)) {
		if (!v->boundary()) continue; // 只考虑边界点

		processedPoints++;
		if (processedPoints % 1000 == 0) {
			std::cout << "已处理 " << processedPoints << "/" << boundaryPointCount << " 个边界点" << std::endl;
		}

		// 找到这个点所在的NURBS边界曲线
		bool foundCurve = false;
		double minDistance = std::numeric_limits<double>::max();
		int meshPointId = v->id();

		// 计算点到所有曲线上采样点的最小距离
		for (auto qbm_boundarySample : m_boundarySamples)
		{
			for (const auto& sample : qbm_boundarySample) {
				// 使用安全的距离计算
				CPoint diff = v->point() - sample.point;
				double distance = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
				minDistance = std::min(minDistance, distance);
			}
		}

		matchedPoints++;
		pointDistances.push_back(std::make_pair(meshPointId, minDistance));
		maxMinDistance = std::max(maxMinDistance, minDistance);

	}

	std::cout << "处理完成，匹配到曲线的边界点数量: " << matchedPoints << "/" << boundaryPointCount << std::endl;
	std::cout << "Hausdorff距离: " << maxMinDistance << std::endl;

	//Step 4 & 5: 返回Hausdorff距离并排序距离列表
	if (!pointDistances.empty()) {
		std::cout << "开始排序距离列表..." << std::endl;
		std::sort(pointDistances.begin(), pointDistances.end(),
			[](const std::pair<int, double>& a, const std::pair<int, double>& b) {
				return a.second > b.second; // 从大到小排序
			});

		std::cout << "排序完成，共 " << pointDistances.size() << " 个点" << std::endl;
		if (!pointDistances.empty()) {
			std::cout << "最大距离点: ID=" << pointDistances[0].first
				<< ", 距离=" << pointDistances[0].second << std::endl;
		}
	}
	else {
		std::cout << "警告: 没有找到匹配的边界点，无法计算Hausdorff距离" << std::endl;
	}

	sortedDistances = pointDistances;
	std::cout << "=== Hausdorff距离计算完成 ===" << std::endl;
	return maxMinDistance;
}

double CCG_QMSLib::CCG_QMS_model::computeHausdorffDistance_allNURBSBoundaryCurve2(M* objTriMesh, double sampleStep, std::vector<std::pair<int, double>>& sortedDistances)
{
	/*std::cout << "=== start computing Hausdorff diatance ===" << std::endl;
	std::cout << " # NURBS boundary curve: " << nbcs().size() << std::endl;
	std::cout << "sampleStep: " << sampleStep << std::endl;*/

	//clear and intialize sampling points vector
	m_boundarySamples.clear();
	m_boundarySamples.resize(nbcs().size());

	int totalSamples = 0; //record the total number of sampling points

	//step1：Sampling points are taken on each NURBS boundary curve
	for (size_t curveIdx = 0; curveIdx < nbcs().size(); ++curveIdx) {
		const NURBS_boundaryCurve& curve = nbcs()[curveIdx];
		int patchId = curve.patchId;
		int curveType = curve.curveType;
	
		std::shared_ptr<CCG_QMS_bsplineSurf> bsurf = nullptr;
		for (auto& surf : bsplineSurfs()) {
			if (surf->id() == patchId) {
				bsurf = surf;
				break;
			}
		}

		if (!bsurf) {
			std::cerr << "  no find surface ID: " << patchId << "，neglect curve" << std::endl;
			continue;
		}

		//Determine the range of sampling parameters based on the type of boundary curve
		double fixedParam = 0.0;
		double startParam = 0.0;
		double endParam = 1.0;
		bool isUFixed = false;//Is the "u" parameter fixed？

		//Set the parameter ranges for different curve types
		switch (curveType) {
		case 1: // u = u_min
			isUFixed = true;
			fixedParam = bsurf->knotsU()[bsurf->degree()[0]]; // u_min
			startParam = bsurf->knotsV()[bsurf->degree()[1]]; // v_min
			endParam = bsurf->knotsV()[bsurf->knotsV().size() - bsurf->degree()[1] - 1]; // v_max
			break;
		case 2: // v = v_min
			isUFixed = false;
			fixedParam = bsurf->knotsV()[bsurf->degree()[1]]; // v_min
			startParam = bsurf->knotsU()[bsurf->degree()[0]]; // u_min
			endParam = bsurf->knotsU()[bsurf->knotsU().size() - bsurf->degree()[0] - 1]; // u_max
			break;
		case 3: // u = u_max
			isUFixed = true;
			fixedParam = bsurf->knotsU()[bsurf->knotsU().size() - bsurf->degree()[0] - 1]; // u_max
			startParam = bsurf->knotsV()[bsurf->degree()[1]]; // v_min
			endParam = bsurf->knotsV()[bsurf->knotsV().size() - bsurf->degree()[1] - 1]; // v_max
			break;
		case 4: // v = v_max
			isUFixed = false;
			fixedParam = bsurf->knotsV()[bsurf->knotsV().size() - bsurf->degree()[1] - 1]; // v_max
			startParam = bsurf->knotsU()[bsurf->degree()[0]]; // u_min
			endParam = bsurf->knotsU()[bsurf->knotsU().size() - bsurf->degree()[0] - 1]; // u_max
			break;
		default:
			std::cerr << "  unkown curve type: " << curveType << "，neglect curve" << std::endl;
			continue;
		}

		//the sampling points on the curve
		int actualSamples = 0;
		for (double t = startParam; t <= endParam; t += sampleStep) {
			BoundaryCurveSample sample;
			if (isUFixed) {
				sample.uvParam = CPoint2(fixedParam, t);
			}
			else {
				sample.uvParam = CPoint2(t, fixedParam);
			}

			try {
				sample.point = bsurf->DeBoorAlgorithmSurface(
					bsurf->cpts(),
					bsurf->knotsU(),
					bsurf->knotsV(),
					bsurf->degree()[0],
					bsurf->degree()[1],
					sample.uvParam[0],
					sample.uvParam[1]
				);
				sample.curveIndex = curveIdx;
				m_boundarySamples[curveIdx].push_back(sample);
				actualSamples++;
			}
			catch (const std::exception& e) {
				std::cerr << "  error in computing sampling points： " << e.what() << std::endl;
				continue;
			}
		}

		// Make sure that the end point is sampled
		if (actualSamples > 0 && std::abs(m_boundarySamples[curveIdx].back().uvParam[isUFixed ? 1 : 0] - endParam) > 1e-6) {
			try {
				BoundaryCurveSample sample;
				if (isUFixed) {
					sample.uvParam = CPoint2(fixedParam, endParam);
				}
				else {
					sample.uvParam = CPoint2(endParam, fixedParam);
				}
				sample.point = bsurf->DeBoorAlgorithmSurface(
					bsurf->cpts(),
					bsurf->knotsU(),
					bsurf->knotsV(),
					bsurf->degree()[0],
					bsurf->degree()[1],
					sample.uvParam[0],
					sample.uvParam[1]
				);
				sample.curveIndex = curveIdx;
				m_boundarySamples[curveIdx].push_back(sample);
				actualSamples++;
			}
			catch (const std::exception& e) {
				std::cerr << "  error in computing end sampling points： " << e.what() << std::endl;
			}
		}

		totalSamples += actualSamples;
		//std::cout << "  # actual Samples: " << actualSamples << std::endl;
	}

	//std::cout << "# totalSamples: " << totalSamples << std::endl;

	// step 2：Calculate the minimum distance from the boundary points of the triangular mesh to the NURBS curve.
	// Define the comparison structure
	struct CompareCPoint {
		bool operator()(const CPoint& a, const CPoint& b) const {
			if (a[0] != b[0]) return a[0] < b[0];
			if (a[1] != b[1]) return a[1] < b[1];
			return a[2] < b[2];
		}
	};

	// Using a custom comparator in map
	std::map<CPoint, int, CompareCPoint> meshBoundaryPointsMap;

	// 1. Add all boundary vertices
	for (auto v : It::MVIterator(objTriMesh)) {
		if (v->boundary()) {
			meshBoundaryPointsMap[v->point()] = v->id();
		}
	}

	// 2. Sampling at the boundary
	int edgeSamplesAdded = 0;
	int maxVertexId = objTriMesh->numVertices(); // Obtain the maximum vertex ID for the purpose of allocating a new ID

	//Adding sampling points by linear interpolation along each halfedge of the triangle
	for (auto e : It::MEIterator(objTriMesh)) {
		if (e->boundary()) {
			H* h = objTriMesh->edgeHalfedge(e, 0);
			V* v1 = objTriMesh->halfedgeSource(h);
			V* v2 = objTriMesh->halfedgeTarget(h);
			CPoint p1 = v1->point();
			CPoint p2 = v2->point();
			CPoint edgeVec = p2 - p1;

			// Uniformly interpolate 3 points on each edge (including the 3 internal points between the two endpoints)
			const int numSamplesPerEdge = 5;
			for (int i = 1; i < numSamplesPerEdge; ++i) {
				double t = static_cast<double>(i) / numSamplesPerEdge;
				CPoint samplePoint;
				samplePoint[0] = p1[0] + t * edgeVec[0];
				samplePoint[1] = p1[1] + t * edgeVec[1];
				samplePoint[2] = p1[2] + t * edgeVec[2];
				// Assign new IDs to the sampling points
				/*int sampleId = maxVertexId + edgeSamplesAdded;
				meshBoundaryPointsMap[samplePoint] = sampleId;*/
				//id = halfedge source's id
				meshBoundaryPointsMap[samplePoint] = v1->id();
				edgeSamplesAdded++;
			}
		}
	}
	//std::cout << "#Add sampling points on boundary edge: " << edgeSamplesAdded << std::endl;
	//std::cout << "# sampling points: " << meshBoundaryPointsMap.size() << std::endl;

	double maxMinDistance = 0.0; // Hausdorff distance
	std::vector<std::pair<int, double>> pointDistances;
	std::vector<std::pair<CPoint, CPoint>> triPointToNurbsPoints;
	int processedPoints = 0;

	//Calculate the minimum distance for each mesh boundary vertices (including edge sampling points)
	for (auto& entry : meshBoundaryPointsMap) {
		CPoint meshPoint = entry.first; //mesh boundary vertices's coordinate
		int pointId = entry.second;

		processedPoints++;
		if (processedPoints % 1000 == 0) {
			/*std::cout << "#process:  " << processedPoints << "/" << meshBoundaryPointsMap.size()
				<< std::endl;*/
		}

		double minDistance = std::numeric_limits<double>::max();
		CPoint closestSamplePoint;

		// Calculate the minimum distance to all the NURBS sampling points
		for (const auto& curveSamples : m_boundarySamples) {
			for (const auto& sample : curveSamples) {
				CPoint diff = meshPoint - sample.point;
				double distance = diff.norm(); // 
				if (distance < minDistance) {
					minDistance = distance;
					closestSamplePoint = sample.point;
				}
			}
		}

		pointDistances.emplace_back(pointId, minDistance);
		triPointToNurbsPoints.emplace_back(meshPoint, closestSamplePoint);
		maxMinDistance = std::max(maxMinDistance, minDistance);
	}

	// step 3：sorting
	if (!pointDistances.empty()) {
		//Sort by distance in descending order
		std::sort(pointDistances.begin(), pointDistances.end(),
			[](const auto& a, const auto& b) {
				return a.second > b.second;
			});

		//Reorder the point pairs to match the sorted distances
		std::vector<std::pair<CPoint, CPoint>> sortedPointPairs;
		for (const auto& pd : pointDistances) {
			size_t idx = &pd - &pointDistances[0];
			sortedPointPairs.push_back(triPointToNurbsPoints[idx]);
		}
		triPointToNurbsPoints = sortedPointPairs;
	}

	sortedDistances = pointDistances;
	return maxMinDistance;
}

void CCG_QMSLib::CCG_QMS_model::outputSortedDistances(const std::vector<std::pair<int, double>>& sortedDistances, const std::string& filename)
{
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "can not open file：" << filename << std::endl;
		return;
	}
	// output format：ID distance
	for (const auto& pair : sortedDistances) {
		outFile << std::left << std::setw(10) << pair.first << std::left << std::setw(10) << pair.second << std::endl;
	}
	outFile.close();
}

void CCG_QMSLib::CCG_QMS_model::computeRealSamplingPos(int faceId, int firstVId, CPoint2 uv, std::vector<int>& controlIndexs, std::vector<double>& controlWeights)
{
	auto bf = this->bfs()[faceId - 1];
	std::vector<int> degree;
	degree.push_back(bf->degree(0));
	degree.push_back(bf->degree(1));
	std::vector<double> bfw;
	bfw = bezierSurfcontrolPtsWeights(degree, uv, bf->maxU(), bf->maxV());
	if (abs(vectorSum(bfw) - 1.0) > 1e-8)
	{
		std::cout << "basis controlWeights sum: " << vectorSum(bfw) << std::endl;
	}
	/*
	* Calculate the control points and corresponding weights of the association. 
	* Here, the local parameterization of the first point (0,0) needs to be considered. 
	* Before that, the control points in the Bezier surface should be sorted counterclockwise according to the firstVId.
	*/
	BE_linkingCtlPoints tempV;
	for (int i = 0; i < bf->linkCptPWs().size(); i++)
	{
		for (auto cwIt : bf->linkCptPWs()[i].linkingCtlPId_weights)
		{
			std::pair<int, double> tempCw = cwIt;
			//std::cout << " before weight: " << tempCw.second << std::endl;
			//std::cout << " coff weight: " << bfw[i] << std::endl;
			tempCw.second *= bfw[i];
			if (tempV.linkingCtlPId_weights.find(tempCw.first) == tempV.linkingCtlPId_weights.end())
			{
				tempV.linkingCtlPId_weights.insert(tempCw);
			}
			else
			{
				//std::cout << "****************" << std::endl;
				std::map<int, double>::iterator tIt = tempV.linkingCtlPId_weights.find(tempCw.first);
				tIt->second += tempCw.second;
			}
			/*check tempv update*/
			/*std::cout << "++++" << std::endl;
			for (auto it : tempV.linkingCtlPId_weights)
			{
				std::cout << " index:" << std::setw(8) << it.first << " weight:" << std::setw(8) << it.second << std::endl;
			}*/
		}
	}

	for (auto it : tempV.linkingCtlPId_weights)
	{
		controlIndexs.push_back(it.first);
		controlWeights.push_back(it.second);
		//std::cout << " index:" << std::setw(8) << it.first << " weight:" << std::setw(8) << it.second << std::endl;
	}

	double error = 1e-8;
	if (abs(vectorSum(controlWeights) - 1.0) > error)
	{
		for (auto it : tempV.linkingCtlPId_weights)
		{
			std::cout << " index:" << std::setw(8) << it.first << " weight:" << std::setw(8) << it.second << std::endl;
		}
		std::cout << " #bezier ctl pts: " << bfw.size() << std::endl;
		for (auto _bfw : bfw)
		{
			std::cout << "Bezier control point weight:  " << _bfw << std::endl;
		}

		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			std::cout << " ---bezier control point: " << i << std::endl;
			for (auto cwIt : bf->linkCptPWs()[i].linkingCtlPId_weights)
			{
				std::cout << " index:" << std::setw(8) << cwIt.first << " weight:" << std::setw(8) << cwIt.second << std::endl;
			}
		}

		std::cout << "faceId: " << faceId << "  firstVid: " << firstVId << std::endl;
		std::cout << "bf->maxU(): " << bf->maxU() << " bf->maxV(): " << bf->maxV() << std::endl;
		std::cout << "uv: " << uv[0] << " " << uv[1] << std::endl;
		std::cout << "controlWeights sum: " << vectorSum(controlWeights) << std::endl;
	}
	/*else
	{
		std::cout << "controlWeights sum: " << vectorSum(controlWeights) << std::endl;
	}*/
}

void CCG_QMSLib::CCG_QMS_model::obtainBezierCptsBy_CWs_scale(M* pMesh, double scale)
{
	/*将所有beizier控制点初始化为{0，0，0}*/
	for (auto bf : this->bfs())
	{
		CPoint tempP;
		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			bf->cpts()[i] = tempP;
		}
	}
	/*通过初始控制点和权重计算beizier控制点*/
	for (auto bf : this->bfs())
	{
		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			CPoint tempP;
			for (auto cwIt : bf->linkCptPWs()[i].linkingCtlPId_weights)
			{
				tempP += ((pMesh->idVertex(cwIt.first)->point()) * cwIt.second);
			}
			bf->cpts()[i] = tempP / scale;
			//bf->cpts()[i] = bf->cpts()[i] / scale;
		}
	}
}

void CCG_QMSLib::CCG_QMS_model::create_bsurBoundarys_general(M* pMesh)
{
	std::vector<H*> tempMarkHs;
	/* travese all b spline surfaces*/
	for (auto bf : this->bsplineSurfs())
	{
		/*case 1: torus type,no boundary*/
		if (bf->closedAlongU() && bf->closedAlongV()) continue;

		/* create boundary by bf's qbb_quadIds*/
		/*The first boundary, v = v_min*/
		{
			NUBE_bsurfBoudnary* firstBoundary = new NUBE_bsurfBoudnary();
			firstBoundary->bsurfId() = bf->id();
			firstBoundary->localId() = 2;
			/*case 2: U is not closed, V is not closed or U is closed, V is not closed*/
			//if ((!bf->closedAlongU() && !bf->closedAlongV()) || (bf->closedAlongU() && !bf->closedAlongV()))
			{
				int firstVId = bf->quadIds()[0];
				int endVId = bf->quadIds()[1];
				V* firstV = pMesh->idVertex(firstVId);
				H* firstH = NULL;
				if (bf->closedAlongU())
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeEdge(vh)->boundary() && pMesh->halfedgeFace(vh)->f_spline_patchId() == firstBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}
				else
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeFace(vh)->f_spline_patchId() == firstBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}
				if (firstH == NULL)
				{
					std::cout << "error!!!, there is no finding the first halfedges on first boundary!" << std::endl;
					return;
				}
				firstBoundary->halfedges().push_back(firstH);
				while (pMesh->halfedgeTarget(firstH)->id() != endVId)
				{
					firstH = pMesh->halfedgeNext(pMesh->halfedgeSym(pMesh->halfedgeNext(firstH)));
					firstBoundary->halfedges().push_back(firstH);
				}
			}
			bf->boundarys().push_back(firstBoundary);
		}

		/*The second boundary, u = u_max*/
		{
			NUBE_bsurfBoudnary* secondBoundary = new NUBE_bsurfBoudnary();
			secondBoundary->bsurfId() = bf->id();
			secondBoundary->localId() = 3;
			/*case 3: U is not closed, V is not closed or U is not closed, V is closed*/
			//if ((!bf->closedAlongU() && !bf->closedAlongV()) || (!bf->closedAlongU() && bf->closedAlongV()))
			{
				int firstVId = bf->quadIds()[1];
				int endVId = bf->quadIds()[2];
				V* firstV = pMesh->idVertex(firstVId);
				H* firstH = NULL;
				if (bf->closedAlongV())
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeEdge(vh)->boundary() && pMesh->halfedgeFace(vh)->f_spline_patchId() == secondBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}
				else
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeFace(vh)->f_spline_patchId() == secondBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}

				if (firstH == NULL)
				{
					std::cout << "error!!!, there is no finding the first halfedges on second boundary!" << std::endl;
					return;
				}
				secondBoundary->halfedges().push_back(firstH);
				while (pMesh->halfedgeTarget(firstH)->id() != endVId)
				{
					firstH = pMesh->halfedgeNext(pMesh->halfedgeSym(pMesh->halfedgeNext(firstH)));
					secondBoundary->halfedges().push_back(firstH);
				}
			}
			bf->boundarys().push_back(secondBoundary);
		}

		/*The third boundary, v = v_max*/
		{
			NUBE_bsurfBoudnary* thirdBoundary = new NUBE_bsurfBoudnary();
			thirdBoundary->bsurfId() = bf->id();
			thirdBoundary->localId() = 4;
			/*case 2: U is not closed, V is not closed or U is closed, V is not closed*/
			//if ((!bf->closedAlongU() && !bf->closedAlongV()) || (bf->closedAlongU() && !bf->closedAlongV()))
			{
				int firstVId = bf->quadIds()[2];
				int endVId = bf->quadIds()[3];
				V* firstV = pMesh->idVertex(firstVId);
				H* firstH = NULL;

				if (bf->closedAlongU())
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeEdge(vh)->boundary() && pMesh->halfedgeFace(vh)->f_spline_patchId() == thirdBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}
				else
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeFace(vh)->f_spline_patchId() == thirdBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}

				if (firstH == NULL)
				{
					std::cout << "error!!!, there is no finding the first halfedges on third boundary!" << std::endl;
					return;
				}
				thirdBoundary->halfedges().push_back(firstH);
				while (pMesh->halfedgeTarget(firstH)->id() != endVId)
				{
					firstH = pMesh->halfedgeNext(pMesh->halfedgeSym(pMesh->halfedgeNext(firstH)));
					thirdBoundary->halfedges().push_back(firstH);
				}
			}
			bf->boundarys().push_back(thirdBoundary);
		}

		/*The fourth boundary, v = v_min*/
		{
			NUBE_bsurfBoudnary* fourthBoundary = new NUBE_bsurfBoudnary();
			fourthBoundary->bsurfId() = bf->id();
			fourthBoundary->localId() = 1;
			/*case 3: U is not closed, V is not closed or U is not closed, V is closed*/
			//if ((!bf->closedAlongU() && !bf->closedAlongV()) || (!bf->closedAlongU() && bf->closedAlongV()))
			{
				int firstVId = bf->quadIds()[3];
				int endVId = bf->quadIds()[0];
				V* firstV = pMesh->idVertex(firstVId);
				H* firstH = NULL;
				if (bf->closedAlongV())
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeEdge(vh)->boundary() && pMesh->halfedgeFace(vh)->f_spline_patchId() == fourthBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}
				else
				{
					for (auto vh : It::VCcwOutHEIterator(pMesh, firstV))
					{
						if (pMesh->halfedgeFace(vh)->f_spline_patchId() == fourthBoundary->bsurfId())
						{
							if (vh->visit()) continue;
							firstH = vh;
							vh->visit() = true;
							tempMarkHs.push_back(vh);
							break;
						}
					}
				}

				if (firstH == NULL)
				{
					std::cout << "error!!!, there is no finding the first halfedges on second boundary!" << std::endl;
					return;
				}
				fourthBoundary->halfedges().push_back(firstH);
				while (pMesh->halfedgeTarget(firstH)->id() != endVId)
				{
					firstH = pMesh->halfedgeNext(pMesh->halfedgeSym(pMesh->halfedgeNext(firstH)));
					fourthBoundary->halfedges().push_back(firstH);
				}
			}
			bf->boundarys().push_back(fourthBoundary);
		}
	}

	for (auto h : tempMarkHs)
	{
		h->visit() = false;
	}

	/*create segments for every boundary*/
	for (auto bf : this->bsplineSurfs())
	{
		for (auto bf_boundary : bf->boundarys())
		{
			if (bf_boundary->halfedges().size() == 0) continue;
			H* firstH = bf_boundary->halfedges()[0];
			std::vector<std::pair<int, int>> segments2HesIndexs;
			int count = 1;
			int firstSegHeindex = 1;
			for (int i = 0; i < bf_boundary->halfedges().size(); i++)
			{
				if (pMesh->idVertex(bf_boundary->halfedges()[i]->target()->id())->ifCorner())
				{
					NUBE_bsurfSegment* bfBSeg = new NUBE_bsurfSegment;
					bfBSeg->bsurfBoundaryLoacalId() = bf_boundary->localId();
					bfBSeg->bsurfId() = bf_boundary->bsurfId();
					bfBSeg->id() = count;
					bfBSeg->firstHalfedgeIndex() = firstSegHeindex;
					bfBSeg->endHalfedgeIndex() = i + 1;
					bf_boundary->segements().push_back(bfBSeg);
					count++;
					firstSegHeindex = i + 2;
				}
			}
		}
	}

	/*obtaining symmetry segment for segment on boundry*/
	{
		for (auto bf : this->bsplineSurfs())
		{
			for (auto bf_boundary : bf->boundarys())
			{
				for (auto seg : bf_boundary->segements())
				{
					H* segFirstH = bf_boundary->halfedges()[seg->firstHalfedgeIndex() - 1];
					H* segFirsHSym = pMesh->halfedgeSym(segFirstH);
					if (segFirsHSym != NULL)
					{
						F* segFirsHSymF = pMesh->halfedgeFace(segFirsHSym);
						for (auto _bf_boudnary : this->bsplineSurfs()[segFirsHSymF->f_spline_patchId() - 1]->boundarys())
						{
							for (auto _seg : _bf_boudnary->segements())
							{
								if (_bf_boudnary->halfedges()[_seg->firstHalfedgeIndex() - 1]->source()->id() == bf_boundary->halfedges()[seg->endHalfedgeIndex() - 1]->target()->id() &&
									_bf_boudnary->halfedges()[_seg->endHalfedgeIndex() - 1]->target()->id() == bf_boundary->halfedges()[seg->firstHalfedgeIndex() - 1]->source()->id())
								{
									seg->symBsurfSeg() = _seg;
									break;
								}
							}
							if (seg->symBsurfSeg() != NULL) break;
						}
					}
				}
			}
		}
	}

	/*check the size of boundary on b spline surface*/
	/*for (auto bf : this->bsplineSurfs())
	{
		std::cout << " bf->id(): " << bf->id() << "  bf->boundarys().size(): " << bf->boundarys().size()<< std::endl;
		for (auto bfB : bf->boundarys())
		{
			std::cout << " ------bfB->localId(): " << bfB->localId() << "  bfB->segements().size(): " << bfB->segements().size() << std::endl;
			for (auto seg : bfB->segements())
			{
				if (seg->symBsurfSeg() != NULL)
				{
					std::cout << "---------------------- seg->symBsurfSeg()->bsurfId()" << seg->symBsurfSeg()->bsurfId() << std::endl;
				}
				else
				{
					std::cout << " seg->bsurfId(): " << seg->bsurfId() << " seg->bsurfBoundaryLoacalId(): " << seg->bsurfBoundaryLoacalId() << " seg->id(): " << seg->id() << std::endl;
				}
			}
		}
	}*/
}

void CCG_QMSLib::CCG_QMS_model::createTopology_bspline_TMesh(M* pMesh)
{
	for (auto v : It::MVIterator(pMesh))
	{
		v->ifVisit() = false;
	}

	/*Create the T (polygon) mesh corresponding to the B-spline surface patch*/
	for (auto bp : this->bsplineSurfs())
	{
		if (bp->closedAlongU() || bp->closedAlongV()) continue;//The periodic B-spline surface does not have a T (polygon) mesh topology corresponding to it, and thus requires special treatment.

		for (auto bpQuadId : bp->quadIds())
		{
			V* bpQuadIdV = pMesh->idVertex(bpQuadId);
			if (!bpQuadIdV->ifVisit())
			{
				bpQuadIdV->ifVisit() = true;
				//create vertex
				V* v = this->bsplineSurfacesTopology_quad().createVertex(bpQuadIdV->id());
				v->point() = bpQuadIdV->point();
			}

		}
	}

	//Create the T (polygon) mesh corresponding to the B-spline surface patch
	for (auto bp : this->bsplineSurfs())
	{
		if (bp->closedAlongU() || bp->closedAlongV()) continue;//The periodic B-spline surface does not have a T (polygon) mesh topology corresponding to it, and thus requires special treatment.
		//create face
		std::vector<V*> vs;
		for (auto bpB : bp->boundarys())
		{
			for (auto bpBS : bpB->segements())
			{
				int vid = bpB->halfedges()[bpBS->firstHalfedgeIndex() - 1]->source()->id();
				vs.push_back((this->bsplineSurfacesTopology_quad().idVertex(vid)));
			}
		}
		this->bsplineSurfacesTopology_quad().createFace(vs, bp->id());
	}

	for (auto v : It::MVIterator(pMesh))
	{
		v->ifVisit() = false;
	}
}

void CCG_QMSLib::CCG_QMS_model::outputTopology_bspline_TMesh_obj(M* pMesh, const char* output)
{
	std::ofstream outFile(output);
	if (!outFile.is_open()) {
		std::cerr << "Failed to open file: " << output << std::endl;
		return;
	}

	for (auto v : It::MVIterator(pMesh))
	{
		v->ifVisit() = false;
	}

	std::map<int, int> idMaps;
	int vCount = 1;

	/*Create the T (polygon) mesh corresponding to the B-spline surface patch*/
	for (auto bp : this->bsplineSurfs())
	{
		if (bp->closedAlongU() || bp->closedAlongV()) continue;//The periodic B-spline surface does not have a T (polygon) mesh topology corresponding to it, and thus requires special treatment.

		for (auto bpQuadId : bp->quadIds())
		{
			V* bpQuadIdV = pMesh->idVertex(bpQuadId);
			if (!bpQuadIdV->ifVisit())
			{
				bpQuadIdV->ifVisit() = true;
				//output v
				//V* v = this->bsplineSurfacesTopology_quad().createVertex(bpQuadIdV->id());
				idMaps.insert(std::pair<int, int>(bpQuadIdV->id(), vCount));
				vCount++;
				outFile << "v " << bpQuadIdV->point()[0] << " " << bpQuadIdV->point()[1] << " " << bpQuadIdV->point()[2] << "\n";
				//outFile << "v " << bpQuadIdV->point()[0] << " " << bpQuadIdV->point()[1] << " " << 0.0 << "\n";
			}

		}

	}

	//create faces
	for (auto bp : this->bsplineSurfs())
	{
		if (bp->closedAlongU() || bp->closedAlongV()) continue;//The periodic B-spline surface does not have a T (polygon) mesh topology corresponding to it, and thus requires special treatment.
		outFile << "f";
		std::vector<V*> vs;
		for (auto bpB : bp->boundarys())
		{
			for (auto bpBS : bpB->segements())
			{
				int vid = bpB->halfedges()[bpBS->firstHalfedgeIndex() - 1]->source()->id();
				//vs.push_back((this->bsplineSurfacesTopology_quad().idVertex(vid)));
				outFile << " " << idMaps.at(vid);
			}
		}
		outFile << "\n";
	}
	//this->bsplineSurfacesTopology_quad().createFace(vs, bp->id());
	for (auto v : It::MVIterator(pMesh))
	{
		v->ifVisit() = false;
	}
}

double CCG_QMSLib::CCG_QMS_model::computeHausdorffDistance_allNURBSBoundaryCurve3(M* objTriMesh, int _numSamplesPerEdge, std::vector<std::pair<int, double>>& sortedDistances)
{
	/*std::cout << "=== 开始计算Hausdorff距离 ===" << std::endl;
	std::cout << " NURBS边界曲线数量: " << nbcs().size() << std::endl;
	std::cout << "采样步长: " << sampleStep << std::endl;*/

	//---计算三角网格边界点到NURBS曲线的最小距离----
	//std::cout << "开始计算三角网格边界点到NURBS曲线的最小距离..." << std::endl;

	// 定义比较结构体
	struct CompareCPoint {
		bool operator()(const CPoint& a, const CPoint& b) const {
			if (a[0] != b[0]) return a[0] < b[0];
			if (a[1] != b[1]) return a[1] < b[1];
			return a[2] < b[2];
		}
	};

	// 使用自定义比较器的map
	std::map<CPoint, int, CompareCPoint> meshBoundaryPointsMap;


	// 收集所有网格边界点和边采样点
	//std::map<CPoint,int> meshBoundaryPointsMap; // 存储边界点ID和坐标, CPoint> meshBoundaryPointsMap; 
	//std::vector<CPoint> meshBoundaryPoints;

	// 1. 添加所有边界顶点
	for (auto v : It::MVIterator(objTriMesh)) {
		if (v->boundary()) {
			meshBoundaryPointsMap[v->point()] = v->id();
		}
	}

	// 2. 在边界边上采样
	//std::cout << "开始在三角网格边界边上采样..." << std::endl;
	int edgeSamplesAdded = 0;
	int maxVertexId = objTriMesh->numVertices(); // 获取最大顶点ID，用于分配新ID

	//通过步长添加采样点
	//for (auto e : NUBELib::It::MEIterator(objTriMesh)) {
	//	if (e->boundary()) {
	//		H* h = objTriMesh->edgeHalfedge(e, 0);
	//		V* v1 = objTriMesh->halfedgeSource(h);
	//		V* v2 = objTriMesh->halfedgeTarget(h);
	//		CPoint p1 = v1->point();
	//		CPoint p2 = v2->point();
	//		CPoint edgeVec = p2 - p1;
	//		double edgeLength = sqrt(edgeVec[0] * edgeVec[0] + edgeVec[1] * edgeVec[1] + edgeVec[2] * edgeVec[2]);

	//		// 计算这条边上需要的采样点数
	//		int samplesOnEdge = static_cast<int>(edgeLength / 0.0001);
	//		if (samplesOnEdge < 1) samplesOnEdge = 1;

	//		std::cout << "  边界边 " << v1->id() << "-" << v2->id()
	//			<< ", 长度: " << edgeLength
	//			<< ", 采样点数: " << samplesOnEdge << std::endl;

	//		// 在边上均匀采样
	//		for (int i = 1; i < samplesOnEdge; ++i) {
	//			double t = static_cast<double>(i) / samplesOnEdge;
	//			CPoint samplePoint;
	//			samplePoint[0] = p1[0] + t * edgeVec[0];
	//			samplePoint[1] = p1[1] + t * edgeVec[1];
	//			samplePoint[2] = p1[2] + t * edgeVec[2];

	//			// 为采样点分配新ID（从maxVertexId开始）
	//			int sampleId = maxVertexId + edgeSamplesAdded;
	//			meshBoundaryPointsMap[samplePoint] = sampleId;
	//			edgeSamplesAdded++;
	//		}
	//	}
	//}

	//每条三角形半边线性插值添加采样点
	for (auto e : It::MEIterator(objTriMesh)) {
		if (e->boundary()) {
			H* h = objTriMesh->edgeHalfedge(e, 0);
			V* v1 = objTriMesh->halfedgeSource(h);
			V* v2 = objTriMesh->halfedgeTarget(h);
			CPoint p1 = v1->point();
			CPoint p2 = v2->point();
			CPoint edgeVec = p2 - p1;

			// 在每个边上均匀插值 _numSamplesPerEdge 个点（包括两个端点之间的两个内点）
			const int numSamplesPerEdge = _numSamplesPerEdge + 2; // 包括两个端点，实际新增 _numSamplesPerEdge 个点
			for (int i = 1; i < numSamplesPerEdge; ++i) {
				double t = static_cast<double>(i) / numSamplesPerEdge;
				CPoint samplePoint;
				samplePoint[0] = p1[0] + t * edgeVec[0];
				samplePoint[1] = p1[1] + t * edgeVec[1];
				samplePoint[2] = p1[2] + t * edgeVec[2];
				// 为采样点分配新ID
				/*int sampleId = maxVertexId + edgeSamplesAdded;
				meshBoundaryPointsMap[samplePoint] = sampleId;*/
				//id=半边的source的id
				meshBoundaryPointsMap[samplePoint] = v1->id();
				edgeSamplesAdded++;
			}
		}
	}
	//std::cout << "在边界边上添加了 " << edgeSamplesAdded << " 个采样点" << std::endl;
	//std::cout << "三角网格边界总采样点数: " << meshBoundaryPointsMap.size() << std::endl;

	double maxMinDistance = 0.0; // Hausdorff距离
	std::vector<std::pair<int, double>> pointDistances;
	std::vector<std::pair<CPoint, CPoint>> triPointToNurbsPoints;
	int processedPoints = 0;

	// 为每个网格边界点(包括边采样点)计算最小距离
	for (auto& entry : meshBoundaryPointsMap) {
		CPoint meshPoint = entry.first; // 网格边界点坐标
		int pointId = entry.second;

		processedPoints++;
		if (processedPoints % 1000 == 0) {
			/*std::cout << "已处理 " << processedPoints << "/" << meshBoundaryPointsMap.size()
				<< " 个边界点" << std::endl;*/
		}

		double minDistance = std::numeric_limits<double>::max();
		CPoint closestSamplePoint;

		// 计算到所有NURBS采样点的最小距离
		//for (const auto& curveSamples : qbm_boundarySamples) {
		//	for (const auto& sample : curveSamples) {
		//		CPoint diff = meshPoint - sample.point;
		//		double distance = diff.norm(); // 
		//		if (distance < minDistance) {
		//			minDistance = distance;
		//			closestSamplePoint = sample.point;
		//		}
		//	}
		//}
		for (auto bf : this->bsplineSurfs())
		{
			for (auto bfBoundary : bf->boundarys())
			{
				if (!bfBoundary->boundary()) continue;
				for (auto bfBoundarySeg : bfBoundary->segements())
				{
					for (auto sample : bfBoundarySeg->samplings()) {
						CPoint diff = meshPoint - sample->pos();
						double distance = diff.norm(); // 
						if (distance < minDistance) {
							minDistance = distance;
							closestSamplePoint = sample->pos();
						}
					}
				}
			}
		}

		pointDistances.emplace_back(pointId, minDistance);
		triPointToNurbsPoints.emplace_back(meshPoint, closestSamplePoint);
		maxMinDistance = std::max(maxMinDistance, minDistance);
	}

	// 第三步：排序并记录前五大距离
	if (!pointDistances.empty()) {
		// 按距离降序排序
		std::sort(pointDistances.begin(), pointDistances.end(),
			[](const auto& a, const auto& b) {
				return a.second > b.second;
			});

		// 重新排序点对以匹配排序后的距离
		std::vector<std::pair<CPoint, CPoint>> sortedPointPairs;
		for (const auto& pd : pointDistances) {
			size_t idx = &pd - &pointDistances[0];
			sortedPointPairs.push_back(triPointToNurbsPoints[idx]);
		}
		triPointToNurbsPoints = sortedPointPairs;

		//// 存储前五大距离点对
		//topDistancePoints.clear();
		//int topCount = std::min(5, (int)pointDistances.size()); // 改为记录前5大距离
		//topDistancePoints.resize(topCount);

		//for (int i = 0; i < topCount; i++) {
		//	topDistancePoints[i] = triPointToNurbsPoints[i];

		//	// 输出前五大距离信息
		//	/*std::cout << "第 " << i + 1 << " 大距离: " << pointDistances[i].second
		//		<< ", 网格点ID: " << pointDistances[i].first << std::endl;
		//	std::cout << "  网格点坐标: (" << topDistancePoints[i].first[0] << ", "
		//		<< topDistancePoints[i].first[1] << ", " << topDistancePoints[i].first[2] << ")" << std::endl;
		//	std::cout << "  NURBS采样点坐标: (" << topDistancePoints[i].second[0] << ", "
		//		<< topDistancePoints[i].second[1] << ", " << topDistancePoints[i].second[2] << ")" << std::endl;*/
		//}
	}

	sortedDistances = pointDistances;
	//std::cout << "=== Hausdorff距离计算完成 ===" << std::endl;
	//std::cout << "最大Hausdorff距离: " << maxMinDistance << std::endl;

	/*computing spline surface boundary's hausdorff distance*/
	for (auto surf : this->bsplineSurfs()) {
		if (!surf) continue;
		// 这块曲面上的所有边界 (通常有 4 条)
		const std::vector<NUBE_bsurfBoudnary*>& boundaries = surf->boundarys();
		if (boundaries.empty()) {
			// 如果没有边界，就跳过
			//std::cout << "曲面 ID=" << surf->id() << " 没有边界，跳过 Hausdorff 计算。" << std::endl;
			continue;
		}

		// 记录这块曲面上的所有边界hausdorff距离的最大值
		double surfhausdorff = 0.0;
		// 遍历这块曲面的每条边界
		for (auto bnd : boundaries) {
			if (!bnd || !bnd->boundary()) continue;
			// 在该边界参数域上按照 sampleStep 做均匀采样
			std::vector<CPoint> boundarySamplePoints;

			// 对这条边界上的每个采样点，计算它到“网格边界点集”中所有点的最小距离
			double hausdorff = 0.0;
			for (auto bndSeg : bnd->segements())
			{
				for (auto sampP : bndSeg->samplings())
				{
					double minDist = std::numeric_limits<double>::infinity();

					// 暴力遍历 meshBoundaryPointsMap 中的所有点
					for (auto& kv : meshBoundaryPointsMap) {
						const CPoint& meshP = kv.first;
						CPoint diff = sampP->pos() - meshP;
						double d = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
						if (d < minDist) {
							minDist = d;
						}
					}
					sampP->diffG0() = minDist;//distance to all boundary vertex or edges from tri mesh
					// 更新这条边界的 Hausdorff 距离 = max(当前最大, 该采样点的最小距离)
					if (minDist > hausdorff) {
						hausdorff = minDist;
					}
				}
			}
			// 把计算得到的 Hausdorff 距离写入边界对象
			bnd->hausdorffDist() = hausdorff;
			/*std::cout << "  B surf ID=" << surf->id()
				<< "B surf boundary localId=" << bnd->localId()
				<< " Hausdorff distance: " << hausdorff << std::endl;*/

				// 更新这个面片的最大 Hausdorff 距离
			if (hausdorff > surfhausdorff) {
				surfhausdorff = hausdorff;
			}
		}
		surf->maxHausdorff() = surfhausdorff;
	}
	this->maxHausdorff() = maxMinDistance;
	return maxMinDistance;
}

void CCG_QMSLib::CCG_QMS_model::computeGlobalSurfaceMetrics()
{
	// 先取出所有 B 样条曲面指针
	const std::vector< std::shared_ptr<CCG_QMS_bsplineSurf>>& allSurfs = this->bsplineSurfs();

	// 初始化：maxRatio 取极小值，minDist 取极大值，hausdorff 取极小值
	double maxRatio = std::numeric_limits<double>::lowest();   // 更安全地用 lowest() 而不是 0，以防所有 ratio 都是负数（一般不会）
	double minDist = std::numeric_limits<double>::infinity();
	double maxHausdorff = std::numeric_limits<double>::lowest();

	// 遍历所有曲面片
	for (auto surf : allSurfs) {
		if (!surf) continue;

		// 取出这个曲面片的 boundaryEdgeRatio() 和 minCtrlPtDist()
		double ratio = surf->boundaryEdgeRatio();
		double dist = surf->minCtrlPtDist();
		double hausdorff = surf->maxHausdorff();

		// 更新最大 ratio
		if (ratio > maxRatio) {
			maxRatio = ratio;
		}
		// 更新最小 dist
		if (dist < minDist) {
			minDist = dist;
		}
		// 更新最大 hausdorff
		if (hausdorff > maxHausdorff) {
			maxHausdorff = hausdorff;
		}
	}

	// 如果所有曲面片的 minCtrlPtDist 都是 +inf（即未被更新过），说明它们都没有有效值，
	// 此时把 minDist 置为 0.0
	if (minDist == std::numeric_limits<double>::infinity()) {
		minDist = 0.0;
	}
	// 如果所有曲面片的 ratio 都小于等于 lowest()（这几乎不可能，但以防万一），把 maxRatio 置为 0.0
	if (maxRatio == std::numeric_limits<double>::lowest()) {
		maxRatio = 0.0;
	}
	if (maxHausdorff == std::numeric_limits<double>::lowest()) {
		maxHausdorff = 0.0;
	}

	this->boundaryEdgeRatio() = maxRatio;
	this->minCtrlPtDist() = minDist;
	//this->maxHausdorff() = maxHausdorff;
}

void CCG_QMSLib::CCG_QMS_model::outputMetricsToFile(const std::string& filename)
{
	std::ofstream ofs(filename);
	if (!ofs.is_open()) {
		std::cerr << "not open report.txt，please check path!" << std::endl;
	}

	// 为了让输出表格更整洁，可以先输出一个标题行和分割线
	ofs << "=====================================================\n";
	ofs << "             Spline Surface Quality Report             \n";
	ofs << "=====================================================\n\n";

	// 辅助函数：将指针转换为字符串，若为 nullptr 则返回 "NULL"
	auto ptrToString = [&](const void* ptr) {
		if (ptr == nullptr) {
			return std::string("NULL");
		}
		std::ostringstream oss;
		oss << ptr;
		return oss.str();
	};

	// 设置输出对齐：左对齐标签，右对齐数值或地址
	const int LABEL_WIDTH = 20; // 标签最长约 20 字符，留点空余
	const int VALUE_WIDTH = 30; // 数值/地址区域宽度

	// 1. qbs_diffG0_sampling
	/*ofs << std::left << std::setw(LABEL_WIDTH) << "G0："
		<< std::right << std::setw(VALUE_WIDTH) << std::fixed << std::setprecision(16) << qbs_diffG0_sampling->diffG0() << "\n";*/
	if (m_diffG0_sampling != NULL)
	{
		ofs << std::left << std::setw(LABEL_WIDTH) << "G0："
			<< std::right << std::setw(VALUE_WIDTH) << m_diffG0_sampling->diffG0() << "\n";
	}
	else
	{
		ofs << std::left << std::setw(LABEL_WIDTH) << "G0："
			<< std::right << std::setw(VALUE_WIDTH) << " This is a single patch! " << "\n";
	}

	// 2. qbs_diffG1_sampling
	/*ofs << std::left << std::setw(LABEL_WIDTH) << "G1："
		<< std::right << std::setw(VALUE_WIDTH) << std::fixed << std::setprecision(16) <<qbs_diffG1_sampling->diffG1() << "\n";*/
	if (m_diffG1_sampling != NULL)
	{
		ofs << std::left << std::setw(LABEL_WIDTH) << "G1："
			<< std::right << std::setw(VALUE_WIDTH) << m_diffG1_sampling->diffG1() << "\n";
	}
	else
	{
		ofs << std::left << std::setw(LABEL_WIDTH) << "G1："
			<< std::right << std::setw(VALUE_WIDTH) << " This is a single patch! " << "\n";
	}

	// 3. qbs_boundaryEdgeRatio
	/*ofs << std::left << std::setw(LABEL_WIDTH) << "maxBoundaryEdgeRatio："
		<< std::right << std::setw(VALUE_WIDTH) << std::fixed << std::setprecision(16)
		<< qbs_boundaryEdgeRatio << "\n";*/
	ofs << std::left << std::setw(LABEL_WIDTH) << "maxBoundaryEdgeRatio："
		<< std::right << std::setw(VALUE_WIDTH)
		<< m_boundaryEdgeRatio << "\n";

	// 4. qbs_minCtrlPtDist
	/*ofs << std::left << std::setw(LABEL_WIDTH) << "minCtrlPtDist："
		<< std::right << std::setw(VALUE_WIDTH) << std::fixed << std::setprecision(16)
		<< qbs_minCtrlPtDist << "\n";*/
	ofs << std::left << std::setw(LABEL_WIDTH) << "minCtrlPtDist："
		<< std::right << std::setw(VALUE_WIDTH)
		<< m_minCtrlPtDist << "\n";

	// 5. qbs_maxHausdorff
	/*ofs << std::left << std::setw(LABEL_WIDTH) << "boundaryHausdorff："
		<< std::right << std::setw(VALUE_WIDTH) << std::fixed << std::setprecision(16)
		<< qbs_maxHausdorff << "\n";*/
	ofs << std::left << std::setw(LABEL_WIDTH) << "boundaryHausdorff："
		<< std::right << std::setw(VALUE_WIDTH)
		<< m_maxHausdorff << "\n";

	ofs << "\n=====================================================\n";
	ofs << "                     End of Report                  \n";
	ofs << "=====================================================\n";

	ofs.close();
}

void CCG_QMSLib::CCG_QMS_model::computeSamplingPtUV(M* pMesh)
{
	for (int i = 0; i < this->samplings().size(); i++)
	{
		std::shared_ptr<QB_sampling> sp = this->samplings()[i];
		F* f = pMesh->idFace(sp->fId());
		H* fh = pMesh->faceHalfedge(f);
		E* fhE = pMesh->halfedgeEdge(fh);
		/*指向四边形网格中第一个点的半边应该是fh->nextHalfedge*/
		H* fhn = pMesh->halfedgeNext(fh);
		E* fhnE = pMesh->halfedgeEdge(fhn);
		CPoint2 cornerA, cornerB, cornerC, cornerD;
		cornerA[0] = 0.0, cornerA[1] = 0.0;
		cornerB[0] = fhE->knotInterval(), cornerB[1] = 0.0;
		cornerC[0] = fhE->knotInterval(), cornerC[1] = fhnE->knotInterval();
		cornerD[0] = 0.0, cornerD[1] = fhnE->knotInterval();
		sp->uv() = cornerA * sp->ws()[0] + cornerB * sp->ws()[1] + cornerC * sp->ws()[2] + cornerD * sp->ws()[3];

		///*check 参数的计算*/
		//std::cout << "------sp->fId(): " << sp->fId() << std::endl;
		//std::cout << "sampling " << i << " weights: ";
		//for (int ii = 0; ii < 4; ii++)
		//{
		//	std::cout << " " << sp->ws()[ii];
		//}
		//std::cout << std::endl;
		//std::cout << "u: " << sp->uv()[0] << " v: " << sp->uv()[1] << std::endl;
		//std::cout << "-----------------------------------------------------------------" << std::endl;

		/*交换uv*/
		/*double temp_u = sp->uv()[0];
		sp->uv()[0] = sp->uv()[1];
		sp->uv()[1] = temp_u;*/

		if (sp->uv()[0] > 1.0 || sp->uv()[1] > 1.0 || sp->uv()[0] < 0.0 || sp->uv()[1] < 0.0)
		{
			sp->mark() = false;
			/*std::cout << "sampling " << i << " weights: ";
			for (int ii = 0; ii < 4; ii++)
			{
				std::cout << " " << sp->ws()[ii];
			}
			std::cout << std::endl;
			std::cout << "u: " << sp->uv()[0] << " v: " << sp->uv()[1] << std::endl;*/
		}
	}
}

void CCG_QMSLib::CCG_QMS_model::computeBezierSufIntegral_mine(std::shared_ptr<CCG_QMS_bezierSurf> bf, std::map<int, BE_linkingCtlPoints>& ctl_ws)
{
	/*基函数0，1，2阶导数的表示*/
	struct powerBaseFun
	{
		int pbf_min_deg;//最小阶数
		int pbf_max_deg;//最大阶数
		std::vector<double> pbf_coeffs;//不同的阶数对应的系数
		int pbf_order;//导数0，1，2
		void show()
		{
			std::cout << "powerBaseFun degree: ";
			for (int ii = pbf_min_deg; ii <= pbf_max_deg; ii++)
			{
				std::cout << " " << ii;
			}
			std::cout << std::endl;
			std::cout << "powerBaseFun coefficient: ";
			for (auto pbf_coeff : pbf_coeffs)
			{
				std::cout << " " << pbf_coeff;
			}
			std::cout << std::endl;
		}
	};

	/*根据bf的信息创建所有基函数的0阶导数，以幂基的形式表示*/
	std::vector<powerBaseFun> bf_basisFun0Diffs;
	for (int i = 0; i <= bf->degree(0); i++)
	{
		powerBaseFun bf_basisFun0Diff;
		bf_basisFun0Diff.pbf_min_deg = i;//给当前基函数的最小阶数赋值
		bf_basisFun0Diff.pbf_max_deg = bf->degree(0);//给当前基函数的最大阶数赋值
		/*计算当前基函数用幂基形式表示时不同阶数对应的系数*/
		for (int ii = bf_basisFun0Diff.pbf_min_deg; ii <= bf_basisFun0Diff.pbf_max_deg; ii++)
		{
			double pbf_coeff = pow(-1, ii + i) * binomial(bf->degree(0), ii) * binomial(ii, i);
			bf_basisFun0Diff.pbf_coeffs.push_back(pbf_coeff);
		}
		bf_basisFun0Diff.pbf_order = 0;
		bf_basisFun0Diffs.push_back(bf_basisFun0Diff);
		/*check bf_basisFun0Diff infor*/
		/*std::cout << "------Basis function 0 differentional------" << std::endl;
		std::cout << "powerBaseFun degree: ";
		for (int ii = bf_basisFun0Diff.pbf_min_deg; ii <= bf_basisFun0Diff.pbf_max_deg; ii++)
		{
			std::cout << " " << ii;
		}
		std::cout << std::endl;
		std::cout << "powerBaseFun coefficient: ";
		for (auto pbf_coeff : bf_basisFun0Diff.pbf_coeffs)
		{
			std::cout <<" " << pbf_coeff;
		}
		std::cout << std::endl;*/
	}

	/*根据bf0阶导数的信息创建所有基函数的1阶导数，以幂基的形式表示*/
	std::vector<powerBaseFun> bf_basisFun1Diffs;
	for (auto bf_basisFun0Diff : bf_basisFun0Diffs)
	{
		powerBaseFun bf_basisFun1Diff;
		/*根据1阶导数与0阶导数的关系计算1阶导数中所有幂基的最小次数*/
		bf_basisFun1Diff.pbf_min_deg = (bf_basisFun0Diff.pbf_min_deg - 1) < 0 ? 0 : (bf_basisFun0Diff.pbf_min_deg - 1);
		/*根据1阶导数与0阶导数的关系计算1阶导数中所有幂基的最大次数*/
		bf_basisFun1Diff.pbf_max_deg = bf_basisFun0Diff.pbf_max_deg - 1;
		for (int ii = bf_basisFun1Diff.pbf_min_deg; ii <= bf_basisFun1Diff.pbf_max_deg; ii++)
		{
			/*1阶导数ii次幂基对应的系数计算*/
			int index = ii + 1 - bf_basisFun0Diff.pbf_min_deg;
			double pbf_coeff = bf_basisFun0Diff.pbf_coeffs[index] * (ii + 1);
			bf_basisFun1Diff.pbf_coeffs.push_back(pbf_coeff);
		}
		bf_basisFun1Diff.pbf_order = 1;
		bf_basisFun1Diffs.push_back(bf_basisFun1Diff);
		/*check bf_basisFun1Diff infor*/
		/*std::cout << "------Basis function 1 differentional------" << std::endl;
		std::cout << "powerBaseFun degree: ";
		for (int ii = bf_basisFun1Diff.pbf_min_deg; ii <= bf_basisFun1Diff.pbf_max_deg; ii++)
		{
			std::cout << " " << ii;
		}
		std::cout << std::endl;
		std::cout << "powerBaseFun coefficient: ";
		for (auto pbf_coeff : bf_basisFun1Diff.pbf_coeffs)
		{
			std::cout << " " << pbf_coeff;
		}
		std::cout << std::endl;*/
	}

	/*根据bf1阶导数的信息创建所有基函数的2阶导数，以幂基的形式表示*/
	std::vector<powerBaseFun> bf_basisFun2Diffs;
	for (auto bf_basisFun1Diff : bf_basisFun1Diffs)
	{
		powerBaseFun bf_basisFun2Diff;
		/*根据1阶导数与0阶导数的关系计算1阶导数中所有幂基的最小次数*/
		bf_basisFun2Diff.pbf_min_deg = (bf_basisFun1Diff.pbf_min_deg - 1) < 0 ? 0 : (bf_basisFun1Diff.pbf_min_deg - 1);
		/*根据1阶导数与0阶导数的关系计算1阶导数中所有幂基的最大次数*/
		bf_basisFun2Diff.pbf_max_deg = bf_basisFun1Diff.pbf_max_deg - 1;
		for (int ii = bf_basisFun2Diff.pbf_min_deg; ii <= bf_basisFun2Diff.pbf_max_deg; ii++)
		{
			/*1阶导数ii次幂基对应的系数计算*/
			int index = ii + 1 - bf_basisFun1Diff.pbf_min_deg;
			double pbf_coeff = bf_basisFun1Diff.pbf_coeffs[index] * (ii + 1);
			bf_basisFun2Diff.pbf_coeffs.push_back(pbf_coeff);
		}
		bf_basisFun2Diff.pbf_order = 2;
		bf_basisFun2Diffs.push_back(bf_basisFun2Diff);
		/*check bf_basisFun1Diff infor*/
		/*std::cout << "------Basis function 2 differentional------" << std::endl;
		std::cout << "powerBaseFun degree: ";
		for (int ii = bf_basisFun2Diff.pbf_min_deg; ii <= bf_basisFun2Diff.pbf_max_deg; ii++)
		{
			std::cout << " " << ii;
		}
		std::cout << std::endl;
		std::cout << "powerBaseFun coefficient: ";
		for (auto pbf_coeff : bf_basisFun2Diff.pbf_coeffs)
		{
			std::cout << " " << pbf_coeff;
		}
		std::cout << std::endl;*/
	}

	/*
	* 得到bernstein基函数的各阶导数后，就可以计算积分了
	* 求积分要分为两个方向，这里默认两个方向的阶数是一致的
	*/
	std::vector< std::vector<double>> total_ctl_integrals;
	/*check total_ctl_integrals的初始化*/
	/*for (auto total_ctl_integral : total_ctl_integrals)
	{
		std::cout << total_ctl_integral << std::endl;
	}*/
	/*对当前bezier曲面上的每个控制点分别求偏导*/
	for (int i = 0; i <= bf->degree(1); i++)
	{
		/*获取被积函数v方向的一部分（总共是两部分的乘积）*/
		powerBaseFun v2_part1_diff = bf_basisFun2Diffs[i];//v方向的二阶导
		powerBaseFun v1_part1_diff = bf_basisFun1Diffs[i];//v方向的一阶导
		powerBaseFun v0_part1_diff = bf_basisFun0Diffs[i];//v方向的0阶导
		for (int j = 0; j <= bf->degree(0); j++)
		{
			//std::cout << "----------------differential: " << j << "," << i << "----------------" << std::endl;
			/*ji对应在当前控制点的索引是"i*bf->degree(0)+j"*/
			/*
			* 对应在当前bezier曲面每个控制点上的积分系数，积分系数与每个控制点相乘得到的是薄板能量关于ji控制点的偏导结果
			* 注意这里的顺序要跟当前曲面片上控制点顺序保持一致(下面循环中的顺序为(jj,ii))
			*/
			/*获取被积函数u方向的一部分（总共是两部分的乘积）*/
			powerBaseFun u2_part1_diff = bf_basisFun2Diffs[j];//u方向的二阶导
			powerBaseFun u1_part1_diff = bf_basisFun1Diffs[j];//u方向的一阶导
			powerBaseFun u0_part1_diff = bf_basisFun0Diffs[j];//u方向的0阶导

			std::vector<double> ctl_integrals;
			/*计算每个控制点上的积分系数*/
			for (int ii = 0; ii <= bf->degree(1); ii++)
			{
				/*获取被积函数v方向的第二部分（总共是两部分的乘积）*/
				powerBaseFun v2_part2_diff = bf_basisFun2Diffs[ii];//v方向的二阶导
				powerBaseFun v1_part2_diff = bf_basisFun1Diffs[ii];//v方向的一阶导
				powerBaseFun v0_part2_diff = bf_basisFun0Diffs[ii];//v方向的0阶导
				for (int jj = 0; jj <= bf->degree(0); jj++)
				{
					/*初始化当前积分项系数*/
					double ctl_integral = 0.0;
					/*获取被积函数u方向的第二部分（总共是两部分的乘积）*/
					powerBaseFun u2_part2_diff = bf_basisFun2Diffs[jj];//u方向的二阶导
					powerBaseFun u1_part2_diff = bf_basisFun1Diffs[jj];//u方向的一阶导
					powerBaseFun u0_part2_diff = bf_basisFun0Diffs[jj];//u方向的0阶导
					/*积分系数分为三部分的和*/
					/*1. 关于u方向的二阶导积分*v方向的0阶导积分系数*/
					{
						//std::cout << "---u2v0---" << std::endl;
						/*声明两个方向被积函数*/
						powerBaseFun u2_diff, v0_diff;
						/*----------u2_diff的计算开始------------------*/
						/*当前被积函数的幂基表示的最小次数是两部分最小次数的和*/
						u2_diff.pbf_min_deg = u2_part1_diff.pbf_min_deg + u2_part2_diff.pbf_min_deg;
						/*当前被积函数的幂基表示的最大次数是两部分最大次数的和*/
						u2_diff.pbf_max_deg = u2_part1_diff.pbf_max_deg + u2_part2_diff.pbf_max_deg;
						/*根据获得的最小次数和最大次数来初始化此时被积函数系数*/
						for (int iii = u2_diff.pbf_min_deg; iii <= u2_diff.pbf_max_deg; iii++)
						{
							u2_diff.pbf_coeffs.push_back(0.0);
						}
						/*计算系数*/
						for (int iii = 0; iii < u2_part1_diff.pbf_coeffs.size(); iii++)
						{
							for (int jjj = 0; jjj < u2_part2_diff.pbf_coeffs.size(); jjj++)
							{
								/*第一部分第iii个幂基函数次数+第二部分第jjj个幂基函数次数=对应被积函数次数*/
								int u2_diff_degree = (iii + u2_part1_diff.pbf_min_deg) + (jjj + u2_part2_diff.pbf_min_deg);
								/*计算当前对应被积函数次数所贡献的系数*/
								double u2_diff_degree_coeff = u2_part1_diff.pbf_coeffs[iii] * u2_part2_diff.pbf_coeffs[jjj];
								/*将系数添加到被积函数系数对应的位置上*/
								u2_diff.pbf_coeffs[u2_diff_degree - u2_diff.pbf_min_deg] += u2_diff_degree_coeff;
							}
						}
						/*----------u2_diff的计算完成------------------*/
						/*check 被积函数u2_diff的计算*/
						/*std::cout << " integral function part 1: " << std::endl;
						u2_part1_diff.show();
						std::cout << " integral function part 2: " << std::endl;
						u2_part2_diff.show();
						std::cout << " integral function: " << std::endl;
						u2_diff.show();*/
						/*对被积函数u2_diff进行积分*/
						double u2_diff_integral = 0.0;
						for (int iii = 0; iii < u2_diff.pbf_coeffs.size(); iii++)
						{
							int coeff_degree = iii + u2_diff.pbf_min_deg;
							u2_diff_integral += (1.0 / (coeff_degree + 1) * u2_diff.pbf_coeffs[iii]);
						}

						/*----------v0_diff的计算开始------------------*/
						/*当前被积函数的幂基表示的最小次数是两部分最小次数的和*/
						v0_diff.pbf_min_deg = v0_part1_diff.pbf_min_deg + v0_part2_diff.pbf_min_deg;
						/*当前被积函数的幂基表示的最大次数是两部分最大次数的和*/
						v0_diff.pbf_max_deg = v0_part1_diff.pbf_max_deg + v0_part2_diff.pbf_max_deg;
						/*根据获得的最小次数和最大次数来初始化此时被积函数系数*/
						for (int iii = v0_diff.pbf_min_deg; iii <= v0_diff.pbf_max_deg; iii++)
						{
							v0_diff.pbf_coeffs.push_back(0.0);
						}
						/*计算系数*/
						for (int iii = 0; iii < v0_part1_diff.pbf_coeffs.size(); iii++)
						{
							for (int jjj = 0; jjj < v0_part2_diff.pbf_coeffs.size(); jjj++)
							{
								/*第一部分第iii个幂基函数次数+第二部分第jjj个幂基函数次数=对应被积函数次数*/
								int v0_diff_degree = (iii + v0_part1_diff.pbf_min_deg) + (jjj + v0_part2_diff.pbf_min_deg);
								/*计算当前对应被积函数次数所贡献的系数*/
								double v0_diff_degree_coeff = v0_part1_diff.pbf_coeffs[iii] * v0_part2_diff.pbf_coeffs[jjj];
								/*将系数添加到被积函数系数对应的位置上*/
								v0_diff.pbf_coeffs[v0_diff_degree - v0_diff.pbf_min_deg] += v0_diff_degree_coeff;
							}
						}
						/*----------v0_diff的计算完成------------------*/
						/*check 被积函数v0_diff的计算*/
						/*std::cout << " integral function part 1: " << std::endl;
						v0_part1_diff.show();
						std::cout << " integral function part 2: " << std::endl;
						v0_part2_diff.show();
						std::cout << " integral function: " << std::endl;
						v0_diff.show();*/

						/*对被积函数v0_diff进行积分*/
						double v0_diff_integral = 0.0;
						for (int iii = 0; iii < v0_diff.pbf_coeffs.size(); iii++)
						{
							int coeff_degree = iii + v0_diff.pbf_min_deg;
							v0_diff_integral += (1.0 / (coeff_degree + 1) * v0_diff.pbf_coeffs[iii]);
						}
						/*将两个方向的积分相乘作为当前控制点权重的一部分*/
						ctl_integral += (u2_diff_integral * v0_diff_integral);
					}
					/*2. 关于u方向的一阶导积分*v方向的1阶导积分系数*/
					{
						//std::cout << "---u1v1---" << std::endl;
						/*声明两个方向被积函数*/
						powerBaseFun u1_diff, v1_diff;
						/*----------u1_diff的计算开始------------------*/
						/*当前被积函数的幂基表示的最小次数是两部分最小次数的和*/
						u1_diff.pbf_min_deg = u1_part1_diff.pbf_min_deg + u1_part2_diff.pbf_min_deg;
						/*当前被积函数的幂基表示的最大次数是两部分最大次数的和*/
						u1_diff.pbf_max_deg = u1_part1_diff.pbf_max_deg + u1_part2_diff.pbf_max_deg;
						/*根据获得的最小次数和最大次数来初始化此时被积函数系数*/
						for (int iii = u1_diff.pbf_min_deg; iii <= u1_diff.pbf_max_deg; iii++)
						{
							u1_diff.pbf_coeffs.push_back(0.0);
						}
						/*计算系数*/
						for (int iii = 0; iii < u1_part1_diff.pbf_coeffs.size(); iii++)
						{
							for (int jjj = 0; jjj < u1_part2_diff.pbf_coeffs.size(); jjj++)
							{
								/*第一部分第iii个幂基函数次数+第二部分第jjj个幂基函数次数=对应被积函数次数*/
								int u1_diff_degree = (iii + u1_part1_diff.pbf_min_deg) + (jjj + u1_part2_diff.pbf_min_deg);
								/*计算当前对应被积函数次数所贡献的系数*/
								double u1_diff_degree_coeff = u1_part1_diff.pbf_coeffs[iii] * u1_part2_diff.pbf_coeffs[jjj];
								/*将系数添加到被积函数系数对应的位置上*/
								u1_diff.pbf_coeffs[u1_diff_degree - u1_diff.pbf_min_deg] += u1_diff_degree_coeff;
							}
						}
						/*----------u1_diff的计算完成------------------*/
						/*check 被积函数u1_diff的计算*/
						/*std::cout << " integral function part 1: " << std::endl;
						u1_part1_diff.show();
						std::cout << " integral function part 2: " << std::endl;
						u1_part2_diff.show();
						std::cout << " integral function: " << std::endl;
						u1_diff.show();*/
						/*对被积函数u2_diff进行积分*/
						double u1_diff_integral = 0.0;
						for (int iii = 0; iii < u1_diff.pbf_coeffs.size(); iii++)
						{
							int coeff_degree = iii + u1_diff.pbf_min_deg;
							u1_diff_integral += (1.0 / (coeff_degree + 1) * u1_diff.pbf_coeffs[iii]);
						}

						/*----------v1_diff的计算开始------------------*/
						/*当前被积函数的幂基表示的最小次数是两部分最小次数的和*/
						v1_diff.pbf_min_deg = v1_part1_diff.pbf_min_deg + v1_part2_diff.pbf_min_deg;
						/*当前被积函数的幂基表示的最大次数是两部分最大次数的和*/
						v1_diff.pbf_max_deg = v1_part1_diff.pbf_max_deg + v1_part2_diff.pbf_max_deg;
						/*根据获得的最小次数和最大次数来初始化此时被积函数系数*/
						for (int iii = v1_diff.pbf_min_deg; iii <= v1_diff.pbf_max_deg; iii++)
						{
							v1_diff.pbf_coeffs.push_back(0.0);
						}
						/*计算系数*/
						for (int iii = 0; iii < v1_part1_diff.pbf_coeffs.size(); iii++)
						{
							for (int jjj = 0; jjj < v1_part2_diff.pbf_coeffs.size(); jjj++)
							{
								/*第一部分第iii个幂基函数次数+第二部分第jjj个幂基函数次数=对应被积函数次数*/
								int v1_diff_degree = (iii + v1_part1_diff.pbf_min_deg) + (jjj + v1_part2_diff.pbf_min_deg);
								/*计算当前对应被积函数次数所贡献的系数*/
								double v1_diff_degree_coeff = v1_part1_diff.pbf_coeffs[iii] * v1_part2_diff.pbf_coeffs[jjj];
								/*将系数添加到被积函数系数对应的位置上*/
								v1_diff.pbf_coeffs[v1_diff_degree - v1_diff.pbf_min_deg] += v1_diff_degree_coeff;
							}
						}
						/*----------v1_diff的计算完成------------------*/
						/*check 被积函数v0_diff的计算*/
						/*std::cout << " integral function part 1: " << std::endl;
						v1_part1_diff.show();
						std::cout << " integral function part 2: " << std::endl;
						v1_part2_diff.show();
						std::cout << " integral function: " << std::endl;
						v1_diff.show();*/

						/*对被积函数v0_diff进行积分*/
						double v1_diff_integral = 0.0;
						for (int iii = 0; iii < v1_diff.pbf_coeffs.size(); iii++)
						{
							int coeff_degree = iii + v1_diff.pbf_min_deg;
							v1_diff_integral += (1.0 / (coeff_degree + 1) * v1_diff.pbf_coeffs[iii]);
						}
						/*将两个方向的积分相乘作为当前控制点权重的一部分*/
						ctl_integral += (u1_diff_integral * v1_diff_integral);
					}
					/*3. 关于u方向的0阶导积分*v方向的二阶导积分系数*/
					{
						//std::cout << "---u2v0---" << std::endl;
						/*声明两个方向被积函数*/
						powerBaseFun u0_diff, v2_diff;
						/*----------u2_diff的计算开始------------------*/
						/*当前被积函数的幂基表示的最小次数是两部分最小次数的和*/
						u0_diff.pbf_min_deg = u0_part1_diff.pbf_min_deg + u0_part2_diff.pbf_min_deg;
						/*当前被积函数的幂基表示的最大次数是两部分最大次数的和*/
						u0_diff.pbf_max_deg = u0_part1_diff.pbf_max_deg + u0_part2_diff.pbf_max_deg;
						/*根据获得的最小次数和最大次数来初始化此时被积函数系数*/
						for (int iii = u0_diff.pbf_min_deg; iii <= u0_diff.pbf_max_deg; iii++)
						{
							u0_diff.pbf_coeffs.push_back(0.0);
						}
						/*计算系数*/
						for (int iii = 0; iii < u0_part1_diff.pbf_coeffs.size(); iii++)
						{
							for (int jjj = 0; jjj < u0_part2_diff.pbf_coeffs.size(); jjj++)
							{
								/*第一部分第iii个幂基函数次数+第二部分第jjj个幂基函数次数=对应被积函数次数*/
								int u0_diff_degree = (iii + u0_part1_diff.pbf_min_deg) + (jjj + u0_part2_diff.pbf_min_deg);
								/*计算当前对应被积函数次数所贡献的系数*/
								double u0_diff_degree_coeff = u0_part1_diff.pbf_coeffs[iii] * u0_part2_diff.pbf_coeffs[jjj];
								/*将系数添加到被积函数系数对应的位置上*/
								u0_diff.pbf_coeffs[u0_diff_degree - u0_diff.pbf_min_deg] += u0_diff_degree_coeff;
							}
						}
						/*----------u0_diff的计算完成------------------*/
						/*check 被积函数u0_diff的计算*/
						/*std::cout << " integral function part 1: " << std::endl;
						u0_part1_diff.show();
						std::cout << " integral function part 2: " << std::endl;
						u0_part2_diff.show();
						std::cout << " integral function: " << std::endl;
						u0_diff.show();*/
						/*对被积函数u0_diff进行积分*/
						double u0_diff_integral = 0.0;
						for (int iii = 0; iii < u0_diff.pbf_coeffs.size(); iii++)
						{
							int coeff_degree = iii + u0_diff.pbf_min_deg;
							u0_diff_integral += (1.0 / (coeff_degree + 1) * u0_diff.pbf_coeffs[iii]);
						}

						/*----------v2_diff的计算开始------------------*/
						/*当前被积函数的幂基表示的最小次数是两部分最小次数的和*/
						v2_diff.pbf_min_deg = v2_part1_diff.pbf_min_deg + v2_part2_diff.pbf_min_deg;
						/*当前被积函数的幂基表示的最大次数是两部分最大次数的和*/
						v2_diff.pbf_max_deg = v2_part1_diff.pbf_max_deg + v2_part2_diff.pbf_max_deg;
						/*根据获得的最小次数和最大次数来初始化此时被积函数系数*/
						for (int iii = v2_diff.pbf_min_deg; iii <= v2_diff.pbf_max_deg; iii++)
						{
							v2_diff.pbf_coeffs.push_back(0.0);
						}
						/*计算系数*/
						for (int iii = 0; iii < v2_part1_diff.pbf_coeffs.size(); iii++)
						{
							for (int jjj = 0; jjj < v2_part2_diff.pbf_coeffs.size(); jjj++)
							{
								/*第一部分第iii个幂基函数次数+第二部分第jjj个幂基函数次数=对应被积函数次数*/
								int v2_diff_degree = (iii + v2_part1_diff.pbf_min_deg) + (jjj + v2_part2_diff.pbf_min_deg);
								/*计算当前对应被积函数次数所贡献的系数*/
								double v2_diff_degree_coeff = v2_part1_diff.pbf_coeffs[iii] * v2_part2_diff.pbf_coeffs[jjj];
								/*将系数添加到被积函数系数对应的位置上*/
								v2_diff.pbf_coeffs[v2_diff_degree - v2_diff.pbf_min_deg] += v2_diff_degree_coeff;
							}
						}
						/*----------v2_diff的计算完成------------------*/
						/*check 被积函数v2_diff的计算*/
						/*std::cout << " integral function part 1: " << std::endl;
						v2_part1_diff.show();
						std::cout << " integral function part 2: " << std::endl;
						v2_part2_diff.show();
						std::cout << " integral function: " << std::endl;
						v2_diff.show();*/

						/*对被积函数v2_diff进行积分*/
						double v2_diff_integral = 0.0;
						for (int iii = 0; iii < v2_diff.pbf_coeffs.size(); iii++)
						{
							int coeff_degree = iii + v2_diff.pbf_min_deg;
							v2_diff_integral += (1.0 / (coeff_degree + 1) * v2_diff.pbf_coeffs[iii]);
						}
						/*将两个方向的积分相乘作为当前控制点权重的一部分*/
						ctl_integral += (u0_diff_integral * v2_diff_integral);
					}
					/*注意这里的积分系数都乘以了2.0*/
					ctl_integrals.push_back(ctl_integral * 2.0);
					//std::cout << "--------------weight: " << jj << "," << ii << "-----------------" << ii * (bf->degree(0)+1) + jj << std::endl;
				}
			}
			/*check 积分项系数*/
			/*std::cout << "------ctl " << j << "," << i << " coeff------" << std::endl;
			for (auto ctl_integral : ctl_integrals)
			{
				std::cout << ctl_integral << std::endl;
			}
			std::cout << "-----------------------------------------------------------------------------------------" << std::endl;*/
			total_ctl_integrals.push_back(ctl_integrals);
		}
	}

	/*获取当前bezier曲面片上薄板能量关于一些网格点偏导在网格点上的权重*/
	/*1.找出关联的网格点id*/
	std::vector<int> bf_link_meshVIds;
	for (auto bf_ctl_link : bf->linkCptPWs())
	{
		for (auto bf_ctl_link_ws : bf_ctl_link.linkingCtlPId_weights)
		{
			int bf_link_meshVId = bf_ctl_link_ws.first;
			bool find_bf_link_meshVId = false;
			for (auto t_bf_link_meshVId : bf_link_meshVIds)
			{
				if (t_bf_link_meshVId == bf_link_meshVId)
				{
					find_bf_link_meshVId = true;
					break;
				}
			}
			if (!find_bf_link_meshVId)
			{
				bf_link_meshVIds.push_back(bf_link_meshVId);
			}
		}

	}
	/*check 找到关联的控制点个数*/
	//std::cout << "bf_link_meshVIds.size(): " << bf_link_meshVIds.size() << std::endl;
	/*2.计算每个网格点的偏导关于控制点上的权重*/
	std::vector<std::vector<double>> bf_link_meshVId_weights(bf_link_meshVIds.size(), std::vector<double>(bf->cpts().size(), 0.0));
	for (int j = 0; j < bf_link_meshVIds.size(); j++)
	{
		int bf_link_meshVId = bf_link_meshVIds[j];
		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			for (auto bf_ctl_link_ws : bf->linkCptPWs()[i].linkingCtlPId_weights)
			{
				int t_bf_link_meshVId = bf_ctl_link_ws.first;
				if (t_bf_link_meshVId == bf_link_meshVId)
				{
					double t_bf_link_meshVId_w = bf_ctl_link_ws.second;
					for (int ii = 0; ii < total_ctl_integrals[i].size(); ii++)
					{
						bf_link_meshVId_weights[j][ii] += (total_ctl_integrals[i][ii] * t_bf_link_meshVId_w);
					}
				}
			}
		}
	}
	/*check 每个网格点的偏导关于控制点上的权重*/
	/*std::cout << "bf_link_meshVId_weights.size(): " << bf_link_meshVId_weights.size() << std::endl;
	for (auto bf_link_meshVId_weight : bf_link_meshVId_weights)
	{
		for (auto bf_link_meshVId_w : bf_link_meshVId_weight)
		{
			std::cout << bf_link_meshVId_w << std::endl;
		}
		std::cout << "---------------------------------------------------------------" << std::endl;
	}*/
	/*3.计算每个网格点的偏导关于网格点上的权重*/
	for (int j = 0; j < bf_link_meshVIds.size(); j++)
	{
		BE_linkingCtlPoints tempV;
		for (int i = 0; i < bf->linkCptPWs().size(); i++)
		{
			for (auto cwIt : bf->linkCptPWs()[i].linkingCtlPId_weights)
			{
				std::pair<int, double> tempCw = cwIt;
				tempCw.second *= bf_link_meshVId_weights[j][i];
				if (tempV.linkingCtlPId_weights.find(tempCw.first) == tempV.linkingCtlPId_weights.end())
				{
					tempV.linkingCtlPId_weights.insert(tempCw);
				}
				else
				{
					std::map<int, double>::iterator tIt = tempV.linkingCtlPId_weights.find(tempCw.first);
					tIt->second += tempCw.second;
				}
			}
		}
		std::pair<int, BE_linkingCtlPoints> tempV_ws(bf_link_meshVIds[j], tempV);
		ctl_ws.insert(tempV_ws);
	}
	/*check ctl_ws*/
	/*for (auto ctl_w : ctl_ws)
	{
		std::cout << "------------" << ctl_w.first << "--------------" << std::endl;
		for (auto ctl_w_t : ctl_w.second.linkingCtlPId_weights)
		{
			std::cout <<"index: " << ctl_w_t.first <<" weight: "<<ctl_w_t.second << std::endl;
		}
	}*/
}

void CCG_QMSLib::CCG_QMS_model::approximateSurface_optimize_check(M* pMesh, double approxi_thread, double approxi_smooth_weight)
{
	using namespace Eigen;

	// ============================================================
	// 1. 预处理：建立内部控制点索引
	// ============================================================
	// Map: 原始控制点ID (0-based) -> 矩阵列索引 (0-based)
	std::map<int, int> innerCtlPtsIndexes;
	std::vector<int> boundaryVs;
	int countInnerCtlPts = 0;
	int nCtrlPts = (int)this->controlPoints().size();

	for (int i = 0; i < nCtrlPts; i++) {
		// pMesh->idVertex 通常是 1-based
		V* v = pMesh->idVertex(i + 1);
		if (v->boundary() || v->feature() || v->fixed()) {
			boundaryVs.push_back(i);
		}
		else {
			innerCtlPtsIndexes.insert(std::pair<int, int>(i, countInnerCtlPts++));
		}
	}

	int nUnknowns = countInnerCtlPts; // 待求解的变量数
	// 约束数量 = 采样点数 + 平滑约束数(控制点数)
	int nConstraints = (int)this->samplings().size() + nCtrlPts;

	// ============================================================
	// 2. 构建稀疏矩阵 A (只构建一次，XYZ通用)
	// ============================================================
	std::vector<Triplet<double>> tripletlist;
	// 存储边界点（已知量）对右端项的影响
	std::vector<Triplet<double>> tripletlist_boundarys;

	// 预分配内存以加速
	tripletlist.reserve(this->samplings().size() * 4 + nCtrlPts * 5);

	// --- (A) 采样点约束 ---
	for (int i = 0; i < this->samplings().size(); i++) {
		if (!(this->samplings()[i]->mark())) continue;

		std::vector<int> index;
		std::vector<double> weight;
		// 假设此处为串行，避免多线程内存竞争
		computeRealSamplingPos(this->samplings()[i]->fId(), this->samplings()[i]->vId(), this->samplings()[i]->uv(), index, weight);

		for (size_t j = 0; j < index.size(); j++) {
			int vIdx_0based = index[j] - 1;
			V* v = pMesh->idVertex(index[j]); // idVertex is 1-based

			if (v->boundary() || v->feature() || v->fixed()) {
				// 已知项：放入 boundary 列表，后续移至等式右边
				tripletlist_boundarys.push_back(Triplet<double>(i, vIdx_0based, weight[j]));
			}
			else {
				// 未知项：放入 A 矩阵
				tripletlist.push_back(Triplet<double>(i, innerCtlPtsIndexes.at(vIdx_0based), weight[j]));
			}
		}
	}

	// --- (B) Thin-Plate 平滑能量约束 ---
	double E_thinPlate_weight = approxi_smooth_weight;
	int startRow = (int)this->samplings().size(); // 平滑约束从这一行开始

	for (auto bf : this->bfs()) {
		std::map<int, BE_linkingCtlPoints> ctl_ws;
		computeBezierSufIntegral_mine(bf, ctl_ws);

		for (auto ctl_w : ctl_ws) {
			for (auto ctl_w_t : ctl_w.second.linkingCtlPId_weights) {
				// 行号：startRow + ctl_w.first (偏移)
				int row = startRow + ctl_w.first - 1;
				int col_0based = ctl_w_t.first - 1;
				double val = ctl_w_t.second * E_thinPlate_weight;

				V* v = pMesh->idVertex(col_0based + 1);

				if (v->boundary() || v->feature() || v->fixed()) {
					tripletlist_boundarys.push_back(Triplet<double>(row, col_0based, val));
				}
				else {
					tripletlist.push_back(Triplet<double>(row, innerCtlPtsIndexes.at(col_0based), val));
				}
			}
		}
	}

	// ============================================================
	// 3. 矩阵组装与分析 (最耗时步骤之一，只需做一次)
	// ============================================================
	SparseMatrix<double> A(nConstraints, nUnknowns);
	A.setFromTriplets(tripletlist.begin(), tripletlist.end());

	// 初始化 LSCG 求解器
	LeastSquaresConjugateGradient<SparseMatrix<double>> lscg;
	lscg.setTolerance(approxi_thread);
	// lscg.setMaxIterations(200); // 可选：限制迭代次数以保证速度
	lscg.compute(A); // 分析矩阵结构，分配内存

	if (lscg.info() != Success) {
		// std::cerr << "Matrix analysis failed!" << std::endl;
		return;
	}

	// ============================================================
	// 4. 预计算边界条件对 RHS 的贡献
	// ============================================================
	// 边界点是固定的，它们对 Ax=b 的贡献可以移到右边： A_inner * x = b - A_boundary * x_fixed
	std::vector<VectorXd> rhs_boundary_fix(3);
	for (int d = 0; d < 3; d++) rhs_boundary_fix[d] = VectorXd::Zero(nConstraints);

	for (const auto& t : tripletlist_boundarys) {
		int vIdx_1based = t.col() + 1; // Triplet col 存的是 0-based index
		V* v = pMesh->idVertex(vIdx_1based);
		for (int d = 0; d < 3; d++) {
			rhs_boundary_fix[d][t.row()] += v->point()[d] * t.value();
		}
	}

	// ============================================================
	// 5. 并行求解 XYZ 三个维度 (核心优化)
	// ============================================================
	std::vector<VectorXd> results(3);

	// 开启 OpenMP 并行，同时求解 X, Y, Z
#pragma omp parallel for
	for (int dim = 0; dim < 3; dim++)
	{
		// --- 5.1 构建右端项 b ---
		VectorXd b = VectorXd::Zero(nConstraints);

		// 采样点部分：目标是 pos[dim]
		for (int j = 0; j < this->samplings().size(); j++) {
			b[j] = this->samplings()[j]->pos()[dim] - rhs_boundary_fix[dim][j];
		}
		// 平滑部分：目标是 0 (最小化能量)
		for (int j = 0; j < nCtrlPts; j++) {
			b[startRow + j] = 0.0 - rhs_boundary_fix[dim][startRow + j];
		}

		// --- 5.2 构建初始猜测 x0 (Warm Start) ---
		// 使用当前的控制点位置作为初值
		VectorXd x0(nUnknowns);

		// 【修改部分】兼容 C++11/C++14 的写法
		for (auto const& entry : innerCtlPtsIndexes) {
			int originalIdx = entry.first;  // Map 的 Key
			int solverIdx = entry.second;   // Map 的 Value

			x0[solverIdx] = this->controlPoints()[originalIdx][dim];
		}

		// --- 5.3 带初值的求解 ---
		// solveWithGuess 会利用 x0 加速收敛
		results[dim] = lscg.solveWithGuess(b, x0);

		// Debug:
		// if (lscg.info() != Success) { /* handle error */ }
	}

	// ============================================================
	// 6. 更新控制点坐标
	// ============================================================
	for (int i = 0; i < nCtrlPts; i++) {
		// 只有内部非固定点才更新
		if (innerCtlPtsIndexes.find(i) != innerCtlPtsIndexes.end()) {
			int idx = innerCtlPtsIndexes[i];
			for (int dim = 0; dim < 3; dim++) {
				this->controlPoints()[i][dim] = results[dim][idx];
			}
		}
		else {
			// 边界/固定点保持原样 (这部分其实已经在 computeRealSamplingPos 中处理了，
			// 但为了保险起见，或者如果需要显式重置为 mesh 坐标，可在此处处理)
			V* v = pMesh->idVertex(i + 1);
			for (int dim = 0; dim < 3; dim++) {
				this->controlPoints()[i][dim] = v->point()[dim];
			}
		}
	}
}

void CCG_QMSLib::CCG_QMS_model::approximateSurface_optimize_check2(M* pMesh, double approxi_thread, double approxapproxi_smooth_weight)
{
	using namespace Eigen;
	/*
	* 根据采样点列方程
	* 这里采用稀疏矩阵的形式存储
	* 稀疏矩阵的维数为 （采样点个数+网格点个数）*网格点个数
	*/
	SparseMatrix<double> A(this->samplings().size() + this->controlPoints().size(), this->controlPoints().size());
	std::vector<Eigen::Triplet<double>> tripletlist;
	/*获取采样点对应的稀疏矩阵*/
	for (int i = 0; i < this->samplings().size(); i++)
	{
		if (!(this->samplings()[i]->mark()))
		{
			//std::cout << "-------bad sampling!" << std::endl;
			continue;
		}
		//std::vector<double> uv = { this->samplings()[i]->uv()[0],this->samplings()[i]->uv()[1] };
		std::vector<int> index;
		std::vector<double> weight;
		//std::cout << "-----------sampling pts " << i << "---------------------" << std::endl;
		computeRealSamplingPos(this->samplings()[i]->fId(), this->samplings()[i]->vId(), this->samplings()[i]->uv(), index, weight);
		for (int j = 0; j < index.size(); j++)
		{
			tripletlist.push_back(Triplet<double>(i, index[j] - 1, weight[j]));
		}
	}
	//A.setFromTriplets(tripletlist.begin(), tripletlist.end());

	/*加上thinplate energy*/
	/*薄板能量对每个网格点求偏导，让偏导等于0*/
	int count = this->samplings().size();
	/*存放薄板能量关于所有网格点的偏导*/
	std::vector<Triplet<double>> E_thinPlates;
	double E_thinPlates_weight = approxapproxi_smooth_weight;
	/*遍历所有Bezier面片*/
	for (auto bf : this->bfs())
	{
		/*获得与当前bezier曲面积分的计算有关的初始四边形网格点的id,以及薄板能量关于网格点偏导对应在一系列网格点的权重*/
		std::map<int, BE_linkingCtlPoints> ctl_ws;
		computeBezierSufIntegral_mine(bf, ctl_ws);
		/*check ctl_ws*/
		/*for (auto ctl_w : ctl_ws)
		{
			std::cout << "------------" << ctl_w.first << "--------------" << std::endl;
			for (auto ctl_w_t : ctl_w.second.linkingCtlPId_weights)
			{
				std::cout << "index: " << ctl_w_t.first << " weight: " << ctl_w_t.second << std::endl;
			}
		}*/
		for (auto ctl_w : ctl_ws)
		{
			for (auto ctl_w_t : ctl_w.second.linkingCtlPId_weights)
			{
				E_thinPlates.push_back(Triplet<double>(count + ctl_w.first - 1, ctl_w_t.first - 1, ctl_w_t.second * E_thinPlates_weight));
			}
		}
	}
	/*将薄板能量约束加到最小二乘约束的矩阵中*/
	for (auto e_thinPlate : E_thinPlates)
	{
		tripletlist.push_back(e_thinPlate);
	}
	/*从三元组中获得方程组的系数矩阵，是个稀疏矩阵*/
	A.setFromTriplets(tripletlist.begin(), tripletlist.end());

	/*三个分量分别进行求解*/
	for (int i = 0; i < 3; i++)
	{
		VectorXd b = VectorXd::Constant(this->samplings().size() + this->controlPoints().size(), 0.0);
		for (int j = 0; j < this->samplings().size(); j++)
		{
			b[j] = this->samplings()[j]->pos()[i];
		}

		/*求解*/
		A.makeCompressed();
		LeastSquaresConjugateGradient<SparseMatrix<double>> Solver_sparse;
		Solver_sparse.setTolerance(approxi_thread);
		Solver_sparse.compute(A);
		VectorXd r = Solver_sparse.solve(b);

		for (int j = 0; j < r.size(); j++)
		{
			this->controlPoints()[j][i] = r[j];
		}

	}
}

/*-------------------function for CCG_QMS_model end------------------------------------*/


/*-------------------function for CCG_QMS_bezierSurf start------------------------------------*/
void CCG_QMSLib::CCG_QMS_bezierSurf::raiseDegree3to5(M* pMesh)
{
	/*
	* The definition of the matrix needs to be correlated with the distance between the corresponding knots.
	* For now, it is assumed that the distance between the knots is 1.0.
	*/
	/*
	* Define the matrices corresponding to the third to fourth order increments
	*/
	Eigen::MatrixXd raiseMatrix3To4(5, 4);
	raiseMatrix3To4 << 1, 0, 0, 0,
		1.0 / 4, 3.0 / 4, 0, 0,
		0, 1.0 / 2, 1.0 / 2, 0,
		0, 0, 3.0 / 4, 1.0 / 4,
		0, 0, 0, 1;
	//std::cout << "raiseMatrix3To4:\n" << raiseMatrix3To4 << std::endl;

	/*
	* Define the matrices corresponding to the fourth to fifth order increments
	*/
	Eigen::MatrixXd raiseMatrix4To5(6, 5);
	raiseMatrix4To5 << 1, 0, 0, 0, 0,
		1.0 / 5, 4.0 / 5, 0, 0, 0,
		0, 2.0 / 5, 3.0 / 5, 0, 0,
		0, 0, 3.0 / 5, 2.0 / 5, 0,
		0, 0, 0, 4.0 / 5, 1.0 / 5,
		0, 0, 0, 0, 1;
	//std::cout << "raiseMatrix4To5:\n" << raiseMatrix4To5 << std::endl;

	/*Define the matrices corresponding to the third to fifth order increments*/
	Eigen::MatrixXd raiseMatrix3To5 = raiseMatrix4To5 * raiseMatrix3To4;
	//std::cout << "raiseMatrix3To5:\n" << raiseMatrix3To5 << std::endl;

	/*Update knot vector*/
	//...
	//std::cout << "be_cpts: " << std::endl;
	//outputCpt(be_cpts);

	/*-----------------------------------------------------------------------------------------------------------*/
	/*U then V*/
	//Update control points in one direction
	std::vector<MeshLib::CPoint> n_cpt1(24);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				n_cpt1[i * 6 + j] += (m_cpts[i + k * 4] * raiseMatrix3To5(j, k));
			}
		}
	}
	//std::cout << "n_cpt1: " << std::endl;
	//outputCpt(n_cpt1);
	/*Update mapping between control points and quad mesh vertex in one direction*/
	std::vector<BE_linkingCtlPoints> n_cw1(24);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			BE_linkingCtlPoints tempCW;
			for (int k = 0; k < 4; k++)
			{
				//n_cpt1[i * 6 + j] += (be_cpts[i + k * 4] * raiseMatrix3To5(j, k));
				for (auto cwIt : m_linkCtlPWs[i + k * 4].linkingCtlPId_weights)
				{
					std::pair<int, double> tempV = cwIt;
					tempV.second *= (raiseMatrix3To5(j, k));
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
			}
			n_cw1[i * 6 + j] = tempCW;
		}
	}
	//std::cout << "n_cw1: " << std::endl;
	//outputCpt(pMesh,n_cw1);
	//std::cout << "111" << std::endl;
	//Update control points in another direction
	std::vector<MeshLib::CPoint> n_cpt2(36);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				n_cpt2[i * 6 + j] += (n_cpt1[k * 6 + i] * raiseMatrix3To5(j, k));
			}
		}
	}
	//std::cout << "n_cpt2: " << std::endl;
	//outputCpt(n_cpt2);
	/*Update mapping between control points and quad mesh vertex in another direction*/
	std::vector<BE_linkingCtlPoints> n_cw2(36);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			BE_linkingCtlPoints tempCW;
			for (int k = 0; k < 4; k++)
			{
				for (auto cwIt : n_cw1[k * 6 + i].linkingCtlPId_weights)
				{
					std::pair<int, double> tempV = cwIt;
					tempV.second *= (raiseMatrix3To5(j, k));
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
			}
			n_cw2[i * 6 + j] = tempCW;
		}
	}
	/*-----------------------------------------------------------------------------------------------------------*/

	//std::cout << "222" << std::endl;
	/*Update control points*/
	m_cpts.resize(36);
	m_cpts = n_cpt2;
	//std::cout << "be_cpts: " << std::endl;
	//outputCpt(be_cpts);

	/*Update mapping control points and quad mesh vertex*/
	m_linkCtlPWs.resize(36);
	m_linkCtlPWs = n_cw2;
	//std::cout << "n_cw2: " << std::endl;
	//outputCpt(pMesh,n_cw2);

	m_singular_linkCtlPWs.resize(36);
}

/*-------------------function for CCG_QMS_bezierSurf end------------------------------------*/

/*-------------------function for CCG_QMS_bsplineSurf start------------------------------------*/
void CCG_QMSLib::CCG_QMS_bsplineSurf::RemoveCurveKnot_UVConsistency2(double error_threshold)
{
	//error_threshold = getTOL(error_threshold);
	/*std::cout << "Uqbb_cpts:" << std::endl;
	for (int ii = 0; ii < qbb_cpts.size(); ii++)
	{
		for (int jj = 0; jj < qbb_cpts[0].size(); jj++) show(qbb_cpts[ii][jj]);
		std::cout << std::endl;
	}*/
	//U direction
	std::vector<double> qbb_knotsU_temp = m_knotsU;
	//std::vector<std::vector<CPoint>> qbb_cpts_temp = qbb_cpts;

	/*cfw 2025/3/1 modified, At this point, the dimension of the first dimension of the two-dimensional control points is the number of control points in the U direction.*/
	std::vector<std::vector<CPoint>> qbb_cptsU;
	for (int i = 0; i < m_cpts[0].size(); i++)
	{
		std::vector<CPoint> qbb_tempV;
		for (int j = 0; j < m_cpts.size(); j++)
		{
			qbb_tempV.push_back(m_cpts[j][i]);
		}
		qbb_cptsU.push_back(qbb_tempV);
	}
	std::vector<std::vector<CPoint>> qbb_cpts_temp = qbb_cptsU;
	std::vector<std::vector<double>> knot_and_redundancy;//knot index and its repetition degree
	//compute knot index and its repetition degree
	for (int aa = 0, bb = 0; aa < m_knotsU.size() - m_degree[0] - 1; aa++)
	{
		int count = 0;
		if (m_knotsU[aa] != m_knotsU[0])
		{
			while (m_knotsU[aa] == m_knotsU[aa + 1])
			{
				count++;
				aa++;
			}
			count++;
			if (bb >= knot_and_redundancy.size())
			{
				knot_and_redundancy.resize(bb + 1);
			}
			knot_and_redundancy[bb].push_back(aa);
			knot_and_redundancy[bb].push_back(count);
			bb++;
		}
	}

	int remove_count = 0;
	for (int index1 = 0; index1 < knot_and_redundancy.size(); index1++)//Remove the repetitive knot in U, multiplicity is knot_and_redundancy.size()
	{
		for (int index2 = 0; index2 < knot_and_redundancy[index1][1]; index2++)
		{
			int flag = 1;
			for (int index = 0; index < qbb_cpts_temp.size(); index++)//Remove the knots along the U direction corresponding to each V.
			{
				int r = knot_and_redundancy[index1][0] - remove_count;  //r is the index of the deleted knot
				int s = knot_and_redundancy[index1][1] - index2;  //s is the multiplicity now
				double u = qbb_knotsU_temp[r];
				/* 1.order p
				* 2. knot vector array U
				* 3. Control point array Pw
				* 4. Delete knot u with duplicate degree of s and their indices r
				* 5. Delete times num */
				int i; i = r - m_degree[0];  //Here we assume that all nodes can be deleted.
				int j; j = r - s;
				int x = 0; //When removing a knot, it is used to record the index of the newly calculated temporary control point array, facilitating the use of the corresponding new control points during the iteration process.
				int removeflag = 0;
				std::vector<MeshLib::CPoint> newpoint_i;
				std::vector<MeshLib::CPoint> newpoint_j;
				std::vector<MeshLib::CPoint> newpoint_i_j;
				newpoint_i.push_back(qbb_cpts_temp[index][i - 1]);
				newpoint_j.push_back(qbb_cpts_temp[index][j + 1]);
				while (j - i > 0)
				{
					double quan_i = (u - qbb_knotsU_temp[i]) / (qbb_knotsU_temp[i + m_degree[0] + 1] - qbb_knotsU_temp[i]);
					double quan_j = (u - qbb_knotsU_temp[j]) / (qbb_knotsU_temp[j + m_degree[0] + 1] - qbb_knotsU_temp[j]);
					MeshLib::CPoint temp1 = qbb_cpts_temp[index][i];
					MeshLib::CPoint temp2 = newpoint_i[x];
					MeshLib::CPoint temp3 = qbb_cpts_temp[index][j];
					MeshLib::CPoint temp4 = newpoint_j[x];
					newpoint_i.push_back((temp1 - temp2 * (1 - quan_i)) / quan_i);
					newpoint_j.push_back((temp3 - temp4 * quan_j) / (1 - quan_j));

					/*if (!(((temp3 - temp4 * quan_j) / (1 - quan_j) - (temp1 - temp2 * (1 - quan_i)) / quan_i).norm() < error_threshold))
					{
						newpoint_i_j.push_back((temp3 - temp4 * quan_j) / (1 - quan_j));
					}*/
					if (j - i >= 2)
					{
						newpoint_i_j.insert(newpoint_i_j.begin() + x, (temp1 - temp2 * (1 - quan_i)) / quan_i);
						/*newpoint_i_j.push_back((temp3 - temp4 * quan_j) / (1 - quan_j));*/
						newpoint_i_j.insert(newpoint_i_j.begin() + x + 1, (temp3 - temp4 * quan_j) / (1 - quan_j));
					}
					else if (j - i == 1)
					{
						newpoint_i_j.insert(newpoint_i_j.begin() + x, (((temp1 - temp2 * (1 - quan_i)) / quan_i) + (temp3 - temp4 * quan_j) / (1 - quan_j)) / 2);
					}
					x++;
					i++;
					j--;
				}
				if (j - i < 0)
				{
					//std::cout << "case 1 computation error:" << (newpoint_i[x] - newpoint_j[x]).norm() << std::endl;
					if ((newpoint_i[x] - newpoint_j[x]).norm() < error_threshold)//allowable error
					{
						/*std::cout << "can be removed by case 1" << std::endl;*/
						removeflag = 1;
					}
				}
				else if (j - i == 0) {
					MeshLib::CPoint temp = qbb_cpts_temp[index][(i + j) / 2];
					double quan_ij = (u - qbb_knotsU_temp[i]) / (qbb_knotsU_temp[i + m_degree[0] + 1] - qbb_knotsU_temp[i]);
					//std::cout << "case 2 computation error:" << (temp - newpoint_j[x] * quan_ij - newpoint_i[x] * (1 - quan_ij)).norm() << std::endl;
					if ((temp - newpoint_j[x] * quan_ij - newpoint_i[x] * (1 - quan_ij)).norm() < error_threshold)//allowable error
					{
						/*std::cout << "can be removed by case 2" << std::endl;*/
						removeflag = 1;
					}
				}
				if (removeflag == 1)
				{
					int insert_position = r - m_degree[0];
					int erase_start = r - m_degree[0];
					int erase_end = r - s + 1;
					qbb_cpts_temp[index].erase(qbb_cpts_temp[index].begin() + erase_start, qbb_cpts_temp[index].begin() + erase_end);
					qbb_cpts_temp[index].insert(qbb_cpts_temp[index].begin() + insert_position, newpoint_i_j.begin(), newpoint_i_j.end());
				}
				else {
					//std::cout << "can't be removed" << std::endl;
					flag = 0;
					break;
				}
			}
			if (flag == 0)
			{
				qbb_knotsU_temp = m_knotsU;
				qbb_cpts_temp = qbb_cptsU;
				break;
			}
			else if (flag == 1)
			{
				int r = knot_and_redundancy[index1][0] - remove_count;
				qbb_knotsU_temp.erase(qbb_knotsU_temp.begin() + r);
				m_knotsU = qbb_knotsU_temp;

				qbb_cptsU = qbb_cpts_temp;
				/*for (auto ku : qbb_cptsU)
				{
					for (int i = 0;i < ku.size();i++)
					{
						std::cout << ku[i][0] << " " << ku[i][1] << " " << ku[i][2] << std::endl;
					}
				}*/
				/*for (int i = 0;i < qbb_cptsU[0].size();i++)
				{
					std::cout << qbb_cptsU[0][i][0] << " " << qbb_cptsU[0][i][1] << " " << qbb_cptsU[0][i][2] << std::endl;
				}*/
				remove_count++;
			}
		}
	}

	//transposition
	std::vector<std::vector<CPoint>> qbb_cpts_temp1;
	for (int i = 0; i < qbb_cptsU[0].size(); i++)
	{
		std::vector<CPoint> qbb_temp1;
		for (int j = 0; j < qbb_cptsU.size(); j++)
		{
			qbb_temp1.push_back(qbb_cptsU[j][i]);
		}
		qbb_cpts_temp1.push_back(qbb_temp1);
	}
	m_cpts = qbb_cpts_temp1;


	/*V direction*/
	//Transposing a two-dimensional array is the same as handling the U-direction method.
	//Transpose a two-dimensional array
	/*std::vector<std::vector<CPoint>> qbb_cptsV;
	for (int i = 0; i < qbb_cpts[0].size(); i++)
	{
		std::vector<CPoint> qbb_tempV;
		for (int j = 0; j < qbb_cpts.size(); j++)
		{
			qbb_tempV.push_back(qbb_cpts[j][i]);
		}
		qbb_cptsV.push_back(qbb_tempV);
	}*/
	std::vector<std::vector<CPoint>> qbb_cptsV = m_cpts;
	std::vector<double> qbb_knotsV_temp = m_knotsV;
	//std::vector<std::vector<CPoint>> qbb_cptsV_temp = qbb_cptsV;
	/*cfw 2025/3/1 modified,At this point, the dimension of the first dimension of the two-dimensional control points is the number of control points in the U direction.*/
	std::vector<std::vector<CPoint>> qbb_cptsV_temp = m_cpts;

	knot_and_redundancy.clear();//knot index and its repetition degree
	for (int aa = 0, bb = 0; aa < m_knotsV.size() - m_degree[1] - 1; aa++)
	{
		int count = 0;
		if (m_knotsV[aa] != m_knotsV[0])
		{
			while (m_knotsV[aa] == m_knotsV[aa + 1])
			{
				count++;
				aa++;
			}
			count++;
			if (bb >= knot_and_redundancy.size())
			{
				knot_and_redundancy.resize(bb + 1);
			}
			knot_and_redundancy[bb].push_back(aa);
			knot_and_redundancy[bb].push_back(count);
			bb++;
		}
	}
	
	remove_count = 0;
	for (int index1 = 0; index1 < knot_and_redundancy.size(); index1++)//Remove the repetitive knot in U, multiplicity is knot_and_redundancy.size()
	{
		for (int index2 = 0; index2 < knot_and_redundancy[index1][1]; index2++)
		{
			int flag = 1;
			for (int index = 0; index < qbb_cptsV_temp.size(); index++)//Remove the knots along the V direction corresponding to each U.
			{
				int r = knot_and_redundancy[index1][0] - remove_count;  //r is the index of the deleted knot
				int s = knot_and_redundancy[index1][1] - index2;  //s is the multiplicity now
				double u = qbb_knotsV_temp[r];
				/* 1.order p
				* 2. knot vector array U
				* 3. Control point array Pw
				* 4. Delete knot u with duplicate degree of s and their indices r
				* 5. Delete times num */
				int i; i = r - m_degree[1];  //Here we assume that all nodes can be deleted.
				int j; j = r - s;
				int x = 0; //When removing a knot, it is used to record the index of the newly calculated temporary control point array, facilitating the use of the corresponding new control points during the iteration process.
				int removeflag = 0;
				std::vector<MeshLib::CPoint> newpoint_i;
				std::vector<MeshLib::CPoint> newpoint_j;
				std::vector<MeshLib::CPoint> newpoint_i_j;
				newpoint_i.push_back(qbb_cptsV_temp[index][i - 1]);
				newpoint_j.push_back(qbb_cptsV_temp[index][j + 1]);
				while (j - i > 0)
				{
					double quan_i = (u - qbb_knotsV_temp[i]) / (qbb_knotsV_temp[i + m_degree[1] + 1] - qbb_knotsV_temp[i]);
					double quan_j = (u - qbb_knotsV_temp[j]) / (qbb_knotsV_temp[j + m_degree[1] + 1] - qbb_knotsV_temp[j]);
					MeshLib::CPoint temp1 = qbb_cptsV_temp[index][i];
					MeshLib::CPoint temp2 = newpoint_i[x];
					MeshLib::CPoint temp3 = qbb_cptsV_temp[index][j];
					MeshLib::CPoint temp4 = newpoint_j[x];
					newpoint_i.push_back((temp1 - temp2 * (1 - quan_i)) / quan_i);
					newpoint_j.push_back((temp3 - temp4 * quan_j) / (1 - quan_j));

					/*if (!(((temp3 - temp4 * quan_j) / (1 - quan_j) - (temp1 - temp2 * (1 - quan_i)) / quan_i).norm() < error_threshold))
					{
						newpoint_i_j.push_back((temp3 - temp4 * quan_j) / (1 - quan_j));
					}*/
					if (j - i >= 2)
					{
						newpoint_i_j.insert(newpoint_i_j.begin() + x, (temp1 - temp2 * (1 - quan_i)) / quan_i);
						/*newpoint_i_j.push_back((temp3 - temp4 * quan_j) / (1 - quan_j));*/
						newpoint_i_j.insert(newpoint_i_j.begin() + x + 1, (temp3 - temp4 * quan_j) / (1 - quan_j));
					}
					else if (j - i == 1)
					{
						newpoint_i_j.insert(newpoint_i_j.begin() + x, (((temp1 - temp2 * (1 - quan_i)) / quan_i) + (temp3 - temp4 * quan_j) / (1 - quan_j)) / 2);
					}
					x++;
					i++;
					j--;
				}
				if (j - i < 0)
				{
					/*std::cout << "case 1 computation error::" << (newpoint_i[x] - newpoint_j[x]).norm() << std::endl;*/
					if ((newpoint_i[x] - newpoint_j[x]).norm() < error_threshold)//allowable error
					{
						//std::cout << "can be removed by  case 1"" << std::endl;
						removeflag = 1;
					}
				}
				else if (j - i == 0) {
					MeshLib::CPoint temp = qbb_cptsV_temp[index][(i + j) / 2];
					double quan_ij = (u - qbb_knotsV_temp[i]) / (qbb_knotsV_temp[i + m_degree[1] + 1] - qbb_knotsV_temp[i]);
					//std::cout << "case 2 computation error:" << (temp - newpoint_j[x] * quan_ij - newpoint_i[x] * (1 - quan_ij)).norm() << std::endl;
					if ((temp - newpoint_j[x] * quan_ij - newpoint_i[x] * (1 - quan_ij)).norm() < error_threshold)//allowable error
					{
						//std::cout << "can be removed by case 2" << std::endl;
						removeflag = 1;
					}
				}
				if (removeflag == 1)
				{
					int insert_position = r - m_degree[1];
					int erase_start = r - m_degree[1];
					int erase_end = r - s + 1;
					qbb_cptsV_temp[index].erase(qbb_cptsV_temp[index].begin() + erase_start, qbb_cptsV_temp[index].begin() + erase_end);
					qbb_cptsV_temp[index].insert(qbb_cptsV_temp[index].begin() + insert_position, newpoint_i_j.begin(), newpoint_i_j.end());
				}
				else {
					//std::cout << "can't be removed" << std::endl;
					flag = 0;
					break;
				}
			}
			if (flag == 0)
			{
				qbb_knotsV_temp = m_knotsV;
				qbb_cptsV_temp = qbb_cptsV;
				break;
			}
			else if (flag == 1)
			{
				int r = knot_and_redundancy[index1][0] - remove_count;
				qbb_knotsV_temp.erase(qbb_knotsV_temp.begin() + r);
				m_knotsV = qbb_knotsV_temp;
				qbb_cptsV = qbb_cptsV_temp;
				remove_count++;
			}
		}
	}
	////transposition
	//std::vector<std::vector<CPoint>> qbb_cpts_temp1;
	//for (int i = 0; i < qbb_cptsV[0].size(); i++)
	//{
	//	std::vector<CPoint> qbb_temp1;
	//	for (int j = 0; j < qbb_cptsV.size(); j++)
	//	{
	//		qbb_temp1.push_back(qbb_cptsV[j][i]);
	//	}
	//	qbb_cpts_temp1.push_back(qbb_temp1);
	//}
	//qbb_cpts = qbb_cpts_temp1;
	m_cpts = qbb_cptsV;
}

void CCG_QMSLib::CCG_QMS_bsplineSurf::computeUniqueKnotNum()
{
	//U direction
	std::vector<std::vector<double>> knot_and_redundancy;//knot index and its repetition degree
	//compute knot index and its repetition degree
	for (int aa = 0, bb = 0; aa < m_knotsU.size() - m_degree[0] - 1; aa++)
	{
		int count = 0;
		if (m_knotsU[aa] != m_knotsU[0])
		{
			while (m_knotsU[aa] == m_knotsU[aa + 1])
			{
				count++;
				aa++;
			}
			count++;
			if (bb >= knot_and_redundancy.size())
			{
				knot_and_redundancy.resize(bb + 1);
			}
			knot_and_redundancy[bb].push_back(aa);
			knot_and_redundancy[bb].push_back(count);
			bb++;
		}
	}
	/*Consider the first value of the U-direction knot vector*/
	{
		m_knotsU_unique_num.insert(std::pair<double, int>(m_knotsU.front(), m_degree[0] + 1));
	}
	/*Consider the vector of the middle knot in the U direction*/
	for (auto kr : knot_and_redundancy)
	{
		m_knotsU_unique_num.insert(std::pair<double, int>(m_knotsU[kr[0]], kr[1]));
	}
	/*Consider the last value of the U-direction knot vector*/
	{
		m_knotsU_unique_num.insert(std::pair<double, int>(m_knotsU.back(), m_degree[0] + 1));
	}
	
	/*V direction*/
	knot_and_redundancy.clear();//knot index and its repetition degree
	for (int aa = 0, bb = 0; aa < m_knotsV.size() - m_degree[1] - 1; aa++)
	{
		int count = 0;
		if (m_knotsV[aa] != m_knotsV[0])
		{
			while (m_knotsV[aa] == m_knotsV[aa + 1])
			{
				count++;
				aa++;
			}
			count++;
			if (bb >= knot_and_redundancy.size())
			{
				knot_and_redundancy.resize(bb + 1);
			}
			knot_and_redundancy[bb].push_back(aa);
			knot_and_redundancy[bb].push_back(count);
			bb++;
		}
	}
	/*Consider the first value of the V-direction knot vector*/
	{
		m_knotsV_unique_num.insert(std::pair<double, int>(m_knotsV.front(), m_degree[1] + 1));
	}
	/*Consider the vector of the middle knot in the V direction*/
	for (auto kr : knot_and_redundancy)
	{
		m_knotsV_unique_num.insert(std::pair<double, int>(m_knotsV[kr[0]], kr[1]));
	}
	/*Consider the last value of the V-direction knot vector*/
	{
		m_knotsV_unique_num.insert(std::pair<double, int>(m_knotsV.back(), m_degree[1] + 1));
	}
}

CPoint CCG_QMSLib::CCG_QMS_bsplineSurf::DeBoorAlgorithmSurface(const std::vector<std::vector<CPoint>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, double u, double v)
{
	//Search for the knot intervals of u and v
	int i = FindSpan1(degreeU, knotVectorU, u);
	int j = FindSpan1(degreeV, knotVectorV, v);

	//Recursion along the U direction
	std::vector<CPoint> temp(degreeV + 1);
	for (int col = 0; col <= degreeV; ++col) {
		std::vector<CPoint> column(degreeU + 1);
		for (int k = 0; k <= degreeU; ++k) {
			column[k] = controlPoints[i - degreeU + k][j - degreeV + col];
		}
		//Deboor algorithm along the U direction
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

	//Recursion along the V direction
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

/*-------------------function for CCG_QMS_bsplineSurf end------------------------------------*/
