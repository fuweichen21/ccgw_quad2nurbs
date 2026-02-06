/*
*	- (c++) Yiming Zhu
*	- 2019/09/10
*	- Algorithm : A mesh iterator base "c++ 11"
*/

#ifndef _MyLib_MESH_ITERATOR_H_
#define _MyLib_MESH_ITERATOR_H_

#include <list>
#include <vector>
#include <iterator>
#include <algorithm>
/*Need link BaseMesh.h address!*/
//#include <PreMeshLib/core/Mesh/BaseMesh.h>
//#include"MeshLib\core\Mesh\BaseMesh.h"//FuweiChen modified 2021/11/27
#include"../3rdParty/MeshLib/core/Mesh/BaseMesh.h"


namespace MyLib 
{
	template<typename M>
	struct My_MeshIterator 
	{
	private:
		typedef typename M::V	V;
		typedef typename M::E	E;
		typedef typename M::F	F;
		typedef typename M::HE	HE;

	public:
		/*
		* Iterator to access the out halfedge of a vertex in counter-clockwise direction.
		* Avalible only in manifold surface.
		* It is faster than the VoutHEIterator, so try to use this in manifold case.
		* \param pV the pointer to the vertex
		*/
		class VCcwOutHEIterator : public std::iterator<std::forward_iterator_tag, HE*> {
		public:
			VCcwOutHEIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostClwOutHalfEdge(pV)) {};
			VCcwOutHEIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VCcwOutHEIterator& operator++() {
				_pHE = (_pHE == _pMesh->vertexMostCcwOutHalfEdge(_pV)) ? NULL : _pMesh->vertexNextCcwOutHalfEdge(_pHE);
				return *this;
			};

			bool operator==(const VCcwOutHEIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VCcwOutHEIterator& otherIter) { return _pHE != otherIter._pHE; }
			HE* operator*() { return _pHE; }
			HE* value() { return _pHE; }

			VCcwOutHEIterator begin() { return VCcwOutHEIterator(_pMesh, _pV); }
			VCcwOutHEIterator end() { return  VCcwOutHEIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
		};

		/*
		* Iterator to access the out halfedge of a vertex in clockwise direction.
		* Avalible only in manifold surface.
		* It is faster than the VoutHEIterator, so try to use this in manifold case.
		* \param pV the pointer to the vertex
		*/
		class VClwInHEIterator : public std::iterator<std::forward_iterator_tag, HE*> {
		public:
			VClwInHEIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostCcwInHalfEdge(pV)) {};
			VClwInHEIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VClwInHEIterator& operator++() {
				_pHE = (_pHE == _pMesh->vertexMostClwInHalfEdge(_pV)) ? _pHE = NULL : _pMesh->vertexNextClwInHalfEdge(_pHE);
				return *this;
			};

			bool operator==(const VClwInHEIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VClwInHEIterator& otherIter) { return _pHE != otherIter._pHE; }
			HE* operator*() { return _pHE; }
			HE* value() { return _pHE; }

			VClwInHEIterator begin() { return VClwInHEIterator(_pMesh, _pV); }
			VClwInHEIterator end() { return  VClwInHEIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
		};

		class VCcwFIterator : public std::iterator<std::forward_iterator_tag, F*> {
		public:
			VCcwFIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostClwOutHalfEdge(pV)) {};
			VCcwFIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VCcwFIterator& operator++() {
				_pHE = (_pHE == _pMesh->vertexMostCcwOutHalfEdge(_pV)) ? NULL : _pMesh->vertexNextCcwOutHalfEdge(_pHE);
				return *this;
			};

			bool operator==(const VCcwFIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VCcwFIterator& otherIter) { return _pHE != otherIter._pHE; }
			F* operator*() { return _pMesh->halfedgeFace(_pHE); }
			F* value() { return _pMesh->halfedgeFace(_pHE); }

			VCcwFIterator begin() { return VCcwFIterator(_pMesh, _pV); }
			VCcwFIterator end() { return  VCcwFIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
		};

		class VClwFIterator : public std::iterator<std::forward_iterator_tag, F*> {
		public:
			VClwFIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostCcwInHalfEdge(pV)) {};
			VClwFIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VClwFIterator& operator++() {
				_pHE = (_pHE == _pMesh->vertexMostClwInHalfEdge(_pV)) ? _pHE = NULL : _pMesh->vertexNextClwInHalfEdge(_pHE);
				return *this;
			};

			bool operator==(const VClwFIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VClwFIterator& otherIter) { return _pHE != otherIter._pHE; }
			F* operator*() { return _pMesh->halfedgeFace(_pHE); }
			F* value() { return _pMesh->halfedgeFace(_pHE); }

			VClwFIterator begin() { return VClwFIterator(_pMesh, _pV); }
			VClwFIterator end() { return  VClwFIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
		};

		class VCcwEIterator : public std::iterator<std::forward_iterator_tag, E*> {
		public:
			VCcwEIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostClwOutHalfEdge(pV)) {};
			VCcwEIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VCcwEIterator& operator++() {
				if (_pMesh->isBoundary(_pV)) {
					if (_pHE == _pMesh->vertexMostCcwOutHalfEdge(_pV)) {
						_pHE = (HE*)_pHE->he_prev();
						reachBoundary = true;
					}
					else if (reachBoundary) {
						_pHE = NULL;
					}
					else {
						_pHE = _pMesh->vertexNextCcwOutHalfEdge(_pHE);
					}
				}
				else {
					_pHE = (_pHE == _pMesh->vertexMostCcwOutHalfEdge(_pV)) ? _pHE = NULL : _pMesh->vertexNextCcwOutHalfEdge(_pHE);
				}
				return *this;
			};

			bool operator==(const VCcwEIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VCcwEIterator& otherIter) { return _pHE != otherIter._pHE; }
			E* operator*() { return _pMesh->halfedgeEdge(_pHE); }
			E* value() { return _pMesh->halfedgeEdge(_pHE); }

			VCcwEIterator begin() { return VCcwEIterator(_pMesh, _pV); }
			VCcwEIterator end() { return  VCcwEIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
			bool reachBoundary = false;
		};

		class VClwEIterator : public std::iterator<std::forward_iterator_tag, E*> {
		public:
			VClwEIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostCcwInHalfEdge(pV)) {};
			VClwEIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VClwEIterator& operator++() {
				if (_pMesh->isBoundary(_pV)) {
					if (_pHE == _pMesh->vertexMostClwInHalfEdge(_pV)) {
						_pHE = (HE*)_pHE->he_next();
						reachBoundary = true;
					}
					else if (reachBoundary) {
						_pHE = NULL;
					}
					else {
						_pHE = _pMesh->vertexNextClwInHalfEdge(_pHE);
					}
				}
				else {
					_pHE = (_pHE == _pMesh->vertexMostClwInHalfEdge(_pV)) ? _pHE = NULL : _pMesh->vertexNextClwInHalfEdge(_pHE);
				}
				return *this;
			};

			bool operator==(const VClwEIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VClwEIterator& otherIter) { return _pHE != otherIter._pHE; }
			E* operator*() { return _pMesh->halfedgeEdge(_pHE); }
			E* value() { return _pMesh->halfedgeEdge(_pHE); }

			VClwEIterator begin() { return VClwEIterator(_pMesh, _pV); }
			VClwEIterator end() { return  VClwEIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
			bool reachBoundary = false;
		};

		class VCcwVIterator : public std::iterator<std::forward_iterator_tag, V*> {
		public:
			VCcwVIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostClwOutHalfEdge(pV)) {};
			VCcwVIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VCcwVIterator& operator++() {
				if (_pMesh->isBoundary(_pV)) {
					if (_pHE == _pMesh->vertexMostCcwOutHalfEdge(_pV)) {
						_pHE = (HE*)_pHE->he_prev();
						reachBoundary = true;
					}
					else if (reachBoundary) {
						_pHE = NULL;
					}
					else {
						_pHE = _pMesh->vertexNextCcwOutHalfEdge(_pHE);
					}
				}
				else {
					_pHE = (_pHE == _pMesh->vertexMostCcwOutHalfEdge(_pV)) ? _pHE = NULL : _pMesh->vertexNextCcwOutHalfEdge(_pHE);
				}
				return *this;
			};

			bool operator==(const VCcwVIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VCcwVIterator& otherIter) { return _pHE != otherIter._pHE; }
			V* operator*() {
				if (reachBoundary)
					return _pMesh->halfedgeSource(_pHE);
				else
					return _pMesh->halfedgeTarget(_pHE);
			}
			V* value() {
				if (reachBoundary)
					return _pMesh->halfedgeSource(_pHE);
				else
					return _pMesh->halfedgeTarget(_pHE);
			}

			VCcwVIterator begin() { return VCcwVIterator(_pMesh, _pV); }
			VCcwVIterator end() { return  VCcwVIterator(_pMesh, _pV, NULL); }

		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
			bool reachBoundary = false;
		};

		class VClwVIterator : public std::iterator<std::forward_iterator_tag, V*> {
		public:
			VClwVIterator(M* m, V* pV) : _pMesh(m), _pV(pV), _pHE(_pMesh->vertexMostCcwInHalfEdge(pV)) {};
			VClwVIterator(M* m, V* pV, HE* pHE) : _pMesh(m), _pV(pV), _pHE(pHE) {};

			VClwVIterator& operator++() {
				if (_pMesh->isBoundary(_pV)) {
					if (_pHE == _pMesh->vertexMostClwInHalfEdge(_pV)) {
						_pHE = (HE*)_pHE->he_next();
						reachBoundary = true;
					}
					else if (reachBoundary) {
						_pHE = NULL;
					}
					else {
						_pHE = _pMesh->vertexNextClwInHalfEdge(_pHE);
					}
				}
				else {
					_pHE = (_pHE == _pMesh->vertexMostClwInHalfEdge(_pV)) ? _pHE = NULL : _pMesh->vertexNextClwInHalfEdge(_pHE);
				}
				return *this;
			};

			bool operator==(const VClwVIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const VClwVIterator& otherIter) { return _pHE != otherIter._pHE; }
			V* operator*() {
				if (reachBoundary)
					return _pMesh->halfedgeTarget(_pHE);
				else
					return _pMesh->halfedgeSource(_pHE);
			}
			V* value() {
				if (reachBoundary)
					return _pMesh->halfedgeTarget(_pHE);
				else
					return _pMesh->halfedgeSource(_pHE);
			}

			VClwVIterator begin() { return VClwVIterator(_pMesh, _pV); }
			VClwVIterator end() { return  VClwVIterator(_pMesh, _pV, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			V* _pV;
			HE* _pHE;
			bool reachBoundary = false;
		};

		class FHEIterator : public std::iterator<std::forward_iterator_tag, HE*> {
		public:
			FHEIterator(M* m, F* pF) : _pMesh(m), _pF(pF), _pHE((HE*)pF->halfedge()) {};
			FHEIterator(M* m, F* pF, HE* pHE) : _pMesh(m), _pF(pF), _pHE(pHE) {};

			FHEIterator& operator++()
			{
				_pHE = (HE*)_pHE->he_next();
				if (_pHE == (HE*)_pF->halfedge()) _pHE = NULL;
				return *this;
			};

			bool operator==(const FHEIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const FHEIterator& otherIter) { return _pHE != otherIter._pHE; }
			HE* operator*() { return _pHE; }
			HE* value() { return _pHE; }

			FHEIterator begin() { return FHEIterator(_pMesh, _pF); }
			FHEIterator end() { return FHEIterator(_pMesh, _pF, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			HE* _pHE;
			F* _pF;
		};

		class FEIterator : public std::iterator<std::forward_iterator_tag, E*> {
		public:
			FEIterator(M* m, F* pF) : _pMesh(m), _pF(pF), _pHE((HE*)pF->halfedge()) {};
			FEIterator(M* m, F* pF, HE* pHE) : _pMesh(m), _pF(pF), _pHE(pHE) {};

			FEIterator& operator++()
			{
				_pHE = (HE*)_pHE->he_next();
				if (_pHE == (HE*)_pF->halfedge()) _pHE = NULL;
				return *this;
			};

			bool operator==(const FEIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const FEIterator& otherIter) { return _pHE != otherIter._pHE; }
			E* operator*() { return (E*)_pHE->edge(); }
			E* value() { return (E*)_pHE->edge(); }

			FEIterator begin() { return FEIterator(_pMesh, _pF); }
			FEIterator end() { return FEIterator(_pMesh, _pF, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			HE* _pHE;
			F* _pF;
		};

		class FVIterator : public std::iterator<std::forward_iterator_tag, V*> {
		public:
			FVIterator(M* m, F* pF) : _pMesh(m), _pF(pF), _pHE(_pMesh->faceHalfedge(pF)) {};
			FVIterator(M* m, F* pF, HE* pHE) : _pMesh(m), _pF(pF), _pHE(pHE) {};

			FVIterator& operator++()
			{
				_pHE = (HE*)_pHE->he_next();
				if (_pHE == (HE*)_pF->halfedge()) _pHE = NULL;
				return *this;
			};

			bool operator==(const FVIterator& otherIter) { return _pHE == otherIter._pHE; }
			bool operator!=(const FVIterator& otherIter) { return _pHE != otherIter._pHE; }
			V* operator*() { return (V*)_pHE->vertex(); }
			V* value() { return (V*)_pHE->vertex(); }

			FVIterator begin() { return FVIterator(_pMesh, _pF); }
			FVIterator end() { return FVIterator(_pMesh, _pF, NULL); }

			HE* get() { return _pHE; }
		private:
			M* _pMesh;
			HE* _pHE;
			F* _pF;
		};

		class MVIterator : public std::iterator<std::forward_iterator_tag, V*> {
		public:
			MVIterator(M* pM) : _pM(pM), _iter(pM->vertices().begin()) {};

			MVIterator& operator++() { ++_iter; return *this; };

			bool operator==(const MVIterator& otherIter) { return _iter == otherIter._iter; }
			bool operator!=(const MVIterator& otherIter) { return _iter != otherIter._iter; }
			V* operator*() { return *_iter; }
			V* value() { return *_iter; }

			MVIterator begin() { return MVIterator(_pM); }
			MVIterator end() { return MVIterator(_pM, _pM->vertices().end()); }

		private:
			MVIterator(M* pM, typename std::list<V*>::iterator iter) : _pM(pM), _iter(iter) {};
			typename std::list<V*>::iterator _iter;
			M* _pM;
		};

		class MFIterator : public std::iterator<std::forward_iterator_tag, F*> {
		public:
			MFIterator(M* pM) : _pM(pM), _iter(pM->faces().begin()) {};

			MFIterator& operator++() { ++_iter; return *this; };

			bool operator==(const MFIterator& otherIter) { return _iter == otherIter._iter; }
			bool operator!=(const MFIterator& otherIter) { return _iter != otherIter._iter; }
			F* operator*() { return *_iter; }
			F* value() { return *_iter; }

			MFIterator begin() { return MFIterator(_pM); }
			MFIterator end() { return MFIterator(_pM, _pM->faces().end()); }

		private:
			MFIterator(M* pM, typename std::list<F*>::iterator iter) : _pM(pM), _iter(iter) {};
			typename std::list<F*>::iterator _iter;
			M* _pM;
		};

		class MEIterator : public std::iterator<std::forward_iterator_tag, E*> {
		public:
			MEIterator(M* pM) : _pM(pM), _iter(pM->edges().begin()) {};

			MEIterator& operator++() { ++_iter; return *this; };

			bool operator==(const MEIterator& otherIter) { return _iter == otherIter._iter; }
			bool operator!=(const MEIterator& otherIter) { return _iter != otherIter._iter; }
			E* operator*() { return *_iter; }
			E* value() { return *_iter; }

			MEIterator begin() { return MEIterator(_pM); }
			MEIterator end() { return MEIterator(_pM, _pM->edges().end()); }

		private:
			MEIterator(M* pM, typename std::list<E*>::iterator iter) : _pM(pM), _iter(iter) {};
			typename std::list<E*>::iterator _iter;
			M* _pM;
		};

	};

}
#endif //_MyLib_MESH_ITERATOR_H_ defined