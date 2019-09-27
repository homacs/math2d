/*
 * MultiGraph.h
 *
 *  Created on: 20 Jun 2019
 *      Author: homac
 */

#ifndef MATH2D_TRIANGULATOR_H_
#define MATH2D_TRIANGULATOR_H_


#include <iostream>

#include <string>
#include <stdio.h>

#include <assert.h>
#include <stdlib.h>

#include <glm/gtx/vector_angle.hpp>
using namespace glm;


#include <math2d/MatrixMxM.h>

#include <limits>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_set>
#include <algorithm>
using namespace std;


#include <math2d/line.h>
#include <math2d_config.h>

#include <math2d/Graph2DPlanar.h>

struct Test_Triangulator;

namespace math2d {

// TODO: Triangulator consider using integer only
// TODO: performance: floating point causes to many inaccuracies,
//       which have to be solved through branches and usage of double.
//       It slows down comparison between vectors and
//       it adds additional nodes due to inaccuracies.
class Triangulator : public Graph2DPlanar {
	typedef Graph2DPlanar super;
	friend struct ::Test_Triangulator;

public:

	/** Winding order, used when generating edges of polygons. */
	enum WindingOrder {
		/** faces have clockwise winding order */
		WO_CW,
		/** faces have counter-clockwise winding order */
		WO_CCW
	};

	/** rule used to interpret type of area (e.g. even-odd or non-zero)*/
	enum AreaRule {
		/** faces have even weight, and holes have odd weight */
		AR_EVEN_ODD,
		/** faces have non-zero weight, and holes have a weight of zero */
		AR_NON_ZERO,
	};

	/** type of area (visible face or hole/void)*/
	enum AreaType {
		/** area represents a hole or void */
		AREA_HOLE = false,
		/** area represents a face */
		AREA_FACE = true
	};


	struct Area;

	/** We extend edges by area adjacency information. */
	struct AreaEdge : super::Edge {
		/** area to the left according to direction a->b */
		Area* area_left;
		/** area to the right according to direction a->b */
		Area* area_right;

		AreaEdge() : super::Edge() {
			area_left = NULL;
			area_right = NULL;
		}

		AreaEdge(Node* a, Node* b, int weight) : super::Edge() {
			init(a,b,weight);
		}
		virtual ~AreaEdge() {}


		virtual void init(Node* _a, Node* _b, int _weight) {
			super::Edge::init(_a, _b, _weight);
			area_left = NULL;
			area_right = NULL;
		}

		bool isSplitEdge() {
			assert(area_left != NULL && area_right != NULL);
			return area_left == area_right;
		}


		void addToAreaLeft(Area* area) {
			assert(area_left == NULL);
			area_left = area;
			area->add(this);
		}
		void addToAreaRight(Area* area) {
			assert(area_right == NULL);
			area_right = area;
			area->add(this);
		}


		void applyWindingOrder(AreaRule rule, WindingOrder winding) {
			bool isUpOrRight = weight > 0;

			AreaType leftAreaType = area_left->type(rule);

			// apply winding order
			if (leftAreaType == AREA_FACE) {
				if (winding == WO_CCW) {
					if (!isUpOrRight) reverse();
				} else {
					// winding = WO_CW
					if (isUpOrRight) reverse();
				}
			} else {
				// type == AREA_HOLE
				if (winding == WO_CCW) {
					if (isUpOrRight) reverse();
				} else {
					// winding = WO_CW
					if (!isUpOrRight) reverse();
				}
			}

		}


		Node* getStartByWeight() {
			return (weight>0) ? a : b;
		}

		Area* getAreaByWeight(bool left) {
			// return left_area if up and left or !up and !left
			bool up = (weight > 0);
			left = (left == up);
			return (left) ? area_left : area_right;
		}

		Area* getAreaLeftByStartNode(Node* n) {
			if (n == a) {
				return area_left;
			} else {
				assert(n == b);
				return area_right;
			}
		}
		Area* getAreaRightByStartNode(Node* n) {
			if (n == a) {
				return area_right;
			} else {
				assert(n == b);
				return area_left;
			}
		}
	};

	/** overriding createEdge() to replace with AreaEdges */
	virtual Edge* createEdge(Node* a, Node* b, int weight) {
		return new AreaEdge(a,b,weight);
	}

	struct AreaEdges : vector<AreaEdge*>  {
		void removeAll(AreaEdges& edges) {
			set<AreaEdge*> filter;
			for (iterator it = begin(); it != end(); ) {
				if (filter.find(*it) != filter.end()) {
					it = erase(it);
				} else {
					it++;
				}
			}

		}
	};

	typedef set<Area*> Neighbours;

	struct Area {
		AreaEdges edges;
		int weight;



#ifndef NDEBUG
		static int total_areas_in_system;
#endif

	public:

		Area(int _weight) {
			this->weight = _weight;
#ifndef NDEBUG
			total_areas_in_system++;
#endif
		}

		~Area() {
#ifndef NDEBUG
			total_areas_in_system--;
#endif
		}

		void reset() {
			edges.clear();
		}

		void init(int _weight) {
			this->weight = _weight;
		}

		void add(AreaEdge* e) {
			assert(e->area_left == this || e->area_right == this);
			edges.push_back(e);
		}

		AreaType type(AreaRule rule) {
			switch(rule) {
			case AR_EVEN_ODD:
				return AreaType((weight & 0x1f) > 0);
			case AR_NON_ZERO:
				return AreaType(weight != 0);
			default:
				assert(0 != "missing case in switch");
				return AreaType(-1);
			}
		}



	};




	struct Areas : set<Area*> {
	private:



	public:
		Area* AREA_VOID;

		Areas() {

			AREA_VOID = new Area(0);
			init();
		}

		void init() {
		}

		~Areas() {
			reset();
			delete AREA_VOID;
		}


		bool exists(Area* area) {
			if (area == AREA_VOID) return true;
			iterator it = find(area);
			return it != end();
		}


		void reset() {
			for (iterator it = begin(); it != end(); it++) {
				delete (*it);
			}
			clear();
			AREA_VOID->reset();
			init();
		}

		Area* create(Area* prev, AreaEdge* left) {
			assert(prev != NULL);
			assert(left != NULL);
			int weight = prev->weight;
			weight += left->weight;
			Area* area = new Area(weight);
			emplace(area);
			return area;
		}

		void remove(Area* area) {
			iterator it = find(area);
			assert(it != end());
			delete (*it);
			erase(it);
		}

		/** Merges edges of B into A and deletes B.
		 * Updates area assignment of edges.
		 * */
		void mergeEdges(Area* A, Area* B) {

			for (AreaEdges::iterator it = B->edges.begin(); it != B->edges.end(); it++) {
				AreaEdge* e = (*it);
				if (e->area_left == B) {
					e->area_left = A;
				} else {
					assert(e->area_right == B);
					e->area_right = A;
				}
			}
			A->edges.insert(A->edges.end(), B->edges.begin(), B->edges.end());

			remove(B);
		}




		/**
		 * Merge neighbours and remove common edges.
		 */
		void mergeNeighbours(Area* A, Area* B) {
			// remove all splitting edges between A and B
			for (AreaEdges::iterator it = A->edges.begin(); it != A->edges.end(); ) {
				AreaEdge* e = (*it);
				if (e->area_left == A && e->area_right != B)
				{
					// keep
					it++;
				}
				else if (e->area_right == A && e->area_left != B)
				{
					// keep
					it++;
				}
				else
				{
					it = A->edges.erase(it);
				}
			}

			// add edges of B
			for (AreaEdges::iterator it = B->edges.begin(); it != B->edges.end(); it++) {
				AreaEdge* e = (*it);
				if (e->area_left == B && e->area_right != A)
				{
					e->area_left = A;
					A->edges.push_back(e);
				}
				else if (e->area_right == B && e->area_left != A)
				{
					e->area_right = A;
					A->edges.push_back(e);
				}
			}

		}

	};


	typedef map<double, Node*> NodeDistances;



	/*
	 * This is just a mask used in reinterpret casts of EdgeRefs::iterator-s.
	 * There is no actual instance of this structure.
	 */
	struct AreaEdgeRef {
		/** reference on connected node. */
		Node* to;

		/** reference on the edge data */
		AreaEdge* edge;

		/** orientation of b-a */
		double orientation;

	private:
		AreaEdgeRef();

	};

	typedef vector<AreaEdgeRef*> AreaEdgeRefs;

	enum Status {
		STATUS_INIT,
		STATUS_AREAS_FOUND,
		STATUS_AREAS_MERGED,
		STATUS_AREAS_MONOTONIZED,
	};

private:

	AreaRule areaRule;
	WindingOrder winding;

	Status status;

	Areas areas;

	/** reference on super::edges of type Edges stored in base class but we interpret them as AreaEdges */
	AreaEdges& edges;





public:

	Triangulator() : Triangulator(WindingOrder::WO_CCW, AreaRule::AR_NON_ZERO, false, 0) {
	}
	Triangulator(WindingOrder winding, AreaRule areaRule)
	: Triangulator(winding, areaRule, false, 0)
	{
	}
	Triangulator(WindingOrder winding, AreaRule areaRule, bool keepAdjacentAreas, int min_nodes) : super(min_nodes, keepAdjacentAreas), edges(forced_cast<AreaEdges>(super::edges)) {
		assert(assertEdgeRefCompliance());
		this->status = STATUS_INIT;
		this->winding = winding;
		this->areaRule = areaRule;
	}


	virtual ~Triangulator() {
		reset();
	}

	void reset() {
		super::reset();
		areas.reset();
		this->status = STATUS_INIT;
	}

	/** advance to required status */
	void advance(Status requiredStatus) {
		while (status != requiredStatus) {
			switch(status) {
			case STATUS_INIT:
				findAreas(); // -> STATUS_AREAS_FOUND
				break;
			case STATUS_AREAS_FOUND:
				mergeAreasOfSameType(); // -> STATUS_AREAS_MERGED
				break;
			case STATUS_AREAS_MERGED:
				monotonize_X(); // -> STATUS_AREAS_MONOTONIZED
				break;
			case STATUS_AREAS_MONOTONIZED:
				assert(false);
				break;
			}
		}
	}

public:

#ifndef NDEBUG

	bool assertEdgeRefCompliance() {
		// this test fails, if EdgeRef != AreaEdgeRef

		EdgeRef* e = NULL;
		AreaEdgeRef* ae = (AreaEdgeRef*)(e);
		// check size
		assert(sizeof(EdgeRef) == sizeof(AreaEdgeRef));
		// check member's addresses
		assert(&e->to == &ae->to);
		assert((void*)&e->edge == (void*)&ae->edge);
		assert(&e->orientation == &ae->orientation);
		return true;
	}
#endif


	template <bool HORIZONTAL>
	struct AELEntry {
		// TODO: move to appropriate place
		/** coordinate x or y according to sorting order (horizontal or vertical). */
		float coord;
		/** scan line end coordinate. The scan line height, where this entry will be discarded */
		float scan_end;
		/** edge orientation against Y-axis (horizontal) or X-axis (vertical) */
		double orientation;
		/** edge */
		AreaEdge* e;
		/** node */
		Node* n;


		AELEntry(float coord, double orientation) {
			this->coord = coord;
			this->orientation = orientation;
#ifndef NDEBUG
			scan_end = coord;
			e = 0;
			n = 0;
#endif


		}

		AELEntry(AreaEdge* e) {
			init(e);
		}

		void init(AreaEdge* _e) {
			e = _e;
			n = 0;
			Node* start = HORIZONTAL ? e->a : (e->a->x() < e->b->x()) ? e->a : e->b;
			Node* end = e->opposite(start);
			coord = HORIZONTAL ? start->x() : start->y();
			scan_end = HORIZONTAL ? end->y() : end->x();
			vec2 b_a = normalize(e->b->pos - e->a->pos);
			orientation = HORIZONTAL ? pseudo_orientation_e2(b_a) : pseudo_orientation_e1(b_a);
		}

		AELEntry(Node* n) {
			init(n);
		}

		void init(Node* _n) {
			e = 0;
			n = _n;
			coord = HORIZONTAL ? n->pos.x : n->pos.y;
			scan_end = HORIZONTAL ? n->pos.y : n->pos.x;
			orientation = HORIZONTAL ? 0 : PSEUDO_PI_HALF;
		}


		bool isNode() {
			return n;
		}

	};


	/** defines less operator for AEL entries a, b :
	 * [a->x, a->orientation] < [b->x, b->orientation]
	 */
	template <bool HORIZONTAL>
	struct AELEntry_less {
		bool operator() (AELEntry<HORIZONTAL> * const a, AELEntry<HORIZONTAL> * const b) {
			double diff = a->coord - b->coord;
			if (diff == 0) {
				diff = a->orientation - b->orientation;
			}
			return diff < 0;
		}
	};


	/**
	 * A list of AELEntries ordered by AELEntry_less
	 */
	template <bool HORIZONTAL>
	struct T_AEL : set<AELEntry<HORIZONTAL>*, AELEntry_less<HORIZONTAL>> {
		typedef set<AELEntry<HORIZONTAL>*, AELEntry_less<HORIZONTAL>> super;
		AELEntry_less<HORIZONTAL> LT;
		typedef typename super::iterator iterator;


		~T_AEL() {
			reset();
		}

		void reset() {
			for (iterator it = super::begin(); it != super::end(); it++) {
				delete *it;
			}
			super::clear();
		}

		bool isActiveEdge(AreaEdge* edge, float scan_coord) {
			if (HORIZONTAL) return edge->a->y() <= scan_coord && edge->b->y() > scan_coord;
			else return edge->a->x() <= scan_coord && edge->b->x() > scan_coord;
		}

		iterator put (AreaEdge* e) {
			pair<iterator, bool> result = super::insert(new AELEntry<HORIZONTAL>(e));
			assert(result.second);
			return result.first;
		}

		iterator put (Node* n) {
			pair<iterator, bool> result = super::insert(new AELEntry<HORIZONTAL>(n));
			assert(result.second);
			return result.first;
		}

		iterator put_behind (AreaEdge* e, iterator it_prev) {
			iterator hint(it_prev);
			if (hint != super::end()) ++hint;

			return super::emplace_hint(hint, new AELEntry<HORIZONTAL>(e));
		}

		iterator put_behind (Node* n, iterator it_prev) {
			iterator hint(it_prev);
			if (hint != super::end()) ++hint;

			return super::emplace_hint(hint, new AELEntry<HORIZONTAL>(n));
		}

		iterator find(AreaEdge* e) {
			AELEntry<HORIZONTAL> entry(e);
			return super::lower_bound(&entry);
		}

		/** find first entry of the given node */
		iterator find_first(Node* n) {
			float coord = HORIZONTAL ? n->pos.x : n->pos.y;
			double orientation = HORIZONTAL ? -PSEUDO_PI_HALF : 0;
			AELEntry<HORIZONTAL> entry(coord, orientation);
			return super::lower_bound(&entry);
		}

		/** find last entry of the given node */
		iterator find_last(Node* n) {
			float coord = HORIZONTAL ? n->pos.x : n->pos.y;
			double orientation = HORIZONTAL ? +PSEUDO_PI_HALF : +PSEUDO_PI;
			AELEntry<HORIZONTAL> entry(coord, orientation);
			return --super::upper_bound(&entry);
		}

		iterator remove (iterator it) {
			iterator result = super::erase(it);
			delete *it;
			return result;
		}

		iterator prev (AreaEdge* e) {
			return prev(find(e));
		}

		iterator prev (iterator it) {
			if (it == super::begin()) return super::end();
			else return --iterator(it);
		}
		iterator next (iterator it) {
			if (it == super::end()) return super::end();
			else return ++iterator(it);
		}
	};


	Areas* findAreas() {
		// TODO: it may be possible to combine findAreas() and monotonize_X().
		typedef T_AEL<true> AEL;

		areas.reset();
		Area* AREA_VOID = areas.AREA_VOID;
		AEL ael;
		AreaEdge* last_horizontal;
		float scan_y;


		// INITIALISATION
		if (nodes.size() == 0) return &areas;
		scan_y = -INFINITY;
		last_horizontal = NULL;

		// go through all nodes, from bottom to top
		for (Nodes::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++ ) {
			Node* n = (*it_n);

			// if new scan line
			if (scan_y < n->pos.y) {
				scan_y = n->pos.y;
				// reset last horizontal
				last_horizontal = NULL;

				//
				// Tasks:
				// 1. update AEL entries,
				// 2. remove edges ending at this scan line from AEL
				//
				AreaEdge* edge = NULL;
				for (AEL::iterator it_ael = ael.begin(); it_ael != ael.end();) {
					AEL::value_type entry = (*it_ael);
					edge = entry->e;

					if (edge->b->pos.y == scan_y) {
						// end node on scan line -> edge will be removed

						// remove ending edges
						it_ael = ael.remove(it_ael);
					} else {
						float x;
						scanline_horizontal_intersection_m_below(scan_y, edge->a->pos, edge->b->pos, x);
						// hit only, if set contains crossing edges
						entry->coord = x;
						assert(it_ael == ael.begin() || ael.LT(*ael.prev(it_ael), entry));

						// advance in list
						it_ael++;
					}
				}

			} // AEL update finished


			//
			// add new edges to AEL and assign areas
			//


			assert(n->grad() >= 2);
			// first edge
			EdgeRefs::iterator _it_er  = n->findEdgeOriented(-PSEUDO_PI_HALF);
			// edge >= PSEUDO_PI_HALF
			EdgeRefs::iterator _it_er_end = n->findEdgeOriented(+PSEUDO_PI_HALF);

			AreaEdgeRefs::iterator& it_er = forced_cast<AreaEdgeRefs::iterator>(_it_er);
			AreaEdgeRefs::iterator& it_er_end = forced_cast<AreaEdgeRefs::iterator>(_it_er_end);

			AreaEdgeRefs& n_edgerefs = forced_cast<AreaEdgeRefs>(n->edgerefs);

			if (it_er_end == n_edgerefs.end() || (*it_er_end)->orientation > PSEUDO_PI_HALF) {
				it_er_end--;
			}
			if (it_er <= it_er_end) {

				//
				// process first edge of the node, going upwards
				//

				AEL::iterator it_ael_prev = ael.end();
				Area* area_prev = NULL;
				if (it_er == it_er_end) {
					// this edge either points to the right or even downward
					// leave to post processing
				} else {
					if ((*it_er)->edge == last_horizontal) {
						// edge pointing left
						// --> was initialised by last node on same Y
						// skip to next edge
						it_er++;
					}

					if (it_er != it_er_end) {
						Area* area_left = NULL;
						// first edge actually pointing upward
						AEL::iterator it_ael_cur = ael.put((*it_er)->edge);
						it_ael_prev = ael.prev(it_ael_cur);

						if (it_ael_prev != ael.end()) {
							area_left = (*it_ael_prev)->e->area_right;
						} else {
							area_left = AREA_VOID;
						}
						(*it_ael_cur)->e->addToAreaLeft(area_left);
						area_prev = area_left;
						it_ael_prev = it_ael_cur;
						it_er++;
					}
				}

				//
				// process all edges between first and last in upward direction
				//
				for (;it_er != it_er_end; it_er++) {
					// add new areas between all edges pointing upward

					assert(it_ael_prev != ael.end());

					// first edge actually pointing upward
					AEL::iterator it_ael_cur = ael.put_behind((*it_er)->edge, it_ael_prev);

					// preconditions
					assert(it_ael_prev == ael.prev(it_ael_cur));
					assert((*it_ael_prev)->e->area_right == NULL);

					Area* area_left = areas.create(area_prev, (*it_ael_prev)->e);
					(*it_ael_prev)->e->addToAreaRight(area_left);
					(*it_er)->edge->addToAreaLeft(area_left);

					area_prev = area_left;
					it_ael_prev = it_ael_cur;
				}


				//
				// processing last edge (up or right)
				//

				assert(it_er == it_er_end);
				if ((*it_er)->edge == last_horizontal) {
					// first == last and horizontal to left
					// This can happen if this a maximum of
					// a polygon (i.e. no upward pointing edge)

					// end of a chain of horizontal edges
					last_horizontal = NULL;
					// horizontal edge to left has already been initialised by previous node
					// ==> nothing to do
				}
				else
				{
					// at least one edge pointing upward or right

					AEL::iterator it_ael_cur;

					if ((*it_er)->orientation == PSEUDO_PI_HALF) {
						// pointing horizontal right -> don't add to AEL, but determine its theoretical position to get its predecessor
						it_ael_cur = ael.find((*it_er)->edge);

						last_horizontal = (*it_er)->edge;
					} else {
						// still pointing upward

						// it_ael_prev may or may not be set -> use it anyway as hint
						it_ael_cur = ael.put_behind((*it_er)->edge, it_ael_prev);
						// end of horizontal chain
						last_horizontal = NULL;
					}

					// make sure it_ael_prev is set
					if (it_ael_prev == ael.end()) {
						it_ael_prev = ael.prev(it_ael_cur);
					}

					Area* area_left = NULL;
					if (it_ael_cur == ael.begin()) {
						area_left = AREA_VOID;
					} else if ((*it_ael_prev)->e->area_right == NULL) {
						// it_ael_prev  is an edge of this node going up
						area_left = areas.create(area_prev, (*it_ael_prev)->e);
						(*it_ael_prev)->e->addToAreaRight(area_left);
					} else {
						// it_ael_prev  is an edge of another node
						area_left = (*it_ael_prev)->e->area_right;
					}
					assert(area_left != NULL);
					(*it_er)->edge->addToAreaLeft(area_left);

					// To determine the right area of the last edge, we have to
					// find the next edge of this node in clockwise order.
					// This implicitly handles split nodes too.
					Area* area_right = NULL;
					EdgeRefs::iterator __it_er_next = n->clockwise_next(forced_cast<EdgeRefs::iterator>(it_er));
					AreaEdgeRefs::iterator it_er_next = forced_cast<AreaEdgeRefs::iterator>(__it_er_next);
					assert(it_er_next != n_edgerefs.end());
					assert(it_er_next != it_er);
					if ((*it_er_next)->to == (*it_er_next)->edge->a) {
						// from n downward or to the left pointing edge
						area_right = (*it_er_next)->edge->area_right;
					} else {
						// upward pointing edge
						area_right = (*it_er_next)->edge->area_left;
					}
					assert(area_right != NULL);
					assert((*it_er)->edge->area_right == NULL);
					(*it_er)->edge->addToAreaRight(area_right);
				}

			}


			//
			// handle merge nodes
			//
			// find the two nodes left of 0 and right of 0.
			// For n to be a merge node, both edges have to
			// be going downwards (e->b == n).
			//
			// Merge areas if left is going down or horizontally
			// left and right is going down.
			{
				EdgeRefs::iterator it_er = n->findEdgeOriented(0);
				if (it_er == n->edgerefs.end()) it_er = n->edgerefs.begin();

				assert(it_er != n->edgerefs.end());
				AreaEdge* right = (AreaEdge*)((*it_er)->edge);

				it_er = n->clockwise_prev(it_er);
				assert(it_er != n->edgerefs.end());
				AreaEdge* left = (AreaEdge*)((*it_er)->edge);
				assert(left != right);
				if (left->b == n) {
					// left edge is downward or horizontally left pointing edge

					if (right->b == n) {
						// right edge is downward pointing edge
						if (left->area_left != right->area_right) {
							// TODO: find out, why some areas are already merged
							assert(left->area_left != NULL);
							assert(right->area_right != NULL);
							areas.mergeEdges(left->area_left, right->area_right);
							assert(left->area_left == right->area_right);
						}
					}

				}
			}
		} // nodes finished

		status = STATUS_AREAS_FOUND;

		return &areas;
	}


	AreaEdges::iterator removeAreaEdge(AreaEdges::iterator it_ae) {
		Edges::iterator it_e = forced_cast<Edges::iterator>(it_ae);
		it_e = super::removeEdge(it_e);
		it_ae = forced_cast<AreaEdges::iterator>(it_e);
		return it_ae;
	}

	/** Does not consider neighbourship in any way! */
	int mergeAreasOfSameType() {

		// merge areas of the same type and cleanup internal structures
		int merged = 0;

		for (AreaEdges::iterator it_e = edges.begin(); it_e != edges.end();) {
			AreaEdge* edge = *it_e;

			Area* A = edge->area_left;
			Area* B = edge->area_right;

			if (!(areas.exists(A) && areas.exists(B))) {
				// at least one of the areas does not exist anymore
				// => areas have been merged already
				it_e = removeAreaEdge(it_e);
			} else if (A->type(areaRule) == B->type(areaRule)) {
				// Found adjacent areas A and B of same type

				// Merge B into A where A == VOID or A < B and B != VOID
				if ((B == areas.AREA_VOID) || (A != areas.AREA_VOID && A > B)) {
					// swap A <-> B
					Area* tmp = A;
					A = B;
					B = tmp;
				}

				// both areas exist: we have not merged yet
				areas.mergeNeighbours(A, B);
				areas.remove(B);

				it_e = removeAreaEdge(it_e);
			} else {
				// keep edge
				it_e++;
			}

		}

		super::removeEmptyNodes();

		status = STATUS_AREAS_MERGED;
		return merged;
	}

	static double vec2_cmp_x(const vec2& a, const vec2& b) {
		double result = a.x - b.x;
		if (result == 0) result = a.y - b.y;
		return result;
	}

	struct NodeLessX {
		bool operator () (Node* const a, Node* const b) {
			return vec2_cmp_x(a->pos, b->pos) < 0;
		}
	};


	AreaEdge* createSplitEdge(Node* m, Node* n, Area* area) {
		AreaEdge* edge = (AreaEdge*)createEdge(n, m, 1);
		// add edge to nodes and edge list
		insertEdge(edge);
		// assign to area (both sides)
		edge->area_left = area;
		edge->area_right = area;
		area->add(edge);
		return edge;
	}

	AreaEdges::iterator findAreaEdge(AreaEdge* edge) {
		Edges::iterator it = super::findEdge(edge);
		return forced_cast<AreaEdges::iterator>(it);
	}

	void changeSplit(AreaEdge* split, Node* a, Node* b) {
		assert(Node_lt(a, b));

		edges.erase(findAreaEdge(split));
		if (split->a == a) {
			assert(split->b != b);
			split->b->removeEdge(split);
			split->b = b;
			b->addEdge(split);
		} else {
			assert(split->b == b);
			assert(split->a != a);
			split->a->removeEdge(split);
			split->a = a;
			a->addEdge(split);
		}
		edges.insert(findAreaEdge(split),split);
	}

	bool isSplitEdge(AreaEdge* e) {
		return e->area_left == e->area_right;
	}

	void updateAEL(T_AEL<true>& ael, Nodes::iterator it_n_begin, float scan_y, bool downSplits) {
		typedef T_AEL<true> AEL;
		// update coords of AEL entries
		for (AEL::iterator it_ael = ael.begin(); it_ael != ael.end(); it_ael++ ) {
			AEL::value_type entry = *it_ael;
			if (!entry->isNode()) {
				// update coord of active edges
				AreaEdge* edge = entry->e;
				scanline_horizontal_intersection_m_below(scan_y, edge->a->pos, edge->b->pos, entry->coord);
				if (it_ael != ael.begin()) {
					AEL::iterator it_prev = ael.prev(it_ael);
					// find previous edge
					while (it_prev != ael.end() && (*it_prev)->isNode()) it_prev = ael.prev(it_prev);
					if (it_prev != ael.end()
						&& (*it_prev)->coord >= (*it_ael)->coord
						&& edge->b != (*it_prev)->e->b
						&& edge->a != (*it_prev)->e->b
						&& edge->b != (*it_prev)->e->a
					)
					{

						assert(edge->a != (*it_prev)->n && edge->a != (*it_prev)->e->a);

						// found intersection between split and another edge
						AEL::iterator it_s;
						AreaEdge* split;
						AreaEdge* neighbour;
						if (isSplitEdge(edge)) {
							// split edge is current
							it_s = it_ael;
							split = (*it_s)->e;
							neighbour = (*it_prev)->e;
							if (downSplits) {
								Node* m = neighbour->a->x() > neighbour->b->x() ? neighbour->a : neighbour->b;
								// split right down
								changeSplit(split, m, split->b);
							} else {
								Node* m = neighbour->a->x() < neighbour->b->x() ? neighbour->a : neighbour->b;
								// split left up
								changeSplit(split, split->a, m);
							}
							assert(split->a->x() > split->b->x());
						} else {
							// split edge is previous
							it_s = it_prev;
							split = (*it_s)->e;
							AreaEdge* neighbour = (*it_ael)->e;
							assert(isSplitEdge(split));
							if (downSplits) {
								// split left down
								Node* m = neighbour->a->x() < neighbour->b->x() ? neighbour->a : neighbour->b;
								changeSplit(split, m, split->b);
							} else {
								// split right up
								Node* m = neighbour->a->x() > neighbour->b->x() ? neighbour->a : neighbour->b;
								changeSplit(split, split->a, m);
							}
							assert(split->a->x() < split->b->x());
						}
						assert(!isSplitEdge(neighbour));
						// update ael entry for changed split
						scanline_horizontal_intersection_m_below(scan_y, split->a->pos, split->b->pos, (*it_s)->coord);
					}
				}
			}
		}

		// TODO: this can be moved into the loop above if we remove only previous entries
		//       and update it_ael like this it_ael = ++ael.remove(it_ael_prev);

		//
		// remove entries below scan_y
		//
		for (AEL::iterator it_ael = ael.begin(); it_ael != ael.end(); ) {

			if ((*it_ael)->scan_end <= scan_y) {
				it_ael = ael.remove(it_ael);
			} else {
				it_ael++;
			}
		}

#ifndef NDEBUG
		// check consistency of AEL
		for (AEL::iterator it_ael = ael.begin(); it_ael != ael.end(); it_ael++ ) {
			// hit only, if ael still contains crossing edges or nodes
			assert(it_ael == ael.begin() || ael.LT(*ael.prev(it_ael),*it_ael));
		}
#endif



		//
		// add all edges of nodes on current scan line, going up
		//

		// this iterator is used for ael insert hints and points always to the last entry added
		AEL::iterator it_ael_prev = ael.begin();
		for (Nodes::iterator it_n2 = it_n_begin; it_n2 != nodes.end(); it_n2++) {
			Node* n2 = *it_n2;
			if (n2->pos.y > scan_y) {
				break;
			} else {
				// first edge
				EdgeRefs::iterator it_er  = n2->findEdgeOriented(-PSEUDO_PI_HALF);
				if (it_er != n2->edgerefs.end() && (*it_er)->orientation == -PSEUDO_PI_HALF) it_er++;
				// edge >= PSEUDO_PI_HALF
				EdgeRefs::iterator it_er_end = n2->findEdgeOriented(+PSEUDO_PI_HALF);

				if (it_er == it_er_end) {
					// add node without any edges going up
					it_ael_prev = ael.put_behind(n2, it_ael_prev);
					// we don't need to check for intersections with split edges,
					// because a split edge would have had an intersection with
					// at least one of the edges ending in this node earlier.
					// TODO: test case
				}
				else
				{
					for (;it_er != it_er_end; it_er++) {
						AreaEdgeRef* edgeref = (AreaEdgeRef*)(*it_er);
						it_ael_prev = ael.put_behind(edgeref->edge, it_ael_prev);
					}
				}
				// check for intersections between horizontal edge and split
				if (it_er_end != n2->edgerefs.end() && (*it_er_end)->orientation == +PSEUDO_PI_HALF) {
					// last edge is a horizontal edge to the right
					AreaEdge* horizontal = (AreaEdge*)((*it_er_end)->edge);
					// all edges in ael between a and b do intersect with this horizontal edge.
					AEL::iterator it_ael = it_ael_prev;
					AEL::iterator it_ael_test_end = ael.end();
					AEL::iterator it_ael_end = ael.find_first(horizontal->b);
					// step over horizontal->a but not beyond it_ael_end
					if (it_ael != it_ael_end) it_ael++;
					while (it_ael != it_ael_end) {
						AreaEdge* split = (*it_ael)->e;
						// this has to be a split
						assert(isSplitEdge(split));
						if (split->a->x() < split->b->x()) {
							// split is previous
							if (downSplits) {
								// split left down
								changeSplit(split, horizontal->a, split->b);
							} else {
								// split right up
								changeSplit(split, split->a, horizontal->b);
							}
						} else {
							// split is current
							assert(split->a->x() > split->b->x());
							if (downSplits) {
								// split right down
								changeSplit(split, horizontal->b, split->b);
							} else {
								// split left up
								changeSplit(split, split->a, horizontal->a);
							}
						}
						// split now ends on current scan_y
						it_ael = ael.remove(it_ael);
					}
				}
			}

		}
	}



#ifndef NDEBUG
	static struct VECTOR_MOCK {
		void push_back(AreaEdge* e){}
	}vector_mock;

	template <class DEBUG_SPLITS = VECTOR_MOCK>
	void monotonize_X(DEBUG_SPLITS& debug_splits = vector_mock)
#else
	/** creates weak X-monotony */
	void monotonize_X ()
#endif
	{

		typedef T_AEL<true> AEL;

		if (nodes.size() == 0) return;

		AEL ael;
		float scan_y = -INFINITY;

		for (Nodes::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++) {
			Node* n = *it_n;


			if (scan_y < n->pos.y) {
				scan_y = n->pos.y;
				updateAEL(ael, it_n, scan_y, false);
			}

			// find split or merge edges of n


			//
			// split node: two neighbouring edges going right and common FACE to the left of n
			//

			// edge with orientation >= 0   or   no edge
			EdgeRefs::iterator it_er_1  = n->findEdgeOriented(0);
			// neighbouring edge in counter clockwise direction
			EdgeRefs::iterator it_er_2 = n->clockwise_prev(it_er_1);
			// all nodes have at least 2 edges, thus end()-1 is an edge
			assert (it_er_2 != n->edgerefs.end());

			// Note: To consider vertical edges, we will also accept
			//       first edge up and second right
			if (it_er_1 != n->edgerefs.end() && (*it_er_1)->orientation >= 0 && (*it_er_2)->orientation > 0 && (*it_er_2)->orientation < PSEUDO_PI) {
				// both edges point to the right
				assert((*it_er_2)->orientation <= PSEUDO_PI);
				assert((*it_er_1)->orientation >= 0 && (*it_er_1)->orientation < PSEUDO_PI);
				assert((*it_er_1)->orientation < (*it_er_2)->orientation);

				AreaEdge* edge_1 = (AreaEdge*)(*it_er_1)->edge;
				AreaEdge* edge_2 = (AreaEdge*)(*it_er_2)->edge;
				Area* area = edge_1->getAreaLeftByStartNode(n);
				assert(edge_2->getAreaRightByStartNode(n) == area);
				if (area->type(areaRule) == AREA_FACE)
				{
					// split node found!

					//
					// find a node m to the left of n to create split edge
					//
					AEL::iterator it_ael = ael.find_first(n);
					it_ael--; // closest ael entry to the left
					AEL::value_type entry = *it_ael;
					Node* m;
					bool addToAel = false;
					if (entry->isNode()) {
						// don't add horizontal edge to ael
						m = entry->n;
					} else {
						AreaEdge* edge = entry->e;
						if (edge->a->y() == n->y() && edge->a->x() < n->x()) {
							// prefer horizontal edge
							m = edge->a;
						} else if (edge->b->x() < n->x()) {
							m = edge->b;
							addToAel = true;
						} else {
							// down split
							m = edge->a;
						}
					}
					// edge(n,m) is supposed to point to the left
					assert(n->x() >= m->x());
					//
					// create split edge
					//

					AreaEdge* split = createSplitEdge(n, m, area);

					// add edge to AEL
					if (addToAel) {
						ael.put_behind(split, it_ael);
					}
#ifndef NDEBUG
					debug_splits.push_back(split);
#endif
				}
			}

			//
			// merge node: neighbouring edges going left with FACE to the right
			//

			// edge with orientation >= -PI
			it_er_1  = n->edgerefs.begin();
			// neighbouring edge in counter clockwise direction
			it_er_2 = n->edgerefs.end()-1;
			// all nodes have at least 2 edges, thus end()-1 is an edge
			assert (it_er_1 != it_er_2);
			// Note: To consider vertical edges, we will also accept
			//       first edge left and second up
			if ((*it_er_1)->orientation < 0 && (*it_er_2)->orientation <= 0) {
				// both edges point to the left

				// (1) predecessor of (2)
				assert((*it_er_1)->orientation < (*it_er_2)->orientation);

				AreaEdge* edge_1 = (AreaEdge*)(*it_er_1)->edge;
				AreaEdge* edge_2 = (AreaEdge*)(*it_er_2)->edge;
				Area* area = edge_1->getAreaLeftByStartNode(n);

				// neighbouring edges always border same area
				assert(edge_2->getAreaRightByStartNode(n) == area);
				if (area->type(areaRule) == AREA_FACE)
				{
					// split node found!

					//
					// find a node m to the right of n to create split edge
					//
					AEL::iterator it_ael = ael.find_last(n);
					it_ael++; // closest ael entry to the right
					AEL::value_type entry = *it_ael;
					Node* m;
					bool addToAel = false;
					if (entry->isNode()) {
						// don't add horizontal edge to ael
						m = entry->n;
					} else {
						AreaEdge* edge = entry->e;
						if (edge->a->y() == n->y() && n->x() < edge->a->x()) {
							// prefer horizontal edge
							m = edge->a;
						} else if (n->x() < edge->b->x()) {
							// up split
							addToAel = true;
							m = edge->b;
						} else {
							// down split
							m = edge->a;
						}
					}
					// edge(n,m) is supposed to point to the right
					assert(n->x() <= m->x());

					//
					// create split edge
					//
					AreaEdge* split = (AreaEdge*)createSplitEdge(n, m, area);
					// add edge to AEL
					if (addToAel) {
						ael.put_behind(split, --it_ael);
					}
#ifndef NDEBUG
					debug_splits.push_back(split);
#endif
				}

			}

		}

		// check for intersections between down splits and regular edges
		ael.reset();
		scan_y = -INFINITY;
		for (Nodes::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++) {
			Node* n = *it_n;
			if (scan_y < n->pos.y) {
				scan_y = n->pos.y;
				updateAEL(ael, it_n, scan_y, true);
			}
		}

		status = STATUS_AREAS_MONOTONIZED;
	}

	struct Polygon : vector<vec2> {
		typedef vector<vec2> super;

		bool closed;
		WindingOrder winding;

		Polygon() : Polygon(WO_CCW){}

		Polygon(WindingOrder winding) {
			reset(winding);
		}

		void reset(WindingOrder winding) {
			clear();
			closed = false;
			this->winding = winding;
		}

		void setClosed(bool closed) {
			this->closed = closed;
		}

		void swap (vec2& a, vec2& b) {
			vec2 tmp = a;
			a = b;
			b = tmp;
		}


		void reverseWinding() {
			for (unsigned i = 0; i < size()/2; i++) {
				vec2& a = (*this)[i];
				vec2& b = (*this)[size()-1-i];
				swap(a,b);
			}
			winding = winding == WO_CW ? WO_CCW : WO_CW;
		}

		void push_back(Edge* e) {
			assert(!closed);

			if (size() <= 2) {
				if (size() < 2) {
					// store first two nodes
					super::push_back(e->a->pos);
					super::push_back(e->b->pos);
					return;
				} else /* size() == 2 */ {

					// determine correct order of first two nodes
					iterator it_0 = begin();
					iterator it_1 = ++begin();
					if (e->a->pos == *it_0) {
						swap(*it_0, *it_1);
					} else if (e->b->pos == *it_0) {
						swap(*it_0, *it_1);
					}
				}
			}

			Node* n;
			if (e->b->pos == back()) {
				n = e->a;
			} else {
				assert(e->a->pos == back());
				n = e->b;
			}
			if (n->pos == front()) closed = true;
			else super::push_back(n->pos);
		}
	};

	template<class OUTPUT>
	void triangulate(OUTPUT& out) {

		advance(STATUS_AREAS_MONOTONIZED);


		SimpleTriangulation T;
		typedef vector<Polygon*> Polygons;
		Polygons polygons;
		for (Areas::iterator it_a = areas.begin(); it_a != areas.end(); it_a++) {
			Area* area = *it_a;
			if (area->type(areaRule) == AREA_FACE) {
				polygonize(area, polygons, winding);
				triangulate(polygons, out);
			}
		}
	}

	void triangulate(vector<Polygon*> polygons, vector<vec2>& output) {
		for (vector<Polygon*>::iterator it_p = polygons.begin(); it_p != polygons.end(); it_p++) {
			triangulate(*(*it_p), output);
		}
	}


	struct SimpleTriangulation {

		struct SortLessI {
			// array with vertices
			vec2* v;
			SortLessI(vec2* _v) {
				this->v = _v;
			}
			bool operator () (unsigned i_a, unsigned i_b) {
				return less(v[i_a], v[i_b]);
			}
		};

		WindingOrder winding;

		// index of first of lower side vertices (smallest x)
		unsigned split;
		// number of vertices of the polygon
		unsigned size;

		/** array of vertices, containing:
		 * - [0,split[    : sequence of upper side vertices
		 * - [split,size[ : sequence of lower side vertices
		*/
		vec2* v;


		SimpleTriangulation()
		{
			winding = WO_CCW;
			split = -1;
			size = 0;
			v = NULL;
		}

		~SimpleTriangulation() {
			reset();
		}

		void init(Polygon& polygon) {
			reset();
			winding = polygon.winding;
			size = polygon.size();
			v = new vec2[size];

			// find extrema (x->y ordered)
			// index of smallest v
			unsigned smallest = 0;
			// index of largest v
			unsigned largest = 0;
			for (unsigned i = 0; i < polygon.size(); i++) {
				if (less(polygon[i], polygon[smallest])) smallest = i;
				if (less(polygon[largest], polygon[i])) largest = i;
			}

			// copy vertices int v and initialise this->split
			// to achieve [largest:smallest[ [smallest:largest[
			vec2* data = polygon.data();
			unsigned upper_start = largest;
			unsigned lower_start = smallest;
			if (winding == WO_CW) {
				upper_start = (smallest < size-1)? smallest+1 : 0;
				lower_start = (largest < size-1)? largest+1 : 0;
			}
			if (upper_start < lower_start) {
				split = (lower_start-upper_start);
				unsigned chunk_size = polygon.size()-upper_start;
				// copy data[largest:smallest[ and data[smallest:size[ to v
				vec2* to = v;
				memcpy(to, data + upper_start, chunk_size * sizeof(vec2));
				// append data[0:largest[ to v
				to += chunk_size;
				chunk_size = upper_start;
				memcpy(to, data, chunk_size * sizeof(vec2));
			} else /* lower_start < upper_start */ {
				split = size - upper_start + lower_start;
				unsigned chunk_size = (polygon.size()-upper_start);
				// copy data[largest:size[ to v
				vec2* to = v;
				memcpy(to, data + upper_start, chunk_size * sizeof(vec2));
				// append data[0:smallest[ and [smallest:largest[ to v
				to += chunk_size;
				chunk_size = upper_start;
				memcpy(to, data, chunk_size * sizeof(vec2));
			}
		}



		void reset() {
			if (v) {
				delete [] v;
				v = NULL;
			}
		}


		bool differentSides(unsigned a, unsigned b) {
			return (isLower(a) != isLower(b));
		}

		unsigned Top(vector<unsigned>& S) {
			return S.back();
		}

		void Push(vector<unsigned>& S, unsigned p) {
			S.push_back(p);
		}

		unsigned Pop(vector<unsigned>& S) {
			unsigned result = S.back();
			S.pop_back();
			return result;
		}

		template<class OUTPUT>
		void triangulate(Polygon& polygon, OUTPUT& out) {
			init(polygon);
			// stack with indices to vertices in v
			vector<unsigned> S;
			// list of indices into v, sorted by v[p[i]] < v[p[i+1]]
			vector<unsigned> p;

			// fill p with indices
			//prepare_strict(p);

			prepare_nonstrict(p);


			Push(S, p[0]);
			Push(S, p[1]);
			unsigned k = size-1;
			for (unsigned i = 2; i < k; i++) {
				if (differentSides(p[i], Top(S))) {
					unsigned p_prev = S[0]; // first may or may not be on different sides
					for (unsigned j = 1; j < S.size(); j++) {
						unsigned p_j = S[j];
						// all (p_j , p[i]) are on different sides
						triangle(p_prev, p_j, p[i], out);
						p_prev = p_j;
					}
					S.clear();
					assert(p_prev == p[i-1]);
					// these entries are on different sides,
					// because p[i-1] was Top(S) when we
					// entered this branch.
					Push(S, p[i-1]);
					Push(S, p[i]);
				}
				else /* p[i] and Top(S) are on the same side */
				{
					// all entries in S [1:n] are on the same side as p[i]
					// S[0] may or may not be on the same side with p[i]
					unsigned point = Pop(S);
					while (S.size() && try_triangle(Top(S), point, p[i], out)) {
						point = Pop(S);
					}
					Push(S, point);
					Push(S,p[i]);
				}
			}
			if (S.size() >= 2) {
				for (unsigned i = 1; i < S.size(); i++) {
					triangle(S[i-1], S[i], p[k], out);
				}
			}
			reset();
		}

		void prepare_strict(vector<unsigned>& p) {
			p.reserve(this->size);
			for (unsigned i = 0; i < this->size; i++) {
				p.push_back(i);
			}

			// sort indices by v[a] < v[b]
			SortLessI cmp(v);
			sort(p.begin(), p.end(), cmp);
		}


		void prepare_nonstrict(vector<unsigned>& p) {
			// sort indices by v[a] < v[b]
			SortLessI less(v);

			if (winding == WO_CCW) {
				unsigned upper = split-1; // smallest upper
				unsigned lower = split;   // smallest lower
				unsigned i = 0;
				for (; i < size && upper >= 0 && lower < size; i++) {
					// consistency checks
					assert(isUpper(upper) && isLower(lower));

					if (less(lower, upper)) {
						p.push_back(lower);
						lower++;
					} else {
						p.push_back(upper);
						upper--;
					}
				}
				if (upper >= 0) {
					for (; i < size; i++) {
						assert(upper >= 0);
						p.push_back(upper);
						upper--;
					}
				}
				else {
					for (; i < size; i++) {
						assert(lower < size);
						p.push_back(lower);
						lower++;
					}
				}
			} else {
				unsigned upper = 0;        // smallest upper
				unsigned lower = size-1;   // smallest lower
				unsigned i = 0;
				for (; i < size && upper < split && lower >= split; i++) {
					// consistency checks
					assert(isUpper(upper) && isLower(lower));
					if (less(lower, upper)) {
						p.push_back(lower);
						lower--;
					} else {
						p.push_back(upper);
						upper++;
					}
				}
				if (upper < split) {
					for (; i < size; i++) {
						assert(upper < split);
						p.push_back(upper);
						upper++;
					}
				}
				else {
					for (; i < size; i++) {
						assert(lower >= split);
						p.push_back(lower);
						lower--;
					}
				}
			}


		}


		bool isUpper(unsigned i) {
			return i < split;
		}

		bool isLower(unsigned i) {
			return !isUpper(i);
		}


		static bool less(const vec2& a, const vec2& b) {
			return a.x < b.x || (a.x == b.x && a.y < b.y);
		}

#ifndef NDEBUG
		bool lessByOrder(unsigned i, unsigned j) {
			if (isUpper(i) && isUpper(j)) {
				if (winding == WO_CCW) {
					return j < i;
				} else {
					return i < j;
				}
			} else if (isLower(i) && isLower(j)) {
				if (winding == WO_CW) {
					return j < i;
				} else {
					return i < j;
				}
			} else {
				return less(v[i],v[j]);
			}
		}
#endif

		/**
		 * Tests if this triangle describes a face,
		 * inside of the polygon.
		 *
		 * Vertices have to be naturally ordered.
		 * b decides where the inside of the polygon lies
		 * example: if b is at the upper side, inside is below
		 */
		bool isFace(unsigned a, unsigned b, unsigned c) {
			assert(lessByOrder(a, b) && lessByOrder(b,c));
			double result;
			if (isLower(b))
			{
				result = cross_z(v[c]-v[b], v[a]-v[b]);
			}
			else
			{
				result = cross_z(v[a]-v[b], v[c]-v[b]);
			}
			return result >= 0;
		}

		/**
		 * returns true, if a->b has to be reversed according
		 * to winding order.
		 */
		bool reverse(unsigned a, unsigned b) {
			assert(lessByOrder(a, b));
			bool upper = isUpper(b);
			bool cw = (winding == WO_CW);
			bool reverse = upper ^ cw; // upper && !cw || !upper && cw
			return reverse;
		}

		/**
		 * output triangle only if this set of vertices creates a face.
		 * @return whether triangle was emitted
		 */
		template<class OUTPUT>
		bool try_triangle(unsigned a, unsigned b, unsigned c, OUTPUT& out) {
			assert(lessByOrder(a, b) && lessByOrder(b, c));
			if (isFace(a,b,c)) {
				triangle(a,b,c, out);
				return true;
			}
			return false;
		}

		/**
		 * Emits a triangle.
		 * vertices will be reversed if necessary according to winding order.
		 */
		template<class OUTPUT>
		void triangle(unsigned a, unsigned b, unsigned c, OUTPUT& out) {
			vec2& v_a = v[a];
			vec2& v_b = v[b];
			vec2& v_c = v[c];

			assert(lessByOrder(a, b) && lessByOrder(b, c));
			assert(isFace(a,b,c));

			if (reverse(a,b)) {
				out.push_back(v_c);
				out.push_back(v_b);
				out.push_back(v_a);
			} else {
				out.push_back(v_a);
				out.push_back(v_b);
				out.push_back(v_c);
			}
		}


	};


	/**
	 * Simple triangulation of an X-monotone polygon.
	 * Requirements:
	 * 1. findAreas()
	 * 2. mergeAreasOfSameType()
	 * 3. monotonize_X()
	 */
	void triangulate(Polygon& polygon, vector<vec2>& output) {
		SimpleTriangulation T;

		T.triangulate(polygon, output);
	}





	template<class POLYGON = Polygon>
	void polygonizeFaces(vector<POLYGON*>& polygons, WindingOrder winding) {
		for (Areas::iterator it = areas.begin(); it != areas.end(); it++) {
			Area* area = *it;
			if (area->type(areaRule) == AREA_FACE) {
				polygonize<POLYGON>(area, polygons, winding);
			}
		}
	}


	/**
	 * Requires: findAreas(), monotonize_X()
	 * Returns list of polygons for given area.
	 * Each polygon is ordered list of edges.
	 */
	template<class POLYGON = Polygon>
	void polygonize(Area* area, vector<POLYGON*>& polygons, WindingOrder winding) {

		AreaEdges splits;
		AreaEdges work;
		AreaEdge* e = *area->edges.begin();
		work.push_back(e);
		if (isSplitEdge(e)) {
			splits.push_back(e);
		}


		while (work.size()) {
			AreaEdge* edge = work.back(); work.pop_back();
			Node* start;
			Node* prev_node;
			Node* node;
			if (edge->area_left == area) {
				if (isSplitEdge(edge)) {
					// reset area, to indicate, that we
					// have used this split edge left area
					edge->area_left = NULL;
				}
				if (winding == WO_CCW) start = edge->a;
				else start = edge->b;
			} else {
				// not a split node
				assert(edge->area_right == area);
				if (winding == WO_CCW) start = edge->b;
				else start = edge->a;
			}

			POLYGON* polygon = new POLYGON(winding);
			polygons.push_back(polygon);
			polygon->push_back(edge);

			prev_node = start;
			node = edge->opposite(prev_node);
			AreaEdge* prev_edge = edge;
			do {

				EdgeRefs::iterator it_er_prev = node->findEdgeRef(prev_node);
				assert(it_er_prev != node->edgerefs.end());
				EdgeRefs::iterator it_er_cur;
				if (winding == WO_CCW) it_er_cur = node->clockwise_next(it_er_prev);
				else it_er_cur = node->clockwise_prev(it_er_prev);
				edge = (AreaEdge*)((*it_er_cur)->edge);
				assert(edge->area_left == area || edge->area_right == area);
				assert(edge != prev_edge);
				polygon->push_back(edge);
				if (isSplitEdge(edge)) {
					if (edge->a == node) {
						if (winding == WO_CCW) {
							edge->area_left = NULL;
						} else {
							edge->area_right = NULL;
						}
					} else {
						if (winding == WO_CCW) {
							edge->area_right = NULL;
						} else {
							edge->area_left = NULL;
						}
					}
					splits.push_back(edge);
					work.push_back(edge);
				}

				prev_edge = edge;
				prev_node = node;
				node = edge->opposite(prev_node);
			} while (node != start);
		}

		// reset area_left for all split edges
		for (AreaEdges::iterator it = splits.begin(); it != splits.end(); it++) {
			AreaEdge* e = *it;
			e->area_left = area;
			e->area_right = area;
		}

	}

};



} /* namespace math */

#endif /* MATH2D_TRIANGULATOR_H_ */
