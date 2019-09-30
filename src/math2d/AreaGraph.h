/*
 * MultiGraph.h
 *
 *  Created on: 20 Jun 2019
 *      Author: homac
 */

#ifndef MATH2D_AREAGRAPH_H_
#define MATH2D_AREAGRAPH_H_

#include "math2d-config.h"

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
#include <math2d-config.h>

#include <math2d/Graph2DPlanar.h>


struct Test_AreaGraph;

namespace math2d {

// TODO: AreaGraph consider using adjacency matrix
// TODO: AreaGraph consider using integer only
// TODO: performance: floating point causes to many inaccuracies,
//       which have to be solved through branches.
//       It also slows down comparison between vectors.
//       Floating point also adds too many nodes due to inaccuracy.
// TODO: performance: Use ActiveEdgeList to reduce number of intersection tests
class AreaGraph : public Graph2DPlanar {
	typedef Graph2DPlanar super;
	friend struct ::Test_AreaGraph;

public:


	enum WindingOrder {
		WO_CW,
		WO_CCW
	};

	enum AreaRule {
		AR_EVEN_ODD,
		AR_NON_ZERO,
	};

	enum AreaType {
		AREA_HOLE = false,
		AREA_FACE = true
	};


	struct Area;

	/** we extend edges by area adjacency */
	struct AreaEdge : super::Edge {

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

		/** area to the left according to direction a->b */
		Area* area_left;
		void addToAreaLeft(Area* area) {
			assert(area_left == NULL);
			area_left = area;
			area->add(this);
		}
		/** area to the right according to direction a->b */
		Area* area_right;
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
		Neighbours neighbours;

		Area* parent;
		AreaEdges outline;


#ifndef NDEBUG
		static int total_areas_in_system;
#endif

	public:

		Area(int _weight) {
			this->weight = _weight;
			this->parent = NULL;
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
			for (Neighbours::iterator it_n = neighbours.begin(); it_n != neighbours.end(); it_n++) {
				Area* neighbour = *it_n;
				neighbour->neighbours.erase(this);
			}
			neighbours.clear();
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
		 * ATTENTION: Updates area assignment of edges, but not neighbourship of areas.
		 * @see updateNeighbourship()
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



		int mergeSameNeighbours(AreaRule areaRule, AreaEdges& removed) {
			// requires neighbourship info
			assert(AREA_VOID->neighbours.size() > 0);


			int merged = 0;
			Area* area = AREA_VOID;
			AreaType type = area->type(areaRule);

			vector<Area*> merge;
			for (Neighbours::iterator it_n = area->neighbours.begin(); it_n != area->neighbours.end(); it_n++) {
				Area* neighbour = *it_n;
				if (neighbour->type(areaRule) == type) {
					merge.push_back(neighbour);
				}
			}

			for (vector<Area*>::iterator it = merge.begin(); it != merge.end(); it++) {
				Area* neighbour = *it;
				mergeNeighbours(area, neighbour, removed);
				remove(neighbour);
				merged++;
			}
			merge.clear();


			for (Areas::iterator it_a = begin(); it_a != end(); it_a++) {
				Area* area = (*it_a);
				AreaType type = area->type(areaRule);
				for (Neighbours::iterator it_n = area->neighbours.begin(); it_n != area->neighbours.end(); it_n++) {
					Area* neighbour = *it_n;
					if (neighbour->type(areaRule) == type) {
						merge.push_back(neighbour);
					}
				}

				if (merge.size()) {
					for (vector<Area*>::iterator it = merge.begin(); it != merge.end(); it++) {
						Area* neighbour = *it;
						mergeNeighbours(area, neighbour, removed);
						remove(neighbour);
						merged++;
					}
					merge.clear();
					// since we have altered the set, we have to reset the iterator
					it_a = find(area);
					assert(it_a != end());
				}

			}
			return merged;
		}


		/**
		 * Merges disjunct edges of B into A, and removes all splitting edges between them.
		 * Requires neighbourship info.
		 * Returns a list of removed edges.
		 */
		void mergeNeighbours(Area* A, Area* B, AreaEdges& removed) {
			// requires neighbourship info
			assert(A->neighbours.size() > 0);

			AreaEdges merged;

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
					removed.push_back(e);
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


			mergeNeighbourship(A, B);
		}

		/**
		 * Merge neighbours and remove common edges.
		 * Does not consider neighbourship!
		 */
		void mergeNeighbours(Area* A, Area* B) {
			// don't use with neighbourship info!
			assert(A->neighbours.size() == 0);
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


		/** Turns all neighbours of B into neighbours of A */
		void mergeNeighbourship(Area* A, Area* B) {
			for (Neighbours::iterator it_n = B->neighbours.begin(); it_n != B->neighbours.end(); it_n++) {
				Area* neighbour = *it_n;
				if (neighbour != A) {
					neighbour->neighbours.erase(B);
					neighbour->neighbours.emplace(A);
				}
			}
			A->neighbours.erase(B);
			B->neighbours.erase(A);
			A->neighbours.insert(B->neighbours.begin(), B->neighbours.end());
		}



	};



private:

	AreaRule areaRule;
	WindingOrder winding;


	Areas areas;

	/** reference on super::edges of type Edges stored in base class but we interpret them as AreaEdges */
	AreaEdges& edges;

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


	AreaGraph() : AreaGraph(WindingOrder::WO_CCW, AreaRule::AR_NON_ZERO, false, 0) {
	}
	AreaGraph(WindingOrder winding, AreaRule areaRule)
	: AreaGraph(winding, areaRule, false, 0)
	{
	}
	AreaGraph(WindingOrder winding, AreaRule areaRule, bool keepAdjacentAreas, int min_nodes) : super(min_nodes, keepAdjacentAreas), edges(forced_cast<AreaEdges>(super::edges)) {
		assert(assertEdgeRefCompliance());
		this->winding = winding;
		this->areaRule = areaRule;
	}


	virtual ~AreaGraph() {
		reset();
	}

	void reset() {
		super::reset();
		areas.reset();
	}


public:

	struct AELEntry {
		// TODO: move to appropriate place
		float x;
		double orientation;
		AreaEdge* e;


		AELEntry(AreaEdge* e) {
			init(e);
		}

		void init(AreaEdge* _e) {
			e = _e;
			x = e->a->pos.x;
			vec2 b_a = normalize(e->b->pos - e->a->pos);
			orientation = pseudo_orientation_e2(b_a);
		}


	};


	struct AELEntry_less {
		bool operator() (AELEntry * const a, AELEntry * const b) {
			double diff = a->x - b->x;
			if (diff == 0) {
				diff = a->orientation - b->orientation;
			}
			return diff < 0;
		}
	};


	/**
	 * A list of AELEntries e ordered by [e->x, e->orientation]
	 */
	struct AEL : set<AELEntry*, AELEntry_less> {
		typedef set<AELEntry*, AELEntry_less> super;
		key_compare LT;


		iterator put (AreaEdge* e) {
			pair<iterator, bool> result = super::insert(new AELEntry(e));
			assert(result.second);
			return result.first;
		}

		iterator put_behind (AreaEdge* e, iterator it_prev) {
			iterator hint(it_prev);
			if (hint != end()) ++hint;

			return super::emplace_hint(hint, new AELEntry(e));
		}

		iterator find(AreaEdge* e) {
			AELEntry entry(e);
			return super::lower_bound(&entry);
		}

		iterator remove (iterator it) {
			iterator result = erase(it);
			delete *it;
			return result;
		}

		iterator prev (AreaEdge* e) {
			return prev(find(e));
		}

		iterator prev (iterator it) {
			if (it == begin()) return end();
			else return --iterator(it);
		}
	};



	Areas* findAreas() {
		areas.reset();
		Area* AREA_VOID = areas.AREA_VOID;
		AEL ael;
		AreaEdge* last_horizontal;
		float scan_y;


		// INITIALISATION
		if (nodes.size() == 0) return &areas;
		scan_y = (*nodes.begin())->pos.y-1;
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
				// 3. merge areas ending at this scan line, when conflicting
				//


				AreaEdge* prev = NULL;
				AreaEdge* edge = NULL;
				for (AEL::iterator it_ael = ael.begin(); it_ael != ael.end();) {
					AELEntry* entry = (*it_ael);
					edge = entry->e;

					if (edge->b->pos.y == scan_y) {
						// end node on scan line -> edge will be removed

						// merge areas at convex nodes
						if (prev == NULL || prev->b != edge->b) {
							// First edge of at least one edges ending in the same 
							// node n_b on height scan_y.

							// For edges ending in the same node (to be removed
							// further down here) we have to check now for a
							// special case. Two of them may build the border
							// of a concave area above and to the left and right
							// of this node. Up to this point, it was unknown
							// that the areas to the left and right will finally
							// connect in this point and have been treated as
							// separate areas instead.
							//
							// Thus, the two edges coming from the lower
							// left and the lower right have been
							// assigned different areas. If they are not separated
							// by another edge between them, those areas have to
							// be merged.

							// Check leftmost and rightmost outer edges of b pointing
							// downward. If they are neighbouring by angle
							// (left->right, clockwise), then they may have been
							// assigned to different areas for the outer area.
							Node* n_b = edge->b;
							// Find next edge in clockwise orientation
							AreaEdge* right = (AreaEdge*)n_b->getCWNeighbour(edge);
							if (right->b == edge->b) {
								// either both pointing down or neighbour is horizontal left
								if (right->area_right == NULL) {
									// this is a horizontal edge to the left of this node
									// its area is not yet assigned
								} else if (edge->area_left != right->area_right) {
									// these two edges point downward
									areas.mergeEdges(edge->area_left, right->area_right);
								}
							}
						}

						// remove edges ending in this node and update iterator
						it_ael = ael.remove(it_ael);
					} else {
						float x;
						scanline_horizontal_intersection_m_below(scan_y, edge->a->pos, edge->b->pos, x);
						// hit only, if set contains crossing edges
						entry->x = x;
						assert(it_ael == ael.begin() || ael.LT(*ael.prev(it_ael), entry));

						// advance in list
						it_ael++;
					}
					prev = edge;
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
			if (it_er > it_er_end) {
				// all edges pointing downward
				// --> nothing to do
				continue;
			}

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
				// find the next edge of this node in clockwise order
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

		} // nodes finished



		return &areas;
	}

	/**
	 * set neighbourship for each area based on its edges
	 */
	void determineNeighbourship() {

		for (AreaEdges::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			AreaEdge* e = *it_e;
			e->area_left->neighbours.emplace(e->area_right);
			e->area_right->neighbours.emplace(e->area_left);
		}
	}

	/**
	 * Requires determineNeighbourship().
	 *
	 * Use this method only, if you need neighbourship information anyway.
	 */
	int mergeAreasOfSameType_() {
		// merge areas of the same type and cleanup internal structures
		AreaEdges removed;
		int merged = areas.mergeSameNeighbours(areaRule, removed);
		for (AreaEdges::iterator it_e = removed.begin(); it_e != removed.end(); it_e++) {
			removeEdge((*it_e));
		}
		removed.clear();
		super::removeEmptyNodes();
		return merged;
	}


	AreaEdges::iterator removeAreaEdge(AreaEdges::iterator it_ae) {
		Edges::iterator it_e = forced_cast<Edges::iterator>(it_ae);
		it_e = super::removeEdge(it_e);
		it_ae = forced_cast<AreaEdges::iterator>(it_e);
		return it_ae;
	}

	/** Does not consider neighbourship in any way! */
	int mergeAreasOfSameType2() {
		// don't use it, if you need neighbourship info!
		assert(areas.AREA_VOID->neighbours.size() == 0);


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

		return merged;
	}


	void applyWindingOrder() {
		// TODO: test
		// Apply winding order for all edges
		for (AreaEdges::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			AreaEdge* edge = *it_e;
			edge->applyWindingOrder(areaRule, winding);
		}
	}

	/**
	 * Turns network of adjacent areas into a tree of nested areas.
	 *
	 * Precondition: Areas of the same type must have been merged
	 *               beforehand (mergeAreasOfSameType().
	 *
	 * Each may have nested areas, which have a different type as
	 * the parent area. All adjacent areas are guaranteed to have
	 * different types.
	 */
	int determineOutlinesAndNesting() {

		// TODO: test or delete

		int merged = 0;
		// For the remaining areas, the following statements are true:
		// 1. Any neighbour of an area has a different type as the area itself.
		// 2. Areas of the same type can never have a common edge.
		// 3. Areas of the same type can only meet in a common node.
		// 4. Each edge separates two areas of different type.
		// 5. Nodes have a multiple of 2 edges.

		// An area can only have nested areas of different type.
		// Thus, nested areas must be either truly separated, or
		// touch each other in at least one node
		// which is also shared by the parent.
		//
		// Examples:
		// |----*----|
		// |   / \   |
		// |   \ /   |
		// |----*----|
		// |----*----|
		// |   / \___|
		// |   |  ---*
		// |   \ /   |
		// |----*----|
		//





		// To determine a clean nesting relation, two alternatives are possible:
		//
		// A)  All children, which share a node with each other and its
		//     parent, can be merged. --> The set of children of one area
		//     is the set of all closed paths, consisting of edges, shared
		//     with the parent area.
		// B)  All adjacent areas of any node shared between child and parent
		//     becomes a child of the parent. If this new child is of the same
		//     type, it will be merged with the parent.


		// Beginnend bei Void:
		// Für jeden Nachbarbereich, der noch nicht bearbeitet wurde:
		//
		// Suche eine irgendeine Start-Kante, die zu dem Bereich gehört.
		//
		// Durchlaufe Kanten eines Bereichs gemäß WindingOrder und Bereichstyp
		// und suche für jede Kante im Endknoten die jeweils fortsetzende Kante
		// entlang desselben Bereichs.
		//
		//
		// Füge alle jene Kanten zu einer Outline hinzu und merke dir, welche Bereiche
		// auf der anderen Seite liegen in einem set<Area*>.
		//
		// Alle hinzugefügten Kanten, entferne aus den Kanten der Parent-Area.
		//
		// Wenn ein Kantenzug geschlossen wurde, dann markiere alle Bereiche,
		// die nun im set sind, als bereits bearbeitet und merge sie miteinander.
		// Der gemergte Bereich ist nun des neue Kind. Entferne den Parent-Bereich aus
		// der Liste der Nachbarbereiche des Kindes und füge ihn als Parent hinzu.
		// Füge alle Kanten der Outline zu seiner outline hinzu und entferne die
		// outline-kanten aus der edge-list. Übrig bleiben also die Kanten der Kindeskinder.
		//
		// Fahre mit den verbleibenden Nachbarbereichen desselben Parents fort, bis
		// alle Bereiche als bearbeitet markiert wurden.
		//
		// Wenn Parent abgeschlossen, verfahre mit den Kindern auf die gleiche Weise.
		// Beachte, dass die Kantenlisten der Kinder nun nur noch die Kanten mit
		// deren Kindern beinhaltet.


		queue<Area*> work;
		work.push(areas.AREA_VOID);


		while (work.size() > 0) {

			Area* parent = work.front(); work.pop();

			Neighbours done;
			Neighbours children(parent->neighbours);
			Neighbours associated;
			for (Areas::iterator it_n = children.begin(); it_n != children.end(); it_n++) {
				Area* child = *it_n;

				// find an outline and merge all associated children
				if (done.find(child) != done.end()) {
					AreaEdges outline;
					getOutline(parent, child, outline, associated);

					// merge associated into one area (neighbour)
					// and mark associated areas as done
					for (Neighbours::iterator it_an = associated.begin(); it_an != associated.end(); it_an++) {
						areas.mergeEdges(child, *it_an);
						areas.mergeNeighbourship(child, *it_an);
						done.emplace(*it_an);
						merged++;
					}

					// turn neighbour into a child
					child->neighbours.erase(parent);
					child->parent = parent;

					// assign outline
					child->outline = outline;

					// remove outline edges from parent and child
					child->edges.removeAll(outline);
					parent->edges.removeAll(outline);

					outline.clear();

					// mark neighbour as done
					done.emplace(child);
				}
			}

			// empty temporary storage
			done.clear();
			children.clear();
			associated.clear();

			// add children to work queue
			// note that the set of children was reduced during this loop iteration
			for (Neighbours::iterator it_n = parent->neighbours.begin(); it_n != parent->neighbours.end(); it_n++) {
				work.push(*it_n);
			}

			assert(parent->edges.size() == 0);
		}
		return merged;
	}


	/** 
	 * Requires: determineNeighbourship(), mergeAreasOfSameType() and applyWindingOrder()
	 *
	 * Determines the path, which separates the parent area from a set of connected
	 * child areas. This path is the outline of the set of the children, which can
	 * be merged.
	 * 
	 * Outputs: The closed path in 'outline' and the set of 'children' enclosed in
	 *          this path.
	 *
	 * @param parent The parent area.
	 * @param child  The first child area, adjacent to the parent.
	 *               Required to find a start edge with adjacent areas
	 *               parent<->child.
	 * @param outline Output: list of edges of the outline.
	 * @param children Output: Set of children enclosed by the outline.
	 */
	void getOutline(Area* parent, Area* child, AreaEdges& outline, Neighbours& children) {
		// TODO: test! or delete
		//
		// find any first edge
		//
		AreaEdge* edge = NULL;
		// area is in application given edge direction on the left side
		bool left = false;
		for (AreaEdges::iterator it_e = parent->edges.begin(); it_e != parent->edges.end(); it_e++) {
			edge = *it_e;
			if (edge->area_left == parent && edge->area_right == child) {
				left = (edge->weight > 0);
				break;
			} else if (edge->area_right == parent && edge->area_left == child) {
				left = (edge->weight < 0);
				break;
#ifndef NDEBUG
			} else {
				edge = NULL;
#endif
			}
		}
		assert(edge != NULL);

		// now follow along to the next edge

		Node* start = edge->getStartByWeight();
		for (Node* next = edge->opposite(start); next != start; next = edge->opposite(next)) {

			// find next edge from next node

			EdgeRefs::iterator it_ref = next->findEdgeRef(edge->b);
			assert ((*it_ref)->edge == (Edge*)edge);
			EdgeRef* ref = *(next->clockwise_next(it_ref));
			edge = (AreaEdge*)ref->edge;
			if (parent != edge->getAreaByWeight(left)) {
				ref = *(next->clockwise_prev(it_ref));
				edge = (AreaEdge*)ref->edge;
			}
			assert(parent == edge->getAreaByWeight(left));

			// remember other area
			Area* other = edge->getAreaByWeight(!left);
			if (other != child) {
				children.emplace(other);
				child = other;
			}

			// append edge to outline
			outline.push_back(edge);

		}


	}

};



} /* namespace math */

#endif /* MATH2D_AREAGRAPH_H_ */
