/*
 * MultiGraph.h
 *
 *  Created on: 20 Jun 2019
 *      Author: homac
 */

#ifndef MATH2D_GRAPH2DPLANAR_H_
#define MATH2D_GRAPH2DPLANAR_H_

#include "math2d-config.h"
#include <iostream>

#include <string>
#include <stdio.h>

#include <assert.h>
#include <stdlib.h>

#include <glm/gtx/vector_angle.hpp>
using namespace glm;


#include <limits>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <algorithm>
using namespace std;


#include <math2d/line.h>
#include <math2d-config.h>



namespace math2d {





// TODO: Graph2DPlanar remove line-point intersections
// TODO: Graph2DPlanar consider using integer only
// TODO: performance: floating point causes to many inaccuracies,
//       which have to be solved through branches.
//       It also slows down comparison between vectors.
//       Floating point also adds too many nodes due to inaccuracy.
// TODO: performance: Use ActiveEdgeList to reduce number of intersection tests
class Graph2DPlanar {
public:
	static constexpr const float_mantissa_mask_t GRAPH_FLOAT_PRECISION = float_mantissa_mask(6);
	static constexpr const float FLOAT_TOLERANCE = (0.00001f);
	struct Edge;
	typedef std::vector<Edge*> Edges;
	typedef typename Edges::iterator edge_iterator;

	struct EdgeRef;
	typedef std::vector<EdgeRef*> EdgeRefs;


	struct Node;

	typedef std::vector<Node*> Nodes;

	typedef typename Nodes::iterator node_iterator;

	struct Area;


	struct Node {
		vec2 pos;

		EdgeRefs edgerefs;

#ifndef NDEBUG
		static int total_nodes_in_system;
#endif

		Node() {
			init(vec2(0,0));
		}

		Node(const vec2& _pos)	{
			init(_pos);
		}

		void init(const vec2& _pos) {
			this->pos = _pos;
#ifndef NDEBUG
			total_nodes_in_system++;
#endif
		}


		virtual ~Node() {
#ifndef NDEBUG
			total_nodes_in_system--;
#endif
			for (EdgeRefs::iterator e = edgerefs.begin(); e != edgerefs.end(); e++) {
				EdgeRef* edgeref = *e;
				delete edgeref;
			}
		}

		float& x() {return pos.x;}
		float& y() {return pos.y;}

		unsigned grad() {
			return edgerefs.size();
		}

		bool equals(const Node& b) {
			return pos == b.pos;
		}

		void addEdge(Edge* edge) {
			EdgeRef* ref = new EdgeRef(edge->opposite(this), edge);
			EdgeRefs::iterator e = findEdgeLowerBound(ref);
#ifndef NDEBUG
			EdgeRef* foundref = (e == edgerefs.end())?NULL:(*e);
			// hit if edgeref is not unique
			assert(e == edgerefs.end() || foundref->to != ref->to);
#endif
			edgerefs.insert(e, ref);
		}

		void removeEdge(Edge* edge) {
			EdgeRef ref(edge->opposite(this), edge);
			EdgeRefs::iterator e = findEdgeLowerBound(&ref);
			// hit if edgeref was not in list
			assert(e != edgerefs.end() && (*e)->to == ref.to);
			delete (*e);
			edgerefs.erase(e);
		}

		Edge* getCWNeighbour(Edge* edge) {
			EdgeRefs::iterator e = findEdgeRef(edge->opposite(this)); assert((*e)->edge == edge);
			e = clockwise_next(e);
			return (*e)->edge;
		}

		Edge* getCCWNeighbour(Edge* edge) {
			EdgeRefs::iterator e = findEdgeRef(edge->opposite(this));
			e = clockwise_prev(e);
			return (*e)->edge;
		}

		EdgeRefs::iterator clockwise_next(EdgeRefs::iterator it) {
			if (++it == edgerefs.end()) it = edgerefs.begin();
			return it;
		}

		EdgeRefs::iterator clockwise_prev(EdgeRefs::iterator it) {
			if (it == edgerefs.begin()) it = edgerefs.end();
			return --it;
		}

		/**
		 * Returns first edge, which is greater equal to orientation.
		 * warning: requires edges to be unique by angle */
		EdgeRefs::iterator findEdgeOriented(double orientation) {
			// TODO: make sure edges of one node are unique by angle
			EdgeRef ref;
			ref.orientation = orientation;
			// ref.to = NULL; <-- not used
			EdgeRefs::iterator e = lower_bound(edgerefs.begin(), edgerefs.end(), &ref, EdgeRef_oriented_lt);
			return e;
		}

		EdgeRefs::iterator findEdgeLowerBound(EdgeRef* ab) {
			EdgeRefs::iterator e = lower_bound(edgerefs.begin(), edgerefs.end(), ab, EdgeRef_lt);
			return e;
		}
		EdgeRef* findEdgeExact(Node* to) {
			EdgeRefs::iterator e = findEdgeRef(to);
			return *e;
		}
		EdgeRefs::iterator findEdgeRef(Node* to) {
			Edge ab(this, to, 0);
			EdgeRef ref(to, &ab);
			EdgeRefs::iterator e = findEdgeLowerBound(&ref);
			for (; e != edgerefs.end() && (*e)->to != ref.to; e++) {
				if ((*e)->orientation != ref.orientation) {
					e = edgerefs.end();
					break;
				}
			}
			return e;
		}
	};



	/**
	 * An edge is the connection between to points a and b.
	 *
	 * To maintain uniqueness of edges, each edge will be added
	 * to the graph exactly once. To achieve that, the following
	 * constraint has to hold for each edge:
	 * 		a < b
	 *
	 * Given that constraint holds, each edge is added to the list
	 * of unique edges of the graph.
	 *
	 * The edges connected to a node are stored as references to the
	 * graphs unique edges.
	 *
	 * The application given direction of an edge depends on its weight.
	 * - weight  > 0: a -> b
	 * - weight  < 0: b -> a
	 * - weight == 0: no edge here
	 *
	 * When adding multiple edges between the same pair of nodes, each
	 * edge will add to the weight:
	 *   +1 if edge points towards b
	 *   -1 if edge points towards a
	 * The accumulated weight can be used apply the even-odd and non-zero
	 * rule to polygon areas. In both cases, the edge can be entirely
	 * ignored, if its weight is zero.
	 *
	 * However, if adjacent areas of the same kind (face/hole) should be
	 * removed or kept.
	 */
	struct Edge {
		Node* a;
		Node* b;

		/**
		 * Weight is an indicator for the set of edges between the same
		 * pair of nodes when applying the even-odd or non-zero rule.
		 *
		 * Edges may have the exact opposite direction in comparison to a->b.
		 * - Each edge in the same     direction adds +1 to the weight.
		 * - Each edge in the opposite direction adds -1 to weight.
		 *
		 */
		int weight;

		// TODO: do we still need the edge length?
		float length;
#ifndef NDEBUG

		float orientation;

		static int total_edges_in_system;
#endif




		Edge() {
			a = NULL;
			b = NULL;
			weight = 0;
			length = 0;
#ifndef NDEBUG
			orientation = 0;
			total_edges_in_system++;
#endif
		}


		Edge(Node* a, Node* b, int weight) {
			init(a, b, weight);
#ifndef NDEBUG
			total_edges_in_system++;
#endif
		}

		virtual ~Edge(){
#ifndef NDEBUG
			total_edges_in_system--;
#endif
		}


		virtual void init(Node* _a, Node* _b, int _weight) {
			assert((_a==NULL && _b==NULL) || !_a->equals(*_b));
			if (Node_lt(_a,_b)) {
				this->a = _a;
				this->b = _b;
				this->weight = +_weight;
			} else {
				this->a = _b;
				this->b = _a;
				this->weight = -_weight;
			}

			this->length = glm::length(this->b->pos - this->a->pos);
#ifndef NDEBUG
			this->orientation = pseudo_orientation_e2(b->pos - a->pos);
#endif
		}


		Node* opposite(Node* from) {
			Node * to = a;
			if (from == to) {
				to = b;
			}
#ifndef NDEBUG
			else
			{
				// note: this method only works for unique node pointers!
				assert(from == b);
			}
#endif
			return to;
		}


		bool equals(Edge* that) {
			return this->a->pos == that->a->pos && this->b->pos == that->b->pos;
		}

		void reverse() {
			weight = - weight;
		}

	};

	/**
	 * Edge references are stored at each node, ordered
	 * by the angle between the node and the connected
	 * node.
	 * For each edge there exist two edge references -
	 * one for each node on either end.
	 *
	 * EdgeRefs of one node are ordered by their
	 * outgoing direction (node -> to)
	 * (not the actual application level direction).
	 */
	struct EdgeRef {
		/** reference on connected node. */
		Node* to;
		/** reference on the edge data */
		Edge* edge;

		/** orientation of b-a */
		double orientation;

		EdgeRef() {
			this->to = NULL;
			this->edge = NULL;
			this->orientation = 0;
		}

		EdgeRef(Node* _to, Edge* _edge) {
			init(_to, _edge);
		}

		void init(Node* _to, Edge* _edge) {
			this->to = _to;
			this->edge = _edge;
			Node* from = _edge->opposite(_to);
			vec2 dir(_to->pos - from->pos);

			dir = dir / length();
			this->orientation = pseudo_orientation_e2(dir);
		}

		int weight(Node* from) {
			return edge->a == from ? edge->weight : - edge->weight;
		}

		float length() {
			return edge->length;
		}

	};


	static double vec2_cmp(const vec2& a, const vec2& b) {
		double result = a.y - b.y;
		if (result == 0) result = a.x - b.x;
		return result;
	}


	static bool Node_lt( Node * const a, Node* const b) {
		return vec2_cmp(a->pos, b->pos) < 0;
	}

	static bool EdgeRef_lt(EdgeRef* const a, EdgeRef* const b) {
		// order by (1) outgoing direction (clockwise) and then by (2) targeted node
		double lt = a->orientation - b->orientation;
		if (lt == 0) {
			return Node_lt(a->to,b->to);
		}
		return lt < 0.0;
	}

	static bool EdgeRef_oriented_lt(EdgeRef* const a, EdgeRef* const b) {
		// order by outgoing direction (clockwise)
		double lt = a->orientation - b->orientation;
		return lt < 0.0;
	}

	static bool Edge_lt(Edge* const m, Edge* const n) {
		// global edges unique by node pair(m,n)
		double diff = vec2_cmp(m->a->pos, n->a->pos);
		if (diff == 0) {
			diff = vec2_cmp(m->b->pos, n->b->pos);
		}
		return diff < 0;
	}



protected:

	// unique node references
	Nodes nodes;
	// unique edge list
	Edges edges;
	// whether to keep adjacent areas of the same kind
	bool keepAdjacentAreas;


	typedef map<double, Node*> NodeDistances;




public:


	Graph2DPlanar(bool _keepAdjacentAreas) {
		keepAdjacentAreas = _keepAdjacentAreas;
	}
	Graph2DPlanar(int min_nodes, bool keepAdjacentAreas) {
		this->keepAdjacentAreas = keepAdjacentAreas;
		nodes.reserve(min_nodes);
	}


	virtual ~Graph2DPlanar() {
		reset();
	}

	virtual void reset() {
		for (Nodes::iterator it = nodes.begin(); it != nodes.end(); it++) {
			Node* node = *it;
			delete node;
		}
		nodes.clear();
		for (Edges::iterator it = edges.begin(); it != edges.end(); it++) {
			Edge* edge = *it;
			delete edge;
		}
		edges.clear();
	}


	size_t nodes_size() {
		return nodes.size();
	}
	node_iterator nodes_begin() {
		return nodes.begin();
	}
	node_iterator nodes_end() {
		return nodes.end();
	}

	size_t edges_size() {
		return edges.size();
	}

	edge_iterator edges_begin() {
		return edges.begin();
	}
	edge_iterator edges_end() {
		return edges.end();
	}

public:
	/** adds node if it does not exist.
	 * @return Node representing given coordinate.*/
	Node* addNode(const vec2& a) {
		node_iterator it = findNode(a);
		if (it != nodes.end() && (*it)->pos ==  a) {
			// already in list
			return (*it);
		} else {
			// TODO: optimise: use temp Node
			Node* node = new Node(a);
			// insert node
			nodes.insert(it, node);
			return node;
		}
	}



	void addLineSegment(const vec2& a, const vec2& b) {
		// hit if a == b
		// TODO: bezier: this is legal
		assert(a != b);

		Node* origin = addNode(a);
		Node* target = addNode(b);

		// hit if node was not registered or internal error
		assert(origin != NULL && vec2_cmp(origin->pos, a) == 0);
		assert(target != NULL && vec2_cmp(target->pos, b) == 0);

		Edge* edge = createEdge(origin, target, 1);


		NodeDistances crossed;

		Edges::iterator e_ab = findEdge(edge);
		if (e_ab == edges.end() || !edge->equals(*e_ab)) {
			// edge is not unique
			// 1. check for node crossings
			// 2. check for edge crossings
			// TODO: i doubt it is a good idea to have these two conflicting methods
			// Example: Due to floating point inaccuracy several things may happen:
			// (1) An edge might cross two nodes but not overlap an edge between those nodes.
			// (2) An edge can have multiple crossings with edges very close to a node
			//     instead of crossing the node.
			findNodesCrossed(edge, crossed);
			splitCrossedEdges(edge, crossed);
		}
		if (crossed.size() == 0) {
			applyEdge(edge);
		} else {
			NodeDistances::iterator it_nd = crossed.begin();
			// adjust the edge we have
			Node* last = edge->b;
			edge->init(edge->a, it_nd->second, edge->weight);
			applyEdge(edge);
			Node* prev = it_nd->second;

			// add remaining edges
			// make sure, we don't create edges pointing on itself, due to floating point.
			for (++it_nd; it_nd != crossed.end(); it_nd++) {
				if (prev != it_nd->second) {
					applyEdge(createEdge(prev, it_nd->second, edge->weight));
				}
				prev = it_nd->second;
			}
			if (prev != last) {
				applyEdge(createEdge(prev, last, edge->weight));
			}
		}
	}

	/**
	 *
	 * Removes nodes of degree 0 and edges with at least one node of degree 1.
	 *
	 * Dangling nodes may occur, if edges have been remove, since they
	 * had a weight of zero.
	 *
	 * Dangling edges can only occur, if the given paths were inconsistent
	 * (not closed).
	 *
	 * @return Returns true, if any dangling nodes or edges were found and removed.
	 */
	bool cleanup() {
		// this is not at all optimised because we are dealing here with very rare cases


		unordered_set<Edge*> removeEdges;
		for (Edges::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			Edge* e =(*it_e);
			if (e->a->grad() == 1 || e->b->grad() == 1) {
				removeEdges.emplace(e);
				// remove all edges, which will have a node with degree 1 after removal of e
				Node* n =   (e->a->grad() == 2) ? e->a : ((e->b->grad() == 2) ? e->b : 0);
				while (n) {
					EdgeRefs::iterator it_er = n->edgerefs.begin();
					e = (e != (*it_er)->edge)? (*it_er)->edge : (* ++it_er)->edge;
					removeEdges.emplace(e);
					n = e->opposite(n);
					n = n->grad() == 2 ? n : NULL;
				}
			}
		}

		for (unordered_set<Edge*>::iterator it = removeEdges.begin(); it != removeEdges.end(); it++) {
			removeEdge(*it);
		}

		int removedNodes = removeEmptyNodes();

		return (removeEdges.size() + removedNodes) > 0;
	}

	int removeEmptyNodes() {
		vector<Node*> removeNodes;
		for (Nodes::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++) {
			Node* n = (*it_n);
			if (n->grad() == 0) {
				removeNodes.push_back(n);
			}
		}

		for (vector<Node*>::iterator it = removeNodes.begin(); it != removeNodes.end(); it++) {
			removeNode(*it);
		}
		return removeNodes.size();
	}


public:

	/**
	 * Add edge if it is unique, otherwise apply
	 * edge weight to the existing edge.
	 */
	void applyEdge(Edge* edge) {
		// attention: this is the EDGE or the INSERTION point!
		Edges::iterator e_ab = findEdge(edge);


		if (likely(e_ab == edges.end() || !edge->equals(*e_ab))) {
			// new unique edge
			insertEdge(e_ab, edge);
		} else {
			// edge exists
			(*e_ab)->weight += edge->weight;
			delete edge;

			if ((*e_ab)->weight == 0 && !keepAdjacentAreas) {
				// remove edge
				removeEdge(e_ab);
			}
		}
	}
	/** inserts a unique edge at it_ab */
	void insertEdge(Edges::iterator e_ab, Edge* edge) {
		assert(e_ab == edges.end() || !edge->equals(*e_ab));

		// add references to each node
		edge->a->addEdge(edge);
		edge->b->addEdge(edge);

		// add edge
		edges.insert(e_ab, edge);
	}

	/** inserts a unique edge */
	void insertEdge(Edge* edge) {
		Edges::iterator e_ab = findEdge(edge);
		insertEdge(e_ab, edge);
	}


	// for testing purposes only
	Node* getNode(unsigned index) {
		assert(index >= 0 && index < nodes.size());
		return nodes[index];
	}


	node_iterator findNode(const vec2& a) {
		Node n(a);
		return findNode(&n);
	}

	node_iterator findNode(Node* const n) {
		node_iterator it = lower_bound (nodes.begin(), nodes.end(), n, Node_lt);
		return it;
	}

	Edge* findEdge(const vec2 a, const vec2 b) {
		Nodes::iterator it_na = findNode(a);
		if (it_na == nodes.end()) return NULL;
		Nodes::iterator it_nb = findNode(b);
		if (it_nb == nodes.end()) return NULL;

		Edge ab(*it_na,*it_nb,0);
		Edges::iterator it = findEdge(&ab);
		if (it == edges.end())return NULL;
		else return *it;
	}

protected:
	virtual Edge* createEdge(Node* a, Node* b, int weight) {
		return new Edge(a,b,weight);
	}

	virtual EdgeRef* createEdgeRef(Node* n, Edge* edge) {
		return new EdgeRef(n, edge);
	}


protected:

	void removeNode(Node* n) {
		Nodes::iterator it_n = findNode(n);
		delete n;
		nodes.erase(it_n);
	}

	void removeEdge(Edge* e) {
		Edges::iterator it = findEdge(e);
		assert(it != edges.end() && e->equals(*it));
		removeEdge(it);
	}

	Edges::iterator removeEdge(Edges::iterator it_e) {
		Edge* edge = *it_e;
		edge->a->removeEdge(*it_e);
		edge->b->removeEdge(*it_e);
		delete edge;

		return edges.erase(it_e);
	}


	Edges::iterator findEdge(Edge* edge) {
		return lower_bound(edges.begin(), edges.end(), edge, Edge_lt);
	}

	unsigned splitCrossedEdges(Edge* e, NodeDistances& crossed) {
		Edges add;
		Edges remove;
		vec2 p;
		double s;

		// TODO: optimise crossing edge detection
		// Use scan line algorithm with active edge list (AEL).
		// AEL contains edges, sorted by their x coord. for the y coord.
		// of the current sweep line height.
		// Sweep through nodes from bottom to top and for each node:
		// (a) add all edges going up to the active edge list (AEL).
		// (b) update x coordinate of all edges in AEL for y of current node
		//     -> each edge in AEL, which is now in wrong place crosses
		//        all edges between last and new place according to x-order.
		// (c) for each edge going down
		//     -> remove from AEL
		//     -> test for crossings with edges in AEL.



		// Given that for each edge a < b, e can only have
		// crossing edges c with
		//      c->a.y <= e->b.y and e->a.y <= c->b.y
		//      min(c->_.x) <= max(e->_.x) and max(e->_.x) <= max(c->_.x)
		// Edges are ordered ascending by the following
		// attribute priority: a.y -> a.x -> b.y -> b.x.
		// Thus, the upper bound is c->a <= (max(e->_.x), e->b.y)
		// corners of the rectangle covering e

		vec2 lower_left(fmin(e->a->x(), e->b->x()), e->a->y());
		vec2 top_right(fmax(e->a->x(), e->b->x()), e->b->y());

		Node top_right_a  (top_right);
		Node top_right_b  (vec2(top_right.x,  top_right.y+1));
		Edge top_right_edge (&top_right_a, &top_right_b, 0);
		Edges::iterator it_e_end = findEdge(&top_right_edge);
		if (it_e_end != edges.end()) it_e_end++;

		for (Edges::iterator it_e = edges.begin(); it_e != it_e_end; it_e++) {
			Edge* crossing = (*it_e);

			if (crossing->b->pos.y < lower_left.y) {
				// edge entirely below -> no crossing possible
				continue;
			}

			if (line_segment_intersection_highp(e->a->pos, e->b->pos, crossing->a->pos, crossing->b->pos, p, s)) {
				Node* n = addNode(p);
				if (n != e->a && n != e->b) {
					// add n to list of crossed nodes at distance s
					crossed[s] = n;
					if (n != crossing->a && n != crossing->b) {
						// split crossing edge in n
						remove.push_back(crossing);
						Edge* a_n = createEdge(crossing->a, n, crossing->weight);
						Edge* n_b = createEdge(n, crossing->b, crossing->weight);
#ifndef NDEBUG
						Edges::iterator __it__test = findEdge(a_n);
						assert (__it__test == edges.end() || !(*__it__test)->equals(a_n));
						__it__test = findEdge(n_b);
						assert (__it__test == edges.end() || !(*__it__test)->equals(n_b));
#endif
						add.push_back(a_n);
						add.push_back(n_b);
					}
				}
			}
		}

		for (Edges::iterator it_e = remove.begin(); it_e != remove.end(); it_e++) {
			removeEdge(*it_e);
		}
		for (Edges::iterator it_e = add.begin(); it_e != add.end(); it_e++) {
			// Floating point arithmetic may cause two edges, which are very
			// close to each other fall together into the same place, when split apart.
			applyEdge(*it_e);
		}
		return crossed.size();
	}

	/**
	 * Search all nodes crossed by the given edge e and store nodes in 'crossed'
	 * Crossing nodes are stored by distance to start node of e according to the
	 * direction of e. Thus, the edge can be applied to the crossing nodes in the
	 * given order.
	 *
	 * @return number of crossing nodes found
	 * */
	unsigned findNodesCrossed(Edge* e, NodeDistances& crossed) {

		// only nodes in the area described by the rectangle
		// lower_left(min(ax,bx), min(ay,by)) -> upper_right(max(ax,bx), max(ay,by))
		// may possibly be on the edge a->b

		// nodes a,b of edge e are ordered y-major
		// thus, we know ay < by.
		vec2 lower_left(e->a->x(), e->a->y());
		vec2 upper_right(e->b->x(), e->b->y());
		// swap left/right if ax > bx
		if (lower_left.x > upper_right.x) {
			float tmp = lower_left.x;
			lower_left.x = upper_right.x;
			upper_right.x = tmp;
		}

		// nodes are ordered y-major
		Nodes::iterator it_n = findNode(lower_left);
		assert (it_n != nodes.end());
		Nodes::iterator it_n_end = findNode(upper_right);

		// it_n: all nodes with y>=lower and y <=upper
		// starting with lower_left and ending with upper_right
		for (; it_n != it_n_end; it_n++) {
			Node* n = *it_n;
			if (n == e->a || n == e->b) {
				continue;
			}
			// we could filter here by x-range, but line_contains_point
			// does that as well, so we do not.
			//if (lower_left.x <= n->x() && n->x() <= upper_right.x) {
			double distance;
			if (line_contains_point_highp(e->a->pos, e->b->pos, n->pos, distance, FLOAT_TOLERANCE)) {
				if(distance > 0 && distance < 1) {
					// edge e contains node n
					crossed[distance] = n;
				}
			}


		}

		return crossed.size();
	}


};



} /* namespace math */

#endif /* MATH2D_GRAPH2DPLANAR_H_ */
