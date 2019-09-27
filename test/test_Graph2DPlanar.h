/*
 * test_MultiGraph.h
 *
 *  Created on: 21 Jun 2019
 *      Author: homac
 */

#ifndef MATH_TEST_GRAPH2DPLANAR_H_
#define MATH_TEST_GRAPH2DPLANAR_H_

#include "perf_clock.h"

#include <glm/gtc/constants.hpp>
#include <math2d/Graph2DPlanar.h>

using namespace glm;
using namespace math2d;

void list_content(Graph2DPlanar& g) {

	// FILE* out = stdout;
	FILE* out = fopen("/tmp/g.txt", "w");

	Graph2DPlanar::node_iterator n;
	n = g.nodes_begin();

	fprintf(out, "NODES\n");
	for (;n < g.nodes_end(); n++) {
		fprintf(out, "node: %1.0f,%1.0f\n", (*n)->x(), (*n)->y());
		Graph2DPlanar::EdgeRefs& edges = (*n)->edgerefs;
		Graph2DPlanar::EdgeRefs::iterator e;
		for (e = edges.begin(); e != edges.end(); e++) {
			Graph2DPlanar::Node* to = (*e)->to;
			fprintf(out, "\tedge: (angle=%1.1f) to: (%1.0f,%1.0f) (weight=%d)\n", (*e)->orientation, to->pos.x, to->pos.y, (*e)->weight(*n));
		}
	}
	fprintf(out, "EDGES\n");
	Graph2DPlanar::edge_iterator e = g.edges_begin();
	for (; e != g.edges_end(); e++) {
		Graph2DPlanar::Edge* edge = *e;
		fprintf(out, "edge: (%1.0f,%1.0f) %s (%1.0f,%1.0f)  [%d]\n", edge->a->pos.x, edge->a->pos.y, edge->weight>0 ? "->" : "<-", edge->b->pos.x, edge->b->pos.y, edge->weight);
	}


	if (out != stdout) fclose(out);
}


void test_performance_Graph2DPlanar() {
	//
	// response performance
	//
#define PATH_REALISTIC   1
#define PATH_CRISS_CROSS 2
#define PATH_BEST_CASE   3
#define PATH_TYPE         PATH_REALISTIC

#define PROGRESS_REPORT 1
	const int NODES = 80000;
	// A path has n nodes and n edges
	const int EDGES = NODES;

	typedef vector<vec2> Path;

	const float SECS_PER_OP = 0.00002; // secs per operation (arch specific)
	const float approx = (float(NODES) + float(EDGES+1)*float(EDGES)/2.0) * SECS_PER_OP;
	size_t size = sizeof(Graph2DPlanar::Node)*NODES + sizeof(Graph2DPlanar::Edge)*EDGES + sizeof(Graph2DPlanar::EdgeRef)*EDGES*2;
	cout << "memory: " << size/1024 << " KiB" << endl;
	cout << "worst case time: " << approx << " s" << endl;

	Graph2DPlanar g(NODES, true);

	Path path;

	path.reserve(NODES);
#if PATH_TYPE == PATH_REALISTIC


	/**
	 * how many vertical pairs of vertical (or horizontal) paths to generate
	 *   1 for paths, which are supposed to describe approx. one polygon.
	 *   2 for the outline of a glyph (text font)
	 * < 10 for paths created by a human being.
	 * < 100 for paths intentionally crossing itself a lot.
	 * > 100 to test very rare worst case scenarios.
	 */
	const float NUM_POLYGONS = 1;
	const float MAX_L = 100;  // maximum length of edges
	const float MAX_DIST = float(NODES)/(NUM_POLYGONS*2.0) * MAX_L; // maximum extend of generated shape
	const float rand_fac = MAX_L/RAND_MAX;

	cout << "generating 'realistic' path with approx. " << NUM_POLYGONS << " face(s)" << endl;


#elif PATH_TYPE == PATH_CRISS_CROSS
	cout << "generating 'criss cross' path" << endl;
	const float rand_fac = float(NODES)/RAND_MAX;
#else
	cout << "generating 'best case' path" << endl;
	const float rand_fac = float(NODES)/RAND_MAX;
#endif
	vec2 first(0,0);
	vec2 prev = first;
	path.push_back(prev);

	vec2 prev_dir(1,0);
	vec2 first_dir(1,0);
	int seed = clock();
	srand(seed);

	for (int i = 1; i<NODES-1; i++) {
		vec2 v;
		do {
#if PATH_TYPE == PATH_REALISTIC
			v.x = float(rand())-RAND_MAX;
			v.y = float(rand())-RAND_MAX;
			// the further away from first node and the closer we are to the end
			// the more we want to get closer to first.
			float dist = length(prev-first);
			float f_i = float(i)/NODES;
			float f_first_dir = pow(dist / MAX_DIST  * f_i, 0.5);
			// float f_first_dir = (dist / MAX_DIST  + f_i) / 2;
			// keep the general orientation we had before
			float f_prev_dir = (1.f-f_first_dir) * 0.7;
			float f_v_dir = 1.f - f_first_dir - f_prev_dir;

			float len = float(rand())*rand_fac;

			vec2 v_dir = normalize(
					normalize(v)* f_v_dir
					+ prev_dir  * f_prev_dir
					+ first_dir * f_first_dir);
			v = v_dir*len;
			v = v + prev;


#elif PATH_TYPE == PATH_CRISS_CROSS
			// WORST CASE SCENARIO (CRISS CROSS)
			v.x = rand()*rand_fac;
			v.y = rand()*rand_fac;
#else
			// ALL EDGES IN A LINE
			v.x = i;
			v.y = i;
#endif
		} while (v == prev);
		path.push_back(v);

#if PATH_TYPE == PATH_REALISTIC
		first_dir = normalize(first - v);
		prev_dir = normalize(v-prev);
#endif
		prev = v;
	}
	path.push_back(first);



	cout << "start measurement" << endl;

 	perf_time_t end;
	perf_time_t start;
	double secs;
#if PROGRESS_REPORT
	cout << "adding nodes and edges" << endl;
#endif
	perf_clock(start);

	for (unsigned i = 0; i<path.size(); i++) {
		vec2 v = path[i];
		g.addNode(v);
	}

	int MAX_EDGES = NODES;
#if PROGRESS_REPORT
	double last_secs = 0;
	int percent_1 = fmax(MAX_EDGES/100,1);
	fflush(stdout);
#endif

	Graph2DPlanar::Node* a = (*g.findNode(path[0]));
	Graph2DPlanar::Node* b;

	for (int i = 1; i<MAX_EDGES; i++) {
		b = (*g.findNode(path[i]));
		g.addLineSegment(a->pos,b->pos);
		a = b;
#if PROGRESS_REPORT
		if (i%percent_1 == 0) {
			perf_clock(end);
			secs = (end-start).secs_nsecs();
			printf("%02d%%: %03.3f s (+%03.3f)\r", i/percent_1, secs, secs-last_secs);
			fflush(stdout);
			last_secs = secs;
		}
#endif
	}

	g.cleanup();


	perf_clock(end);
	secs = (end-start).secs_nsecs();
	printf("total: %03.3f s\n", secs);


//	list_content(g);

#undef PROGRESS_REPORT
#undef PATH_TYPE
#undef PATH_REALISTIC
#undef PATH_CRISS_CROSS
#undef PATH_BEST_CASE
}


int test_Graph2DPlanar_addNode ()
{
	const int NODES = 128;

	Graph2DPlanar g(NODES, false);

	vec2 a(0,1);
	vec2 b(0,2);
	vec2 c(0,3);

	Graph2DPlanar::Node* n;
	n = g.addNode(a);
	assert(g.nodes_size() == 1);
	assert(n != NULL && n->pos == a);
	n = g.addNode(b);
	assert(g.nodes_size() == 2);
	assert(n != NULL && n->pos == b);
	n = g.addNode(c);
	assert(g.nodes_size() == 3);
	assert(n != NULL && n->pos == c);

	Graph2DPlanar::node_iterator it_n;
	it_n = g.findNode(a);
	assert(it_n != g.nodes_end() && (*it_n)->pos == a);
	it_n = g.findNode(b);
	assert(it_n != g.nodes_end() && (*it_n)->pos == b);
	it_n = g.findNode(c);
	assert(it_n != g.nodes_end() && (*it_n)->pos == c);
	return 0;
}


int test_Graph2DPlanar_addLineSegment_separate ()
{
	Graph2DPlanar g(false);

	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(0,2);
	vec2 d(0,3);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	//
	// test non-crossing segments (parallel)
	//

	g.addLineSegment(a,b);
	g.addLineSegment(c,d);

	assert(g.nodes_size() == 4);
	assert(g.edges_size() == 2);



	return 0;
}

int test_Graph2DPlanar_cleanup ()
{
	Graph2DPlanar g(false);

	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(0,2);
	vec2 d(0,3);
	vec2 e(2,2);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);


	g.addLineSegment(a,d); // creates a->b->c->d
	g.addLineSegment(d,e); // creates d->e
	g.addLineSegment(e,a); // closes loop a->b->c->d->e ->a
	// overlap a->d in reverse orientation causes removal between b->c
	g.addLineSegment(c,b);

	assert(g.nodes_size() == 5);
	assert(g.edges_size() == 4);

	g.cleanup();

	assert(g.nodes_size() == 0);
	assert(g.edges_size() == 0);

	return 0;
}

int test_Graph2DPlanar_addLineSegment_duplicate ()
{
	Graph2DPlanar g(false);

	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(0,2);
	vec2 d(0,3);

	g.addNode(a);
	g.addNode(b);

	//
	// test simple crossing
	//

	g.addLineSegment(a,b);
	g.addLineSegment(a,b);

	assert(g.nodes_size() == 2);
	assert(g.edges_size() == 1);



	return 0;
}

int test_Graph2DPlanar_addLineSegment_duplicate_reverse ()
{
	Graph2DPlanar g(false);

	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(0,2);
	vec2 d(0,3);

	g.addNode(a);
	g.addNode(b);

	g.addLineSegment(a,b);
	g.addLineSegment(b,a);

	assert(g.nodes_size() == 2);
	assert(g.edges_size() == 0);

	////////////////////////////////

	Graph2DPlanar f(true);
	f.addNode(a);
	f.addNode(b);

	f.addLineSegment(a,b);
	f.addLineSegment(b,a);

	assert(f.nodes_size() == 2);
	assert(f.edges_size() == 1);
	Graph2DPlanar::Edge* e = (*f.edges_begin());
	assert(e->weight == 0);

	return 0;
}

int test_Graph2DPlanar_addLineSegment_crossing_nodes ()
{
	Graph2DPlanar g(false);

	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(0,2);
	vec2 d(0,3);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	// edge crossing 1 node
	g.addLineSegment(a,c);
	assert(g.edges_size() == 2);
	assert(g.nodes_size() == 4);

	g.reset();

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	// edge crossing 2 nodes
	g.addLineSegment(a,d);

	assert(g.edges_size() == 3);
	assert(g.nodes_size() == 4);

	return 0;
}

int test_Graph2DPlanar_addLineSegment_overlapping_edges ()
{
	Graph2DPlanar g(false);

	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(0,2);
	vec2 d(0,3);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	// 1 edge in the middle
	g.addLineSegment(b,c);
	assert(g.edges_size() == 1);
	assert(g.nodes_size() == 4);

	// edge overlapping previous edge
	g.addLineSegment(a,d);

	assert(g.edges_size() == 3);
	assert(g.nodes_size() == 4);
	Graph2DPlanar::Edge* e;
	e = g.findEdge(a,b);
	assert(e != NULL && e->weight == 1);
	e = g.findEdge(b,c);
	assert(e != NULL && e->weight == 2);
	e = g.findEdge(c,d);
	assert(e != NULL && e->weight == 1);

	////////////////////////////////
	g.reset();
	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	// edge to outer nodes
	g.addLineSegment(a,d);
	assert(g.edges_size() == 3);
	assert(g.nodes_size() == 4);

	// edge overlapping in the middle
	g.addLineSegment(b,c);

	assert(g.edges_size() == 3);
	assert(g.nodes_size() == 4);

	e = g.findEdge(a,b);
	assert(e != NULL && e->weight == 1);
	e = g.findEdge(b,c);
	assert(e != NULL && e->weight == 2);
	e = g.findEdge(c,d);
	assert(e != NULL && e->weight == 1);

	return 0;
}

int test_Graph2DPlanar_addLineSegment_crossing_edges ()
{
	Graph2DPlanar g(false);

	vec2 a(-1,0);
	vec2 b(+1,0);
	vec2 c(0,-1);
	vec2 d(0,+1);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	//
	// test simple crossing
	//

	g.addLineSegment(a,b);
	g.addLineSegment(c,d);

	assert(g.edges_size() == 4);
	assert(g.nodes_size() == 5);

	return 0;
}



int test_Graph2DPlanar_all ()
{

#define ASSERT_MEMORY() 	\
	assert(    Graph2DPlanar::Node::total_nodes_in_system == 0 \
			&& Graph2DPlanar::Edge::total_edges_in_system == 0 \
			)


	test_Graph2DPlanar_addNode();
	ASSERT_MEMORY();

	test_Graph2DPlanar_addLineSegment_separate();
	ASSERT_MEMORY();

	test_Graph2DPlanar_addLineSegment_duplicate();
	ASSERT_MEMORY();

	test_Graph2DPlanar_addLineSegment_duplicate_reverse();
	ASSERT_MEMORY();

	test_Graph2DPlanar_addLineSegment_crossing_nodes();
	ASSERT_MEMORY();

	test_Graph2DPlanar_addLineSegment_overlapping_edges();
	ASSERT_MEMORY();

	test_Graph2DPlanar_addLineSegment_crossing_edges();
	ASSERT_MEMORY();

	test_Graph2DPlanar_cleanup();
	ASSERT_MEMORY();

#undef ASSERT_MEMORY
	return 0;
}




#endif /* MATH_TEST_GRAPH2DPLANAR_H_ */
