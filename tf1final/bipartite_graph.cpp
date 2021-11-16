/*
 BIPARTITE GRAPH PARTIAL COLORING
 COMPILE:
 	g++ -g -Wall -fopenmp -o test bipartite_graph.cpp

	g++ -std=c++14 -fopenmp-simd  -o test bipartite_graph.cpp
 RUN:
 	./test 4 < coPapersDBLP.mtx
*/

#include <iostream>
#include <fstream>
#include <omp.h>
#include "bipartite_graph.h"

using namespace std;
int thread_count;
unsigned int v, e;

//		VERTEX BASED

template<class N, class E>
void CBipartiteGraph<N, E>::color_queue_vertex()
{
	int tam = mleft.size();
	vector<int> forbiddenColors;
	// cambiar el tamaño de los colores prohibidos
	forbiddenColors.resize ( tam, -1);
	unsigned int id_w, id_v;

	#pragma omp parallel for num_threads(thread_count) firstprivate(forbiddenColors)
	for (unsigned int w=0; w<W.size(); w++) {
		find_mleft(W[w]->value, id_w);
		// para cada vértice izquierdo en adj: nets(w)
		for (typename list<Node*>::iterator v=mleft[id_w]->neighbours.begin(); v!=mleft[id_w]->neighbours.end(); v++) {
			// para cada vértice derecho en adj: vtxs(v)
			find_mright((*v)->value, id_v);
			for (typename list<Node*>::iterator u=mright[id_v]->neighbours.begin(); u!=mright[id_v]->neighbours.end(); u++) {
				// marcar el color del nodo (*u) como prohibido al vértice
				if((*u)->color != -1) {
					forbiddenColors[(*u)->color] = W[w]->value;
				}
			}
		}
		int candidate_col = 0;
		// política de coloración first-fit
		while (forbiddenColors[candidate_col] == W[w]->color)
			candidate_col += 1;
		W[w]->color = candidate_col;
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::remove_conflicts_vertex()
{
	unsigned int id_w, id_v;
	#pragma omp parallel for shared(mvertices_need_new_color)
	for (unsigned int w=0; w<W.size(); w++) {
		find_mleft(W[w]->value, id_w);
		for (typename list<Node*>::iterator v=mleft[id_w]->neighbours.begin(); v!=mleft[id_w]->neighbours.end(); v++) {
			find_mright((*v)->value, id_v);
			for (typename list<Node*>::iterator u=mright[id_v]->neighbours.begin(); u!=mright[id_v]->neighbours.end(); u++) {
				if ((*u)->color == W[w]->color && W[w]->value > (*u)->value)
					#pragma omp critical
					mvertices_need_new_color.push_back(W[w]);
					break;
			}
		}
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::greedy_graph_coloring_vertex()
{
	while(!W.empty()) {
		std::cout<<"	coloring\n";
		color_queue_vertex();
		std::cout<<"	remove conflicts\n";
		remove_conflicts_vertex();
		// vértices limpios para colorear
		W.clear();
		for (unsigned int i=0; i<mvertices_need_new_color.size(); i++)
			W.push_back(mvertices_need_new_color[i]);
	}
}

// 		NET BASED

template<class N, class E>
void CBipartiteGraph<N, E>::color_queue_net()
{
	int tam = mleft.size();
	vector<int> forbiddenColors;
	// cambiar el tamaño de los colores prohibidos
	forbiddenColors.resize ( tam, -1);
	unsigned int id_v;
	int candidate_col = 0;

	#pragma omp parallel for num_threads(thread_count) firstprivate(forbiddenColors, candidate_col)
	for (unsigned int v=0; v<mright.size(); v++) {
		candidate_col = 0;
		find_mleft(mright[v]->value, id_v);
		// vtxs(v)
		for (typename list<Node*>::iterator u=mright[id_v]->neighbours.begin(); u!=mright[id_v]->neighbours.end(); u++) {
			// marcar el color del nodo (*u) como prohibido al vértice
			if((*u)->color == -1 || (*u)->color==forbiddenColors[candidate_col]) {
				while (forbiddenColors[candidate_col] == mright[v]->color)
					candidate_col += 1;
				mright[v]->color = candidate_col;
				forbiddenColors[(*u)->color] = mright[v]->value;
			}
		}
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::remove_conflicts_net()
{
	int tam = mleft.size();
	vector<int> forbiddenColors;
	// cambiar el tamaño de los colores prohibidos
	forbiddenColors.resize ( tam, -1);
	unsigned int id_v;

	#pragma omp parallel for num_threads(thread_count) firstprivate(forbiddenColors)
	for (unsigned int v=0; v<mright.size(); v++) {
		find_mleft(mright[v]->value, id_v);
		// vtxs(v)
		for (typename list<Node*>::iterator u=mright[id_v]->neighbours.begin(); u!=mright[id_v]->neighbours.end(); u++) {
			// marcar el color del nodo (*u) como prohibido al vértice
			if((*u)->color != -1 ) {
				if (forbiddenColors[(*u)->color])
					(*u)->color = -1;
				else
					forbiddenColors[(*u)->color] = mright[v]->value;
			}
		}
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::greedy_graph_coloring_net()
{
	std::cout<<"\t  coloring\n";
	color_queue_net();
	std::cout<<"\t  remove conflicts\n";
	remove_conflicts_net();
}

//		ALGORITHMS
template<class N, class E>
void CBipartiteGraph<N, E>::v_v()
{
	while(!W.empty()) {
		std::cout<<"\t  coloring\n";
		color_queue_vertex();
		std::cout<<"\t  remove conflicts\n";
		remove_conflicts_vertex();
		// clear vertices to be colored
		W.clear();
		for (unsigned int i=0; i<mvertices_need_new_color.size(); i++)
			W.push_back(mvertices_need_new_color[i]);
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::v_n1()
{
	int it = 0;
	while(!W.empty()) {
		std::cout<<"\t  coloring\n";
		color_queue_vertex();
		std::cout<<"\t  remove conflicts\n";
		if (it < 1)
			remove_conflicts_net();
		// clear vertices to be colored
		W.clear();
		for (unsigned int i=0; i<mvertices_need_new_color.size(); i++)
			W.push_back(mvertices_need_new_color[i]);
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::v_n2()
{
	int it = 0;
	while(!W.empty()) {
		std::cout<<"\t  coloring\n";
		color_queue_vertex();
		std::cout<<"\t  remove conflicts\n";
		if (it < 2)
			remove_conflicts_net();
		// clear vertices to be colored
		W.clear();
		for (unsigned int i=0; i<mvertices_need_new_color.size(); i++)
			W.push_back(mvertices_need_new_color[i]);
	}
}

template<class N, class E>
void CBipartiteGraph<N, E>::n1_n2()
{
	int it1 = 0;
	int it2 = 0;
	std::cout<<"\t  coloring\n";
	if (it1 < 1)
		color_queue_net();
	std::cout<<"\t  remove conflicts\n";
	if (it2 < 2)
		remove_conflicts_net();
}

template<class N, class E>
void CBipartiteGraph<N, E>::n2_n2()
{
	int it1 = 0;
	int it2 = 0;
	std::cout<<"\t  coloring\n";
	if (it1 < 2)
		color_queue_net();
	std::cout<<"\t  remove conflicts\n";
	if (it2 < 2)
		remove_conflicts_net();
}

int main(int argc, char* argv[])
{
	int a, b;
	double start, finish;
    std::cout << "\t    BIPARTITE GRAPH PARTIAL COLORING" << std::endl;
    std::cout << "\t    --------------------------------" << std::endl << std::endl;
	CBipartiteGraph<int, int> bgraph;

	//thread_count = strtol(argv[1], NULL, 10);
	thread_count = 4;

    std::cout<<"Numero de Vertices: "; std::cin >> v;
    std::cout<<"Numero de Aristas : "; std::cin >> e;

	for(unsigned int i=1; i<=v; i++) {
		bgraph.InsNode(i);
	}
	std::cout<<"	Cargar nodos\n";
	// Numero de nodos y aristas: 5
	// Aristas: 2--1, 1--4, 2--5, 5--4,2--3

	unsigned int id_a;
	unsigned int id_b;
	while(e>0) {
		cin>> a >> b;
		if (bgraph.find_mleft(a, id_a) && bgraph.find_mright(b, id_b)) {
			bgraph.mleft[id_a]->InsEdge(bgraph.mright[id_b]);
		}
		else if (bgraph.find_mright(a, id_a) && bgraph.find_mleft(b, id_b)) {
			bgraph.mright[id_a]->InsEdge(bgraph.mleft[id_b]);
		}
		e--;
	}
	std::cout<<"	Cargar aristas\n";

	bgraph.printGraph();
	bgraph.first_queue();

	ofstream result;
	result.open ("result_bipartite.txt", std::ios_base::app);
	result << "\n Numero de threads: " << thread_count << " - Dataset: coPapersDBLP \n"; // test

	 //cout<<"> Vertex Based: \n";
	 //start = omp_get_wtime();
     //bgraph.greedy_graph_coloring_vertex();
     //finish = omp_get_wtime();
     //std::cout<<"Elapsed time: "<<finish-start<<" seconds\n";
	 //result << "> Vertex Based: "<< finish-start << "\n";

	// cout<<"> Net Based: \n";
	// start = omp_get_wtime();
    // bgraph.greedy_graph_coloring_net();
    // finish = omp_get_wtime();
    // std::cout<<"Elapsed time: "<<finish-start<<" seconds\n";
	// result << "> Net Based: "<< finish-start << "\n";

	// ---------------------------------     ALGORITHMS
	std::cout << "\n\t    --------- ALGORITMOS ----------\n";
	cout<<"\n\t > V-V \n";
	start = omp_get_wtime();
    bgraph.v_v();
    finish = omp_get_wtime();
    std::cout<<"\t  Elapsed time: "<<finish-start<<" seconds\n";
	result << "> V-V: "<< finish-start << "\n";

	bgraph.printGraph();

	cout<<"\n\t > V-N1 \n";
	start = omp_get_wtime();
    bgraph.v_n1();
    finish = omp_get_wtime();
    std::cout<<"\t  Elapsed time: "<<finish-start<<" seconds\n";
	result << "> V-N1: "<< finish-start << "\n";

	bgraph.printGraph();

	cout<<"\n\t > V-N2 \n";
	start = omp_get_wtime();
    bgraph.v_n2();
    finish = omp_get_wtime();
    std::cout<<"\t  Elapsed time: "<<finish-start<<" seconds\n";
	result << "> V-N2: "<< finish-start << "\n";

	bgraph.printGraph();

	cout<<"\n\t > N1-N2 \n";
	start = omp_get_wtime();
    bgraph.n1_n2();
    finish = omp_get_wtime();
    std::cout<<"\t  Elapsed time: "<<finish-start<<" seconds\n";
	result << "> N1-N2: "<< finish-start << "\n";

	bgraph.printGraph();

	cout<<"\n\t > N2-N2 \n";
	start = omp_get_wtime();
    bgraph.n2_n2();
    finish = omp_get_wtime();
    std::cout<<"\t  Elapsed time: "<<finish-start<<" seconds\n";
	result << "> N2-N2: "<< finish-start << "\n";

	bgraph.printGraph();

    std::cout << "\n\t    ------ COLORES DEL GRAFO-------\n";
    for (int i=0; i<bgraph.mleft.size(); i++)
        std::cout<<bgraph.mleft[i]->color<<" ";
    std::cout<<"\n";
    for (int i=0; i<bgraph.mright.size(); i++)
        std::cout<<bgraph.mright[i]->color<<" ";

	result.close();
	return 0;
}
