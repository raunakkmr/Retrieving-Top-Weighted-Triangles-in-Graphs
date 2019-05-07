#ifndef CLIQUE_ENUMERATOR_H
#define CLIQUE_ENUMERATOR_H
/*

Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

To compile:
"gcc kClist.c -O9 -o kClist".

To execute:
"./kClist k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print the number of k-cliques.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#include "graph.h"

namespace wsdm_2019_graph {

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg ;


typedef struct {

	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	unsigned *ns;//ns[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
	//unsigned *map;//oldID newID correspondance

	unsigned char *lab;//lab[i] label of node i
	unsigned **sub;//sub[l]: nodes in G_l

} specialsparse;


void freespecialsparse(specialsparse *g, unsigned char k){
	unsigned char i;
	free(g->ns);
	for (i=2;i<k+1;i++){
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

specialsparse* readedgelist(vector<full_edge>& edgelist){
	unsigned e1=NLINKS;
	specialsparse *g=(specialsparse*) malloc(sizeof(specialsparse));

	g->n=0;
	g->e=0;
	g->edges=(edge*) malloc(e1*sizeof(edge));
	for (auto& e : edgelist) {
		g->edges[g->e].s = e.src;
		g->edges[g->e].t = e.dst;
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		g->e++;
		if (g->e==e1) {
			e1+=NLINKS;
			g->edges=(edge*) realloc(g->edges,e1*sizeof(edge));
		}
	}
	g->n++;

	g->edges=(edge*) realloc(g->edges,g->e*sizeof(edge));

	return g;
}

void relabel(specialsparse *g){
	unsigned i, source, target, tmp;

	for (i=0;i<g->e;i++) {
		source=g->rank[g->edges[i].s];
		target=g->rank[g->edges[i].t];
		if (source<target){
			tmp=source;
			source=target;
			target=tmp;
		}
		g->edges[i].s=source;
		g->edges[i].t=target;
	}

}

///// CORE ordering /////////////////////

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max){
	unsigned i;
	bheap *heap=(bheap*) malloc(sizeof(bheap));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=(unsigned*) malloc(n_max*sizeof(unsigned));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=(keyvalue*)malloc(n_max*sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap,unsigned i, unsigned j) {
	keyvalue kv_tmp=heap->kv[i];
	unsigned pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

void bubble_up(bheap *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap,keyvalue kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_up(heap,heap->n-1);
}

void update(bheap *heap,unsigned key){
	unsigned i=heap->pt[key];
	if (i!=unsigned(-1)){
		((heap->kv[i]).value)--;
		bubble_up(heap,i);
	}
}

keyvalue popmin(bheap *heap){
	keyvalue min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n,unsigned *v){
	unsigned i;
	keyvalue kv;
	bheap* heap=construct(n);
	for (i=0;i<n;i++){
		kv.key=i;
		kv.value=v[i];
		insert(heap,kv);
	}
	return heap;
}

void freeheap(bheap *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(specialsparse* g){
	unsigned i,j,r=0,n=g->n;
	keyvalue kv;
	bheap *heap;

	unsigned *d0=(unsigned*)calloc(g->n,sizeof(unsigned));
	unsigned *cd0=(unsigned*)malloc((g->n+1)*sizeof(unsigned));
	unsigned *adj0=(unsigned*)malloc(2*g->e*sizeof(unsigned));
	for (i=0;i<g->e;i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}
	cd0[0]=0;
	for (i=1;i<g->n+1;i++) {
		cd0[i]=cd0[i-1]+d0[i-1];
		d0[i-1]=0;
	}
	for (i=0;i<g->e;i++) {
		adj0[ cd0[g->edges[i].s] + d0[ g->edges[i].s ]++ ]=g->edges[i].t;
		adj0[ cd0[g->edges[i].t] + d0[ g->edges[i].t ]++ ]=g->edges[i].s;
	}

	heap=mkheap(n,d0);

	g->rank=(unsigned*)malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++){
		kv=popmin(heap);
		g->rank[kv.key]=n-(++r);
		for (j=cd0[kv.key];j<cd0[kv.key+1];j++){
			update(heap,adj0[j]);
		}
	}
	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
//Building the special graph structure
void mkspecial(specialsparse *g, unsigned char k){
	unsigned i,ns,max;
	unsigned *d,*sub;
	unsigned char *lab;

	d=(unsigned*)calloc(g->n,sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		d[g->edges[i].s]++;
	}

	g->cd=(unsigned*)malloc((g->n+1)*sizeof(unsigned));
	ns=0;
	g->cd[0]=0;
	max=0;
	sub=(unsigned*)malloc(g->n*sizeof(unsigned));
	lab=(unsigned char*)malloc(g->n*sizeof(unsigned char));
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		max=(max>d[i-1])?max:d[i-1];
		sub[ns++]=i-1;
		d[i-1]=0;
		lab[i-1]=k;
	}
	//printf("max degree = %u\n",max);

	g->adj=(unsigned*)malloc(g->e*sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->adj[ g->cd[g->edges[i].s] + d[ g->edges[i].s ]++ ]=g->edges[i].t;
	}
	free(g->edges);

	g->ns=(unsigned*)malloc((k+1)*sizeof(unsigned));
	g->ns[k]=ns;

	g->d=(unsigned**)malloc((k+1)*sizeof(unsigned*));
	g->sub=(unsigned**)malloc((k+1)*sizeof(unsigned*));
	for (i=2;i<k;i++){
		g->d[i]=(unsigned*)malloc(g->n*sizeof(unsigned));
		g->sub[i]=(unsigned*)malloc(max*sizeof(unsigned));
	}
	g->d[k]=d;
	g->sub[k]=sub;

	g->lab=lab;
}


// TODO: Maybe this is premature optimization to lift these two parameters
// out of function call. Profile again with these within function
map<int, map<int, int>>* _adjmat;
vector<weighted_clique>* _clique_list;
void add_vertex(weighted_clique& C, int u) {
	C.vertices.push_back(u);
	for (int w : C.vertices) {
		C.weight += (*_adjmat)[w][u];
	}
}
void remove_last(weighted_clique& C) {
	int u = C.vertices.back();
	C.vertices.pop_back();
	for (int w : C.vertices) {
		C.weight -= (*_adjmat)[w][u];
	}
}
void kclique(unsigned l, specialsparse *g, unsigned long long *n, weighted_clique& C) {
	unsigned i,j,k,end,u,v,w;

	if(l==2){
		for(i=0; i<g->ns[2]; i++){//list all edges
			u=g->sub[2][i];
			add_vertex(C, u);
			//(*n)+=g->d[2][u];
			end=g->cd[u]+g->d[2][u];
			for (j=g->cd[u];j<end;j++) {
				v=g->adj[j];

				add_vertex(C, v);
				_clique_list->push_back(C);
				remove_last(C);
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
			remove_last(C);
		}

		return;
	}

	for(i=0; i<g->ns[l]; i++){
		u=g->sub[l][i];
		//printf("%u %u\n",i,u);
		g->ns[l-1]=0;
		end=g->cd[u]+g->d[l][u];
		for (j=g->cd[u];j<end;j++){//relabeling nodes and forming U'.
			v=g->adj[j];
			if (g->lab[v]==l){
				g->lab[v]=l-1;
				g->sub[l-1][g->ns[l-1]++]=v;
				g->d[l-1][v]=0;//new degrees
			}
		}
		for (j=0;j<g->ns[l-1];j++){//reodering adjacency list and computing new degrees
			v=g->sub[l-1][j];
			end=g->cd[v]+g->d[l][v];
			for (k=g->cd[v];k<end;k++){
				w=g->adj[k];
				if (g->lab[w]==l-1){
					g->d[l-1][v]++;
				}
				else{
					g->adj[k--]=g->adj[--end];
					g->adj[end]=w;
				}
			}
		}

		add_vertex(C, u);
		kclique(l-1, g, n, C);
		remove_last(C);

		for (j=0;j<g->ns[l-1];j++){//restoring labels
			v=g->sub[l-1][j];
			g->lab[v]=l;
		}

	}
}


set<weighted_clique> find_cliques(Graph& G, int k, bool print_diagnostic = false){
	specialsparse* g;
	unsigned long long n;

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	map<int, map<int, int>> adjmat;
	vector<full_edge> edgelist;
	for (int i = 0; i < (int) G.size(); i++) {
		for (auto e : G[i]) {
			if (i < e.dst) {
				edgelist.push_back({i, e.dst, e.wt});
				adjmat[i][e.dst] = e.wt;
				adjmat[e.dst][i] = e.wt;
			}
		}
	}

	g=readedgelist(edgelist);
	if (print_diagnostic) {
		printf("Number of nodes = %u\n",g->n);
		printf("Number of edges = %u\n",g->e);

		t2=time(NULL);
		printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
		t1=t2;

		printf("Building the graph structure\n");
	}
	ord_core(g);
	relabel(g);

	mkspecial(g,k);

	if (print_diagnostic) {
		printf("Number of nodes = %u\n",g->n);
		printf("Number of edges = %u\n",g->e);

		t2=time(NULL);
		printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
		t1=t2;

		printf("Iterate over all cliques\n");
	}
	vector<weighted_clique> retval;
	_clique_list = &retval;
	_adjmat = &adjmat;
	weighted_clique C;
	n=0;
	kclique(k, g, &n, C);
	_clique_list = nullptr;
	_adjmat = nullptr;

	if (print_diagnostic) {
		printf("Number of %u-cliques: %llu\n",k,n);

		t2=time(NULL);
		printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
		t1=t2;
	}
	freespecialsparse(g,k);

	if (print_diagnostic) {
		printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));
	}
	return set<weighted_clique>(retval.begin(), retval.end());
}

}

#endif