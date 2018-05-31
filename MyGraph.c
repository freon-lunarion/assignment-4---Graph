#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

// A vertex is a 2D point
typedef struct Vertex { 
    int x; // x-coordinate
    int y; // y-coordinate
} Vertex;
// each edge is a pair of vertices (end-points)
typedef struct Edge {
    Vertex *p1; // first end point
    Vertex *p2; // second end point
} Edge;

typedef struct EdgeNode {
    Vertex *v;
    float dist;
    struct EdgeNode *next;
    struct EdgeNode *prev;
} EdgeNode;

typedef struct VertexNode {
    Vertex *v;
    int visited;
    struct VertexNode *next;
    EdgeNode *edges;
    EdgeNode *lastEdge;
} VertexNode;

typedef struct GraphRep *Graph;
typedef struct GraphRep { // graph header
    VertexNode *vertices; // an array of vertices or a linked list of vertices  
    int nV; // #vertices
    int nE; // #edges
} GraphRep;

// QUEUE-STACK Hybrid
typedef struct QNode {
    Vertex *v;
    struct QNode *next;
    struct QNode *prev;
} QNode;

typedef struct QStack {
    QNode *first;
    QNode *last;
    int n;
} QStack;

QStack *newQStack() {
    QStack *QH;
    QH = malloc(sizeof(QStack));
    QH->first == NULL;
    QH->last == NULL;
    QH->n = 0;

    return QH;
}

void pushQStack(QStack *QH, Vertex *V) {
    QNode *N;
    N = malloc(sizeof(QNode));
    N->v = V;
    if (QH->n == 0) {
        QH->first = N;
        QH->last = N;
        QH->n = 1;
        return;
    }
    QH->last->next = N;
    N->prev = QH->last;
    QH->last = N;
    QH->n = QH->n +1;
}

Vertex *popQueue(QStack *QH) {
    QNode *N;
    N = malloc(sizeof(QNode));

    N = QH->first;
    QH->n = QH->n -1;
    QH->first = N->next;
    return N->v;
}

Vertex *popStack(QStack *QH) {
    QNode *N;
    N = malloc(sizeof(QNode));

    N = QH->last;
    QH->n = QH->n -1;
    QH->last = N->prev;
    return N->v;
}


// end of QUEUE
// PRIORITY QUEUE
typedef struct PQNode { 
  // each node stores the priority (key), name, execution time,
  //  release time and deadline of one task
  int key; //key of this item
  Vertex *v;
  struct PQNode *parent;
  struct PQNode *right;
  struct PQNode *left;
} PQNode;


typedef struct PQ{ //this is heap header
  int  size;      // count of items in the heap
  PQNode *lastNode; // last node pointer 
  PQNode *root; // pointer to the root of heap
} PQ;

PQ *newPQ() {
    PQ *pq;
    pq = malloc(sizeof(PQ));
    assert(pq != NULL);
    pq->size = 0;
    pq->lastNode = NULL;
    pq->root = NULL;
    return pq;
}

PQNode *newPQNode(int key, Vertex *vertex, PQNode *left, PQNode *right, PQNode *parent) {
    PQNode *new;
    new =  malloc(sizeof(PQNode));
    assert(new != NULL);
    new->key = key;
    new->v   = vertex;

    new->left   = left; // left child
    new->right  = right; // righ child
    new->parent = parent; // parent
    return new;
}

PQNode *mostLeftPQ(PQNode *N) {
  while (N->left != NULL) {
    N = N->left;
  }

  return N;
}


void moveUpPQ(PQ *h, PQNode *n){
    PQNode *p = n->parent;
    if (p->parent) {
        if (p == p->parent->right){
            p->parent->right = n;
        } else if (p == p->parent->left) {
            p->parent->left = n;
        }
    } else {
        h->root = n;
    }
    n->parent = p->parent;
    int nside;
    PQNode *c;
    if (n == p->right){
        nside = 1;
        c = p->left;
    } else{
        nside = 0;
        c = p->right;
    }

    p->left = n->left;
    if (p->left) {
        p->left->parent = p;
    }

    p->right = n->right;
    if (p->right) {
        p->right->parent = p;
    }

    if (nside == 0){
        n->left = p;
    } else {
        n->right = p;
    }

    p->parent = n;

    if (nside == 1){
        n->left = c;
    } else {
        n->right = c;
    }

    if (c) {
        c->parent = n;
    }

    if (n == h->lastNode){
        h->lastNode = p;
    }
}

void replacePQ(PQ *h, PQNode *d, PQNode *s) {
    if (d->parent) {
        if (d == d->parent->right) {
            d->parent->right = s;
        } else {
            d->parent->left = s;
        }   

    } else {
        h->root = s;
    }

    s->parent = d->parent;

    s->left = d->left;
    if (s->left) {
        s->left->parent = s;
    } 

    s->right = d->right;
    if (s->right) {
        s->right->parent = s;
    }
}

void insertPQ(PQ *T, int k, Vertex *n) {
    PQNode *newNode,*pNode,*temp;
    if (T->size == 0) {
        newNode = newPQNode(k,n,NULL,NULL,NULL);
        T->root = newNode;
        T->size += 1;
        T->lastNode = newNode;      
        return;
    }

    // int cap = capPQ(lvPQ(T->size));
    if (T->size % 2 == 0){ // time complexity = 1+1+1+3+1; O(1)
        pNode = T->lastNode->parent;
        newNode = newPQNode(k,n,NULL,NULL,pNode);
        pNode->right = newNode;
        T->size += 1;
        T->lastNode = newNode;
    } else {
        pNode = T->lastNode;
        while (T->root != temp){
            if (pNode->parent->right == pNode) {
                pNode = pNode->parent;
            } else {
                break;
            }
        }

        pNode = pNode->parent->right;
        pNode = mostLeftPQ(pNode);
        newNode = newPQNode(k,n,NULL,NULL,pNode);
        pNode->left = newNode;
        T->size +=1;
        T->lastNode = newNode;
        if (newNode->key < pNode->key){
            T->lastNode = pNode;
        }  
    }

    while (newNode->parent != NULL && newNode->key < newNode->parent->key) {
        moveUpPQ(T,newNode);
    } 
}

PQNode* removePQ(PQ *T) {
    PQNode *firstNode, *cur, *result;
    firstNode = T->root;
    result = newPQNode(firstNode->key, firstNode->v, NULL, NULL, NULL);

    if (!firstNode->parent && !firstNode->left && !firstNode->right){
        T->root = NULL;
        T->lastNode = NULL;
        T->size = 0;
        return firstNode;
    }

    cur = T->lastNode;
    while (cur->parent && cur == cur->parent->left) {
        cur = cur->parent;
    } 

    if (cur->parent){
        cur = cur->parent->left;
        while (cur->right) {
            cur = cur->right;
        }
    } else {
        while (cur->right) {
            cur = cur->right;
        }
    }

    // disconect lastNode
    if (T->lastNode == T->lastNode->parent->right) {
        T->lastNode->parent->right = NULL;
    } else {
        T->lastNode->parent->left = NULL;
    }

    if(firstNode == T->lastNode) {
        T->lastNode = cur;
    } else {
        PQNode *lastNode;
        lastNode = T->lastNode;
        replacePQ(T,firstNode, lastNode);
        if (firstNode != cur) {
            T->lastNode = cur;
        }

        if (lastNode->parent && lastNode->key < lastNode->parent->key){
            do {
                moveUpPQ(T,lastNode);
            } while (lastNode->parent && lastNode->key < lastNode->parent->key);
        } else {
            while (lastNode->left || lastNode->right) {
                PQNode *lNode, *rNode;
                lNode = lastNode->left;
                rNode = lastNode->right;
                if (lNode && rNode ){
                    if (lastNode->key < lNode->key && lastNode < rNode->key){
                        break;
                    }      

                    if (lNode->key < rNode->key) {
                        moveUpPQ(T,lNode);            
                    } else {
                        moveUpPQ(T,rNode);
                    }

                } else if(lNode) {
                    PQNode *lNode;
                    lNode = lastNode->left;
                    if (lastNode->key < lNode->key){
                        break;
                    }
                    moveUpPQ(T,lNode);

                }
            }
        }
    }
    T->size -= 1;
    return result;
}
// end of PRIORITY QUEUE
   
float calc_dist_v(Vertex *V1, Vertex *V2) {
    int x,y,x_sq,y_sq;

    x = V1->x - V2->x;
    y = V1->y - V2->y;

    x_sq = pow(x, 2);
    y_sq = pow(y, 2);

    return sqrt(x_sq + y_sq);
}

void reset_vertices(Graph g){
    VertexNode *cV;
    cV = malloc(sizeof(VertexNode));
    // Reset visiting 
    cV = g->vertices;
    while (cV ) {
        cV->visited = 0;
        cV = cV->next;
    }
}

Vertex *CreateVertex(int x, int y) {
    Vertex *v;
    v = malloc(sizeof(Vertex));
    v->x = x;
    v->y = y;

    return v;
}

VertexNode *CreateVertexNode(Vertex *v) {
    VertexNode *VN;
    VN = malloc(sizeof(VertexNode));

    assert (VN != NULL);
    VN->v = v;
    VN->visited = 0;
    VN->edges = NULL;
    VN->lastEdge=NULL;
    VN->next = NULL;

    return VN;
}

EdgeNode *CreateEdgeNode(Vertex *origin, Vertex *v) {
    EdgeNode *EN;
    EN = malloc(sizeof(EdgeNode));
    assert(EN != NULL);
    EN->v = v;
    EN->next = NULL;
    EN->prev = NULL;
    EN->dist = calc_dist_v(origin,v);
    return EN;
}

// Add the time complexity analysis of CreateEmptyGraph() here
Graph CreateEmptyGraph()
{
    Graph G;
    G = malloc(sizeof(Graph));
    G->vertices = NULL;
    G->nV = 0;
    G->nE = 0;
    return G;
}
// Add the time complexity analysis of InsertEdge() here
int InsertEdge(Graph g, Edge *e)
{
    Vertex *V1,*V2;
    VertexNode *vn1, *vn2;
    EdgeNode *ed1, *ed2;
    
    V1 = CreateVertex(e->p1->x,e->p1->y);
    V2 = CreateVertex(e->p2->x,e->p2->y);
    
    if (g->nV == 0) { 
        // if the size of graph is 0
        vn1 = CreateVertexNode(V1);
        ed1 = CreateEdgeNode(V1, V2);
        vn1->edges = ed1;
        vn1->lastEdge = ed1;

        vn2 = CreateVertexNode(V2);
        ed2 = CreateEdgeNode(V2, V1);
        vn2->edges = ed2;
        vn2->lastEdge = ed2;
        
        vn1->next = vn2;
        g->vertices = vn1;

        g->nV = 2;
        g->nE = 1;
        return 1;
    } 
    VertexNode *vc;
    vc = malloc(sizeof(VertexNode));
    EdgeNode *ec;
    ec = malloc(sizeof(EdgeNode));
    vn1 = NULL;
    vn2 = NULL;

    // check already on graph or not
    vc = g->vertices;

    while (1) {
        if (vc->v->x == V1->x && vc->v->y == V1->y) {
            ec = vc->edges;
            // printf("E(%d,%d)\n", ec->v->x,ec->v->y);

            while (ec->next != NULL) {
                if (ec->v->x == V2->x && ec->v->y == V2->y) {
                	// we can assume p2 already on the graph because P1 have edge with P2
                    return 0;
                }
                ec = ec->next;
            }
            // p1 on the graph, but no edge to p2 ( either p2 on graph or not)
            vn1 = vc;
            if (vn2 != NULL){
                break;
            }
        }

        // p2 on the graph
        if (vc->v->x == V2->x && vc->v->y == V2->y){
            vn2 = vc;
            if (vn1 != NULL){
                break;
            }
        }
        if (vc->next){
            vc = vc->next;
            /* code */
        } else {
            break;
        }
    }

    if (vn1 == NULL && vn2 == NULL){
    // both of the vertices not on graph
        vn1 = CreateVertexNode(V1);
        ed1 = CreateEdgeNode(V1, V2);
        vn1->edges = ed1;
        vn1->lastEdge = ed1;
        
        vn2 = CreateVertexNode(V2);
        ed2 = CreateEdgeNode(V2, V1);
        vn2->edges = ed2;
        vn2->lastEdge = ed2;

        g->nV += 2;
        vn1->next = vn2;
        vc->next = vn1;
        
    } else if (vn1 != NULL && vn2 == NULL) {
    // p1 on graph, but p2
        ed1 = CreateEdgeNode(vn1->v, V2);
        
        vn1->lastEdge->next = ed1;
        ed1->prev = vn1->lastEdge;
        vn1->lastEdge = ed1;
        
        vn2 = CreateVertexNode(V2);
        ed2 = CreateEdgeNode(V2, vn1->v);
        vn2->edges = ed2;
        vn2->lastEdge = ed2;

        g->nV += 1;
        vc->next = vn2;

    } else if (vn1 == NULL && vn2 != NULL) {
    // p2 on graph, but p1
        vn1 = CreateVertexNode(V1);
        ed1 = CreateEdgeNode(V1, vn2->v);
        vn1->edges = ed1;
        vn1->lastEdge = ed1;

        ed2 = CreateEdgeNode(vn2->v, V1);
		ed2->prev = vn2->lastEdge;
		vn2->lastEdge->next = ed2;
        vn2->lastEdge = ed2;


        g->nV += 1;
        vc->next = vn1;

    } else {
    // both p1 and p2 on the graph, but not have connection
        ed1 = CreateEdgeNode(vn1->v,vn2->v);
  	  ed1->prev = vn1->lastEdge;      
		vn1->lastEdge->next = ed1;
		vn1->lastEdge = ed1;

        ed2 = CreateEdgeNode(vn2->v,vn1->v);
        ed2->prev = vn2->lastEdge;
        vn2->lastEdge->next = ed2;
        vn2->lastEdge = ed2;
    }
    
    g->nE += 1;

    return 1;
}

// Add the time complexity analysis of DeleteEdge() here
void DeleteEdge(Graph g, Edge *e)
{
    
}

QStack *getReachableVertices(Graph g, Vertex *v) {
    
    QStack *queue, *result;
    result = newQStack();
    queue = newQStack();

    VertexNode *sV, *cV;
    cV = malloc(sizeof(VertexNode));
    sV = malloc(sizeof(VertexNode));
    sV = g->vertices;
    while (sV) {
        if (sV->v->x == v->x && sV->v->y == v->y){
            break;
        }
        sV = sV->next;
    }

    if (sV == NULL){ // didn't found the vertex on the graph
        return result;
    }
    
    Vertex *tV;
    tV = malloc(sizeof(Vertex));

    EdgeNode *EN;
    EN = malloc(sizeof(EdgeNode));

    pushQStack(queue,sV->v);
    while (queue->n != 0) {
        tV = popQueue(queue);
        cV = g->vertices;
        while (cV != NULL) {
            if (cV->v == tV && cV->visited == 0){
                pushQStack(result,tV);
                cV->visited = 1;
                EN = cV->edges;
                while (EN != NULL) {
                    pushQStack(queue,EN->v);
                    EN = EN->next;
                }
            }
            cV = cV->next;
        }
    }

    return result;

}

// Add the time complexity analysis of ReachableVertices() here
void ReachableVertices(Graph g, Vertex *v)
{
    Vertex *tV;
    tV = malloc(sizeof(Vertex));

    QStack *result;
    result = getReachableVertices(g,v);
    
    while (result->n != 0) {
        tV = popQueue(result);
        if (tV->x == v->x && tV->y == v->y)
        {
            continue;
        }
        printf("(%d,%d)",tV->x, tV->y);

        if (result->n != 0){
            printf(", ");
        } else {
            printf("\n");  
        }
    }

    reset_vertices(g);

}

void *dijkstra(VertexNode *origin, VertexNode *target) {

}
// Add the time complexity analysis of ShortestPath() here
void ShortestPath(Graph g, Vertex *u, Vertex *v)
{
    VertexNode *cur, *origin, *target;
    QStack *vertices;
    vertices = getReachableVertices(g, u);


    
}

// Add the time complexity analysis of FreeGraph() here
void FreeGraph(Graph g)
{

}

// Add the time complexity analysis of ShowGraph() here
// O(|V|.|E|); V = vertices, E = edges
void ShowGraph(Graph g)
{
    QStack *queue, *result;
    queue = newQStack();
    Vertex *v;
    v = malloc(sizeof(Vertex));

    VertexNode *VN, *cV;
    cV = malloc(sizeof(VertexNode));
    // VN = g->vertices;
    EdgeNode *EN;
    EN = malloc(sizeof(EdgeNode)); 
    pushQStack(queue,g->vertices->v);
    while (queue->n != 0) {
        v = popQueue(queue);
        cV = g->vertices;
        while (cV != NULL) { // O(|V|.|E|)

            if (cV->v == v && cV->visited == 0){
                cV->visited = 1;
                EN = cV->edges;
                
                while (EN != NULL) { // O(|E|)
                    printf("(%d,%d),(%d,%d), ", v->x,v->y,EN->v->x,EN->v->y);
                    pushQStack(queue,EN->v);
                    EN = EN->next;
                }
            }
            cV = cV->next;
        }

        cV = g->vertices;
        if (queue->n == 0 ){
            while(cV) { // O(|V|)
                if (cV->visited == 0) {
                    pushQStack(queue,cV->v);
                    break;
                }
                cV = cV->next;
            }
        }

    }
    
    reset_vertices(g);

    printf("\n");

}

int main() //sample main for testing 
{ 
    Graph g1;
    Edge *e_ptr; 
    Vertex *v1, *v2;

    // Create an empty graph g1;
    g1=CreateEmptyGraph();

    // Create first connected component 
    // Insert edge (0,0)-(0,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=0;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (0,0)-(5,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=5;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // // Insert edge (0, 10)-(10, 10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=10;
    v2->x=10;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (0,10)-(5,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=10;
    v2->x=5;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (0,0)-(5,4)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=5;
    v2->y=4;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (5, 4)-(10, 4)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=5;
    v1->y=4;
    v2->x=10;
    v2->y=4;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (5,6)-(10,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=5;
    v1->y=6;
    v2->x=10;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (10,10)-(10,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=10;
    v1->y=10;
    v2->x=10;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (10, 6)-(10, 4)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=10;
    v1->y=6;
    v2->x=10;
    v2->y=4;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Create second connected component
    // Insert edge (20,4)-(20,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=20;
    v1->y=4;
    v2->x=20;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (20,10)-(30,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=20;
    v1->y=10;
    v2->x=30;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    free(e_ptr);
    free(v1);
    free(v2);

    // Insert edge (25,5)-(30,10)   
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=25;
    v1->y=5;
    v2->x=30;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n"); 
    free(e_ptr);
    free(v1);
    free(v2);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    v2->x=0;
    v2->y=10;
    ReachableVertices(g1, v2);
    //Display graph g1
    // ShowGraph(g1);

    // Find the shortest path between (0,0) and (10,6) 
    // v1=(Vertex*) malloc(sizeof(Vertex));
    // assert(v1 != NULL);
    // v2=(Vertex *) malloc(sizeof(Vertex));
    // assert(v2 != NULL);
    // v1->x=0;
    // v1->y=0;
    // v2->x=10;
    // v2->y=6;
    // ShortestPath(g1, v1, v2);
    // free(v1);
    // free(v2);
      
    // // Delete edge (0,0)-(5, 6)
    // e_ptr = (Edge*) malloc(sizeof(Edge));
    // assert(e_ptr != NULL);
    // v1=(Vertex*) malloc(sizeof(Vertex));
    // assert(v1 != NULL);
    // v2=(Vertex *) malloc(sizeof(Vertex));
    // assert(v2 != NULL);
    // v1->x=0;
    // v1->y=0;
    // v2->x=5;
    // v2->y=6;
    // e_ptr->p1=v1;
    // e_ptr->p2=v2;    
    // DeleteEdge(g1, e_ptr);
    // free(e_ptr);
    // free(v1);
    // free(v2);
         
    // // Display graph g1
    // ShowGraph(g1);

    // // Find the shortest path between (0,0) and (10,6) 
    // v1=(Vertex*) malloc(sizeof(Vertex));
    // assert(v1 != NULL);
    // v2=(Vertex *) malloc(sizeof(Vertex));
    // assert(v2 != NULL);
    // v1->x=0;
    // v1->y=0;
    // v2->x=10;
    // v2->y=6; 
    // ShortestPath(g1, v1, v2);
    // free(v1);
    // free(v2);

    // // Find the shortest path between (0,0) and (25,5)
    // v1=(Vertex*) malloc(sizeof(Vertex));
    // assert(v1 != NULL);
    // v2=(Vertex *) malloc(sizeof(Vertex));
    // assert(v2 != NULL);
    // v1->x=0;
    // v1->y=0;
    // v2->x=25;
    // v2->y=5;
    // ShortestPath(g1, v1, v2);
    // free(v1);
    // free(v2);   

    // // Find reachable vertices of (0,0)
    // v1=(Vertex*) malloc(sizeof(Vertex));
    // assert(v1 != NULL);
    // v1->x=0;
    // v1->y=0;
    // ReachableVertices(g1, v1);
    // free(v1);

    // // Find reachable vertices of (20,4)
    // v1=(Vertex*) malloc(sizeof(Vertex));
    // assert(v1 != NULL);
    // v1->x=20;
    // v1->y=4;
    // ReachableVertices(g1, v1);
    // free(v1);

    // // Free graph g1
    // FreeGraph(g1);

    // DEBUG
    // VertexNode *VN;
    // printf("V%d, E%d\n", g1->nV, g1->nE);
    // VN = malloc(sizeof(VertexNode));
    // VN = g1->vertices;
    // while (VN != NULL) {
    //     EdgeNode *EC;
    //     printf("(%d,%d)->", VN->v->x,VN->v->y);
    //     EC = malloc(sizeof(EdgeNode));
    //     EC = VN->edges;
    //     while(EC != NULL) {
    //         printf("(%d,%d = %f) ", EC->v->x,EC->v->y, EC->dist );
    //         EC = EC->next;
    //     }
    //     printf("\n");
    //     VN = VN->next;
    // }


    return 0; 
}
