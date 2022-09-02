#define TAG 42
#define POS(x, y, side) ((x) + (y) * (side))

void plassman_p_processors(int pid, int p, MPI_Comm * comm_workers);
void plassman_v_processors(int pid);
int find_my_color(customSet * psetcolors, int nodes);
void print_adjacency_matrix(int pid, int * madjacency, int my_nodes, int nodes) ;
void error(char * msg);
void read_graph_to_adjacency_matrix(FILE * f, int ** pmatrix, int * pnodes);
