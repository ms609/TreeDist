void report_calloc_error() {
  Rprintf("Error allocating memory in calloc");  
}

static R_NativePrimitiveArgType order_edges_number_nodes_t[] = {
  INTSXP, INTSXP, INTSXP
};
void order_edges_number_nodes(int *parent, int *child, const int *n_edge)
{
  int i, q_pos = 0, o_node, next_node;
  const int n_node = *n_edge / 2;
  const int n_allnodes = *n_edge + 1L, root_node = n_node + 2L;
  int * start_p = calloc(*n_edge, sizeof(int));
  int * start_c = calloc(*n_edge, sizeof(int));
  int * child_l  = calloc(n_node, sizeof(int));
  int * child_r  = calloc(n_node, sizeof(int));
  int * queue_p  = calloc(n_node, sizeof(int));
  int * queue_c  = calloc(n_node, sizeof(int));
  // TODO check that calloc has returned a non-null pointer; clean-up and exit if calloc has failed
  for (i = 0; i < *n_edge; i++) {
    // Initialize
    start_p[i] = parent[i];
    start_c[i] = child[i];
    q_pos = parent[i] - root_node;
    if (child_l[q_pos]) {
      if (child_l[q_pos] < child[i] && child[i] > root_node) {
        child_r[q_pos] = child[i];
      } else {
        child_r[q_pos] = child_l[q_pos];
        child_l[q_pos] = child[i];
      }
    } else {
      child_l[q_pos] = child[i];
    }
  }
  o_node = root_node;
  q_pos = 0;
  for (i = 0; i < *n_edge; i++) {
    if (o_node < root_node) { // We've just reached a tip
      parent[i] = queue_p[--q_pos];
      child[i] = queue_c[q_pos];
      o_node = child[i];
    } else { // We're at an internal node
      parent[i] = o_node;
      child[i]  = child_l[o_node - root_node];
      queue_p[q_pos] = o_node;
      queue_c[q_pos++] = child_r[o_node - root_node];
      o_node = child_l[o_node - root_node];
    }
  }
  free(start_p);
  free(start_c);
  free(child_l);
  free(child_r);
  free(queue_p);
  free(queue_c);
  
  // Now number nodes:
  int * renumber = calloc(n_allnodes, sizeof(int));
  if (renumber != NULL) {
    next_node = root_node;
    for (i = 0; i < n_allnodes; i++) renumber[i] = i + 1;
    for (i = 0; i < *n_edge; i++) {
      if (child[i] > root_node) renumber[child[i]-1] = ++(next_node);
    }
    for (i = 0; i < *n_edge; i++) {
      parent[i] = renumber[parent[i]-1L];
      child[i] = renumber[child[i]-1L];
    }
    free(renumber);
  } else {
    report_calloc_error();
  }
}

SEXP RENUMBER(SEXP par, SEXP chi, SEXP ned){   
    int *data, *nr=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   
    double *pvtmp;  
    SEXP DAT, pars, pvec, pscore, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, 4L));
    PROTECT(pars = allocVector(INTSXP, *nr));
    PROTECT(pscore = allocVector(REALSXP, 1L));
    PROTECT(DAT = allocMatrix(INTSXP, nr[0], m));
    PROTECT(pvec = allocVector(REALSXP, m));
    pvtmp = REAL(pvec);
    data = INTEGER(DAT);
    for(i=0; i<m; i++) pvtmp[i] = 0.0;
    for(i=0; i<*nr; i++) INTEGER(pars)[i] = 0L;
    REAL(pscore)[0]=0.0;
    for(i=0; i<(*nr * n); i++)data[i] = INTEGER(dat)[i];
    
    fitch8(data, nr, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), REAL(weight), pvtmp, REAL(pscore));
    
    SET_VECTOR_ELT(RESULT, 0, pscore);
    SET_VECTOR_ELT(RESULT, 1, pars);
    SET_VECTOR_ELT(RESULT, 2, DAT);
    SET_VECTOR_ELT(RESULT, 3, pvec);
    UNPROTECT(5);
    return(RESULT); 
}
