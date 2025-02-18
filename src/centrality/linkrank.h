#ifndef IGRAPH_LINKRANK_H
#define IGRAPH_LINKRANK_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

__BEGIN_DECLS

igraph_error_t igraph_linkrank(const igraph_t *graph,
                               const igraph_vector_t *pagerank,
                               const igraph_vector_t *weights,
                               igraph_matrix_t *linkrank);

igraph_error_t igraph_linkrank_basic(const igraph_t *graph,
                                     igraph_matrix_t *linkrank,
                                     igraph_vector_t *weights,
                                     igraph_real_t damping,
                                     igraph_neimode_t mode);

igraph_error_t igraph_personalized_linkrank(const igraph_t *graph,
                                            igraph_matrix_t *linkrank,
                                            igraph_vector_t *weights,
                                            igraph_vector_t *reset,
                                            igraph_real_t damping,
                                            igraph_neimode_t mode);

__END_DECLS

#endif /* IGRAPH_LINKRANK_H */
