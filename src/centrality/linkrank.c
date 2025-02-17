#include "igraph.h"

igraph_error_t igraph_linkrank(const igraph_t *graph,
                               const igraph_vector_t *pagerank,
                               const igraph_vector_t *weights,
                               igraph_matrix_t *linkrank) {
    long int i, j, no_of_nodes = igraph_vcount(graph);
    igraph_real_t total_strength;
    igraph_vector_t strength;

    IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);
    igraph_vector_fill(&strength, 0.0);

    //Computes strength of each vertex (equal to the sum of outgoing edge weightss)
    for (i = 0; i < no_of_nodes; i++) {
        igraph_es_t es;
        igraph_eit_t eit;
        igraph_real_t sum = 0.0;

        IGRAPH_CHECK(igraph_incident(graph, &es, i, IGRAPH_OUT));
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
        
        while (!IGRAPH_EIT_END(eit)) {
            long int eid = IGRAPH_EIT_GET(eit);
          //outgoing weights sum
            sum += VECTOR(*weights)[eid]; 
            IGRAPH_EIT_NEXT(eit);
        }
        
        IGRAPH_EIT_DESTROY(&eit);
        IGRAPH_ES_DESTROY(&es);
        //Stregnth of the ith node
        VECTOR(strength)[i] = sum; 
    }
    //Compute LinkRank
    //LinkRank(e) = PageRank(v) * weight(e) / strength(v)
    for (i = 0; i < no_of_nodes; i++) {
        igraph_es_t es;
        igraph_eit_t eit;

        IGRAPH_CHECK(igraph_incident(graph, &es, i, IGRAPH_OUT));
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));

        while (!IGRAPH_EIT_END(eit)) {
            long int eid = IGRAPH_EIT_GET(eit);
            long int target;
            igraph_edge(graph, eid, NULL, &target);

            if (VECTOR(strength)[i] > 0) {
                MATRIX(*linkrank, i, target) = 
                    VECTOR(*pagerank)[i] * VECTOR(*weights)[eid] / VECTOR(strength)[i];
            } else {
                MATRIX(*linkrank, i, target) = 0.0;
            }

            IGRAPH_EIT_NEXT(eit);
        }

        IGRAPH_EIT_DESTROY(&eit);
        IGRAPH_ES_DESTROY(&es);
    }

    igraph_vector_destroy(&strength);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_linkrank_basic(
    const igraph_t *graph,
    igraph_matrix_t *linkrank,
    igraph_vector_t *weights,
    igraph_real_t damping,
    igraph_neimode_t mode) {
    
    igraph_vector_t pagerank;
    IGRAPH_VECTOR_INIT_FINALLY(&pagerank, igraph_vcount(graph));

    IGRAPH_CHECK(igraph_pagerank(graph, IGRAPH_PAGERANK_ALGO_PRPACK, &pagerank, NULL,
                                 igraph_vss_all(), mode, damping, weights, NULL));

    IGRAPH_CHECK(igraph_linkrank(graph, &pagerank, weights, linkrank));

    igraph_vector_destroy(&pagerank);
    IGRAPH_FINALLY_CLEAN(1);
    
    return IGRAPH_SUCCESS;
}

//Basic linkrank implementation 
igraph_error_t igraph_linkrank_basic(
    const igraph_t *graph,
    igraph_matrix_t *linkrank,
    igraph_vector_t *weights,
    igraph_real_t damping,
    igraph_neimode_t mode) {
    
    igraph_vector_t pagerank;
    IGRAPH_VECTOR_INIT_FINALLY(&pagerank, igraph_vcount(graph));

    IGRAPH_CHECK(igraph_pagerank(graph, IGRAPH_PAGERANK_ALGO_PRPACK, &pagerank, NULL,
                                 igraph_vss_all(), mode, damping, weights, NULL));

    IGRAPH_CHECK(igraph_linkrank(graph, &pagerank, weights, linkrank));

    igraph_vector_destroy(&pagerank);
    IGRAPH_FINALLY_CLEAN(1);
    
    return IGRAPH_SUCCESS;
}

//Personalized linkrank 
igraph_error_t igraph_personalized_linkrank(
    const igraph_t *graph,
    igraph_matrix_t *linkrank,
    igraph_vector_t *weights,
    igraph_vector_t *reset,
    igraph_real_t damping,
    igraph_neimode_t mode) {

    igraph_vector_t pagerank;
    IGRAPH_VECTOR_INIT_FINALLY(&pagerank, igraph_vcount(graph));

    IGRAPH_CHECK(igraph_personalized_pagerank(graph, IGRAPH_PAGERANK_ALGO_PRPACK, 
                                              &pagerank, NULL, igraph_vss_all(), 
                                              mode, damping, reset, weights, NULL));

    IGRAPH_CHECK(igraph_linkrank(graph, &pagerank, weights, linkrank));

    igraph_vector_destroy(&pagerank);
    IGRAPH_FINALLY_CLEAN(1);
    
    return IGRAPH_SUCCESS;
}

