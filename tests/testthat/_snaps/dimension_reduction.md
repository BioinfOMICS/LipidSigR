# dr_pca scaling & centering parameters

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = NULL, centering = TRUE, clustering = "kmeans",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `dr_pca()`:
      ! scaling must be a logical value.

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = "yes", clustering = "kmeans",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `dr_pca()`:
      ! centering must be a logical value.

# dr_pca clustering & related parameters

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "k-means",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "kmeans", "kmedoids", "hclustering", "dbscan", "group_info"

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "kmeans",
        cluster_num = 0, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! cluster_num must be a numeric value between 2 and 22.

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "kmedoids",
        cluster_num = 2, kmedoids_metric = "invalid", distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! kmedoids_metric must be one of "euclidean" or "manhattan".

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "hclustering",
        cluster_num = 2, kmedoids_metric = NULL, distfun = "invalid", hclustfun = "ward.D",
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! distfun must be one of "pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski".

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "hclustering",
        cluster_num = 2, kmedoids_metric = NULL, distfun = "maximum", hclustfun = "invalid",
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! hclustfun must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid".

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "dbscan",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = 100, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! eps must be must be numeric value >= 0.

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "dbscan",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = 0, minPts = -1, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! minPts must be must be numeric value >= 0.

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "hclustering",
        cluster_num = 24, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! cluster_num must be a numeric value between 2 and 23.

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "group_info",
        cluster_num = 24, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! cluster_num must be the number of groups in group_info

---

    Code
      suppressWarnings(dr_pca(profiling_se, scaling = TRUE, centering = TRUE,
        clustering = "group_info", cluster_num = 24, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2),
        plot_topN = 10))
    Condition
      Error in `.dr_check()`:
      ! clustering select "group_info" in processed_se must provide group_info

# dr_pca plot_topN parameters

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "kmeans",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 0))
    Condition
      Error in `dr_pca()`:
      ! plot_topN must be a numeric value less than the number of lipids.

# dr_pca handles insufficient number of significant lipids condition

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "kmeans",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.deSig_data()`:
      ! Insufficient number of significant lipids.

---

    Code
      suppressWarnings(dr_pca(deSp_se, scaling = TRUE, centering = TRUE, clustering = "kmeans",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10))
    Condition
      Error in `.deSig_data()`:
      ! Insufficient number of significant lipids.

# dr_tsne can return error message for incorrect input

    Code
      suppressWarnings(dr_tsne(profiling_se, pca = NULL, perplexity = 5, max_iter = 500,
        clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `dr_tsne()`:
      ! pca must be a logical value.

---

    Code
      suppressWarnings(dr_tsne(profiling_se, pca = TRUE, perplexity = -1, max_iter = 500,
        clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `dr_tsne()`:
      ! perplexity cannot be a negative value.

---

    Code
      suppressWarnings(dr_tsne(profiling_se, pca = TRUE, perplexity = 5, max_iter = -
        1, clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `dr_tsne()`:
      ! max_iter cannot be a negative value.

# dr_umap returns error message for incorrect input

    Code
      suppressWarnings(dr_umap(profiling_se, n_neighbors = 0, scaling = TRUE,
        umap_metric = "euclidean", clustering = "kmeans", cluster_num = 2,
        kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `dr_umap()`:
      ! n_neighbors must be >= 2.

---

    Code
      suppressWarnings(dr_umap(profiling_se, n_neighbors = 150, scaling = "YES",
        umap_metric = "euclidean", clustering = "kmeans", cluster_num = 2,
        kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "Z", "scale", "maxabs", "range", "colrange", "none"

---

    Code
      suppressWarnings(dr_umap(profiling_se, n_neighbors = 150, scaling = TRUE,
        umap_metric = "Euclidean", clustering = "kmeans", cluster_num = 2,
        kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "euclidean", "cosine", "manhattan", "hamming", "correlation"

# dr_pcoa dist_method & correction parameters

    Code
      suppressWarnings(dr_pcoa(profiling_se, dist_method = "invalid_method",
        correction = "cailliez", clustering = "kmeans", cluster_num = 2,
        kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "euclidean", "manhattan", "mahalanobis", "gower", "chord", "canberra"

---

    Code
      suppressWarnings(dr_pcoa(profiling_se, dist_method = "euclidean", correction = "invalid",
        clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "cailliez", "lingoes", "none"

# dr_pcoa clustering & related parameters

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "k-means", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "kmeans", "kmedoids", "hclustering", "dbscan", "group_info"

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "kmeans", cluster_num = 0, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! cluster_num must be a numeric value between 2 and 22.

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "kmedoids", cluster_num = 2, kmedoids_metric = "invalid",
        distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! kmedoids_metric must be one of "euclidean" or "manhattan".

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "hclustering", cluster_num = 2, kmedoids_metric = NULL, distfun = "invalid",
        hclustfun = "ward.D", eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! distfun must be one of "pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski".

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "hclustering", cluster_num = 2, kmedoids_metric = NULL, distfun = "maximum",
        hclustfun = "invalid", eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! hclustfun must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid".

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "dbscan", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = 100))
    Condition
      Error in `.dr_check()`:
      ! eps must be must be numeric value >= 0.

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "dbscan", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = 0, minPts = -1))
    Condition
      Error in `.dr_check()`:
      ! minPts must be must be numeric value >= 0.

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "hclustering", cluster_num = 24, kmedoids_metric = NULL,
        distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! cluster_num must be a numeric value between 2 and 23.

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "group_info", cluster_num = 24, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! cluster_num must be the number of groups in group_info

---

    Code
      suppressWarnings(dr_pcoa(profiling_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "group_info", cluster_num = 24, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.dr_check()`:
      ! clustering select "group_info" in processed_se must provide group_info

# dr_pcoa handles insufficient number of significant lipids condition

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.deSig_data()`:
      ! Insufficient number of significant lipids.

---

    Code
      suppressWarnings(dr_pcoa(deSp_se, dist_method = "euclidean", correction = "cailliez",
        clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL, distfun = NULL,
        hclustfun = NULL, eps = NULL, minPts = NULL))
    Condition
      Error in `.deSig_data()`:
      ! Insufficient number of significant lipids.

# dr_plsda returns correct output results

    Code
      dr_plsda(profiling_se, ncomp = 2, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL)
    Condition
      Error in `dr_plsda()`:
      ! Incorrect SummarizedExperiment structure. Please use the output from
                   deSp_twoGroup, deSp_multiGroup, deChar_twoGroup, or deChar_multiGroup function.

# dr_plsda returns error message for incorrect input.

    Code
      suppressWarnings(dr_plsda(deSp_se, ncomp = -1, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL))
    Condition
      Error in `dr_plsda()`:
      ! ncomp cannot be a negative value or NULL.

---

    Code
      suppressWarnings(dr_plsda(deSp_se, ncomp = 2, scaling = 1, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL))
    Condition
      Error in `dr_plsda()`:
      ! scaling must be a logical value.

# dr_oplsda returns correct output results

    Code
      dr_oplsda(profiling_se, ncomp = 2, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL)
    Condition
      Error in `dr_oplsda()`:
      ! Incorrect SummarizedExperiment structure. Please use the output from
                   deSp_twoGroup, deSp_multiGroup, deChar_twoGroup, or deChar_multiGroup function.

# dr_oplsda returns error message for incorrect input.

    Code
      suppressWarnings(dr_oplsda(deSp_se, ncomp = -1, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL))
    Condition
      Error in `dr_oplsda()`:
      ! ncomp cannot be a negative value or NULL.

---

    Code
      suppressWarnings(dr_oplsda(deSp_se, ncomp = 2, scaling = 1, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL))
    Condition
      Error in `dr_oplsda()`:
      ! scaling must be a logical value.

---

    Code
      dr_oplsda(deSp_se_multi, ncomp = 2, scaling = TRUE, clustering = "group_info",
        cluster_num = 3, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL)
    Condition
      Error in `dr_oplsda()`:
      ! OPLS-DA only available for binary classification.

---

    Code
      dr_oplsda(deChar_se, ncomp = 2, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL)
    Output
      No model was built because the first predictive component was already not significant
      OPLS-DA
      23 samples x 13 variables and 1 response
      standard scaling of predictors and response(s)
      Empty 'opls' object
    Condition
      Error in `dr_oplsda()`:
      ! Model was not built successfully, skipping downstream analysis

# dr_splsda returns correct output results

    Code
      dr_splsda(profiling_se, ncomp = 2, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL)
    Condition
      Error in `dr_splsda()`:
      ! Incorrect SummarizedExperiment structure. Please use the output from
                   deSp_twoGroup, deSp_multiGroup, deChar_twoGroup, or deChar_multiGroup function.

# dr_splsda returns error message for incorrect input.

    Code
      suppressWarnings(dr_splsda(deSp_se, ncomp = -1, scaling = TRUE, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL))
    Condition
      Error in `dr_splsda()`:
      ! ncomp cannot be a negative value or NULL.

---

    Code
      suppressWarnings(dr_splsda(deSp_se, ncomp = 2, scaling = 1, clustering = "group_info",
        cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
        eps = NULL, minPts = NULL))
    Condition
      Error in `dr_splsda()`:
      ! scaling must be a logical value.

