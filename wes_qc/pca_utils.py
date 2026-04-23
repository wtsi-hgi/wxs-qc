from wes_qc import hail_patches
from gnomad.sample_qc.ancestry import pc_project

def prune_mt(mt: hl.MatrixTable, ld_prune_args, **kwargs) -> hl.MatrixTable:
    """
    Splits multiallelic sites and runs ld pruning
    Filter to autosomes before LD pruning to decrease sample size - autosomes only wanted for later steps
    :param MatrixTable mt: input MT containing variants to be pruned
    :param dict config:
    :return: Pruned MatrixTable
    :rtype: hl.MatrixTable

    `hl.ld_prune` returns a maximal subset of variants that are nearly uncorrelated within each window.
    Requires the dataset to contain only diploid genotype calls and no multiallelic variants.
    See also: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune
    """

    print("=== Filtering to autosomes")
    mt = mt.filter_rows(mt.locus.in_autosome())
    print("=== Splitting multiallelic sites")
    mt = hail_patches.split_multi_hts(
        mt, recalculate_gq=False
    )  # this shouldn't do anything as only biallelic sites are used
    print("=== Performing LD pruning")
    pruned_ht = hl.ld_prune(mt.GT, **ld_prune_args)
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    return pruned_mt

def run_pc_project(mt_ref, mt_study, pca_components):
    print("=== Running PCA")
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_ref.GT, k=pca_components, compute_loadings=True)
    pca_af_ht = mt_ref.annotate_rows(pca_af=hl.agg.mean(mt_ref.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    print("=== Running PC projection")
    projection_pca_scores = pc_project(mt_study, pca_loadings, loading_location="loadings", af_location="pca_af")
    union_pca_scores = pca_scores.union(projection_pca_scores)

    return union_pca_scores, pca_scores, pca_loadings

