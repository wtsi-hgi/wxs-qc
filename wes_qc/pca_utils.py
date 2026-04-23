from wes_qc import hail_patches
from gnomad.sample_qc.ancestry import pc_project

def prune_mt(mt: hl.MatrixTable, ld_prune_args, **kwargs) -> hl.MatrixTable:
    '''
    Prune variants in linkage disequilibrium (LD) on autosomes.
    '''
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
    '''
    Run PCA on reference data and project study samples onto the PC space.
    '''
    print("=== Running PCA")
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_ref.GT, k=pca_components, compute_loadings=True)
    pca_af_ht = mt_ref.annotate_rows(pca_af=hl.agg.mean(mt_ref.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    print("=== Running PC projection")
    projection_pca_scores = pc_project(mt_study, pca_loadings, loading_location="loadings", af_location="pca_af")
    union_pca_scores = pca_scores.union(projection_pca_scores)

    return union_pca_scores, pca_scores, pca_loadings

