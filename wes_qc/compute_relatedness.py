import hail as hl
from utils.utils import path_spark
from wes_qc.pca_utils import prune_mt, run_pc_project

def run_king(mt: hl.MatrixTable, king_args: dict, prune_args: dict) -> (hl.MatrixTable, hl.MatrixTable):
    '''
    Separate related and unrelated samples using the KING kinship estimator.
    '''
    print("=== LD pruning before KING")
    pruned_mt= prune_mt(mt, prune_args["ld_prune_args"])
    pruned_mt=pruned_mt.checkpoint(path_spark(prune_args["pruned_mt_file"]), overwrite=True)
    print("=== Running KING")
    king_mt = hl.king(pruned_mt.GT)
    king_ht=king_mt.entries()
    king_ht=king_ht.checkpoint(path_spark(king_args["king_output_file"]), overwrite=True)
    related_pairs_ht = king_ht.filter((king_ht.phi > king_args["kinship_threshold"]) & (king_ht.s_1!=king_ht.s))
    print("=== Identifying related samples")
    samples_to_remove=hl.maximal_independent_set(related_pairs_ht.s_1, related_pairs_ht.s, keep=False)
    unrelated_mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]), keep=False)
    related_mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]))
    print("=== LD pruning unrelated samples")
    pruned_unrelated_mt = prune_mt(unrelated_mt, prune_args["ld_prune_args"])
    pruned_unrelated_mt=pruned_unrelated_mt.checkpoint(path_spark(king_args["unrelated_king_file"]), overwrite=True)
    related_mt = related_mt.semi_join_rows(pruned_unrelated_mt.rows())
    related_mt=related_mt.checkpoint(path_spark(king_args["related_king_file"]), overwrite=True)
    return related_mt, pruned_unrelated_mt

def prune_pc_relate(
    mt: hl.MatrixTable, prune_params: dict, king_params: dict, pc_relate_params: dict, dataset: str, **kwargs
) -> (hl.Table, hl.Table):
    '''
    Complete PC-Relate workflow for identifying and removing related samples.
    '''
    #Running KING
    related_mt, unrelated_mt = run_king(mt, king_params, prune_params)
    #Running PCA
    union_pca_scores, pca_scores, pca_loadings, _ = run_pc_project(unrelated_mt, related_mt, pc_relate_params["pca_components"])
    union_pca_scores=union_pca_scores.checkpoint(path_spark(pc_relate_params["scores_file"]), overwrite=True)
    pca_scores=pca_scores.checkpoint(path_spark(pc_relate_params["unrelated_samples_scores_file"]), overwrite=True)
    pca_loadings=pca_loadings.checkpoint(path_spark(pc_relate_params["pca_loadings_file_pc_relate"]), overwrite=True)
    if dataset == "study":
        related_mt = related_mt.drop(#For some reason unrelated_mt doesn't have these in study data
            "AD",
            "DP",
            "GQ",
            "MIN_DP",
            "PGT",
            "PID",
            "PL",
            "PS",
            "SB",
            "RGQ"
        )
    pruned_mt = related_mt.union_cols(unrelated_mt)
    print("=== Calculating relatedness")
    #calculate relatedness using PC relate
    relatedness_ht = hl.pc_relate(pruned_mt.GT, scores_expr=union_pca_scores[pruned_mt.col_key].scores, **pc_relate_params["pc_relate_args"])
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht[pc_relate_params["relatedness_column"]] > pc_relate_params["relatedness_threshold"])
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    return related_samples_to_remove, relatedness_ht

