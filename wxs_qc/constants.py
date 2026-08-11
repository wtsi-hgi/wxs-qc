"""
Set of literals and constants used to filter and annotate variants
"""

from dataclasses import dataclass

VEP_CONSEQUENCE_ALL = """
3_prime_utr_variant
5_prime_utr_variant
coding_sequence_variant
downstream_gene_variant
frameshift_variant
incomplete_terminal_codon_variant
intergenic_variant
intron_variant
mature_mirna_variant
missense_variant
nmd_transcript_variant
non_coding_transcript_exon_variant
non_coding_transcript_variant
protein_altering_variant
splice_acceptor_variant
splice_donor_5th_base_variant
splice_donor_region_variant
splice_donor_variant
splice_polypyrimidine_tract_variant
splice_region_variant
start_retained_variant
stop_retained_variant
synonymous_variant
upstream_gene_variant""".split()

VEP_CONSEQUENCE_CODING = """
frameshift_variant
missense_variant
splice_acceptor_variant
splice_donor_variant
start_lost_variant
start_retained_variant
stop_retained_variant
stop_gained_variant
stop_lost_variant
inframe_deletion_variant
inframe_insertion_variant
synonymous_variant
""".split()

VEP_CONSEQUENCE_NON_CODING = """
3_prime_utr_variant
5_prime_utr_variant
downstream_gene_variant
intergenic_variant
intron_variant
upstream_gene_variant
""".split()


# Default threshold used in PC-Relate by gnomAD
# To define parent-child relationship
# They are highly conservative, so there is no sense in putting it in the config file.
@dataclass(frozen=True)
class PedigreeThresholds:
    kin_max: float = 0.4
    kin_min: float = 0.19
    ibd0_max: float = 0.025
    ibd1_min: float = 0.75
    ibd2_max: float = 0.125


PEDIGREE_THRESHOLDS = PedigreeThresholds()


##################
# Visualizations #
##################

plots_text_size = "18pt"
width = 1200
height = 1000
