def modify_vcf_metadata(metadata: dict, csq_file: Optional[str] = None, header_file: Optional[str] = None) -> dict:
    """
    Modify metadata to include information about the CSQ annotation file and header file used for VCF export.
    """
    if header_file is not None and csq_file is not None:
        metadata=extract_csq_header(header_file, metadata)
        if not metadata["info"].get("CSQ"):
                mt = mt.drop(mt.CSQ)
    if csq_file is not None:
        metadata["info"].update({
            "consequence": {
                "Description": "Most severe consequence from VEP",
                "Number": "A",
                "Type": "String",
            },
            "gene": {
                "Description": "Gene affected by the most severe consequence from VEP",
                "Number": "A",
                "Type": "String",
            },
            "hgnc_id": {
                "Description": "HGNC id of the gene affected by the most severe consequence from VEP",
                "Number": "A",
                "Type": "String",
            }
        })
    return metadata

def extract_csq_header(header_file: str, metadata: dict) -> dict:
    """
    Modify metadata to include the CSQ description used for VCF export.
    """
    info={}
    with open(header_file) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                pattern = re.compile(
                r'ID=(?P<ID>[^,]+),'
                r'Number=(?P<Number>[^,]+),'
                r'Type=(?P<Type>[^,]+),'
                r'Description="(?P<Description>.*)"'
                )
                match = pattern.search(line)
                info = match.groupdict()

    if len(info) == 0:
        print("===No CSQ header found. CSQ won't be added to the VCFs===")
    else:
        metadata["info"]["CSQ"] = {
            "Description": info["Description"],
            "Number": info["Number"],
            "Type": info["Type"],
        }
    return metadata
