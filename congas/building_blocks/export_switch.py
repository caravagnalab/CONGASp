def export_switch(x):
    export = ["CNV_probabilities"]
    if  x._params['equal_mixture_weights'] or x._params['multiome']:
        export.append("mixture_weights")
    if 'data_rna' in x._data:
        if not x._params['equal_mixture_weights']:
            export.append("mixture_weights_rna")
        export.append("segment_factor_rna")
        if x._params['likelihood_rna'] == "NB":
            export.append("NB_size_rna")
        if x._params['likelihood_rna'] in ["N", "G"]:
            export.append("norm_sd_rna")
        else:
            export.append("segment_factor_atac")

    if 'data_atac' in x._data:
        if not x._params['equal_mixture_weights']:
            export.append("mixture_weights_atac")
        if x._params['likelihood_atac'] == "NB":
            export.append("NB_size_atac")
        if x._params['likelihood_atac'] in ["N", "G"]:
            export.append("norm_sd_atac")
        else:
            export.append("segment_factor_atac")


    return export
