import os
import pickle as pkl
import numpy as np
import pandas as pd
import sys

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

def deseq(counts_df, metadata, design_factors):
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=design_factors,
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    stat_res = DeseqStats(dds, inference=inference)

    results = pd.DataFrame(stat_res.results_df)
    results.head()
    results["log10_padj"] = -np.log10(results["padj"])

    return results

if __name__ == "__main__":
    # 引数の読み込み
    # カウント値のデータ
    count_file_path = sys.argv[1] 
    counts_df = pd.read_table(count_file_path, sep="\t", index_col=0, header=0)
    # metadata
    metadata_file_path = sys.argv[2]
    metadata = pd.read_table(metadata_file_path, sep="\t", index_col=0, header=0)

    # metadata中のラベル
    design_factors = sys.argv[3] 

    results = deseq(counts_df, metadata,design_factors)
