import pandas as pd
import os,re,sys
import yaml
from tqdm import tqdm
import sqlite3
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import gc

from concurrent.futures import ProcessPoolExecutor, as_completed

def merge_qc(param):
    """
    Merges all QC results into a single file.
    :param param: Dictionary of parameters loaded from YAML.
    """

    output_dir = param['output_dir']
    step3_dir = os.path.join(output_dir, "step3-compare")
    step2_dir = os.path.join(output_dir, "step2-qc")

    # ✅ Ensure step3-compare directory exists
    os.makedirs(step3_dir, exist_ok=True)

    # ✅ Get all QC files
    files = [f for f in os.listdir(step2_dir) if f.endswith("_qc.tsv")]

    if not files:
        print("⚠ No QC files found! Skipping merging step.")
        return

    # ✅ Read and concatenate all QC files
    df_list = [pd.read_csv(os.path.join(step2_dir, f), sep="\t") for f in files]
    merged_df = pd.concat(df_list)

    # ✅ Save merged QC file
    merged_qc_path = os.path.join(step3_dir, "merged_qc.tsv")
    merged_df.to_csv(merged_qc_path, sep="\t", index=False)
    del merged_df
    gc.collect()

    print(f"✅ Merged QC saved to {merged_qc_path}")
        
        

def compare_all(param):
    output_dir = param['output_dir']
    step3_dir = os.path.join(output_dir,"step3-compare")
    step2_dir = os.path.join(output_dir,"step2-qc")
    if not os.path.exists(step3_dir):
        os.makedirs(step3_dir)
    files = [i for i in os.listdir(step2_dir) if i.endswith("_mc2.tsv")]
    
    df_list = []
    for i in files:
        df = pd.read_csv(os.path.join(step2_dir, i), sep="\t")
        df_list.append(df)
    
    common_peps = set(df_list[0]["PEP.StrippedSequence"])
    for i in tqdm(df_list):
        common_peps.intersection(set(i["PEP.StrippedSequence"]))
    print(f"Common peptides: {len(common_peps)}")
    
    med_list = []
    for df in df_list:
        med = df[df['PEP.StrippedSequence'].isin(common_peps)]['Missed.Cleavage.Ratio'].median()
        med_list.append(med)
        
    sns.histplot(med_list)
    
    # Save the figure to a PNG file
    hist_path = os.path.join(step3_dir, "histogram.png")
    plt.savefig(hist_path, dpi=300, bbox_inches="tight")
    plt.close()
    rows = []

    for df in tqdm(df_list):
        for index, row in df.iterrows():
            condition = df.columns[2].split(".")[0].split(" ")[1]
            mc_pep = row['PEP.StrippedSequence']
            nmc_pep = "NA"
            mc_pep_quant = list(row)[2]
            nmc_pep_quant = row["NMC.PEP.Quantity"]
            mcr = row["Missed.Cleavage.Ratio"]
            r = [condition, mc_pep, nmc_pep, mc_pep_quant, nmc_pep_quant, mcr]
            rows.append(r)
            
    rdf = pd.DataFrame(rows, columns=['Condition', 'MC_PEP', 'NMC_PEP', 'MC_PEP_Quant', 'NMC_PEP_Quant', 'MCR'])
    
    # rdf['Log2MC_PEP_Quant'] = np.log2(rdf['MC_PEP_Quant'])
    # rdf['Log2NMC_PEP_Quant'] = np.log2(rdf['NMC_PEP_Quant'])
    rdf['Log2MC_PEP_Quant'] = np.log2(rdf['MC_PEP_Quant'].replace(0, np.nan))
    rdf['Log2NMC_PEP_Quant'] = np.log2(rdf['NMC_PEP_Quant'].replace(0, np.nan))
    
    rdf2 = rdf.replace([np.inf, -np.inf], np.nan)
    rdf2 = rdf2.dropna()
    
    wide_df = rdf2.pivot(index='MC_PEP', columns='Condition', values='MCR').reset_index()
    wdf = wide_df.dropna().set_index('MC_PEP')
    
    # 3. Save the DataFrames to CSV files.
    rdf.to_csv(os.path.join(step3_dir, "rdf.csv"), index=False)
    rdf2.to_csv(os.path.join(step3_dir, "rdf2.csv"), index=False)
    wide_df.to_csv(os.path.join(step3_dir, "wide_df.csv"), index=False)
    wdf.to_csv(os.path.join(step3_dir, "wdf.csv"))
    del wdf, wide_df
    gc.collect()

    
    g = sns.clustermap(wdf, cmap="vlag", center=0, figsize=(20, 8))
    clustermap_path = os.path.join(step3_dir, "clustermap.png")
    g.savefig(clustermap_path)
    plt.close()

    g2 = sns.clustermap(wdf, cmap="vlag", center=0, figsize=(20, 8), z_score=0, )
    clustermap2_path = os.path.join(step3_dir, "clustermap_zscore.png")
    g2.savefig(clustermap2_path)
    plt.close()
    
    fig,ax = plt.subplots(figsize=(20, 8))
    sns.boxplot(data = rdf2, x= 'Condition', y = 'MCR')
    tmp = plt.xticks(rotation=90)
    boxplot_path = os.path.join(step3_dir, "boxplot.png")
    plt.savefig(boxplot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    
    fig,ax = plt.subplots(figsize=(20, 8))
    sns.violinplot(data = rdf, x= 'Condition', y = 'MCR', hue="Condition")
    tmp = plt.xticks(rotation=90)  
    violinplot_path = os.path.join(step3_dir, "violinplot.png")
    plt.savefig(violinplot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    print(f"All figures and DataFrames have been saved to {output_dir}.")
    
    # for trypsin/p
    for ii, df in tqdm(enumerate(df_list)):
        for index, row in df.iterrows():
            condition = df.columns[2].split(".")[0].split(" ")[1]
            mc_pep = row['PEP.StrippedSequence']
            nmc_pep = "NA"
            mc_pep_quant = list(row)[2]
            nmc_pep_quant = row["NMC.PEP.Quantity"]
            mcr = row["Missed.Cleavage.Ratio"]
            mc_count = row["Missed.Cleavages.Count"]
            if mc_count != 1:
                continue
            try:
                pos,aa = row["Missed.Cleavages.Sites"].split(",")
            except:
                print(ii,index, row["Missed.Cleavages.Sites"])
                sys.exit(1)
            
            pre_aa = mc_pep[int(pos)]
            r = [condition, mc_pep, nmc_pep, mc_pep_quant, nmc_pep_quant, mcr, pre_aa, aa]
            rows.append(r)
            
    rdf = pd.DataFrame(rows, columns=['Condition', 'MC_PEP', 'NMC_PEP', 'MC_PEP_Quant', 'NMC_PEP_Quant', 'MCR', "PRE_AA", "POST_AA"])
    
    rdf2 = rdf[rdf['POST_AA'] != "P"]
    
    rdf.to_csv(os.path.join(step3_dir, "new_rdf.csv"), index=False)
    rdf2.to_csv(os.path.join(step3_dir, "new_rdf2.csv"), index=False)
    del rdf, rdf2
    gc.collect()
    
    d = dict()
    for index,row in tqdm(rdf.iterrows()):
        sample = row['Condition']
        pre_aa = row['PRE_AA']
        post_aa = row['POST_AA']
        if sample not in d:
            d[sample] = dict()
        if pre_aa not in d[sample]:
            d[sample][pre_aa] = dict()
        if post_aa not in d[sample][pre_aa]:
            d[sample][pre_aa][post_aa] = 0
        d[sample][pre_aa][post_aa] += 1
    
    rows = []
    aa_list = sorted(["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"])
    for sample in d:
        for pre_aa in ["K", "R"]:
            x = []
            if pre_aa not in d[sample]:
                x = [0] * 20 
            else:
                for post_aa in aa_list:
                    x.append(d[sample][pre_aa].get(post_aa, 0))
            rows.append([sample, pre_aa] + x)
    
    xdf = pd.DataFrame(rows,columns=["Sample", "PRE_AA"] + aa_list)
    
    xdf.to_csv(os.path.join(step3_dir, "mc_aa_count.csv"), index=False)
    del xdf
    gc.collect()
    g3 = sns.clustermap(xdf.set_index(["Sample","PRE_AA"]), cmap="vlag", center=0, figsize=(20, 8), row_cluster=False, col_cluster=False)
    clustermap3_path = os.path.join(step3_dir, "heatmap_mc_aa_count.png")
    g3.savefig(clustermap3_path)
    plt.close() 