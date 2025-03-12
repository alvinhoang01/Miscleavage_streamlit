import pandas as pd
import os,re,sys
import yaml
from tqdm import tqdm
import sqlite3
import tempfile
import shutil


from concurrent.futures import ProcessPoolExecutor, as_completed



def is_unique_peptide(peptide, pep_map):
    if peptide not in pep_map:
        return 'NA'
    if len(pep_map[peptide]) == 1:
        return True
    else:
        return False

# calculate Missed Cleavage Rate (MCR) for single table (sample)
def check_missed_cleavages(sequence, target_AAs = "KR", pos=1, omit_AAs = "P"):
    target_aa_list = list(target_AAs)
    omit_aa_list = list(omit_AAs)
    sequence = list(sequence)
    mc_list = []
    if pos == 1:
        for i, aa in enumerate(sequence):
            if aa in target_aa_list:
                if i < len(sequence) - 1 and sequence[i + pos] not in omit_aa_list:
                    mc_list.append(i)
    elif pos == -1:
        for i, aa in enumerate(sequence):
            if i == 0 :
                continue
            if aa in target_aa_list:
                if i < len(sequence) and sequence[i + pos] not in omit_aa_list:
                    mc_list.append(i)
                    
    return mc_list

def check_missed_cleavages_for_trypsin(sequence,pep_map):
    
    mc_list = []
    if sequence not in pep_map:
        return mc_list
    
    proteins = pep_map.get(sequence,"")
    proteins = [i.split(":") for i in proteins] 
    sequence = list(sequence)
    for i,aa in enumerate(sequence):
        if i == 0 and sequence[0] in ["K","R"] and proteins[0][2] in ["K","R"]:
            continue
        elif i == len(sequence) - 2 and sequence[-1] in ["K","R"] and sequence[-2] in ["K","R"]:
            continue
        elif aa in ['K',"R"] and i < len(sequence) - 1:
            mc_list.append((i,sequence[i+1]))
    return mc_list

 

def qc_all(param):
    print(f"🚀 Running QC for all samples...")

    # ✅ Use a Temporary Directory for QC Output
    temp_dir = tempfile.mkdtemp()
    step2_qc_dir = os.path.join(temp_dir, "step2-qc")
    os.makedirs(step2_qc_dir, exist_ok=True)
    print(f"📂 QC output directory created: {step2_qc_dir}")

    # ✅ Ensure the required input files exist
    step1_dir = os.path.join(param['output_dir'], "step1-split")
    sqlite_path = os.path.join(param['output_dir'], "peptides.sqlite")

    if not os.path.exists(sqlite_path):
        print("❌ Error: Missing peptides.sqlite. Run 'Prepare Task' first.")
        return None, None

    if not os.path.exists(step1_dir) or not os.listdir(step1_dir):
        print("❌ Error: Missing step1-split folder. Run 'Split Task' first.")
        return None, None

    files = [i for i in os.listdir(step1_dir) if re.search(r".split.tsv", i)]
    workers = int(param['workers'])
    enz = param["enzyme"]

    if not files:
        print("⚠ No split files found! Exiting QC process.")
        return None, None

    # ✅ Run QC in parallel using ProcessPoolExecutor
    print(f"🧪 Processing QC with {workers} workers...")
    futures = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        if enz == "trypsin/p":
            for f in files:
                futures.append(executor.submit(qc_one_trypsinp, os.path.join(step1_dir, f), step2_qc_dir, sqlite_path, enz))
        else:
            for f in files:
                futures.append(executor.submit(qc_one, os.path.join(step1_dir, f), step2_qc_dir, sqlite_path, enz))
        
        for future in as_completed(futures):
            try:
                future.result()  # Get results and handle any errors
            except Exception as e:
                print(f"❌ Error processing a file: {e}")

    # ✅ Zip the QC output folder
    zip_qc_path = os.path.join(temp_dir, "step2-qc.zip")
    shutil.make_archive(zip_qc_path.replace(".zip", ""), 'zip', step2_qc_dir)
    print(f"📁 QC output files zipped at: {zip_qc_path}")

    return zip_qc_path, temp_dir  # ✅ Return zip path & temp folder
        
def calc_quant_for_fragment_pep(df_mc, pep_map, pep_quant_map, sample):
    rows = []
    for index, row in tqdm(df_mc.iterrows()):
        pos, aa = row['Missed.Cleavages.Sites'][0]
        seq = row['PEP.StrippedSequence']
        pep1 = seq[:pos+1]
        pep2 = seq[pos+1:]
        row['PEP.1'] = pep1
        row['PEP.2'] = pep2
        
        # uniqueness of peptides
        pep1_uniq = is_unique_peptide(pep1, pep_map)
        pep2_uniq = is_unique_peptide(pep2, pep_map)
        row['PEP.1.Uniquness'] = pep1_uniq
        row['PEP.2.Uniquness'] = pep2_uniq
        
        # quantity of non-missed cleaved peptide candidates
        pep1_quant = pep_quant_map.get(pep1, 0)
        pep2_quant = pep_quant_map.get(pep2, 0)
        
        row['PEP.1.Quantity'] = pep1_quant
        row['PEP.2.Quantity'] = pep2_quant
        
        # quantity of non-missed cleaved peptide
        values = []
        if pep1_uniq == True and pep1_quant > 0:
            values.append(pep1_quant)
        if pep2_uniq == True and pep2_quant > 0:
            values.append(pep2_quant)
        max_value = max(values) if values else 0
        
        row['NMC.PEP.Quantity'] = max_value
        
        # missed cleaved ratio
        mcr = row[sample] / (max_value  +  row[sample]) if max_value > 0 else 1
        row['Missed.Cleavage.Ratio'] = mcr * 100
        
        rows.append(row)
    df_mc2 = pd.DataFrame(rows)
    return df_mc2


def qc_one_trypsinp(path, output_dir,sqlite_path, enz):
    
    # Connect to the SQLite database
    print(f"Fetch the protein information for each peptide from {sqlite_path}")
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    
    # peptide_list = list(df_nodup['PEP.StrippedSequence'])
    # placeholders = ','.join('?' for _ in peptide_list)
    # query = f"SELECT peptide, protein FROM peptides WHERE peptide IN ({placeholders})"
    
    # Execute the query with the peptide list as parameters
    # cursor.execute(query, peptide_list)
    
    cursor.execute("SELECT peptide, protein FROM peptides")

    rows = cursor.fetchall()

    # Convert the results into a dictionary (if each peptide is unique)
    pep_map = {peptide: protein.split(";") for peptide, protein in tqdm(rows)}

    # Close the connection
    conn.close()
    
    
    df = pd.read_csv(path,sep="\t")
    sample = df.columns[-1]
    # QC 1 identified peptides
    # remove nan values, only keep the peptides identified in this sample
    df_na = df.dropna()
    
    
    peptide_count = df_na['PEP.StrippedSequence'].unique().shape[0]
    protein_count = df_na['PG.ProteinNames'].unique().shape[0]

    protein_human_count = df_na[df_na['PG.ProteinNames'].str.contains('_HUMAN')]['PG.ProteinNames'].unique().shape[0]
    protein_mouse_count = df_na[df_na['PG.ProteinNames'].str.contains('_MOUSE')]['PG.ProteinNames'].unique().shape[0]
    
    # remove duplicates
    df_nodup = df_na.drop_duplicates()
    df_nodup = df_nodup.copy()
    

    df_nodup.loc[:,'Missed.Cleavages.Sites'] = df_nodup['PEP.StrippedSequence'].apply(lambda x: check_missed_cleavages_for_trypsin(x, pep_map))
    df_nodup = df_nodup.copy()
    df_nodup.loc[:,'Missed.Cleavages.Count'] = df_nodup['Missed.Cleavages.Sites'].apply(len)
    
    pep_quant_map = dict()
    for index,row in tqdm(df_nodup.iterrows()):
        key = row['PEP.StrippedSequence']
        value = row[sample]
        if key not in pep_quant_map:
            pep_quant_map[key] = value
        else:
            print(key)
            
    # calcualte ID-based missed cleavage rate
    rows = []
    
    for index,row in df_nodup.iterrows():
        mc_sites = row['Missed.Cleavages.Sites']
        notP = False
        for site, aa in mc_sites:
            if aa not in ["P"]:
                notP = True
                break
        row["Missed.Cleavages.notP"] = notP
        rows.append(row)
        
    df_nodup = pd.DataFrame(rows)
    
    mc_pep_count_withP = df_nodup[(df_nodup['Missed.Cleavages.Count'] > 0)].shape[0]
        
    mc_pep_count = df_nodup[(df_nodup['Missed.Cleavages.Count'] > 0) & (df_nodup["Missed.Cleavages.notP"]==True)].shape[0] 
    # mc_pep_count = 0
    
    mcr_pep = mc_pep_count / peptide_count
    
    

    
    print("Calculating the peptide uniqueness...")
    df_nodup = df_nodup.copy()
    df_nodup.loc[:,'Uniquness'] = [is_unique_peptide(i,pep_map) for i in df_nodup['PEP.StrippedSequence']]
    
    ratio_peptide_uniquness = df_nodup[df_nodup['Uniquness'] == True].shape[0] / df_nodup.shape[0]
    
    ratio_multiproteins_in_group = df_nodup[df_nodup['PG.ProteinNames'].str.contains(';')].shape[0]/df_nodup.shape[0]
    
    ratio_uniquness_mc1_peptide = df_nodup[(df_nodup['Missed.Cleavages.Count']==1)&(df_nodup['Uniquness']==True)& (df_nodup["Missed.Cleavages.notP"]==True)].shape[0] / df_nodup[(df_nodup['Missed.Cleavages.Count']==1) & (df_nodup["Missed.Cleavages.notP"]==True)].shape[0]
    
    sum_mc_pep_quant = df_nodup[(df_nodup['Missed.Cleavages.Count'] > 0) & (df_nodup["Missed.Cleavages.notP"]==True)][sample].sum()
   
    
    sum_pep_quant = df_nodup[sample].sum()
    mcr_pep_quant = sum_mc_pep_quant / sum_pep_quant

    
    df_nodup_uniq = df_nodup[df_nodup['Uniquness'] == True]
    
    df_mc  = df_nodup_uniq[(df_nodup_uniq['Missed.Cleavages.Count'] == 1) ]
    
    ratio_mc1_and_uniq_peptide = df_mc.shape[0] / df_nodup_uniq.shape[0]
    
    MC1_peptide_count = df_nodup_uniq[(df_nodup_uniq['Missed.Cleavages.Count'] == 1) & (df_nodup_uniq["Missed.Cleavages.notP"]==True)].shape[0]
    MC2_peptide_count = df_nodup_uniq[(df_nodup_uniq['Missed.Cleavages.Count'] == 2) & (df_nodup_uniq["Missed.Cleavages.notP"]==True)].shape[0]
    
    
    df_mc2 = calc_quant_for_fragment_pep(df_mc, pep_map, pep_quant_map, sample)
    
    
    
    mc100_pep_count = df_mc2[(df_mc2["Missed.Cleavage.Ratio"]==100)& (df_mc2["Missed.Cleavages.notP"]==True) ].shape[0]/df_mc2.shape[0]
    
    mc100_pep_count_len = 0
    for index,row in df_mc2.iterrows():
        if row["Missed.Cleavage.Ratio"] == 100 and row["Missed.Cleavages.notP"] == True:
            if len(row['PEP.1']) >=7 or len(row['PEP.2']) >= 7:
                mc100_pep_count_len += 1
    mc100_pep_count_len = mc100_pep_count_len / df_mc2.shape[0]
    
    base_name = re.sub('.split.tsv','',os.path.basename(path))
    
    
    
    out_mc2_path = os.path.join(output_dir,base_name + "_mc2.tsv")
    
    rows = []
    for index,row in df_mc2.iterrows():
        sites = row['Missed.Cleavages.Sites']
        sites = ";".join([f'{pos},{aa}' for pos,aa in sites])
        row['Missed.Cleavages.Sites'] = sites
        rows.append(row)
    df_mc2 = pd.DataFrame(rows)
    
    df_mc2.to_csv(out_mc2_path,sep="\t",index=False)
    
    
    print(f"Results have been written to {out_mc2_path}")
    output_path = os.path.join(output_dir,base_name + "_qc.tsv")
    
    results = {
        'sample_name': base_name,
        'peptide_count': peptide_count,
        'protein_count': protein_count,
        'protein_human_count': protein_human_count,
        'protein_mouse_count': protein_mouse_count,
        'mc_pep_count_withP': mc_pep_count_withP,
        'mc_pep_count': mc_pep_count,
        'mcr_pep': mcr_pep,
        'ratio_peptide_uniquness': ratio_peptide_uniquness,
        'ratio_multiproteins_in_group': ratio_multiproteins_in_group,
        'ratio_uniquness_mc1_peptide': ratio_uniquness_mc1_peptide,
        'sum_mc_pep_quant': sum_mc_pep_quant,
        'sum_pep_quant': sum_pep_quant,
        'mcr_pep_quant': mcr_pep_quant,
        'ratio_mc1_and_uniq_peptide': ratio_mc1_and_uniq_peptide,
        'mc1_peptide_count': MC1_peptide_count,
        'mc2_peptide_count': MC2_peptide_count,
        'mc100_pep_count': mc100_pep_count,
        "mc100_pep_count_len": mc100_pep_count_len
    }

    rdf = pd.DataFrame([results])

    
    rdf.to_csv(output_path,sep="\t",index=False)

    print(f"Results have been written to {output_path}")
    
def qc_one(path, output_dir,sqlite_path, enz):
    
    # Connect to the SQLite database
    print(f"Fetch the protein information for each peptide from {sqlite_path}")
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    
    # peptide_list = list(df_nodup['PEP.StrippedSequence'])
    # placeholders = ','.join('?' for _ in peptide_list)
    # query = f"SELECT peptide, protein FROM peptides WHERE peptide IN ({placeholders})"
    
    # Execute the query with the peptide list as parameters
    # cursor.execute(query, peptide_list)
    
    cursor.execute("SELECT peptide, protein FROM peptides")

    rows = cursor.fetchall()

    # Convert the results into a dictionary (if each peptide is unique)
    pep_map = {peptide: protein.split(";") for peptide, protein in tqdm(rows)}

    # Close the connection
    conn.close()
    
    
    df = pd.read_csv(path,sep="\t")
    sample = df.columns[-1]
    # QC 1 identified peptides
    # remove nan values, only keep the peptides identified in this sample
    df_na = df.dropna()
    
    
    peptide_count = df_na['PEP.StrippedSequence'].unique().shape[0]
    protein_count = df_na['PG.ProteinNames'].unique().shape[0]

    protein_human_count = df_na[df_na['PG.ProteinNames'].str.contains('_HUMAN')]['PG.ProteinNames'].unique().shape[0]
    protein_mouse_count = df_na[df_na['PG.ProteinNames'].str.contains('_MOUSE')]['PG.ProteinNames'].unique().shape[0]
    
    # remove duplicates
    df_nodup = df_na.drop_duplicates()
    df_nodup = df_nodup.copy()
    

    df_nodup.loc[:,'Missed.Cleavages.Sites'] = df_nodup['PEP.StrippedSequence'].apply(lambda x: check_missed_cleavages(x, target_AAs="KR", pos=-1, omit_AAs="P"))
    df_nodup = df_nodup.copy()
    df_nodup.loc[:,'Missed.Cleavages.Count'] = df_nodup['Missed.Cleavages.Sites'].apply(len)
    
    pep_quant_map = dict()
    for index,row in tqdm(df_nodup.iterrows()):
        key = row['PEP.StrippedSequence']
        value = row[sample]
        if key not in pep_quant_map:
            pep_quant_map[key] = value
        else:
            print(key)
            
    # calcualte ID-based missed cleavage rate
    mc_pep_count = df_nodup[df_nodup['Missed.Cleavages.Count'] > 0].shape[0] 
    mcr_pep = mc_pep_count / peptide_count
    
    

    
    print("Calculating the peptide uniqueness...")
    df_nodup = df_nodup.copy()
    df_nodup.loc[:,'Uniquness'] = [is_unique_peptide(i,pep_map) for i in df_nodup['PEP.StrippedSequence']]
    
    ratio_peptide_uniquness = df_nodup[df_nodup['Uniquness'] == True].shape[0] / df_nodup.shape[0]
    
    ratio_multiproteins_in_group = df_nodup[df_nodup['PG.ProteinNames'].str.contains(';')].shape[0]/df_nodup.shape[0]
    
    ratio_uniquness_mc1_peptide = df_nodup[(df_nodup['Missed.Cleavages.Count']==1)&(df_nodup['Uniquness']==True)].shape[0] / df_nodup[(df_nodup['Missed.Cleavages.Count']==1)].shape[0]
    
    sum_mc_pep_quant = df_nodup[df_nodup['Missed.Cleavages.Count'] > 0][sample].sum()
    sum_pep_quant = df_nodup[sample].sum()
    mcr_pep_quant = sum_mc_pep_quant / sum_pep_quant

    
    df_nodup_uniq = df_nodup[df_nodup['Uniquness'] == True]
    
    df_mc  = df_nodup_uniq[df_nodup_uniq['Missed.Cleavages.Count'] == 1]
    
    ratio_mc1_and_uniq_peptide = df_mc.shape[0] / df_nodup_uniq.shape[0]
    
    MC1_peptide_count = df_nodup_uniq[df_nodup_uniq['Missed.Cleavages.Count'] == 1].shape[0]
    MC2_peptide_count = df_nodup_uniq[df_nodup_uniq['Missed.Cleavages.Count'] == 2].shape[0]
    
    
    df_mc2 = calc_quant_for_fragment_pep(df_mc, pep_map, pep_quant_map, sample)
    
    
    
    mc100_pep_count = df_mc2[df_mc2["Missed.Cleavage.Ratio"]==100].shape[0]/df_mc2.shape[0]
    
    
    base_name = re.sub('.split.tsv','',os.path.basename(path))
    
    output_path = os.path.join(output_dir,base_name + "_qc.tsv")
    
    out_mc2_path = os.path.join(output_dir,base_name + "_mc2.tsv")
    df_mc2.to_csv(out_mc2_path,sep="\t",index=False)
    print(f"Results have been written to {out_mc2_path}")
    
    results = {
        'sample_name': base_name,
        'peptide_count': peptide_count,
        'protein_count': protein_count,
        'protein_human_count': protein_human_count,
        'protein_mouse_count': protein_mouse_count,
        'mc_pep_count': mc_pep_count,
        'mcr_pep': mcr_pep,
        'ratio_peptide_uniquness': ratio_peptide_uniquness,
        'ratio_multiproteins_in_group': ratio_multiproteins_in_group,
        'ratio_uniquness_mc1_peptide': ratio_uniquness_mc1_peptide,
        'sum_mc_pep_quant': sum_mc_pep_quant,
        'sum_pep_quant': sum_pep_quant,
        'mcr_pep_quant': mcr_pep_quant,
        'ratio_mc1_and_uniq_peptide': ratio_mc1_and_uniq_peptide,
        'mc1_peptide_count': MC1_peptide_count,
        'mc2_peptide_count': MC2_peptide_count,
        'mc100_pep_count': mc100_pep_count
    }

    rdf = pd.DataFrame([results])
    rdf.to_csv(output_path,sep="\t",index=False)

    print(f"Results have been written to {output_path}")

    # peptide_count, protein_count, protein_human_count, protein_mouse_count
    