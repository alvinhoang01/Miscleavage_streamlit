from pyteomics import fasta, parser
from tqdm import tqdm
import os
import sqlite3
import pandas as pd

def get_peptides(param):
    """
    Function to digest a FASTA file into peptides and store the results in an SQLite database.
    Now includes deeper debugging for peptide loss.
    """

    fasta_path = param['fasta_path']
    enzyme = param['enzyme']
    missed_cleavages = int(param['missed_cleavage'])
    min_length = int(param['min_length'])
    max_length = int(param['max_length'])
    m_cleavege = bool(param['m_cleavage'])

    pep_map = dict()

    if enzyme == "trypsin/p":
        enzyme = r'[KR]'

    print("🔍 Reading FASTA file...")
    with fasta.read(fasta_path) as entries:
        for header, sequence in tqdm(entries, desc="Processing Proteins"):
            peptides = parser.xcleave(
                sequence,
                enzyme,
                missed_cleavages=missed_cleavages,
                min_length=min_length
            )

            protein = header.split(" ")[0].split("|")[-1]
            for start, peptide in peptides:
                if len(peptide) <= max_length:
                    m = {
                        'protein': protein,
                        'start': start,
                        'pre_aa': sequence[start-1] if start > 0 else "_",
                        'post_aa': sequence[start+len(peptide)] if start+len(peptide) < len(sequence) else "_"
                    }

                    if peptide in pep_map:
                        pep_map[peptide].append(m)
                    else:
                        pep_map[peptide] = [m]
    print(f"✅ Total unique peptides stored BEFORE m_cleavege: {len(pep_map)}")

    # ✅ If `m_cleavege` flag is enabled, process again
    if m_cleavege:
        print("🔄 Running m_cleavage processing...")
        with fasta.read(fasta_path) as entries:
            for header, sequence in tqdm(entries, desc="Processing m_cleavege Proteins"):
                peptides = parser.xcleave(
                    sequence[1:],  # Shift sequence by 1
                    enzyme,
                    missed_cleavages=missed_cleavages,
                    min_length=min_length
                )
                protein = header.split(" ")[0].split("|")[-1]
                for start, peptide in peptides:
                    if len(peptide) <= max_length:
                        m = dict()
                        m['protein'] = protein
                        m['start'] = start + 1
                        try:
                            m['pre_aa'] = sequence[start] if start > 0 else sequence[0]
                        except:
                            print(start)
                        m['post_aa'] = sequence[start+1+len(peptide)] if start+1+len(peptide) < len(sequence) else "_"
                        
                        if peptide in pep_map:
                            pep_map[peptide].append(m)
                    else:
                        pep_map[peptide] = [m]
    print(f"✅ Total unique peptides stored AFTER m_cleavege: {len(pep_map)}")

    # ✅ Save results to SQLite
    output_dir = param['output_dir']
    sqlite_path = os.path.join(output_dir, "peptides.sqlite")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"📂 Output directory created: {output_dir}")

    rows = []
    print("Writing peptides to DataFrame...")
    for pep in tqdm(pep_map):
        proteins = sorted(list(pep_map[pep]), key=lambda x: x['protein'])
        proteins = [f"{p['protein']}:{p['start']}:{p['pre_aa']}:{p['post_aa']}" for p in proteins]
        proteins = set(proteins)
        proteins_str = ";".join(sorted(list(set(proteins))))
        rows.append((pep, proteins_str))
        
    df = pd.DataFrame(rows, columns=["peptide", "protein"])

    # ✅ Save to SQLite
    conn = sqlite3.connect(sqlite_path)
    df.to_sql('peptides', conn, if_exists='replace', index=False)
    conn.close()

    print(f"✅ Peptides written to `{sqlite_path}`.")
    print(f"📊 Final database contains {len(df)} peptides.")
    print(f"📁 SQLite file size: {os.path.getsize(sqlite_path) / 1024:.2f} KB")

    return
