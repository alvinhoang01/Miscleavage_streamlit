from pyteomics import fasta, parser
from tqdm import tqdm
import os
import sqlite3

def get_peptides(param):
    """
    Function to digest a FASTA file into peptides and store the results in an SQLite database.
    Now optimized to prevent memory crashes by writing peptides directly to SQLite in chunks.
    """

    fasta_path = param['fasta_path']
    enzyme = param['enzyme']
    missed_cleavages = int(param['missed_cleavage'])
    min_length = int(param['min_length'])
    max_length = int(param['max_length'])
    m_cleavege = bool(param['m_cleavage'])

    # ✅ Define output SQLite path
    output_dir = param['output_dir']
    sqlite_path = os.path.join(output_dir, "peptides.sqlite")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"📂 Output directory created: {output_dir}")

    # ✅ Connect to SQLite and create table
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS peptides;")  # Remove existing table
    cursor.execute("CREATE TABLE peptides (peptide TEXT, protein TEXT);")
    conn.commit()

    # ✅ Process FASTA File
    print("🔍 Reading FASTA file and processing proteins...")

    batch_data = []  # Store small batches of peptides before inserting into SQLite
    batch_size = 5000  # ✅ Adjust this to control memory usage

    with fasta.read(fasta_path) as entries:
        for header, sequence in tqdm(entries, desc="Processing Proteins"):
            peptides = parser.xcleave(sequence, enzyme, missed_cleavages=missed_cleavages, min_length=min_length)
            protein = header.split(" ")[0].split("|")[-1]
            
            for start, peptide in peptides:
                if len(peptide) <= max_length:
                    pre_aa = sequence[start - 1] if start > 0 else "_"
                    post_aa = sequence[start + len(peptide)] if start + len(peptide) < len(sequence) else "_"
                    protein_info = f"{protein}:{start}:{pre_aa}:{post_aa}"
                    
                    batch_data.append((peptide, protein_info))

                # ✅ Insert batch into SQLite to free memory
                if len(batch_data) >= batch_size:
                    cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
                    conn.commit()
                    batch_data = []  # Clear batch

    # ✅ Process `m_cleavege` if enabled
    if m_cleavege:
        print("🔄 Running m_cleavage processing...")

        with fasta.read(fasta_path) as entries:
            for header, sequence in tqdm(entries, desc="Processing m_cleavege Proteins"):
                peptides = parser.xcleave(sequence[1:], enzyme, missed_cleavages=missed_cleavages, min_length=min_length)
                protein = header.split(" ")[0].split("|")[-1]

                for start, peptide in peptides:
                    if len(peptide) <= max_length:
                        pre_aa = sequence[start] if start > 0 else sequence[0]
                        post_aa = sequence[start + 1 + len(peptide)] if start + 1 + len(peptide) < len(sequence) else "_"
                        protein_info = f"{protein}:{start+1}:{pre_aa}:{post_aa}"

                        batch_data.append((peptide, protein_info))

                    # ✅ Insert batch into SQLite to free memory
                    if len(batch_data) >= batch_size:
                        cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
                        conn.commit()
                        batch_data = []  # Clear batch

    # ✅ Final batch insert (in case there are leftovers)
    if batch_data:
        cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
        conn.commit()

    conn.close()

    # ✅ Print Summary
    print(f"✅ Peptides written to `{sqlite_path}`.")
    print(f"📁 SQLite file size: {os.path.getsize(sqlite_path) / 1024:.2f} KB")

    return sqlite_path
