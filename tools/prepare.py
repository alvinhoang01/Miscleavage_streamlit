from pyteomics import fasta, parser
from tqdm import tqdm
import os
import sqlite3
import pandas as pd
import tempfile
import shutil

def get_peptides(param):
    """
    Function to digest a FASTA file into peptides and store the results in an SQLite database.
    Uses a temporary directory for storage.
    """

    fasta_path = param['fasta_path']
    enzyme = param['enzyme']
    missed_cleavages = int(param['missed_cleavage'])
    min_length = int(param['min_length'])
    max_length = int(param['max_length'])
    m_cleavege = bool(param['m_cleavage'])

    if enzyme == "trypsin/p":
        enzyme = r'[KR]'

    # ✅ Create a Temporary Directory
    temp_dir = tempfile.mkdtemp()
    sqlite_path = os.path.join(temp_dir, "peptides.sqlite")
    print(f"📂 Using temporary directory: {temp_dir}")

    # ✅ Connect to SQLite and create table
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS peptides;")
    cursor.execute("CREATE TABLE peptides (peptide TEXT, protein TEXT);")
    conn.commit()

    print("🔍 Reading FASTA file and processing proteins...")

    batch_data = []  # ✅ Store small batches before inserting into SQLite
    batch_size = 5000

    # ✅ Process FASTA file
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
                    pre_aa = sequence[start - 1] if start > 0 else "_"
                    post_aa = sequence[start + len(peptide)] if start + len(peptide) < len(sequence) else "_"
                    protein_info = f"{protein}:{start}:{pre_aa}:{post_aa}"

                    batch_data.append((peptide, protein_info))

                # ✅ Insert batch into SQLite to free memory
                if len(batch_data) >= batch_size:
                    cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
                    conn.commit()
                    batch_data = []  # Clear batch

    print(f"✅ Total unique peptides stored BEFORE m_cleavege: {len(batch_data)}")

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
                        pre_aa = sequence[start] if start > 0 else "_"
                        post_aa = sequence[start + 1 + len(peptide)] if start + 1 + len(peptide) < len(sequence) else "_"
                        protein_info = f"{protein}:{start+1}:{pre_aa}:{post_aa}"

                        batch_data.append((peptide, protein_info))

                    # ✅ Insert batch into SQLite to free memory
                    if len(batch_data) >= batch_size:
                        cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
                        conn.commit()
                        batch_data = []  # Clear batch

    print(f"✅ Total unique peptides stored AFTER m_cleavege: {len(batch_data)}")

    # ✅ Final batch insert (in case there are leftovers)
    if batch_data:
        cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
        conn.commit()

    conn.close()

    # ✅ Print Summary
    print(f"✅ Peptides written to `{sqlite_path}`.")
    print(f"📁 SQLite file size: {os.path.getsize(sqlite_path) / 1024:.2f} KB")

    return sqlite_path, temp_dir  # ✅ Return SQLite file path & temp directory
