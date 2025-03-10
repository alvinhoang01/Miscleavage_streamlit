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
    Now includes temporary folder handling for cloud processing.
    Uses batch writing to prevent memory overload.
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
    batch_size = 5000  # ✅ Prevent memory overflow

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
                    m = {
                        'protein': protein,
                        'start': start,
                        'pre_aa': sequence[start-1] if start > 0 else "_",
                        'post_aa': sequence[start+len(peptide)] if start+len(peptide) < len(sequence) else "_"
                    }

                    batch_data.append((peptide, f"{m['protein']}:{m['start']}:{m['pre_aa']}:{m['post_aa']}"))

                # ✅ Insert batch into SQLite to free memory
                if len(batch_data) >= batch_size:
                    cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
                    conn.commit()
                    batch_data = []  # ✅ Clear batch

    print(f"✅ Peptides stored BEFORE m_cleavege: {len(batch_data)} (final batch size)")

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
                        m = {
                            'protein': protein,
                            'start': start + 1,
                            'pre_aa': sequence[start] if start > 0 else "_",
                            'post_aa': sequence[start+1+len(peptide)] if start+1+len(peptide) < len(sequence) else "_"
                        }

                        batch_data.append((peptide, f"{m['protein']}:{m['start']}:{m['pre_aa']}:{m['post_aa']}"))

                    # ✅ Insert batch into SQLite to free memory
                    if len(batch_data) >= batch_size:
                        cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
                        conn.commit()
                        batch_data = []  # ✅ Clear batch

    print(f"✅ Peptides stored AFTER m_cleavege: {len(batch_data)} (final batch size)")

    # ✅ Final batch insert (in case there are leftovers)
    if batch_data:
        cursor.executemany("INSERT INTO peptides VALUES (?, ?);", batch_data)
        conn.commit()

    conn.close()

    print(f"✅ Peptides written to `{sqlite_path}`.")
    print(f"📊 Final database contains peptides.")
    print(f"📁 SQLite file size: {os.path.getsize(sqlite_path) / 1024:.2f} KB")

    return sqlite_path, temp_dir  # ✅ Return SQLite path & temp dir
