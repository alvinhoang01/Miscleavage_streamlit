import os
import re
import pandas as pd
import sqlite3
import tempfile
import shutil
import zipfile

def split_dia(param):
    """
    Function to split DIA search results into separate sample files.
    Supports both TSV and SQLite inputs.
    """

    input_path = param['input_file']
    
    # ✅ Create a Temporary Directory for Outputs
    temp_dir = tempfile.mkdtemp()
    output_dir = os.path.join(temp_dir, "step1-split")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"📂 Output directory created: {output_dir}")

    # ✅ Case 1: Read from TSV File
    if input_path.endswith(".tsv") or input_path.endswith(".csv"):
        print("🔍 Reading DIA results from TSV file...")
        data = pd.read_csv(input_path, sep="\t")

    # ✅ Case 2: Read from SQLite Database
    elif input_path.endswith(".sqlite"):
        print("🔍 Reading DIA results from SQLite database...")
        conn = sqlite3.connect(input_path)
        query = "SELECT * FROM peptides;"  # Adjust query if needed
        data = pd.read_sql(query, conn)
        conn.close()
    
    else:
        print("❌ Unsupported input file format.")
        return None, None

    # ✅ Split table into multiple sample files
    columns = data.columns.values
    samples = [i for i in columns if re.search(r".PEP.Quantity", i)]
    headers = ["PG.ProteinNames", "PEP.StrippedSequence"]

    for sample in samples:
        df = data.loc[:, headers + [sample]]
        base_name = sample.split(".")[0].split(" ")[1]
        output_path = os.path.join(output_dir, base_name + ".split.tsv")
        df.to_csv(output_path, index=False, sep="\t")
        print(f"✔ {sample} saved to {output_path}")

    # ✅ Create a ZIP file for download
    zip_path = os.path.join(temp_dir, "split_results.zip")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(output_dir):
            for file in files:
                zipf.write(os.path.join(root, file), file)

    print(f"📁 Output files zipped at: {zip_path}")

    return output_dir, zip_path  # ✅ Return folder & ZIP file path
