import os
import re
import pandas as pd
import tempfile
import shutil

def split_dia(param):
    """
    Function to split DIA search results into separate sample files.
    :param param: Dictionary of parameters loaded from YAML.
    :return: Path to the zip file containing all split files.
    """

    path = param['input_file']

    # ✅ Check if input file exists
    if not os.path.exists(path):
        print(f"❌ Error: Input file not found at {path}")
        return None, None

    print(f"🔍 Reading DIA results from TSV file: {path}")
    
    try:
        data = pd.read_csv(path, sep="\t")
    except Exception as e:
        print(f"❌ Error reading input file: {e}")
        return None, None

    # ✅ Use Temporary Directory for Cloud Processing
    temp_dir = tempfile.mkdtemp()
    output_dir = os.path.join(temp_dir, "step1-split")
    os.makedirs(output_dir, exist_ok=True)
    print(f"📂 Output directory created: {output_dir}")

    # ✅ Extract and Split Data
    columns = data.columns.values
    samples = [i for i in columns if re.search(r".PEP.Quantity", i)]

    # ✅ Debugging: Print sample columns detected
    print(f"🧪 Detected sample columns: {samples}")

    if not samples:
        print("⚠ No valid sample columns found! Check input file format.")
        return None, None

    headers = ["PG.ProteinNames", "PEP.StrippedSequence"]
    num_files = 0

    for sample in samples:
        try:
            df = data.loc[:, headers + [sample]]
            base_name = sample.split(".")[0].split(" ")[1]
            output_path = os.path.join(output_dir, base_name + ".split.tsv")
            df.to_csv(output_path, index=False, sep="\t")
            print(f"✔ File saved: {output_path}")
            num_files += 1
        except Exception as e:
            print(f"❌ Error processing {sample}: {e}")

    if num_files == 0:
        print("⚠ No split files were created! Exiting...")
        return None, None

    # ✅ Zip the folder
    zip_path = os.path.join(temp_dir, "split_results.zip")
    shutil.make_archive(zip_path.replace(".zip", ""), 'zip', output_dir)

    print(f"📁 Output files zipped at: {zip_path}")

    return zip_path, temp_dir  # ✅ Return zip path & temp folder
