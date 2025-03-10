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
    data = pd.read_csv(path, sep="\t")

    # ✅ Use Temporary Directory for Cloud Processing
    temp_dir = tempfile.mkdtemp()
    output_dir = os.path.join(temp_dir, "step1-split")
    os.makedirs(output_dir, exist_ok=True)

    # ✅ Extract and Split Data
    columns = data.columns.values
    samples = [i for i in columns if re.search(r".PEP.Quantity", i)]
    headers = ["PG.ProteinNames", "PEP.StrippedSequence"]

    for sample in samples:
        df = data.loc[:, headers + [sample]]
        base_name = sample.split(".")[0].split(" ")[1]
        output_path = os.path.join(output_dir, base_name + ".split.tsv")
        df.to_csv(output_path, index=False, sep="\t")
        print(f"{sample} was saved to {output_path}")

    # ✅ Zip the folder
    zip_path = os.path.join(temp_dir, "split_results.zip")
    shutil.make_archive(zip_path.replace(".zip", ""), 'zip', output_dir)

    print(f"📁 Split files compressed into: {zip_path}")
    
    return zip_path, temp_dir  # ✅ Return zip path & temp folder
