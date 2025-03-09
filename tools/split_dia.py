import os
import re
import pandas as pd

def split_dia(param):
    """
    Function to split DIA search results into separate sample files.
    :param param: Dictionary of parameters loaded from YAML.
    :return: List of logs/messages for Streamlit UI.
    """
    path = param['input_file']
    data = pd.read_csv(path, sep="\t")
    output_dir = param['output_dir']
    output_dir = os.path.join(output_dir,"step1-split")
    if not os.path.exists(output_dir):  
        os.makedirs(output_dir)
    # split table to multiple tables (samples)
    columns = data.columns.values
    samples = [i for i in columns if re.search(r".PEP.Quantity",i)]
    headers = ["PG.ProteinNames","PEP.StrippedSequence"]
    for sample in samples:
        df = data.loc[:,headers + [sample]]
        base_name = sample.split(".")[0].split(" ")[1]
        output_path = os.path.join(output_dir,base_name + ".split.tsv")
        df.to_csv(output_path,index=False, sep="\t")
        print(f"{sample} was saved to {output_path}")

    return
