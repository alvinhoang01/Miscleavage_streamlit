import os
import re
import pandas as pd

def split_dia(param):
    """
    Function to split DIA search results into separate sample files.
    :param param: Dictionary of parameters loaded from YAML.
    :return: List of logs/messages for Streamlit UI.
    """

    logs = []  # ✅ Store messages for UI
    path = param['input_file']

    # ✅ Check if input file exists
    if not os.path.exists(path):
        error_msg = f"❌ Error: Input file not found at {path}"
        logs.append(error_msg)
        print(error_msg)
        return logs

    # ✅ Read DIA results
    try:
        data = pd.read_csv(path, sep="\t")
    except Exception as e:
        error_msg = f"❌ Error reading input file: {e}"
        logs.append(error_msg)
        print(error_msg)
        return logs

    # ✅ Ensure output directory exists
    output_dir = os.path.join(param['output_dir'], "step1-split")
    os.makedirs(output_dir, exist_ok=True)
    logs.append(f"📂 Output directory created: {output_dir}")

    # ✅ Extract and Split Data
    columns = data.columns.values
    samples = [i for i in columns if re.search(r".PEP.Quantity", i)]

    if not samples:
        warning_msg = "⚠ No valid sample columns found! Check input file format."
        logs.append(warning_msg)
        print(warning_msg)
        return logs

    headers = ["PG.ProteinNames", "PEP.StrippedSequence"]
    num_files = 0

    for sample in samples:
        try:
            df = data.loc[:, headers + [sample]]
            base_name = sample.split(".")[0].split(" ")[1]
            output_path = os.path.join(output_dir, f"{base_name}.split.tsv")
            df.to_csv(output_path, index=False, sep="\t")
            success_msg = f"✔ {sample} was saved to {output_path}"
            logs.append(success_msg)
            print(success_msg)
            num_files += 1
        except Exception as e:
            error_msg = f"❌ Error processing {sample}: {e}"
            logs.append(error_msg)
            print(error_msg)

    if num_files == 0:
        logs.append("⚠ No split files were created! Exiting...")
        print("⚠ No split files were created! Exiting...")

    return logs  # ✅ Return messages for Streamlit UI
