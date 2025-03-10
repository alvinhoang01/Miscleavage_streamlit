import streamlit as st
import yaml
import os
import io
import shutil
import tempfile
import pandas as pd
from tools.split_dia import split_dia
from tools.prepare import get_peptides
from tools.qc import qc_all
from tools.compare import compare_all, merge_qc

# Function to read the YAML parameter file
def read_param(param_path):
    with open(param_path) as file:
        param = yaml.load(file, Loader=yaml.FullLoader)
    return param

# Function wrappers for tasks
def prepare_param(param):
    return get_peptides(param)

def split_task(param):
    return split_dia(param)

def qc_task(param):
    return qc_all(param)

def compare_task(param):
    compare_all(param)
    merge_qc(param)

# Default YAML file path (users can download and modify this)
DEFAULT_YAML_PATH = "param/mc_parser.yml"

# Function to read the YAML parameter file
def load_yaml(uploaded_file):
    return yaml.load(uploaded_file, Loader=yaml.FullLoader)

# Function to allow users to download the actual default YAML file
def provide_yaml_download():
    default_yaml_path = "param/mc_parser.yml"
    with open(default_yaml_path, "r") as file:
        yaml_content = file.read()
    buffer = io.BytesIO(yaml_content.encode("utf-8"))
    st.download_button(
        label="Download Default YAML",
        data=buffer,
        file_name="default_parameters.yml",
        mime="text/yaml",
    )

# Cache function to load large files efficiently
@st.cache_data
def load_large_file(file_path):
    """Load a sample of the large file to reduce memory usage."""
    return pd.read_csv(file_path, nrows=1000)  # Limit rows to avoid memory crash

# Streamlit UI
def main():
    st.title("Run Tasks")

    # Provide the actual YAML file for download
    st.write("### Download Default YAML File")
    if os.path.exists("param/mc_parser.yml"):
        provide_yaml_download()
    else:
        st.error("⚠ Default YAML file not found!")

    # User uploads YAML file
    uploaded_yaml = st.file_uploader("Upload YAML File", type=["yaml", "yml"])
    if not uploaded_yaml:
        st.warning("⚠ Please upload a YAML file to proceed.")
        return

    # Load uploaded YAML file
    param = load_yaml(uploaded_yaml)
    st.success("✅ YAML file uploaded successfully!")

    # Display YAML parameters
    st.write("### Parameters Preview:")
    st.json(param)

    # Create a temp directory to store files (ensures unique files for each user)
    temp_dir = tempfile.mkdtemp()

    # Input file upload (streaming directly to disk, no RAM usage)
    uploaded_input_file = st.file_uploader("Upload Input File (Up to 1GB)", type=["tsv", "csv", "txt"])
    if uploaded_input_file:
        input_file_path = os.path.join(temp_dir, uploaded_input_file.name)

        with open(input_file_path, "wb") as f:
            shutil.copyfileobj(uploaded_input_file, f)  # Stream file to disk

        param["input_file"] = input_file_path
        st.success(f"✔ Input file saved to: {input_file_path}")

        # Load only a small sample for preview (prevents full file loading in RAM)
        try:
            df_sample = load_large_file(input_file_path)
            st.write("📊 Preview of Uploaded File (Limited to 1000 Rows):")
            st.dataframe(df_sample)
        except Exception as e:
            st.error(f"Error loading sample: {e}")

    # FASTA file upload (streaming directly to disk)
    uploaded_fasta = st.file_uploader("Upload FASTA File (Up to 1GB)", type=["fasta", "fa"])
    if uploaded_fasta:
        fasta_file_path = os.path.join(temp_dir, uploaded_fasta.name)

        with open(fasta_file_path, "wb") as f:
            shutil.copyfileobj(uploaded_fasta, f)  # Stream file to disk

        param["fasta_path"] = fasta_file_path
        st.success(f"✔ FASTA file saved to: {fasta_file_path}")

    # Output directory selection
    output_dir = st.text_input("Output Directory (Must be an Empty Folder)", value="")
    param["output_dir"] = output_dir

    # Task Execution Section
    st.write("## Run Tasks")

    col1, col2 = st.columns(2)

    with col1:
        if st.button("▶ Run Prepare Task"):
            st.write("Running Prepare task...")
            get_peptides(param)
            st.success("✔ Prepare task completed!")

        if st.button("▶ Run Split Task"):
            st.write("Running Split task...")
            split_dia(param)
            st.success("✔ Split task completed!")

    with col2:
        if st.button("▶ Run QC Task"):
            st.write("Running QC task...")
            qc_all(param)
            st.success("✔ QC task completed!")

        if st.button("▶ Run Compare Task"):
            st.write("Running Compare task...")
            compare_all(param)
            merge_qc(param)
            st.success("✔ Compare task completed!")

    # Run Full Pipeline
    st.write("## 🔄 Run Full Pipeline")
    if st.button("🔄 Run All Tasks Sequentially"):
        st.write("🛠 Running Prepare task...")
        get_peptides(param)
        st.success("✔ Prepare completed")

        st.write("🔬 Running Split task...")
        split_dia(param)
        st.success("✔ Split completed")

        st.write("📊 Running QC task...")
        qc_all(param)
        st.success("✔ QC completed")

        st.write("📈 Running Compare task...")
        compare_all(param)
        merge_qc(param)
        st.success("✅ Full pipeline completed!")

    # Cleanup temp directory (automatically deletes all temp files)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)  # Deletes temp folder and all files inside it

if __name__ == "__main__":
    main()
