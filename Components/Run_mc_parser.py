import streamlit as st
import yaml
import os
import io
import pandas as pd
import tempfile
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

# Streamlit UI
def main():
    st.title("Run Tasks")

    # Provide the actual YAML file for download
    st.write("### Download Default YAML File")
    if os.path.exists("param/mc_parser.yml"):
        provide_yaml_download()
    else:
        st.error("⚠ Default YAML file not found!")

    # User uploads their YAML file
    uploaded_yaml = st.file_uploader("Upload YAML File", type=["yaml", "yml"])
    if not uploaded_yaml:
        st.warning("⚠ Please upload a YAML file to proceed.")
        return

    # Load uploaded YAML file (from memory)
    param = load_yaml(uploaded_yaml)
    st.success("✅ YAML file uploaded successfully!")

    # Display YAML parameters
    st.write("### Parameters Preview:")
    st.json(param)

    # Input file upload (store on disk, not in RAM)
    uploaded_input_file = st.file_uploader("Upload Input File (Up to 1GB)", type=["tsv", "csv", "txt"])
    if uploaded_input_file:
        temp_dir = tempfile.gettempdir()  # Get system temp directory
        temp_input_path = os.path.join(temp_dir, uploaded_input_file.name)

        # Write file to disk in chunks (NO RAM overload)
        with open(temp_input_path, "wb") as f:
            for chunk in uploaded_input_file.chunks(1024 * 1024):  # 1MB chunks
                f.write(chunk)

        param["input_file"] = temp_input_path
        st.success(f"✔ Input file saved to: {temp_input_path}")

    # FASTA file upload (store on disk)
    uploaded_fasta = st.file_uploader("Upload FASTA File (Up to 1GB)", type=["fasta", "fa"])
    if uploaded_fasta:
        temp_fasta_path = os.path.join(temp_dir, uploaded_fasta.name)

        # Write file to disk in chunks
        with open(temp_fasta_path, "wb") as f:
            for chunk in uploaded_fasta.chunks(1024 * 1024):  # 1MB chunks
                f.write(chunk)

        param["fasta_path"] = temp_fasta_path
        st.success(f"✔ FASTA file saved to: {temp_fasta_path}")

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

    # Cleanup temp files (after running tasks)
    if "input_file" in param and os.path.exists(param["input_file"]):
        os.remove(param["input_file"])

    if "fasta_path" in param and os.path.exists(param["fasta_path"]):
        os.remove(param["fasta_path"])

if __name__ == "__main__":
    main()
