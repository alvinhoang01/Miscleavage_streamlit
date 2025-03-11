import streamlit as st
import yaml
import os
import io
import tempfile
import shutil
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

# ✅ Initialize session state for tracking files
if "uploaded_files" not in st.session_state:
    st.session_state.uploaded_files = []
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = tempfile.mkdtemp()

# ✅ Function to read YAML parameter file
def load_yaml(uploaded_file):
    return yaml.load(uploaded_file, Loader=yaml.FullLoader)

# ✅ Function to allow users to download the default YAML file
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

# ✅ Function to save uploaded files to disk (avoids RAM issues)
def save_uploaded_file(uploaded_file):
    """Stream uploaded file directly to disk to prevent memory issues."""
    temp_dir = st.session_state.temp_dir  # Use shared temporary directory
    file_path = os.path.join(temp_dir, uploaded_file.name)

    with open(file_path, "wb") as f:
        shutil.copyfileobj(uploaded_file, f)  # Stream file to disk

    st.session_state.uploaded_files.append(file_path)  # Track uploaded files
    return file_path

# ✅ Function to delete tracked files
def cleanup_files():
    """Delete all temporary files after processing."""
    temp_dir = st.session_state.temp_dir
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)  # Delete entire temp folder
        st.session_state.uploaded_files = []  # Reset session tracking
        st.session_state.temp_dir = tempfile.mkdtemp()  # Create a new temp directory
        st.write("🗑 Temporary files cleaned up!")

# ✅ Streamlit UI
def main():
    st.title("Run Tasks")

    # ✅ Provide the default YAML file for download
    st.write("### Download Default YAML File")
    if os.path.exists("param/mc_parser.yml"):
        provide_yaml_download()
    else:
        st.error("⚠ Default YAML file not found!")

    # ✅ User uploads YAML file
    uploaded_yaml = st.file_uploader("Upload YAML File", type=["yaml", "yml"])
    if not uploaded_yaml:
        st.warning("⚠ Please upload a YAML file to proceed.")
        return

    # ✅ Load uploaded YAML file
    param = load_yaml(uploaded_yaml)
    st.success("✅ YAML file uploaded successfully!")

    # ✅ Display YAML parameters
    st.write("### Parameters Preview:")
    st.json(param)

    # ✅ User uploads input file (processed in chunks)
    uploaded_input_file = st.file_uploader("Upload Input File", type=["tsv", "csv", "txt"])
    if uploaded_input_file:
        input_file_path = save_uploaded_file(uploaded_input_file)
        param["input_file"] = input_file_path
        st.success(f"✔ Input file uploaded and saved at: {input_file_path}")

    # ✅ User uploads FASTA file (also saved to disk)
    uploaded_fasta = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
    if uploaded_fasta:
        fasta_path = save_uploaded_file(uploaded_fasta)
        param["fasta_path"] = fasta_path
        st.success(f"✅ FASTA file uploaded and saved at: {fasta_path}")

    # ✅ Task Execution Section
    st.write("## Run Tasks")

    col1, col2 = st.columns(2)

    with col1:
        if st.button("▶ Run Prepare Task"):
            st.write("Running Prepare task...")
            sqlite_path, temp_dir = get_peptides(param)  # ✅ Get SQLite file path
            st.success("✔ Prepare task completed!")

            # ✅ Provide a download button for the Prepare task output
            if os.path.exists(sqlite_path):
                with open(sqlite_path, "rb") as file:
                    st.download_button(
                        label="Download SQLite (Prepare Task)",
                        data=file,
                        file_name="peptides.sqlite",
                        mime="application/x-sqlite3"
                    )
                st.success("📂 SQLite database ready for download!")

        if st.button("▶ Run Split Task"):
            st.write("Running Split task...")
            
            zip_path, temp_dir = split_dia(param)  # ✅ Unpack the tuple properly

            if zip_path and os.path.exists(zip_path):
                st.success("✔ Split task completed!")

                # ✅ Provide a download button for the ZIP file
                with open(zip_path, "rb") as file:
                    st.download_button(
                        label="📥 Download Split Files (ZIP)",
                        data=file,
                        file_name="step1-split.zip",
                        mime="application/zip"
                    )
                st.success("📂 Split results ready for download!")

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

    # ✅ Run Full Pipeline
    st.write("## 🔄 Run All Tasks Sequentially")
    if st.button("🔄 Run Full Pipeline"):
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

        # ✅ **Cleanup files after processing**
        cleanup_files()

# ✅ **Delete files when the session resets**
if st.sidebar.button("🗑 Clear Session & Delete Files"):
    cleanup_files()

if __name__ == "__main__":
    main()