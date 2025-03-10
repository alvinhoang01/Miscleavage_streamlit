import streamlit as st
import yaml
import os
import io
import requests
import re
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

import streamlit as st
import yaml
import os
import io
import requests
import re
from tools.split_dia import split_dia
from tools.prepare import get_peptides
from tools.qc import qc_all
from tools.compare import compare_all, merge_qc

# Initialize session state for tracking files
if "uploaded_files" not in st.session_state:
    st.session_state.uploaded_files = []

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

# Function to convert Google Drive "view" links to direct download links
def convert_drive_link(url):
    match = re.search(r"https://drive\.google\.com/file/d/([^/]+)/view", url)
    if match:
        file_id = match.group(1)
        return f"https://drive.google.com/uc?id={file_id}&export=download"
    return url  # Return the same URL if not a Google Drive link

# Function to download a large input file in chunks
def download_large_file(url):
    url = convert_drive_link(url)  # Convert if it's a Google Drive link
    file_path = os.path.join(os.getcwd(), "downloaded_input_file.tsv")  # Fixed filename

    st.write(f"🔄 Downloading file from {url}...")

    try:
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            with open(file_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):  # 8KB chunks
                    if chunk:
                        f.write(chunk)  # Write chunk to disk (NOT RAM)

        st.session_state.uploaded_files.append(file_path)  # Track file
        return file_path

    except requests.exceptions.RequestException as e:
        st.error(f"❌ Error downloading file: {e}")
        return None

# Function to save the FASTA file to disk and return the path
def save_fasta_file(uploaded_fasta):
    fasta_path = os.path.join(os.getcwd(), uploaded_fasta.name)  # Save in working directory

    with open(fasta_path, "wb") as f:
        f.write(uploaded_fasta.getbuffer())  # Save file correctly

    st.session_state.uploaded_files.append(fasta_path)  # Track file
    return fasta_path

# Function to delete tracked files
def cleanup_files():
    for file_path in st.session_state.uploaded_files:
        if os.path.exists(file_path):
            os.remove(file_path)
            st.write(f"🗑 Deleted: {file_path}")
    st.session_state.uploaded_files = []  # Clear the list after cleanup

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

    # **Input file must be uploaded via Google Drive / AWS S3**
    file_url = st.text_input("Enter Input File URL (Google Drive / AWS S3)")

    input_file_path = None
    if file_url:
        input_file_path = download_large_file(file_url)
        if input_file_path:
            param["input_file"] = input_file_path
            st.success(f"✔ Input file downloaded")

    # **User uploads FASTA file (saved to disk instead of storing as string)**
    fasta_path = None
    uploaded_fasta = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
    if uploaded_fasta:
        fasta_path = save_fasta_file(uploaded_fasta)
        param["fasta_path"] = fasta_path  # Save file path instead of string content
        st.success(f"✅ FASTA file uploaded")

    # Task Execution Section
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
            zip_path, temp_dir = split_dia(param)  # ✅ Get zip path
            if zip_path:
                st.success("✔ Split task completed!")

                # ✅ Provide a download button for the ZIP file
                with open(zip_path, "rb") as file:
                    st.download_button(
                        label="📥 Download Split Files (ZIP)",
                        data=file,
                        file_name="split_results.zip",
                        mime="application/zip"
                    )
                st.success("📂 Split results are ready for download!")

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

        # ✅ **Cleanup files after processing**
        cleanup_files()

# ✅ **Delete files when the session resets**
if st.sidebar.button("🗑 Clear Session & Delete Files"):
    cleanup_files()

if __name__ == "__main__":
    main()

