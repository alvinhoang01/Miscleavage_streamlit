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
# ✅ Ensure temp_dir is always initialized at the start of the session
if "temp_dir" not in st.session_state or not os.path.exists(st.session_state["temp_dir"]):
    st.session_state["temp_dir"] = tempfile.mkdtemp()


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

# ✅ Function to save uploaded files to disk
def save_uploaded_file(uploaded_file):
    """Stream uploaded file directly to disk to prevent memory issues."""
    temp_dir = st.session_state.temp_dir
    file_path = os.path.join(temp_dir, uploaded_file.name)
    
    with open(file_path, "wb") as f:
        shutil.copyfileobj(uploaded_file, f)  # Stream file to disk
    
    return file_path

# ✅ Function to run the full pipeline
def run_full_pipeline(param):
    """Runs the entire pipeline and stores results in `output_folder`."""
    output_dir = os.path.join(st.session_state.temp_dir, "output_folder")
    os.makedirs(output_dir, exist_ok=True)
    
    # Update output directory in param
    param["output_dir"] = output_dir  

    # 🔹 Step 1: Prepare
    st.write("🛠 Running Prepare step...")
    get_peptides(param)
    st.success("✔ Prepare step completed!")

    # 🔹 Step 2: Split
    st.write("🔬 Running Split step...")
    split_dia(param)
    st.success("✔ Split step completed!")

    # 🔹 Step 3: QC
    st.write("📊 Running QC step...")
    qc_all(param)
    st.success("✔ QC step completed!")

    # 🔹 Step 4: Compare
    st.write("📈 Running Compare step...")
    compare_all(param)
    merge_qc(param)
    st.success("✅ Full pipeline completed!")

    # ✅ Zip the output folder
    zip_output_path = os.path.join(st.session_state.temp_dir, "pipeline_results.zip")
    shutil.make_archive(zip_output_path.replace(".zip", ""), 'zip', output_dir)
    
    return zip_output_path

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

    # ✅ User uploads input file
    uploaded_input_file = st.file_uploader("Upload Input File", type=["tsv", "csv", "txt"])
    if uploaded_input_file:
        input_file_path = save_uploaded_file(uploaded_input_file)
        param["input_file"] = input_file_path
        st.success(f"✅ Input file uploaded")

    # ✅ User uploads FASTA file
    uploaded_fasta = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
    if uploaded_fasta:
        fasta_path = save_uploaded_file(uploaded_fasta)
        param["fasta_path"] = fasta_path
        st.success(f"✅ FASTA file uploaded")

    # ✅ Task Execution Section
    st.write("## Calculate Miscleavage Rate")

    if st.button("Run Misclevage Parser"):
        zip_output_path = run_full_pipeline(param)

        # ✅ Provide download button for results
        if os.path.exists(zip_output_path):
            with open(zip_output_path, "rb") as file:
                st.download_button(
                    label="📥 Download Full Pipeline Results",
                    data=file,
                    file_name="MC_parser_output.zip",
                    mime="application/zip"
                )
            
            # ✅ Clean up the temp folder safely and reinitialize it
            shutil.rmtree(st.session_state.temp_dir, ignore_errors=True)
            st.session_state["temp_dir"] = tempfile.mkdtemp()  # Reinitialize temp directory



if __name__ == "__main__":
    main()