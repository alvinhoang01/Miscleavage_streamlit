import streamlit as st
import yaml
import os
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

# Streamlit UI
def main():
    st.title("Run Tasks")

    # File uploader for parameter file
    uploaded_file = st.file_uploader("Upload Parameter File (YAML)", type=["yaml", "yml"])

    if uploaded_file:
        # Save the uploaded file temporarily
        param_path = "temp_parameters.yaml"
        with open(param_path, "wb") as f:
            f.write(uploaded_file.getbuffer())

        # Read parameters
        param = read_param(param_path)
        st.success("✅ Parameter file loaded successfully!")

        # Display parameters
        st.write("### Parameters Preview:")
        st.json(param)

        # Run Individual Tasks
        st.write("## Run Tasks Individually")

        col1, col2 = st.columns(2)

        with col1:
            if st.button("Run Prepare Task"):  
                st.write("Running Prepare task...")
                prepare_param(param)
                st.success("✔ Prepare task completed!")

            if st.button("Run Split Task"):
                st.write("Running Split task...")
                split_task(param)
                st.success("✔ Split task completed!")

        with col2:
            if st.button("Run QC Task"):
                st.write("Running QC task...")
                qc_task(param)
                st.success("✔ QC task completed!")

            if st.button("Run Compare Task"):
                st.write("Running Compare task...")
                compare_task(param)
                st.success("✔ Compare task completed!")

        # Run Full Pipeline
        st.write("## Run Full Pipeline")

        if st.button("Run All Tasks Sequentially"):
            st.write("🔄 Running full pipeline...")

            st.write("🛠 Running Prepare task...")
            prepare_param(param)
            st.success("✔ Prepare completed")

            st.write("🔬 Running Split task...")
            split_task(param)
            st.success("✔ Split completed")

            st.write("📊 Running QC task...")
            qc_task(param)
            st.success("✔ QC completed")

            st.write("📈 Running Compare task...")
            compare_task(param)
            st.success("✅ Full pipeline completed!")


if __name__ == "__main__":
    main()
