import streamlit as st
import os

def main():
    st.markdown(
        """
        <style>
            .title {
                font-family: "Arial", sans-serif;
                color: #008080;
                font-size: 42px;
                font-weight: bold;
            }
            .subtitle {
                font-family: "Arial", sans-serif;
                color: #800000;
                font-size: 22px;
            }
            .main-header {
                font-family: "Arial", sans-serif; 
                font-size: 28px;
                font-weight: bold;
                text-decoration: underline;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("<div class='title'>Missed Cleavage Parser</div>", unsafe_allow_html=True)
    st.write("\n")

    # Display a static image (Optional)
    data_path = os.path.join('data', "mc_parser_logo.png")  # Replace with your logo file
    if os.path.exists(data_path):
        st.image(data_path, use_column_width=True)

    st.write("""
    This tool calculate the miscleavage rate of mass spectrometry-based DIA search results by executing the following tasks:
    - **Prepare Task:** Extracting peptides
    - **Split Task:** Splitting DIA search results
    - **QC Task:** Performing quality control checks
    - **Compare Task:** Comparing search results

    Upload your **YAML parameter file** and run individual or all tasks using the **Run Tasks** page.
    """)

    # Example Exit Button
    left_column2, right_column2 = st.columns([18, 2])
    with right_column2:
        st.write("\n")
        st.write("\n")
        if st.button('Exit'):
            st.error("The app has been stopped. Please close the browser window manually.")
            os._exit(0)
