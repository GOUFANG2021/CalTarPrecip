import streamlit as st
import pandas as pd
import os
import subprocess
import gdown
import datetime

# ======================== DEFINE PATHS ===========================
GITHUB_REPO = "https://github.com/GOUFANG2021/CalTarPrecip/raw/main"

# File download URLs from new GitHub repository
MODEL_PY_URL = f"{GITHUB_REPO}/CaTar_Model.py"
DATA_TEMPLATE_URL = f"{GITHUB_REPO}/Wine%20Data.xlsx"
INDICATOR_IMAGE_URL = f"{GITHUB_REPO}/indicator.png"

# ======================== FUNCTION TO DOWNLOAD FILE FROM GITHUB ===========================
def download_from_github(url, output_path):
    """Download a file from GitHub repository."""
    try:
        gdown.download(url, output_path, quiet=False)
        return f"✅ The data has been downloaded successfully!"
    except Exception as e:
        return f"❌ Failed to download {os.path.basename(output_path)}: {e}"

# ======================== FUNCTION TO RUN MODEL DIRECTLY FROM GITHUB ===========================
def run_model_from_github(model_url, data_path, simulation_id):
    """Download the model from GitHub and execute it with the uploaded data file."""
    model_path = "CaTar_Model.py"

    # Download model file from GitHub
    download_result = download_from_github(model_url, model_path)

    # Explicitly use the correct Python environment
    python_executable = os.environ.get("VIRTUAL_ENV", "/home/adminuser/venv/bin/python3")

    # Run the model
    try:
        process = subprocess.Popen(
            [python_executable, model_path, data_path], 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        output, error = process.communicate()
        
        if error:
            return f"❌ Model execution failed: {error}"
        
        return f"✅ {simulation_id} completed successfully!\n\n{output}"
    except Exception as e:
        return f"❌ Error running model: {e}"

# ======================== STREAMLIT UI ===========================
st.set_page_config(layout="wide")  
st.title("🍷 Calcium Tartrate Precipitation Predictor")

# Create two columns
col1, col2 = st.columns([1, 1])

# Ensure session state variables exist
if "simulation_results" not in st.session_state:
    st.session_state.simulation_results = []  # Store results as a list
if "uploaded_data" not in st.session_state:
    st.session_state.uploaded_data = None
if "simulation_count" not in st.session_state:
    st.session_state.simulation_count = 0  # Count the number of simulations

with col1:
    # STEP 1: DOWNLOAD TEMPLATE
    st.subheader("Step 1: Download the data format and enter your wine information in the input sheet.")
    
    template_path = "Wine Data.xlsx"
    download_result_template = download_from_github(DATA_TEMPLATE_URL, template_path)

    # Provide download button
    if os.path.exists(template_path):
        with open(template_path, "rb") as f:
            st.download_button("📥 Download Data Format", f, file_name="Wine Data.xlsx")

    # STEP 2: UPLOAD MODIFIED WINE DATA
    st.subheader("Step 2: Upload Your Modified Wine Data (Excel)")
    uploaded_file = st.file_uploader("📤 Browse files to upload Your Wine Data (Excel)", type=["xlsx"])

    if uploaded_file:
        st.session_state.uploaded_data = uploaded_file  
        st.success(f"✅ Uploaded: {uploaded_file.name}")
     #   st.info("🔄 Please delete the current data to upload the next data file.")

    # STEP 3: RUN MODEL
    st.subheader("Step 3: Run Model")    
    if st.button("🚀 Run Model"):
        if st.session_state.uploaded_data is None:
            st.error("⚠️ Please upload a wine data file before running the model.")
        else:
            # Show processing message
            st.info("⏳ The simulation can take a few minutes. It will always stop by itself when finished.")

            # Increment simulation count
            st.session_state.simulation_count += 1

            # Generate a unique simulation ID
            simulation_id = f"Simulation {st.session_state.simulation_count} for {st.session_state.uploaded_data.name}"

            # Save uploaded file temporarily
            uploaded_file_path = "Wine Data.xlsx"
            with open(uploaded_file_path, "wb") as f:
                f.write(st.session_state.uploaded_data.getbuffer())

            # Run model from GitHub
            results = run_model_from_github(MODEL_PY_URL, uploaded_file_path, simulation_id)

            # Store results
            st.session_state.simulation_results.append(f"### {simulation_id}\n{results}\n")
            st.success(f"✅ {simulation_id} completed! Check results on the right.")

with col2:
    # DISPLAY RESULTS FOR ALL SIMULATIONS
    st.subheader("📊 Simulation Results")
    if st.session_state.simulation_results:
        for result in st.session_state.simulation_results:
            st.write(result)
             
    # INTERPRETATION SECTION
    st.subheader("📌 Interpretation")
    st.write("""
    It is recommended that wines with a supersaturation ratio in the high-risk range should be treated to prevent calcium tartrate formation. 
    It is possible for medium-risk wines to form calcium tartrate, but most wines in this range will not require treatment.
    """)

    # DISPLAY IMAGE FROM GITHUB
    indicator_path = "indicator.png"
    gdown.download(INDICATOR_IMAGE_URL, indicator_path, quiet=False)
    if os.path.exists(indicator_path):
        st.image(indicator_path, caption=" ")

    # WARNING MESSAGE
    st.warning(
        "⚠️ The model may not find a solution if the input data falls outside the simulation range. "
        "If this occurs, please delete the uploaded Excel file and upload a new one with the same format."
    )
