import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import tempfile
import os

# Page Configuration
st.set_page_config(page_title="Drug Discovery Screening Dashboard", layout="wide")

st.title("💊 Professional Virtual Screening & Drug Profiling Pipeline")
st.markdown("---")
st.markdown("*Client Quality Assurance Mode:* Rigorous Lipinski Rule of 5 validation pipeline.")

# Sidebar Thresholds
st.sidebar.header("Validated Parameters")
max_mw = st.sidebar.number_input("Max Molecular Weight (Da)", value=500.0, step=10.0)
max_logp = st.sidebar.number_input("Max LogP", value=5.0, step=0.1)
max_hbd = st.sidebar.number_input("Max Hydrogen Bond Donors", value=5, step=1)
max_hba = st.sidebar.number_input("Max Hydrogen Bond Acceptors", value=10, step=1)

# File Uploader
uploaded_file = st.file_uploader("Upload Client SDF File (.sdf)", type=["sdf"])

mols = []
parsing_errors = 0

if uploaded_file is not None:
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as tmp:
            tmp.write(uploaded_file.getvalue())
            tmp_path = tmp.name
        
        suppl = Chem.SDMolSupplier(tmp_path, sanitize=True)
        for idx, m in enumerate(suppl):
            if m is not None:
                mols.append(m)
            else:
                parsing_errors += 1
        os.unlink(tmp_path)
    except Exception as e:
        st.error(f"File Reading Error: {e}")

# Fallback mechanism if no file uploaded
if len(mols) == 0:
    if uploaded_file is not None:
        st.warning("⚠️ Uploaded file me valid molecules nahi mile. Default test molecule (Aspirin) load kiya gaya hai.")
    else:
        st.info("💡 Kripya client ki SDF file upload karein. Default benchmark molecule (Aspirin) loaded hai.")
    
    aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    if aspirin:
        mols = [aspirin]

# Processing Pipeline
results = []
passed_count = 0

for idx, mol in enumerate(mols):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    
    r_mw = mw <= max_mw
    r_logp = logp <= max_logp
    r_hbd = hbd <= max_hbd
    r_hba = hba <= max_hba
    
    passed = r_mw and r_logp and r_hbd and r_hba
    if passed:
        passed_count += 1
        
    if mol.HasProp("Name"):
        compound_name = mol.GetProp("Name")
    elif mol.HasProp("_Name"):
        compound_name = mol.GetProp("_Name")
    else:
        compound_name = f"CMPD{idx+1:03d}"
    
    results.append({
        "Compound ID": compound_name,
        "Molecular Weight (g/mol)": round(mw, 2),
        "LogP": round(logp, 2),
        "H-Donors": hbd,
        "H-Acceptors": hba,
        "Lipinski Compliance": "PASS" if passed else "FAIL"
    })

df = pd.DataFrame(results)

st.success(f"Audit Complete | Total Compounds: {len(mols)} | Passed Screening: {passed_count} | Errors: {parsing_errors}")
st.dataframe(df, use_container_width=True)

# Download CSV Report (Directly opens in Excel, 100% error-free)
if not df.empty:
    csv_data = df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Client Screening Report (CSV for Excel)",
        data=csv_data,
        file_name="Client_Drug_Screening_Report.csv",
        mime="text/csv"
    )