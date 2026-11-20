import streamlit as st
import pandas as pd
import io
import tempfile
import os

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem, rdFingerprintGenerator, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

import py3Dmol
from stmol import showmol


# =====================================================
# 1. PAGE CONFIGURATION
# =====================================================

st.set_page_config(
    page_title="Nexus AI Cheminformatics Pro",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 Nexus Advanced Cheminformatics & Virtual Screening")

st.markdown(
    "Molecular property analysis | Lipinski screening | PAINS | "
    "Tanimoto similarity | 2D/3D visualization | "
    "Chemistry Assistant | Excel reports"
)


# =====================================================
# 2. FILE UPLOAD
# =====================================================

st.sidebar.header("📂 Molecular Dataset")

uploaded_file = st.sidebar.file_uploader(
    "Upload CSV, TXT, SMI or SDF",
    type=["csv", "txt", "smi", "sdf"]
)


def load_dataset(file):

    if file is None:
        return pd.DataFrame({
            "Molecule_Name": [
                "Aspirin",
                "Paracetamol",
                "Benzene",
                "ZNZ922-001"
            ],
            "SMILES": [
                "CC(=O)Oc1ccccc1C(=O)O",
                "CC(=O)Nc1ccc(O)cc1",
                "c1ccccc1",
                "FC1(F)CN(CC1)C(=O)c2ccc(cc2)C3CCN(CC3)c4ccc5nc(=O)n(C6CCC(=O)NC6=O)c5c4"
            ]
        })

    ext = file.name.rsplit(".", 1)[-1].lower()

    if ext == "sdf":

        temp_path = None

        try:
            with tempfile.NamedTemporaryFile(
                delete=False, suffix=".sdf"
            ) as temp:

                temp.write(file.getvalue())
                temp_path = temp.name

            supplier = Chem.SDMolSupplier(temp_path)

            records = []

            for i, mol in enumerate(supplier):

                if mol is not None:

                    name = (
                        mol.GetProp("_Name")
                        if mol.HasProp("_Name")
                        else f"Molecule_{i+1}"
                    )

                    records.append({
                        "Molecule_Name": name,
                        "SMILES": Chem.MolToSmiles(mol)
                    })

            return pd.DataFrame(records)

        finally:

            if temp_path and os.path.exists(temp_path):
                os.unlink(temp_path)

    if ext == "csv":

        data = pd.read_csv(file)

    else:

        text = file.getvalue().decode(
            "utf-8-sig", errors="replace"
        )

        records = []

        for line in text.splitlines():

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) >= 2:
                records.append({
                    "SMILES": parts[0],
                    "Molecule_Name": " ".join(parts[1:])
                })

            else:
                records.append({
                    "SMILES": parts[0],
                    "Molecule_Name": f"Molecule_{len(records)+1}"
                })

        data = pd.DataFrame(records)

    if data.empty:
        return data

    data.columns = [str(c).strip() for c in data.columns]

    smiles_column = next(
        (c for c in data.columns
         if c.lower() in ["smiles", "canonical_smiles", "smi"]),
        data.columns[0]
    )

    name_column = next(
        (c for c in data.columns
         if c.lower() in ["molecule_name", "name", "compound", "id"]),
        None
    )

    if name_column is None:

        data["Molecule_Name"] = [
            f"Molecule_{i+1}" for i in range(len(data))
        ]

        name_column = "Molecule_Name"

    data = data.rename(columns={
        smiles_column: "SMILES",
        name_column: "Molecule_Name"
    })

    return data[["Molecule_Name", "SMILES"]].copy()


try:

    df = load_dataset(uploaded_file)

except Exception as e:

    st.error(f"Dataset loading error: {e}")
    st.stop()


if df.empty:

    st.error("No molecular records found.")
    st.stop()


# =====================================================
# 3. VALIDATE MOLECULES
# =====================================================

df["SMILES"] = df["SMILES"].fillna("").astype(str)

df["Molecule_Name"] = (
    df["Molecule_Name"].fillna("Unnamed").astype(str)
)

molecules = []

for smi in df["SMILES"]:

    mol = Chem.MolFromSmiles(smi.strip())
    molecules.append(mol)


df["Structure_Status"] = [
    "Valid" if mol is not None else "Invalid SMILES"
    for mol in molecules
]


# =====================================================
# 4. MORGAN FINGERPRINT GENERATOR
# =====================================================

fp_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)

valid_indices = [
    i for i, mol in enumerate(molecules)
    if mol is not None
]

if not valid_indices:

    st.error("No valid SMILES structures in this dataset.")
    st.stop()


# =====================================================
# 5. PAINS FILTER
# =====================================================

params = FilterCatalogParams()

params.AddCatalog(
    FilterCatalogParams.FilterCatalogs.PAINS
)

pains_catalog = FilterCatalog(params)


# =====================================================
# 6. MOLECULAR PROPERTY ANALYSIS
# =====================================================

results = []
fingerprints = {}

for i, mol in enumerate(molecules):

    if mol is None:

        results.append({
            "Molecular_Weight": None,
            "LogP": None,
            "H_Donors": None,
            "H_Acceptors": None,
            "TPSA": None,
            "Rotatable_Bonds": None,
            "Aromatic_Rings": None,
            "Lipinski_Violations": None,
            "Lipinski_Result": "Invalid SMILES",
            "PAINS_Alert": "N/A",
            "PAINS_Details": "N/A"
        })

        continue

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Descriptors.NumRotatableBonds(mol)
    aromatic = Descriptors.NumAromaticRings(mol)

    violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])

    pains_matches = pains_catalog.GetMatches(mol)

    pains_names = [
        match.GetDescription()
        for match in pains_matches
    ]

    fp = fp_generator.GetFingerprint(mol)
    fingerprints[i] = fp

    results.append({

        "Molecular_Weight": round(mw, 2),
        "LogP": round(logp, 2),
        "H_Donors": hbd,
        "H_Acceptors": hba,
        "TPSA": round(tpsa, 2),
        "Rotatable_Bonds": rotatable,
        "Aromatic_Rings": aromatic,
        "Lipinski_Violations": violations,

        "Lipinski_Result": (
            "Pass (0 violations)"
            if violations == 0
            else (
                "Pass (1 violation)"
                if violations == 1
                else f"Review ({violations} violations)"
            )
        ),

        "PAINS_Alert": (
            "Alert Found" if pains_matches else "Clean"
        ),

        "PAINS_Details": (
            "; ".join(pains_names) if pains_names else "None"
        )
    })


properties = pd.DataFrame(results)

analysis_df = pd.concat(
    [df.reset_index(drop=True), properties],
    axis=1
)


# =====================================================
# 7. REFERENCE MOLECULE SELECTION
# =====================================================

st.header("🔬 Molecular Similarity Analysis")

selected_ref_idx = st.selectbox(
    "Select reference molecule:",
    options=valid_indices,
    format_func=lambda i: (
        f"{df.iloc[i]['Molecule_Name']} (Row {i+1})"
    )
)

reference_fp = fingerprints[selected_ref_idx]

similarity_values = []

for i, mol in enumerate(molecules):

    if mol is None:
        similarity_values.append(None)

    else:

        similarity = DataStructs.TanimotoSimilarity(
            reference_fp,
            fingerprints[i]
        )

        similarity_values.append(round(similarity, 3))


analysis_df["Tanimoto_Similarity"] = similarity_values

st.info(
    f"Reference molecule: "
    f"{df.iloc[selected_ref_idx]['Molecule_Name']}"
)


# =====================================================
# 8. SCREENING SUMMARY
# =====================================================

st.header("📊 Screening Summary")

valid_count = len(valid_indices)
invalid_count = len(df) - valid_count

lipinski_count = int(
    (analysis_df["Lipinski_Violations"] == 0).sum()
)

pains_count = int(
    (analysis_df["PAINS_Alert"] == "Alert Found").sum()
)

pains_percentage = (
    (pains_count / valid_count) * 100
    if valid_count > 0
    else 0
)

c1, c2, c3, c4, c5 = st.columns(5)

c1.metric("Total Compounds", len(df))
c2.metric("Valid Structures", valid_count)
c3.metric("Lipinski Zero-Violation", lipinski_count)
c4.metric("PAINS Alerts", pains_count)
c5.metric("PAINS Alert %", f"{pains_percentage:.2f}%")


# =====================================================
# 9. INTERACTIVE FILTERS
# =====================================================

st.header("🔎 Compound Filtering")

col1, col2 = st.columns(2)

with col1:

    lipinski_filter = st.selectbox(
        "Lipinski filter",
        ["All", "Zero violations", "One or fewer violations"]
    )

with col2:

    pains_filter = st.selectbox(
        "PAINS filter",
        ["All", "Clean only", "Alert only"]
    )


filtered_df = analysis_df.copy()

if lipinski_filter == "Zero violations":

    filtered_df = filtered_df[
        filtered_df["Lipinski_Violations"] == 0
    ]

elif lipinski_filter == "One or fewer violations":

    filtered_df = filtered_df[
        filtered_df["Lipinski_Violations"] <= 1
    ]


if pains_filter == "Clean only":

    filtered_df = filtered_df[
        filtered_df["PAINS_Alert"] == "Clean"
    ]

elif pains_filter == "Alert only":

    filtered_df = filtered_df[
        filtered_df["PAINS_Alert"] == "Alert Found"
    ]


# =====================================================
# 10. SCREENING RESULTS
# =====================================================

st.header("🧪 Molecular Screening Results")

st.dataframe(
    filtered_df,
    use_container_width=True,
    hide_index=True
)


# =====================================================
# 11. SIMILARITY RANKING
# =====================================================

st.header("🏆 Similarity Ranking")

exclude_reference = st.checkbox(
    "Exclude reference molecule from ranking",
    value=True
)

ranking_df = analysis_df.dropna(
    subset=["Tanimoto_Similarity"]
).copy()

if exclude_reference:
    ranking_df = ranking_df[
        ranking_df.index != selected_ref_idx
    ]

ranking_df = ranking_df.sort_values(
    "Tanimoto_Similarity",
    ascending=False
)

st.dataframe(
    ranking_df[
        [
            "Molecule_Name",
            "Molecular_Weight",
            "LogP",
            "Lipinski_Violations",
            "Tanimoto_Similarity"
        ]
    ],
    use_container_width=True,
    hide_index=True
)


# =====================================================
# 12. 2D MOLECULAR STRUCTURE VIEWER
# =====================================================

st.header("🧬 2D Molecular Structure Viewer")

selected_2d_idx = st.selectbox(
    "Select molecule for 2D structure:",
    options=valid_indices,
    format_func=lambda i: (
        f"{df.iloc[i]['Molecule_Name']} (Row {i+1})"
    ),
    key="viewer_2d"
)

mol_2d = molecules[selected_2d_idx]

st.image(
    Draw.MolToImage(mol_2d, size=(500, 400)),
    caption=df.iloc[selected_2d_idx]["Molecule_Name"]
)


# =====================================================
# 13. 3D MOLECULAR STRUCTURE VIEWER
# =====================================================

st.header("🌐 3D Molecular Structure Viewer")

selected_3d_idx = st.selectbox(
    "Select molecule for 3D visualization:",
    options=valid_indices,
    format_func=lambda i: (
        f"{df.iloc[i]['Molecule_Name']} (Row {i+1})"
    ),
    key="viewer_3d"
)

if st.button("Generate 3D Structure"):

    try:

        mol_3d = Chem.AddHs(
            Chem.Mol(molecules[selected_3d_idx])
        )

        status = AllChem.EmbedMolecule(
            mol_3d,
            randomSeed=42
        )

        if status != 0:

            st.error(
                "3D conformer generation failed. "
                "Try another molecule."
            )

        else:

            optimization = AllChem.MMFFOptimizeMolecule(
                mol_3d
            )

            pdb_block = Chem.MolToPDBBlock(mol_3d)

            viewer = py3Dmol.view(
                width=700,
                height=450
            )

            viewer.addModel(pdb_block, "pdb")
            viewer.setStyle({"stick": {}})
            viewer.setBackgroundColor("white")
            viewer.zoomTo()

            showmol(
                viewer,
                height=450,
                width=700
            )

            st.success(
                f"3D structure generated. "
                f"MMFF optimization status: {optimization}"
            )

    except Exception as e:

        st.error(f"3D visualization error: {e}")


# =====================================================
# 14. NEXUS CHEMISTRY ASSISTANT
# =====================================================

st.divider()

st.header("🤖 Nexus Chemistry Assistant")

st.markdown(
    "Generate an automatic scientific interpretation "
    "of your selected molecule using RDKit calculations."
)

assistant_idx = st.selectbox(
    "Choose a molecule for scientific interpretation:",
    options=valid_indices,
    format_func=lambda i: (
        f"{df.iloc[i]['Molecule_Name']} (Row {i+1})"
    ),
    key="assistant_molecule"
)


def generate_chemistry_report(idx):

    mol = molecules[idx]
    row = analysis_df.iloc[idx]

    name = str(row["Molecule_Name"])
    smiles = str(row["SMILES"])

    mw = row["Molecular_Weight"]
    logp = row["LogP"]
    hbd = row["H_Donors"]
    hba = row["H_Acceptors"]
    tpsa = row["TPSA"]
    rotatable = row["Rotatable_Bonds"]
    aromatic = row["Aromatic_Rings"]

    violations = int(row["Lipinski_Violations"])
    pains = row["PAINS_Alert"]
    pains_details = row["PAINS_Details"]
    similarity = row["Tanimoto_Similarity"]

    report = []

    report.append(f"# Nexus Chemistry Assistant Report: {name}")

    report.append("## 1. Molecular Identification")
    report.append(f"- Molecule name: {name}")
    report.append(f"- Canonical SMILES: {Chem.MolToSmiles(mol)}")
    report.append(f"- Input SMILES: {smiles}")
    report.append("- Structure status: Valid according to RDKit parsing.")

    report.append("## 2. Molecular Property Interpretation")

    report.append(
        f"- Molecular weight: {mw} g/mol. "
        "This represents the calculated mass of the molecule."
    )

    if logp < 0:
        logp_text = (
            "The calculated LogP is negative, suggesting "
            "preference for the aqueous phase in the octanol/water model."
        )
    elif logp <= 3:
        logp_text = (
            "The calculated LogP indicates relatively low to moderate "
            "lipophilicity in the octanol/water model."
        )
    elif logp <= 5:
        logp_text = (
            "The calculated LogP indicates relatively high lipophilicity."
        )
    else:
        logp_text = (
            "The calculated LogP is above 5, indicating high lipophilicity "
            "and a Lipinski threshold violation."
        )

    report.append(f"- LogP: {logp}. {logp_text}")

    report.append(
        f"- Hydrogen bond donors: {hbd}. "
        "These are commonly groups such as NH and OH."
    )

    report.append(
        f"- Hydrogen bond acceptors: {hba}. "
        "These are atoms or groups that can accept hydrogen bonds."
    )

    report.append(
        f"- TPSA: {tpsa} Å². "
        "TPSA estimates the molecule's polar surface area."
    )

    report.append(
        f"- Rotatable bonds: {rotatable}. "
        "This gives a simple indication of molecular flexibility."
    )

    report.append(
        f"- Aromatic rings: {aromatic}. "
        "These are rings identified as aromatic by RDKit."
    )

    report.append("## 3. Lipinski Rule of Five Analysis")

    report.append(
        "- Molecular weight limit: 500 g/mol."
    )
    report.append("- LogP limit: 5.")
    report.append("- Hydrogen bond donor limit: 5.")
    report.append("- Hydrogen bond acceptor limit: 10.")

    report.append(
        f"- Calculated violations: {violations}."
    )

    if violations == 0:

        report.append(
            "Interpretation: The molecule satisfies all four "
            "selected Lipinski property limits."
        )

    elif violations == 1:

        report.append(
            "Interpretation: The molecule has one violation. "
            "It passes the commonly used screening convention "
            "of allowing up to one violation."
        )

    else:

        report.append(
            "Interpretation: The molecule has multiple violations. "
            "Review its properties in the context of the intended "
            "compound series and biological target."
        )

    report.append(
        "Note: Lipinski screening is a heuristic and does not "
        "prove oral absorption, drug activity or safety."
    )

    report.append("## 4. PAINS Filter Interpretation")

    if pains == "Clean":

        report.append(
            "No PAINS alert was identified by the RDKit PAINS catalog."
        )

    else:

        report.append(
            f"Potential PAINS alert detected: {pains_details}. "
            "Review the structural alert and consider orthogonal "
            "assays or interference controls."
        )

    report.append(
        "Note: A clean PAINS result does not guarantee that a compound "
        "is free from assay interference."
    )

    report.append("## 5. Molecular Similarity Analysis")

    report.append(
        f"Tanimoto similarity to the selected reference molecule "
        f"({df.iloc[selected_ref_idx]['Molecule_Name']}): {similarity}."
    )

    report.append(
        "This value is based on Morgan fingerprints with radius 2 "
        "and 2048 bits. Similarity depends on the fingerprint settings "
        "and does not establish similar biological activity."
    )

    report.append("## 6. Preliminary Research Suggestions")

    suggestions = []

    if violations > 1:
        suggestions.append(
            "Review molecular weight, lipophilicity and hydrogen-bonding "
            "properties for possible optimization."
        )

    if logp > 5:
        suggestions.append(
            "Investigate whether high lipophilicity could affect "
            "aqueous solubility or nonspecific binding."
        )

    if tpsa > 140:
        suggestions.append(
            "Review the high polar surface area in relation to the "
            "intended permeability and route of administration."
        )

    if rotatable > 10:
        suggestions.append(
            "Review molecular flexibility and whether reducing "
            "rotatable bonds is feasible."
        )

    if pains == "Alert Found":
        suggestions.append(
            "Investigate the PAINS alert and use suitable assay "
            "interference controls."
        )

    if not suggestions:
        suggestions.append(
            "Consider experimental solubility, permeability, "
            "stability and biological activity studies."
        )

    for suggestion in suggestions:
        report.append(f"- {suggestion}")

    report.append("## 7. Scientific Disclaimer")

    report.append(
        "This report is an automated cheminformatics interpretation "
        "of calculated molecular descriptors and structural alerts. "
        "It is not a validated biological activity prediction, "
        "docking result, toxicity prediction or clinical assessment. "
        "Experimental validation is required."
    )

    return "\n\n".join(report)


if st.button("🧠 Generate Chemistry Research Report"):

    report_text = generate_chemistry_report(assistant_idx)

    st.session_state["nexus_assistant_report"] = report_text


if "nexus_assistant_report" in st.session_state:

    st.markdown("---")

    st.subheader("📋 Scientific Interpretation")

    st.markdown(
        st.session_state["nexus_assistant_report"]
    )

    st.download_button(
        label="📥 Download Chemistry Research Report",
        data=st.session_state["nexus_assistant_report"],
        file_name="Nexus_Chemistry_Assistant_Report.md",
        mime="text/markdown"
    )


# =====================================================
# 15. EXCEL EXPORT
# =====================================================

st.header("📥 Export Professional Screening Report")

output = io.BytesIO()

try:

    with pd.ExcelWriter(
        output,
        engine="openpyxl"
    ) as writer:

        analysis_df.to_excel(
            writer,
            index=False,
            sheet_name="Screening_Results"
        )

        ranking_df.to_excel(
            writer,
            index=False,
            sheet_name="Similarity_Ranking"
        )

    st.download_button(
        label="📥 Download Nexus Excel Report",
        data=output.getvalue(),
        file_name="Nexus_Advanced_Screening_Report.xlsx",
        mime=(
            "application/vnd.openxmlformats-officedocument."
            "spreadsheetml.sheet"
        )
    )

except Exception as e:

    st.error(f"Excel export error: {e}")


# =====================================================
# 16. SCIENTIFIC DISCLAIMER
# =====================================================

st.divider()

st.caption(
    "Nexus is a cheminformatics research and preliminary "
    "screening tool. Lipinski filters, PAINS alerts and "
    "molecular similarity do not establish biological "
    "activity, safety, efficacy or clinical suitability. "
    "Experimental validation is essential."
)

st.caption(
    "Nexus Advanced Cheminformatics | Portfolio Edition"
)
