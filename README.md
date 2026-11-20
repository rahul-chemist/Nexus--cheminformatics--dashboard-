🧪 Nexus AI Cheminformatics Pro

Advanced Molecular Analysis & Virtual Screening Dashboard

Nexus is a Python-based cheminformatics application designed to explore molecular properties, perform preliminary compound screening, analyze molecular similarity, and visualize chemical structures in 2D and 3D.

The dashboard combines RDKit molecular calculations with an interactive Streamlit interface to support early-stage molecular analysis and research workflows.

🚀 Features

- Molecular Property Analysis: Molecular Weight (MW), LogP, Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA), TPSA, rotatable bonds, and aromatic rings.
- Lipinski Rule of Five: Evaluate selected molecular properties and identify rule violations.
- PAINS Filtering: Identify potential structural alerts using the RDKit PAINS filter catalog.
- Molecular Similarity: Calculate Tanimoto similarity using Morgan fingerprints.
- Virtual Screening: Filter compounds using molecular property and structural alert criteria.
- Similarity Ranking: Rank compounds based on similarity to a selected reference molecule.
- 2D Molecular Visualization: Generate chemical structure images.
- 3D Molecular Visualization: Generate 3D conformations with hydrogen addition and MMFF optimization.
- Chemistry Assistant: Generate automated scientific interpretations of calculated molecular properties and screening results.
- Excel Export: Export molecular screening data and similarity rankings.
- Molecular Dataset Upload: Support CSV, TXT, SMI, and SDF formats.

🧬 Molecular Analysis Workflow

1. Upload a molecular dataset or use the built-in example molecules.
2. Validate molecular SMILES structures using RDKit.
3. Calculate molecular descriptors.
4. Evaluate Lipinski Rule of Five violations.
5. Identify potential PAINS structural alerts.
6. Generate Morgan fingerprints.
7. Calculate Tanimoto similarity against a selected reference molecule.
8. Filter and rank compounds.
9. Visualize molecular structures in 2D and 3D.
10. Generate a scientific interpretation and export results.

🛠️ Technology Stack

- Python
- Streamlit
- RDKit
- Pandas
- py3Dmol
- stmol
- OpenPyXL

📊 Exported Reports

The application generates an Excel workbook containing:

- Screening_Results: Molecular properties, structural alerts, screening information, and similarity scores.
- Similarity_Ranking: Compounds ranked by Tanimoto similarity.

The Chemistry Assistant report can also be downloaded in Markdown format.

🎯 Project Objective

The objective of Nexus is to build an interactive cheminformatics workflow that helps organize molecular data, evaluate chemical properties, prioritize compounds using computational filters, and support preliminary research decisions.

⚠️ Scientific Limitations

Nexus is a preliminary cheminformatics and screening tool. Molecular descriptors, Lipinski rules, PAINS alerts, and fingerprint similarity do not independently establish biological activity, safety, efficacy, drug-likeness, or clinical suitability.

3D conformer generation is not molecular docking. Automated chemical interpretations are rule-based and are not validated AI predictions. Experimental validation is essential.

👨‍🔬 Author

Rahul Survase

Chemistry Research | Cheminformatics | Molecular Analysis | Computational Drug Discovery

---

Portfolio project developed to explore the integration of chemistry, cheminformatics, and computational molecular analysis.
