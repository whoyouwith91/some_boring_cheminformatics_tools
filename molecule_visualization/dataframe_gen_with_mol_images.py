import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
import io

# Function to read the uploaded file
def load_file(uploaded_file):
    # If the file is a CSV, read normally
    if uploaded_file.name.endswith('.csv'):
        df = pd.read_csv(uploaded_file)
    # If the file is a TXT or SMI, attempt to read it as space-separated columns and name the column 'SMILES'
    elif uploaded_file.name.endswith('.txt') or uploaded_file.name.endswith('.smi'):
        df = pd.read_csv(uploaded_file, sep='\s+', names=['SMILES'])
    return df

# Title of the app
st.title("SMILES to Molecule Converter")

# File uploader widget (CSV, TXT, or SMI file)
uploaded_file = st.file_uploader("Upload CSV, TXT, or SMI file", type=["csv", "txt", "smi"])

if uploaded_file is not None:
    # Load the file and display a preview
    df = load_file(uploaded_file)
    st.write("Uploaded file preview:")
    st.dataframe(df.head())

    # Select the SMILES column (if the file contains multiple columns, like in a CSV)
    smiles_column = st.selectbox("Select the SMILES column", df.columns)

    # Check if the SMILES column has been selected
    if smiles_column:
        st.write(f"SMILES column selected: {smiles_column}")

        # Button to process the file and generate an Excel file
        if st.button("Convert SMILES to Molecules and Save to Excel"):
            # Perform the conversion only after the button is clicked
            st.write("Converting SMILES to Molecules... Please wait.")
            
            # Convert SMILES to RDKit mol objects
            PandasTools.AddMoleculeColumnToFrame(df, smiles_column, 'Molecule')

            # Create a bytes buffer to store the Excel file
            output = io.BytesIO()

            # Save the DataFrame with molecular images to Excel
            PandasTools.SaveXlsxFromFrame(df, output, molCol='Molecule')
            output.seek(0)  # Move the cursor to the beginning of the buffer

            # Inform the user that conversion is complete
            st.success("Conversion Complete! You can now download the Excel file.")

            # Provide a download link for the Excel file
            st.download_button(label="Download Excel file", 
                               data=output, 
                               file_name="molecules_with_images.xlsx", 
                               mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
