import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64

def smiles_to_mols(smiles_list):
    return [Chem.MolFromSmiles(smiles) for smiles in smiles_list if smiles.strip()]

def mol_to_svg(mol):
    img = Draw.MolToImage(mol, size=(600, 600))
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

def main():
    st.title("SMILES to Molecule Viewer")
    
    # User input for SMILES strings
    smiles_input = st.text_area("Enter SMILES strings (one per line):", 
                                "Cc1ccc(c(c1)N\\2C(=O)CS/C2=N\\C(=O)Nc3ccc(cc3F)c4ncn(n4)c5ccc(cc5)OC(F)(F)F)COCC(F)(F)F")
    
    if st.button("Display Molecules"):
        smiles_list = smiles_input.split('\n')
        mols = smiles_to_mols(smiles_list)
        
        if mols:
            st.write(f"Displaying {len(mols)} molecule(s):")
            
            # Create a 3-column layout
            cols = st.columns(3)
            
            for i, mol in enumerate(mols):
                svg = mol_to_svg(mol)
                cols[i % 3].image(f"data:image/png;base64,{svg}", use_column_width=True)
                cols[i % 3].text(f"Molecule {i+1}")
        else:
            st.error("No valid SMILES strings found. Please check your input.")

if __name__ == "__main__":
    main()
