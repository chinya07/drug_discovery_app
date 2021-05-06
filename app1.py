# rdkit.Chem.Lipinski.NumRotatableBonds

"""
Docker Container: https://hub.docker.com/r/continuumio/anaconda3
RDKit Installation: https://www.rdkit.org/docs/Install.html
"""
import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors
from rdkit.Chem.Lipinski import NumRotatableBonds


st.title("Filter FDA Approved Drugs by Lipinski's Rule-of-Five with Streamlit")


@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(
        "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
    ).dropna()
    return df

# Calculate descriptors
def calc_mw(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the molecular weight"""
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)

def calc_logp(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the LogP"""
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)

def calc_NumHDonors(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHDonors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)

def calc_NumHAcceptors(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHAcceptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHAcceptors(mol)

def calc_NumRotatableBonds(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHAcceptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumRotatableBonds(mol)


# Copy the dataset so any changes are not applied to the original cached version
df = download_dataset().copy()
df["MW"] = df.apply(lambda x: calc_mw(x["smiles"]), axis=1)
df["LogP"] = df.apply(lambda x: calc_logp(x["smiles"]), axis=1)
df["NumHDonors"] = df.apply(lambda x: calc_NumHDonors(x["smiles"]), axis=1)
df["NumHAcceptors"] = df.apply(lambda x: calc_NumHAcceptors(x["smiles"]), axis=1)
df["NumRotatableBonds"] = df.apply(lambda x: calc_NumRotatableBonds(x["smiles"]), axis=1)


# Sidebar panel
st.sidebar.header('Set parameters')
st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')
weight_cutoff = st.sidebar.slider(
    label="Molecular weight",
    min_value=0,
    max_value=1000,
    value=500,
    step=10,
)
logp_cutoff = st.sidebar.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=5,
    step=1,
)
NumHDonors_cutoff = st.sidebar.slider(
    label="NumHDonors",
    min_value=0,
    max_value=15,
    value=5,
    step=1,
)
NumHAcceptors_cutoff = st.sidebar.slider(
    label="NumHAcceptors",
    min_value=0,
    max_value=20,
    value=10,
    step=1,
)

df_result = df[df["MW"] < weight_cutoff]
df_result2 = df_result[df_result["LogP"] < logp_cutoff]
df_result3 = df_result2[df_result2["NumHDonors"] < NumHDonors_cutoff]
df_result4 = df_result3[df_result3["NumHAcceptors"] < NumHAcceptors_cutoff]

st.write(df_result4.shape)
st.write(df_result4)


raw_html = mols2grid.display(df_result4,
                            #subset=["Name", "img"],
                            subset=["img", "Name", "MW", "LogP", "NumHDonors", "NumHAcceptors"],
                            mapping={"smiles": "SMILES", "generic_name": "Name"})._repr_html_()
components.html(raw_html, width=900, height=1100, scrolling=False)


st.title("rule of three (RO3) for defining lead-like compounds")
st.sidebar.header('Set additional parameters')
st.sidebar.write('*During drug discovery, lipophilicity and molecular weight are often increased in order to improve the affinity and selectivity of the drug candidate. Hence it is often difficult to maintain drug-likeness (i.e., RO5 compliance) during hit and lead optimization. Hence it has been proposed that members of screening libraries from which hits are discovered should be biased toward lower molecular weight and lipophility so that medicinal chemists will have an easier time in delivering optimized drug development candidates that are also drug-like. Hence the rule of five has been extended to the rule of three (RO3) for defining lead-like compounds.*')
weight_cutoff = st.sidebar.slider(
    label="Molecular weight",
    min_value=0,
    max_value=1000,
    value=300,
    step=10,
)
logp_cutoff = st.sidebar.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=3,
    step=1,
)
NumHDonors_cutoff = st.sidebar.slider(
    label="NumHDonors",
    min_value=0,
    max_value=15,
    value=3,
    step=1,
)
NumHAcceptors_cutoff = st.sidebar.slider(
    label="NumHAcceptors",
    min_value=0,
    max_value=20,
    value=3,
    step=1,
)

NumRotatableBonds_cutoff = st.sidebar.slider(
    label="NumRotatableBonds",
    min_value=0,
    max_value=20,
    value=3,
    step=1,
)

df_llc_result = df[df["MW"] < weight_cutoff]
df_llc_result2 = df_llc_result[df_llc_result["LogP"] < logp_cutoff]
df_llc_result3 = df_llc_result2[df_llc_result2["NumHDonors"] < NumHDonors_cutoff]
df_llc_result4 = df_llc_result3[df_llc_result3["NumHAcceptors"] < NumHAcceptors_cutoff]
df_llc_result5 = df_llc_result4[df_llc_result3["NumRotatableBonds"] < NumRotatableBonds_cutoff]

st.write(df_llc_result5.shape)
st.write(df_llc_result5)


raw_html = mols2grid.display(df_llc_result5,
                            #subset=["Name", "img"],
                            subset=["img", "Name", "MW", "LogP", "NumHDonors", "NumHAcceptors", "NumRotatableBonds"],
                            mapping={"smiles": "SMILES", "generic_name": "Name"})._repr_html_()
components.html(raw_html, width=900, height=1100, scrolling=False)