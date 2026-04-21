
import streamlit as st
from PIL import Image
from services.data_manager import DataManager
from services.uniprot_service import UniprotService
from ui.view_uniprot import render_uniprot_view

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="UniProt - MicroValue", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="UniProt - MicroValue", layout="wide")

@st.cache_resource
def get_uniprot_service():
    """Initializes the UniprotService."""
    return UniprotService()

@st.cache_data
def load_all_datasets():
    """Loads the enzyme and transporter datasets."""

    manager_enzymes = DataManager(file_path="data/enzymes_data.csv", data_type="enzyme")
    df_enzymes = manager_enzymes.load_data()
    
    manager_transporters = DataManager(file_path="data/transporters_data.csv", data_type="transporter")
    df_transporters = manager_transporters.load_data()
    
    return df_enzymes, df_transporters

uniprot_service = get_uniprot_service()
df_enzymes, df_transporters = load_all_datasets()

st.title("UniProt")

if df_enzymes.empty and df_transporters.empty:
    st.error("Data files not found or empty!")
else:
    tab_enz, tab_trans = st.tabs(["Enzymes", "Transporters"])
    
    with tab_enz:
        render_uniprot_view(uniprot_service, df_enzymes, dataset_type="enzymes")
        
    with tab_trans:
        render_uniprot_view(uniprot_service, df_transporters, dataset_type="transporters")