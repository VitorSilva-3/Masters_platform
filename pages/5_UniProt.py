
import streamlit as st
from services.uniprot_service import UniprotService
from ui.view_uniprot import render_uniprot_view
from utils import configure_page, load_core_datasets

configure_page("UniProt")

@st.cache_resource
def get_uniprot_service():
    """Initializes the UniprotService."""

    return UniprotService()

uniprot_service = get_uniprot_service()
df_enzymes, df_transporters = load_core_datasets()

st.title("UniProt")

if df_enzymes.empty and df_transporters.empty:
    st.error("Data files not found or empty!")
else:
    tab_enz, tab_trans = st.tabs(["Enzymes", "Transporters"])
    
    with tab_enz:
        render_uniprot_view(uniprot_service, df_enzymes, dataset_type="enzymes")
        
    with tab_trans:
        render_uniprot_view(uniprot_service, df_transporters, dataset_type="transporters")