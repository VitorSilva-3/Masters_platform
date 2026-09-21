
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
df_enzymes, _= load_core_datasets()

st.title("UniProt")

if df_enzymes.empty:
    st.error("Enzime data file not found or empty.")
else:
    render_uniprot_view(uniprot_service, df_enzymes)