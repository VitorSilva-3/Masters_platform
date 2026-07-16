
import streamlit as st
from config import AppConfig
from services.taxonomy_service import TaxonomyService
from ui.view_data import render_data_view
from utils import configure_page, load_core_datasets

configure_page("Data")

@st.cache_resource
def get_taxonomy_service():
    """Initializes and caches the TaxonomyService instance."""

    return TaxonomyService()

taxonomy_service = get_taxonomy_service()
df_enzymes, df_transporters = load_core_datasets()

st.title("Data")

if df_enzymes.empty and df_transporters.empty:
    st.error("Data files not found or empty!")
else:
    render_data_view(df_enzymes, df_transporters, taxonomy_service)