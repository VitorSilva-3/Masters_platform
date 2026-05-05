
import streamlit as st
from config import AppConfig
from services.taxonomy_service import TaxonomyService
from services.pubmed_service import PubMedService
from ui.view_literature import render_literature_view
from utils import configure_page, load_core_datasets

configure_page("Scientific literature (PubMed)")

@st.cache_resource
def get_services():
    """Initialize and cache both the TaxonomyService and PubMedService instances."""

    tax_service = TaxonomyService(email=AppConfig.EMAIL)
    pub_service = PubMedService(email=AppConfig.EMAIL)
    return tax_service, pub_service

taxonomy_service, pubmed_service = get_services()
df_enzymes, df_transporters = load_core_datasets()

st.title("Scientific literature (PubMed)")

if df_enzymes.empty and df_transporters.empty:
    st.error("Data files not found or empty!")
else:
    render_literature_view(pubmed_service, df_enzymes, df_transporters, taxonomy_service)