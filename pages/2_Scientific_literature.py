
import streamlit as st
from PIL import Image
from config import AppConfig
from services.data_manager import DataManager
from services.taxonomy_service import TaxonomyService
from services.pubmed_service import PubMedService
from ui.view_literature import render_literature_view

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="Literature - MicroValue", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="Literature - MicroValue", layout="wide")

@st.cache_resource
def get_services():
    """Initialize and cache both the TaxonomyService and PubMedService instances."""

    tax_service = TaxonomyService(email=AppConfig.EMAIL)
    pub_service = PubMedService(email=AppConfig.EMAIL)
    return tax_service, pub_service

@st.cache_data
def load_all_datasets():
    """Initialize and cache both datasets."""

    manager_enzymes = DataManager(file_path="data/enzymes_data.csv", data_type="enzyme")
    df_enzymes = manager_enzymes.load_data()
    
    manager_transporters = DataManager(file_path="data/transporters_data.csv", data_type="transporter")
    df_transporters = manager_transporters.load_data()
    
    return df_enzymes, df_transporters

taxonomy_service, pubmed_service = get_services()
df_enzymes, df_transporters = load_all_datasets()

st.title("Scientific literature (PubMed)")

if df_enzymes.empty and df_transporters.empty:
    st.error("Data files not found or empty!")
else:
    render_literature_view(pubmed_service, df_enzymes, df_transporters, taxonomy_service)