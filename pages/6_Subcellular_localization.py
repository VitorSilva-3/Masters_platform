
import streamlit as st
from PIL import Image
from config import AppConfig
from services.data_manager import DataManager
from services.taxonomy_service import TaxonomyService
from ui.view_subcellular_localization import render_localization_view

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="Subcellular localization - MicroValue", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="Subcellular localization - MicroValue", layout="wide")

@st.cache_resource
def get_taxonomy_service():
    """Initialize and cache the TaxonomyService instance."""

    return TaxonomyService(email=AppConfig.EMAIL)

@st.cache_data
def load_all_datasets():
    """Initialize and cache both datasets."""

    manager_enzymes = DataManager(file_path="data/enzymes_data.csv", data_type="enzyme")
    df_enzymes = manager_enzymes.load_data()
    
    manager_transporters = DataManager(file_path="data/transporters_data.csv", data_type="transporter")
    df_transporters = manager_transporters.load_data()
    
    return df_enzymes, df_transporters

taxonomy_service = get_taxonomy_service()
df_enzymes, df_transporters = load_all_datasets()

st.title("Subcellular localization")

if df_enzymes.empty and df_transporters.empty:
    st.error("Data files not found or empty!")
else:
    render_localization_view(df_enzymes, df_transporters, taxonomy_service)