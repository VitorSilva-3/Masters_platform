
import streamlit as st
from PIL import Image
from config import AppConfig
from services.data_manager import DataManager
from services.taxonomy_service import TaxonomyService
from ui.view_taxonomy import render_taxonomy_view

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="Taxonomy - MicroValue", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="Taxonomy - MicroValue", layout="wide")

@st.cache_resource
def get_services():
    return DataManager(), TaxonomyService(email=AppConfig.EMAIL)

@st.cache_data
def get_main_dataframe():
    return get_services()[0].load_data()

data_manager, taxonomy_service = get_services()
df = get_main_dataframe()

st.title("Taxonomy")
if df.empty:
    st.error("Data file not found or empty!")
else:
    render_taxonomy_view(taxonomy_service, df)