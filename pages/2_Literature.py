
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
    return DataManager(), TaxonomyService(email=AppConfig.EMAIL), PubMedService(email=AppConfig.EMAIL)

@st.cache_data
def get_main_dataframe():
    return get_services()[0].load_data()

data_manager, taxonomy_service, pubmed_service = get_services()
df = get_main_dataframe()

st.title("Literature (PubMed)")
if df.empty:
    st.error("Data file not found or empty!")
else:
    render_literature_view(pubmed_service, df, taxonomy_service)