
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
def get_services():
    return DataManager(), UniprotService()

@st.cache_data
def get_main_dataframe():
    return get_services()[0].load_data()

data_manager, uniprot_service = get_services()
df = get_main_dataframe()

st.title("UniProt")
if df.empty:
    st.error("Data file not found or empty!")
else:
    render_uniprot_view(uniprot_service, df)