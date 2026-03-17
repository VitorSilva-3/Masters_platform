
import streamlit as st
from PIL import Image
from services.data_manager import DataManager
from ui.view_kegg import render_kegg_view

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="KEGG - MicroValue", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="KEGG - MicroValue", layout="wide")

@st.cache_resource
def get_data_manager():
    return DataManager()

@st.cache_data
def get_main_dataframe():
    return get_data_manager().load_data()

df = get_main_dataframe()

st.title("KEGG")
if df.empty:
    st.error("Data file not found or empty!")
else:
    render_kegg_view(df)