
import streamlit as st
from PIL import Image

# Import Configuration
from config import AppConfig

# Import Services (Backend)
from services.data_manager import DataManager
from services.taxonomy_service import TaxonomyService
from services.pubmed_service import PubMedService
from services.uniprot_service import UniprotService

# Import UI Components (Frontend)
from ui.tab_data import render_data_tab
from ui.tab_kegg import render_kegg_tab 
from ui.tab_taxonomy import render_taxonomy_tab
from ui.tab_articles import render_articles_tab
from ui.tab_uniprot import render_uniprot_tab

# 1. Page Configuration
try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(
        page_title="MicroValue",
        page_icon=icon_image,
        layout="wide"
    )
except FileNotFoundError:
    st.set_page_config(page_title="MicroValue", layout="wide")

# 2. Initialize Services (Dependency Injection)
# O DataManager já sabe que tem de ir à pasta "data/" por causa da alteração que fizemos nele!
data_manager = DataManager() 
taxonomy_service = TaxonomyService(email=AppConfig.EMAIL)
pubmed_service = PubMedService(email=AppConfig.EMAIL)
uniprot_service = UniprotService()

# 3. Sidebar
st.sidebar.header("Information")
with st.sidebar.expander("About this platform", expanded=True):
    st.markdown(
        "This platform integrates data on agro-industrial byproducts and the enzymes produced by various microalgae, "
        "organised by species. It identifies which byproducts can serve as substrates for each microalga, "
        "indicating those with the greatest growth potential."
    )

# 4. Main Title
st.title("MicroValue")

# 5. Load Data (Using the DataManager)
@st.cache_data
def get_main_dataframe():
    return data_manager.load_data()

df = get_main_dataframe()

if df.empty:
    st.error("Data file not found or empty!")
    st.info("Please run the data builder script locally first: `python build_caches.py`")
    st.stop()

# 6. Tabs Layout
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Data", 
    "Literature (PubMed)", 
    "Taxonomy", 
    "KEGG Enzyme", 
    "UniProt"
])

# --- TAB 1: DATA ---
with tab1:
    render_data_tab(df, taxonomy_service)

# --- TAB 2: ARTICLES ---
with tab2:
    render_articles_tab(pubmed_service, df)

# --- TAB 3: TAXONOMY ---
with tab3:
    render_taxonomy_tab(taxonomy_service, df)

# --- TAB 4: KEGG ---
with tab4:
    render_kegg_tab(df)

# --- TAB 5: UNIPROT ---
with tab5:
    render_uniprot_tab(uniprot_service, df)