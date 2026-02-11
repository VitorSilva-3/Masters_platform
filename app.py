
import streamlit as st
from PIL import Image

# Import Services (Backend)
from services.data_manager import DataManager
#from services.ncbi_service import NCBIService
#from services.kegg_service import KeggService
#from services.uniprot_service import UniprotService

# Import UI Components (Frontend)
from ui.tab_data import render_data_tab
# Import other tabs as you fix them (commented out for now to avoid errors)
# from ui.tab_articles import render_articles_tab
# from ui.tab_taxonomy import render_taxonomy_tab
# from ui.tab_kegg import render_kegg_tab
# from ui.tab_uniprot import render_uniprot_tab

# 1. Page Configuration

icon_image = Image.open("images/logo.jpg")

st.set_page_config(
    page_title="MicroValue ",
    page_icon=icon_image,
    layout="wide"
)

# 2. Initialize Services (Dependency Injection)
# We create the objects once to use throughout the app
data_manager = DataManager()
#ncbi_service = NCBIService(email=AppConfig.EMAIL)
#kegg_service = KeggService()
#uniprot_service = UniprotService()

# 3. Sidebar
st.sidebar.header("Information & filters")
with st.sidebar.expander("About this platform", expanded=True):
    st.markdown(
        "This platform integrates data on agro-industrial byproducts and the enzymes produced by various microalgae, "
        "organised by species. It identifies which byproducts can serve as substrates for each microalga, "
        "indicating those with the greatest growth potential."
    )

# 4. Main Title
st.title("MicroValue")

# 5. Load Data (Using the DataManager)
# Using Streamlit cache to avoid reloading CSV on every interaction
@st.cache_data
def get_main_dataframe():
    return data_manager.load_data()

df = get_main_dataframe()

if df.empty:
    st.error("Data file not found or empty!")
    st.info("Please run the data builder script locally first: `python services/data_manager.py`")
    st.stop()

# 6. Tabs Layout
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Data", 
    "Articles (NCBI)", 
    "Taxonomy", 
    "KEGG Enzyme", 
    "UniProt"
])

# --- TAB 1: DATA ---
with tab1:
    # We pass the loaded dataframe and the manager to the UI component
    render_data_tab(df)

# --- TAB 2: ARTICLES ---
with tab2:
    st.info("Module under construction. Connect `ui/tab_articles.py` to `ncbi_service`.")
    # render_articles_tab(ncbi_service, df)

# --- TAB 3: TAXONOMY ---
with tab3:
    st.info("Module under construction. Connect `ui/tab_taxonomy.py` to `ncbi_service`.")
    # render_taxonomy_tab(ncbi_service, df)

# --- TAB 4: KEGG ---
with tab4:
    st.info("Module under construction. Connect `ui/tab_kegg.py` to `kegg_service`.")
    # render_kegg_tab(kegg_service)

# --- TAB 5: UNIPROT ---
with tab5:
    st.info("Module under construction. Connect `ui/tab_uniprot.py` to `uniprot_service`.")
    # render_uniprot_tab(uniprot_service)