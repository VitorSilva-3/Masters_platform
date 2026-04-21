
import streamlit as st
import time
import json
from PIL import Image
from streamlit_lottie import st_lottie

def load_lottiefile(filepath: str):
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None

LOTTIE_DNA = "images/DNA strand.json" 
LOTTIE_FACTORY = "images/Factory pollution city air and water.json"

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="MicroValue - Home", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="MicroValue - Home", layout="wide")


if 'app_loaded' not in st.session_state:
    st.session_state.app_loaded = False

if not st.session_state.app_loaded:
    hide_sidebar_css = """
        <style>
            [data-testid="stSidebar"] {display: none !important;}
            [data-testid="collapsedControl"] {display: none !important;}
        </style>
    """
    st.markdown(hide_sidebar_css, unsafe_allow_html=True)
    
    st.markdown("<br><br>", unsafe_allow_html=True)
    
    anim_col1, anim_col2 = st.columns(2)
    
    with anim_col1:
        anim1 = load_lottiefile(LOTTIE_DNA)
        if anim1:
            st_lottie(anim1, height=350, key="loading_dna")
            
    with anim_col2:
        anim2 = load_lottiefile(LOTTIE_FACTORY)
        if anim2:
            st_lottie(anim2, height=350, key="loading_factory")
            
    st.markdown("<h2 style='text-align: center; color: #4CAF50;'>Loading MicroValue...</h2>", unsafe_allow_html=True)
    
    bar_col1, bar_col2, bar_col3 = st.columns([1, 2, 1])
    with bar_col2:
        progress_bar = st.progress(0)
        status_text = st.empty()

        loading_steps = [
            "Initializing platform environment...",
            "Loading central databases...",
            "Configuring modules...",
            "Preparing user interface...",
            "Ready."
        ]

        for i, step in enumerate(loading_steps):
            status_text.markdown(f"<p style='text-align: center; color: #888; font-size: 1.1em;'><i>{step}</i></p>", unsafe_allow_html=True)
            progress_bar.progress((i + 1) / len(loading_steps))
            time.sleep(1)

    time.sleep(0.5)
    st.session_state.app_loaded = True 
    st.rerun()

else:
    st.title("Welcome to MicroValue")
    st.markdown("#### *Biotechnological valorization of agro-industrial byproducts through microalgae and cyanobacteria*")
    st.divider()

    st.markdown(
        """
        This platform integrates data on agro-industrial byproducts as well as the enzymes and transporters produced by various microalgae 
        and cyanobacteria organised by species. It identifies which byproducts can serve as substrates for each microalga/cyanobacterium,
        indicating those with the greatest growth potential.
        """
    )

    st.info("👈 **Use the left menu to navigate through the platform's specialized modules.**")

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("### Platform modules")

    with st.container(border=True):
        st.page_link("pages/1_Data.py", label = ":blue[**Data**]", use_container_width = True)
        st.write("Explore the comprehensive datasets cross-referencing species, target sugars, degrading enzymes, and transport mechanisms.")

    with st.container(border=True):
        st.page_link("pages/2_Scientific_literature.py", label = ":blue[**Scientific literature**]", use_container_width = True)
        st.write("Access an automated curation of PubMed publications relevant to the catalogued species and their biochemical repertoire.")
    
    with st.container(border=True):
        st.page_link("pages/3_Taxonomy.py", label = ":blue[**Taxonomy**]", use_container_width = True)
        st.write("Navigate the detailed evolutionary lineage and taxonomic classification of the catalogued microalgae and cyanobacteria.")
    
    with st.container(border=True):
        st.page_link("pages/4_KEGG.py", label = ":blue[**KEGG**]", use_container_width = True)
        st.write("Map metabolic networks and analyze the specific catalytic reactions underlying biomass degradation.")

    with st.container(border=True):
        st.page_link("pages/5_UniProt.py", label = ":blue[**UniProt**]", use_container_width = True)
        st.write("Investigate the general biochemical properties, functional annotations, and structural features of the targeted proteins.")

    with st.container(border=True):
        st.page_link("pages/6_Subcellular_localization.py", label = ":blue[**Subcellular localization**]", use_container_width = True)
        st.write("Examine AI-driven predictions of subcellular compartments utilizing DeepLoc 2.0 and DeepLocPro models.")