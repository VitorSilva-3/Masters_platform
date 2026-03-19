
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
    st.markdown("#### *Biotechnological valorization of agro-industrial byproducts through microalgae*")
    st.divider()

    st.markdown(
        """
        <div style='font-size: 1.1em; line-height: 1.6; color: #444;'>
        This platform integrates data on agro-industrial byproducts and the enzymes produced by various microalgae, 
        organised by species. It identifies which byproducts can serve as substrates for each microalga, 
        indicating those with the greatest growth potential.
        </div>
        <br>
        """, 
        unsafe_allow_html=True
    )

    st.info("👈 **Use the left menu to navigate through the platform's specialized modules.**")

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("### Platform modules")

    col1, col2, col3 = st.columns(3)

    with col1:
        with st.container(border=True):
            st.markdown("### Data")
            st.write("Explore the central database crossing species, target sugars, and enzymes.")

        with st.container(border=True):
            st.markdown("### Literature")
            st.write("Access a collection of publications related to the species catalogued in the platform.")

    with col2:
        with st.container(border=True):
            st.markdown("### KEGG")
            st.write("Analyze catalytic reactions and metabolic pathways from the KEGG database.")

        with st.container(border=True):
            st.markdown("### Taxonomy")
            st.write("Compare the complete evolutionary lineage of microalgae side-by-side.")

    with col3:
        with st.container(border=True):
            st.markdown("### UniProt")
            st.write("Access deep biochemical properties of the enzymes from the UniProt database.")