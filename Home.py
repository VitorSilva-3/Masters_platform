
import streamlit as st
from PIL import Image

try:
    icon_image = Image.open("images/logo.jpg")
    st.set_page_config(page_title="MicroValue - Home", page_icon=icon_image, layout="wide")
except FileNotFoundError:
    st.set_page_config(page_title="MicroValue - Home", layout="wide")

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