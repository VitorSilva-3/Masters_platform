
import streamlit as st
from ui.view_kegg import render_kegg_view
from utils import configure_page, load_core_datasets

configure_page("KEGG")

df_enzymes, _ = load_core_datasets()

st.title("KEGG")

if df_enzymes.empty:
    st.error("Enzymes data file not found or empty!")
else:
    render_kegg_view(df_enzymes)