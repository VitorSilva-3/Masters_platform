
import streamlit as st
import pandas as pd
from config import AppConfig
from services.taxonomy_service import TaxonomyService
from ui.view_taxonomy import render_taxonomy_view
from utils import configure_page, load_core_datasets

configure_page("Taxonomy")

@st.cache_resource
def get_taxonomy_service():
    return TaxonomyService(email=AppConfig.EMAIL)

taxonomy_service = get_taxonomy_service()
df_enzymes, df_transporters = load_core_datasets()

dfs_to_concat = [df for df in [df_enzymes, df_transporters] if not df.empty]
if dfs_to_concat:
    df_combined = pd.concat(dfs_to_concat, ignore_index=True).drop_duplicates(subset=["Specie"])
else:
    df_combined = pd.DataFrame()

st.title("Taxonomy")

if df_combined.empty:
    st.error("Data files not found or empty!")
else:
    render_taxonomy_view(taxonomy_service, df_combined)