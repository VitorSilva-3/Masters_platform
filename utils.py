
import streamlit as st
import os
import logging
import json
import pandas as pd
import streamlit as st
from typing import Dict, Any
from PIL import Image
from config import AppConfig
from Bio import Entrez

logger = logging.getLogger(__name__)

def load_lottiefile(filepath: str):
    """Loads a Lottie animation JSON file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None

def get_group_for_species(species: str, taxonomy_service) -> str:
    """Helper function to find which target taxa a species belongs to."""
    
    lineage = taxonomy_service.fetch_taxonomy_lineage(species)
    if isinstance(lineage, list):
        for node in lineage:
            if node.get("name") in AppConfig.TARGET_TAXA:
                return node.get("name")
        for node in lineage:
            if node.get("rank") == "class":
                return node.get("name")
    return "Other / Unknown"

def add_taxonomic_class_column(df: pd.DataFrame, taxonomy_service, loading_text: str = "Classifying species...") -> pd.DataFrame:
    """
    Checks if the 'Class' column exists in the DataFrame. 
    If not, it maps each 'Specie' to its taxonomic class and adds the column.
    Returns the updated DataFrame.
    """
    if df.empty:
        return df
        
    if "Class" not in df.columns and "Specie" in df.columns:
        with st.spinner(loading_text):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)
            
            cols = df.columns.tolist()
            cols.insert(cols.index("Specie") + 1, cols.pop(cols.index("Class")))
            df = df[cols]
            
    return df

def configure_page(page_name: str):
    """Sets the standard page configuration and logo for all pages."""

    try:
        icon_image = Image.open("images/logo.jpg")
        st.set_page_config(page_title=f"{page_name} - MicroValue", page_icon=icon_image, layout="wide")
    except FileNotFoundError:
        st.set_page_config(page_title=f"{page_name} - MicroValue", layout="wide")

@st.cache_data(show_spinner=False)
def load_core_datasets():
    """Loads and caches the main enzyme and transporter datasets globally."""

    from services.data_manager import DataManager

    manager_enzymes = DataManager(file_path="data/enzymes_data.csv", data_type="enzyme")
    df_enzymes = manager_enzymes.load_data()
    
    manager_transporters = DataManager(file_path="data/transporters_data.csv", data_type="transporter")
    df_transporters = manager_transporters.load_data()
    
    return df_enzymes, df_transporters

def load_json_cache(cache_file: str, service_name: str = "Service") -> Dict[str, Any]:
    """Loads a JSON cache."""

    if os.path.exists(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"[{service_name}] Error loading cache: {e}")
    return {}

def save_json_cache(cache_file: str, data: Dict[str, Any], service_name: str = "Service") -> None:
    """Saves a dictionary to a JSON cache."""

    try:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)
    except Exception as e:
        logger.error(f"[{service_name}] Error saving cache: {e}")

def setup_ncbi_entrez():
    """Configures NCBI Entrez globally with email and optional API Key from secrets."""
    
    try:
        Entrez.email = st.secrets["NCBI_EMAIL"]
    except KeyError:
        logger.error(" NCBI_EMAIL not found in Streamlit secrets. Please set it to a valid email address.")
        return False

    if "NCBI_API_KEY" in st.secrets:
        Entrez.api_key = st.secrets["NCBI_API_KEY"]
        logger.info(" NCBI API Key found. Using it for faster requests.")
    else:
        logger.warning("NCBI API Key not found. The extraction will be slower.")
        
    return True 