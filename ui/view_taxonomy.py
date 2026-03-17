
import streamlit as st
import pandas as pd
from config import AppConfig

def get_group_for_species(species: str, taxonomy_service) -> str:
    """Helper function to map species to their Taxonomic Class."""

    lineage = taxonomy_service.fetch_taxonomy_lineage(species)
    if isinstance(lineage, list):
        for node in lineage:
            if node.get("name") in AppConfig.TARGET_TAXA:
                return node.get("name")
        for node in lineage:
            if node.get("rank") == "class":
                return node.get("name")
    return "Other / Unknown"

def render_taxonomy_view(taxonomy_service, df: pd.DataFrame):
    """Renders the Taxonomy tab, allowing side-by-side comparison of species lineages."""

    st.subheader("Species taxonomy")
    st.markdown("Explore and compare the evolutionary lineage of the species catalogued in the platform.")

    if df.empty:
        st.warning("No data available to display taxonomy.")
        return

    if "Class" not in df.columns:
        with st.spinner("Classifying species for taxonomy..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)

    filtered_df = df.copy()

    col_title, col_btn = st.columns([4, 1])
    with col_title:
        st.markdown("### Search filters")
    with col_btn:
        if st.button("Clear filters", use_container_width=True):
            for key in ["tax_f_class", "tax_f_specie"]:
                st.session_state[key] = []

    col1, col2 = st.columns(2)

    with col1:
        groups = sorted([g for g in df["Class"].unique() if g != "Other / Unknown"])
        if "Other / Unknown" in df["Class"].values:
            groups.append("Other / Unknown")
        
        selected_groups = st.multiselect(
            "Select taxonomic class:", 
            options=groups, 
            default=[], 
            placeholder="Select one or more classes...", 
            key="tax_f_class"
        )
        if selected_groups:
            filtered_df = filtered_df[filtered_df["Class"].isin(selected_groups)]

    with col2:
        available_species = sorted(filtered_df["Specie"].dropna().unique().tolist())
        selected_species_list = st.multiselect(
            "Select species:", 
            options=available_species, 
            default=[], 
            placeholder="Select one or more species...", 
            key="tax_f_specie"
        )

    if selected_species_list:
        st.divider()
        st.markdown("### Taxonomic lineage")
        
        MAX_COLS_PER_ROW = 3 
        
        chunks = [selected_species_list[i:i + MAX_COLS_PER_ROW] for i in range(0, len(selected_species_list), MAX_COLS_PER_ROW)]
        
        for chunk in chunks:
            cols = st.columns(len(chunk)) 
            
            for idx, species in enumerate(chunk):
                with cols[idx]:
                    with st.container(border=True):
                        st.markdown(f"<h4 style='text-align: center; color: #4CAF50; margin-bottom: 0;'><i>{species}</i></h4>", unsafe_allow_html=True)
                        st.divider()
                        
                        with st.spinner(f"Consulting NCBI for {species}..."):
                            lineage_data = taxonomy_service.fetch_taxonomy_lineage(species)

                        if isinstance(lineage_data, list):
                            list_html = "<div style='font-family: sans-serif; line-height: 1.8; font-size: 0.9em;'>"
                            
                            for item in lineage_data:
                                rank = item.get("rank", "no rank").capitalize()
                                if rank == "No rank": rank = "Clade"
                                name = item.get("name", "Unknown")
                                
                                list_html += f"<div><span style='color: #888; display: inline-block; width: 95px; font-weight: 600;'>{rank}</span> <span>{name}</span></div>"
                            
                            list_html += f"<div style='margin-top: 8px; padding-top: 8px; border-top: 1px solid #ddd;'><span style='color: #888; display: inline-block; width: 95px; font-weight: 600;'>Species</span> <span style='color: #4CAF50;'><b><i>{species}</i></b></span></div>"
                            list_html += "</div>"
                            
                            st.markdown(list_html, unsafe_allow_html=True)
                        else:
                            st.error(f"Taxonomy not found.")
            
            st.markdown("<br>", unsafe_allow_html=True)