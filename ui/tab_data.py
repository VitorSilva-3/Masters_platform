
import streamlit as st
import pandas as pd

TARGET_TAXA = [
    "Chlorophyceae", "Trebouxiophyceae", "Prasinophyceae", "Mamiellophyceae",
    "Bacillariophyceae", "Eustigmatophyceae", "Chrysophyceae", "Xanthophyceae",
    "Porphyridiophyceae", "Cyanidiophyceae", "Rhodellophyceae",
    "Prymnesiophyceae", "Euglenida", "Dinophyceae", "Cryptophyceae",
    "Cyanobacteria"
]

def get_group_for_species(species: str, taxonomy_service) -> str:
    """Helper function to find which Target Taxa a species belongs to."""
    
    lineage = taxonomy_service.fetch_taxonomy_lineage(species)
    
    if isinstance(lineage, list):
        for node in lineage:
            if node.get("name") in TARGET_TAXA:
                return node.get("name")
        
        for node in lineage:
            if node.get("rank") == "class":
                return node.get("name")
                
    return "Other / Unknown"

def render_data_tab(df: pd.DataFrame, taxonomy_service):
    """
    Renders the Data tab, displaying the main dataset with advanced filtering.
    """

    st.subheader("Global dataset")

    if df.empty:
        st.warning("No data available.")
        return

    if "Class" not in df.columns:
        with st.spinner("Classifying species..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)

    cols = df.columns.tolist()
    if "Class" in cols:
        cols.insert(cols.index("Specie") + 1, cols.pop(cols.index("Class")))
        df = df[cols]

    filtered_df = df.copy()

    st.markdown("### Filters")
    
    filter_type = st.radio(
        "Taxonomic level:",
        options=["View all", "Filter by class", "Filter by species"],
        horizontal=True
    )

    if filter_type == "Filter by class":
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            groups = sorted([g for g in df["Class"].unique() if g != "Other / Unknown"])
            if "Other / Unknown" in df["Class"].values:
                groups.append("Other / Unknown")
            
            selected_group = st.selectbox("Select class:", options=["All"] + groups)
            if selected_group != "All":
                filtered_df = filtered_df[filtered_df["Class"] == selected_group]

        with col2:
            available_species = sorted(filtered_df["Specie"].dropna().unique().tolist())
            selected_species = st.selectbox("Select species:", options=["All"] + available_species)
            if selected_species != "All":
                filtered_df = filtered_df[filtered_df["Specie"] == selected_species]

        with col3:
            available_enzymes = sorted(filtered_df["Enzyme"].dropna().unique().tolist())
            selected_enzyme = st.selectbox("Select enzyme:", options=["All"] + available_enzymes)
            if selected_enzyme != "All":
                filtered_df = filtered_df[filtered_df["Enzyme"] == selected_enzyme]

        with col4:
            if "Status" in df.columns:
                available_status = sorted(filtered_df["Status"].dropna().unique().tolist())
                selected_status = st.selectbox("Select status:", options=["All"] + available_status)
                if selected_status != "All":
                    filtered_df = filtered_df[filtered_df["Status"] == selected_status]
            else:
                st.selectbox("Select status:", options=["All"], disabled=True)

    elif filter_type == "Filter by species":
        col1, col2, col3 = st.columns(3)
        
        with col1:
            species_list = sorted(df["Specie"].dropna().unique().tolist())
            selected_species = st.selectbox("Select species:", options=["All"] + species_list)
            if selected_species != "All":
                filtered_df = filtered_df[filtered_df["Specie"] == selected_species]

        with col2:
            available_enzymes = sorted(filtered_df["Enzyme"].dropna().unique().tolist())
            selected_enzyme = st.selectbox("Select enzyme:", options=["All"] + available_enzymes)
            if selected_enzyme != "All":
                filtered_df = filtered_df[filtered_df["Enzyme"] == selected_enzyme]

        with col3:
            if "Status" in df.columns:
                available_status = sorted(filtered_df["Status"].dropna().unique().tolist())
                selected_status = st.selectbox("Select status:", options=["All"] + available_status)
                if selected_status != "All":
                    filtered_df = filtered_df[filtered_df["Status"] == selected_status]
            else:
                st.selectbox("Select status:", options=["All"], disabled=True)

    else:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.selectbox("Select taxonomy:", options=["All Taxa"], disabled=True)
            
        with col2:
            available_enzymes = sorted(filtered_df["Enzyme"].dropna().unique().tolist())
            selected_enzyme = st.selectbox("Select enzyme:", options=["All"] + available_enzymes)
            if selected_enzyme != "All":
                filtered_df = filtered_df[filtered_df["Enzyme"] == selected_enzyme]
                
        with col3:
            if "Status" in df.columns:
                available_status = sorted(filtered_df["Status"].dropna().unique().tolist())
                selected_status = st.selectbox("Select status:", options=["All"] + available_status)
                if selected_status != "All":
                    filtered_df = filtered_df[filtered_df["Status"] == selected_status]
            else:
                st.selectbox("Select status:", options=["All"], disabled=True)

    st.divider()
    
    st.markdown(f"**Showing {len(filtered_df)} records** based on your current filters.")
    st.dataframe(filtered_df, use_container_width=True, hide_index=True)