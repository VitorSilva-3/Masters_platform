
import streamlit as st
import pandas as pd

def render_taxonomy_tab(taxonomy_service, df: pd.DataFrame):
    """
    Renders the Taxonomy tab in the UI, displaying the lineage of species.
    """

    st.subheader("Species taxonomy")
    st.markdown("Explore the evolutionary lineage of the species catalogued in the platform.")

    if df.empty:
        st.warning("No data available to display taxonomy.")
        return

    species_list = sorted(df["Specie"].dropna().unique().tolist())

    col1, col2 = st.columns([2, 1]) 
    
    with col1:
        selected_species = st.selectbox(
            "Select a species to view its taxonomic lineage:",
            options=species_list,
            index=0
        )

    if selected_species:
        st.divider()
        st.markdown(f"### Taxonomic lineage for ***{selected_species}***")
        
        with st.spinner("Consulting NCBI Taxonomy..."):
            lineage_data = taxonomy_service.fetch_taxonomy_lineage(selected_species)

        error_messages = ["Taxonomy ID not found.", "No details found.", "Lineage not available"]
        
        if isinstance(lineage_data, list):
            
            names_only = [item.get("name", "Unknown") for item in lineage_data]
            formatted_lineage = " ➔ ".join(f"**{name}**" for name in names_only)
            st.info(formatted_lineage)
            
            st.markdown("#### Tree structure")
            tree_html = "<div style='font-family: monospace; line-height: 1.8;'>"
            
            for i, item in enumerate(lineage_data):
                indent = "&nbsp;" * (i * 6) 
                prefix = "↳ " if i > 0 else "🌳 "
                
                rank = item.get("rank", "no rank").capitalize()
                
                if rank == "No rank":
                    rank = "Clade"
                    
                name = item.get("name", "Unknown")
                
                tree_html += f"{indent}{prefix}<span style='color: #888; font-size: 0.85em;'>[{rank}]</span> <b>{name}</b><br>"
            
            final_indent = "&nbsp;" * (len(lineage_data) * 6)
            tree_html += f"{final_indent}↳ <span style='color: #888; font-size: 0.85em;'>[Species]</span> <span style='color: #4CAF50;'><b><i>{selected_species}</i></b></span><br>"
            tree_html += "</div>"
            
            st.markdown(tree_html, unsafe_allow_html=True)
            
        else:
            st.error(f"Could not retrieve taxonomy for *{selected_species}*.")
            st.caption(f"Reason: {lineage_data}")