
import streamlit as st
import pandas as pd

def render_uniprot_view(uniprot_service, df: pd.DataFrame):
    """
    Renders the UniProt tab, showing general biochemical properties and species-specific links.
    """

    st.subheader("UniProt database")
    st.markdown("Explore biochemical properties and direct species-specific database entries.")

    if df.empty:
        st.warning("No data available.")
        return

    unique_enzymes = df[['Enzyme', 'EC number']].dropna().drop_duplicates().sort_values(by="Enzyme")
    
    enzyme_dict = {
        f"{row['Enzyme']} (EC {row['EC number']})": (row['Enzyme'], row['EC number']) 
        for _, row in unique_enzymes.iterrows()
    }

    selected_label = st.selectbox(
        "Select an enzyme to analyze:",
        options=list(enzyme_dict.keys()),
        index=None,
        placeholder="Choose an enzyme..."
    )

    if selected_label:
        st.divider()
        enzyme_name, ec_number = enzyme_dict[selected_label]
        
        with st.spinner("Consulting UniProt general database..."):
            data = uniprot_service.fetch_enzyme_data(enzyme_name, ec_number)
        
        if not data or 'general_info' not in data:
            st.info(f"No detailed information found in UniProt for {enzyme_name} (EC {ec_number}).")
        else:
            info = data['general_info']
            
            col_title, col_btn = st.columns([4, 1])
            with col_title:
                st.markdown(f"### {enzyme_name}") 
            with col_btn:
                if 'uniprot_link' in data:
                    st.link_button("View general entry", data['uniprot_link'], use_container_width=True)
            
            st.markdown(f"*(General profile summarized from the top **{data.get('total_entries_analyzed', 0)}** annotated entries globally)*")
            
            def display_list(title, items):
                if items:
                    st.markdown(f"#### {title}")
                    for item in items:
                        st.markdown(f"- {item}")
                    st.write("") 
            
            with st.container(border=True):
                c1, c2 = st.columns(2)
                with c1:
                    display_list("Biological function", info.get('functions', []))
                    display_list("Pathways", info.get('pathways', []))
                    display_list("Subunit structure", info.get('subunit_structure', []))
                    
                with c2:
                    display_list("Catalytic activity", info.get('catalytic_activities', []))
                    display_list("Cofactors", info.get('cofactors', []))

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### Specie UniProt entries")
        st.markdown(f"Direct database access for **{enzyme_name}** in the specific species catalogued in the dataset.")

        associated_species = df[df['EC number'] == ec_number]['Specie'].dropna().unique()
        
        link_data = []
        
        with st.spinner("Fetching specific UniProt entries..."):
            for sp in associated_species:
                link_info = uniprot_service.fetch_species_link(ec_number, sp)
                link_data.append({
                    "Species": sp, 
                    "Status": "Found" if link_info["status"] == "Found" else "Not mapped",
                    "Accession ID": link_info["accession"],
                    "Database Link": link_info["url"]
                })

        if link_data:
            link_df = pd.DataFrame(link_data)
            st.dataframe(
                link_df, 
                use_container_width=True, 
                hide_index=True,
                column_config={
                    "Species": st.column_config.TextColumn("Microalgae Species", width="medium"),
                    "Status": st.column_config.TextColumn("Status", width="small"),
                    "Accession ID": st.column_config.TextColumn("Accession ID", width="small"),
                    "Database Link": st.column_config.LinkColumn("Database Link", display_text="Open in UniProt ↗", width="medium")
                }
            )