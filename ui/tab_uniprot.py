
import streamlit as st
import pandas as pd

def render_uniprot_tab(uniprot_service, df: pd.DataFrame):
    """
    Renders the UniProt tab, showing general biochemical properties of the enzymes.
    """

    st.subheader("UniProt database")
    st.markdown("Explore general biological and chemical properties of the enzymes catalogued in the platform.")

    if df.empty:
        st.warning("No data available.")
        return

    unique_enzymes = df[['Enzyme', 'EC number']].drop_duplicates().sort_values(by="Enzyme")
    
    enzyme_dict = {
        f"{row['Enzyme']} (EC {row['EC number']})": (row['Enzyme'], row['EC number']) 
        for _, row in unique_enzymes.iterrows()
    }

    selected_label = st.selectbox(
        "Select an enzyme to view its general profile:",
        options=list(enzyme_dict.keys()),
        index=0
    )

    if selected_label:
        st.divider()
        enzyme_name, ec_number = enzyme_dict[selected_label]
        
        with st.spinner("Consulting UniProt database..."):
            data = uniprot_service.fetch_enzyme_data(enzyme_name, ec_number)
        
        if not data or 'general_info' not in data:
            st.info(f"No detailed information found in UniProt for {enzyme_name} (EC {ec_number}).")
        else:
            info = data['general_info']
            
            col1, col2 = st.columns([3, 1])
            with col1:
                st.markdown(f"### {info.get('protein_name', enzyme_name)}")
                st.caption(f"**Target enzyme:** {enzyme_name} | **EC number:** {ec_number}")
            with col2:
                if 'uniprot_link' in data:
                    st.link_button("View on UniProt", data['uniprot_link'])
            
            st.markdown(f"*(Profile summarized from the top **{data.get('total_entries_analyzed', 0)}** UniProt entries)*")
            st.markdown("---")
            
            def display_list(title, items):
                if items:
                    st.markdown(f"#### {title}")
                    for item in items:
                        st.markdown(f"- {item}")
                    st.write("") 
            
            c1, c2 = st.columns(2)
            
            with c1:
                display_list("Biological function", info.get('functions', []))
                display_list("Pathways", info.get('pathways', []))
                display_list("Subunit structure", info.get('subunit_structure', []))
                
            with c2:
                display_list("Catalytic activity", info.get('catalytic_activities', []))
                display_list("Cofactors", info.get('cofactors', []))