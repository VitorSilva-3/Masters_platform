
import streamlit as st
import pandas as pd

def render_uniprot_view(uniprot_service, df_enzymes: pd.DataFrame, df_transporters: pd.DataFrame):
    """Renders the UniProt page, showing general biochemical properties in separate tabs."""

    st.markdown("Explore general biochemical properties, functional annotations, and pathways catalogued in **UniProt**.")

    if df_enzymes.empty and df_transporters.empty:
        st.warning("No data available.")
        return

    def render_protein_details(protein_name, identifier, id_type):
        with st.spinner(f"Consulting UniProt general database for {protein_name}..."):
            data = uniprot_service.fetch_protein_data(protein_name, identifier, id_type)
        
        display_id = f"EC {identifier}" if id_type == "EC" else f"Target: {identifier}"

        if not data or 'general_info' not in data:
            st.info(f"No detailed information found in UniProt for {protein_name} ({display_id}).")
        else:
            info = data['general_info']
            
            col_title, col_btn = st.columns([4, 1])
            with col_title:
                st.markdown(f"### {protein_name.capitalize()}") 
            with col_btn:
                if 'uniprot_link' in data:
                    st.link_button("View on UniProt", data['uniprot_link'], use_container_width=True)
            
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

    tab_enzymes, tab_transporters = st.tabs(["Enzymes", "Transporters"])

    with tab_enzymes:
        if df_enzymes.empty:
            st.info("No enzyme data available.")
        else:
            df_e = df_enzymes.copy()
            df_e['Display_Label'] = df_e['Enzyme'] + " (EC " + df_e['EC number'].astype(str) + ")"
            
            unique_e = df_e[['Enzyme', 'EC number', 'Display_Label']].dropna().drop_duplicates().sort_values(by='Display_Label')
            dict_e = {row['Display_Label']: (row['Enzyme'], row['EC number']) for _, row in unique_e.iterrows()}
            
            sel_e = st.selectbox(
                "Select enzyme:",
                options=list(dict_e.keys()),
                index=None,
                placeholder="Search enzyme...",
                key="sel_enzyme"  
            )
            
            if sel_e:
                st.divider()
                p_name, p_id = dict_e[sel_e]
                render_protein_details(p_name, p_id, "EC")

    with tab_transporters:
        if df_transporters.empty:
            st.info("No transporter data available.")
        else:
            df_t = df_transporters.copy()
            df_t['Display_Label'] = df_t['Transporter'].str.capitalize() + " (Target: " + df_t['Target sugar'].astype(str) + ")"
            
            unique_t = df_t[['Transporter', 'Target sugar', 'Display_Label']].dropna().drop_duplicates().sort_values(by='Display_Label')
            dict_t = {row['Display_Label']: (row['Transporter'], row['Target sugar']) for _, row in unique_t.iterrows()}
            
            sel_t = st.selectbox(
                "Select transporter:",
                options=list(dict_t.keys()),
                index=None,
                placeholder="Search transporter...",
                key="sel_transporter" 
            )
            
            if sel_t:
                st.divider()
                p_name, p_id = dict_t[sel_t]
                render_protein_details(p_name, p_id, "NAME")