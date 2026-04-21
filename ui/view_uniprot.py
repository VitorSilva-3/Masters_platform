
import streamlit as st
import pandas as pd

def render_uniprot_view(uniprot_service, df: pd.DataFrame, dataset_type: str = "enzymes"):
    """Renders the UniProt page, showing general biochemical properties."""

    st.subheader("UniProt database")
    st.markdown("Explore general biochemical properties, functional annotations, and pathways catalogued in UniProt.")

    if df.empty:
        st.warning("No data available.")
        return

    item_col = "Enzyme" if dataset_type == "enzymes" else "Transporter"
    id_col = "EC number" if dataset_type == "enzymes" else "TC number"
    id_type = "EC" if dataset_type == "enzymes" else "TC"

    unique_items = df[[item_col, id_col]].dropna().drop_duplicates().sort_values(by=item_col)
    
    protein_dict = {
        f"{row[item_col]} ({id_type} {row[id_col]})": (row[item_col], row[id_col]) 
        for _, row in unique_items.iterrows()
    }

    label_noun = "an enzyme" if dataset_type == "enzymes" else "a transporter"
    selected_label = st.selectbox(
        f"Select {label_noun} to analyze:",
        options=list(protein_dict.keys()),
        index=None,
        placeholder=f"Choose {label_noun}..."
    )

    if selected_label:
        st.divider()
        protein_name, identifier = protein_dict[selected_label]
        
        with st.spinner("Consulting UniProt general database..."):
            data = uniprot_service.fetch_protein_data(protein_name, identifier, id_type=id_type)
        
        if not data or 'general_info' not in data:
            st.info(f"No detailed information found in UniProt for {protein_name} ({id_type} {identifier}).")
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