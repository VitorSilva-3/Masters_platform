
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

def render_literature_tab(df: pd.DataFrame, item_col: str, id_col: str, tab_id: str, pubmed_service, taxonomy_service):
    """Renders the literature search for a specific domain (Enzymes or Transporters)."""
    
    if df.empty:
        st.info("No data available.")
        return

    if "Class" not in df.columns:
        with st.spinner("Classifying species for literature..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)

    key_class = f"lit_class_{tab_id}"
    key_spec = f"lit_spec_{tab_id}"

    for k in [key_class, key_spec]:
        if k not in st.session_state:
            st.session_state[k] = []

    def clear_lit_filters():
        st.session_state[key_class] = []
        st.session_state[key_spec] = []

    def get_filtered_df(exclude_col=None):
        temp_df = df.copy()
        if exclude_col != "Class" and st.session_state[key_class]:
            temp_df = temp_df[temp_df["Class"].isin(st.session_state[key_class])]
        if exclude_col != "Specie" and st.session_state[key_spec]:
            temp_df = temp_df[temp_df["Specie"].isin(st.session_state[key_spec])]
        return temp_df

    def get_safe_options(col_name, state_key):
        opts = get_filtered_df(col_name)[col_name].dropna().unique().tolist()
        for sel in st.session_state[state_key]:
            if sel not in opts:
                opts.append(sel)
        return sorted(opts)

    col_title, col_btn = st.columns([4, 1])
    with col_title:
        st.markdown("### Search filters")
    with col_btn:
        st.button("Clear filters", key=f"btn_clear_{tab_id}", use_container_width=True, on_click=clear_lit_filters)

    col1, col2 = st.columns(2)

    with col1:
        opts_class = get_safe_options("Class", key_class)
        groups = [g for g in opts_class if g != "Other / Unknown"]
        if "Other / Unknown" in opts_class:
            groups.append("Other / Unknown")
        
        st.multiselect(
            "Select taxonomic class:", 
            options=groups,
            placeholder="Select one or more classes...", 
            key=key_class
        )

    with col2:
        opts_specie = get_safe_options("Specie", key_spec)
        st.multiselect(
            "Select species:", 
            options=opts_specie,  
            placeholder="Select one or more species...", 
            key=key_spec
        )

    selected_species_list = st.session_state[key_spec]

    if selected_species_list:
        st.divider()
        st.markdown("### Search results")
        
        all_articles = []
        seen_pmids = set()
        
        for current_species in selected_species_list:
            species_df = df[df['Specie'] == current_species]
            
            unique_items = species_df[[item_col, id_col]].dropna().drop_duplicates()
            items_to_search = [(row[item_col], row[id_col]) for _, row in unique_items.iterrows()]

            with st.spinner(f"Fetching literature for {current_species}..."):
                for p_name, p_id in items_to_search:
                    arts = pubmed_service.search_articles(current_species, p_name, p_id)
                    for art in arts:
                        pmid = art.get('pmid')
                        if not pmid or pmid not in seen_pmids:
                            all_articles.append(art)
                            if pmid:
                                seen_pmids.add(pmid)

        if not all_articles:
            st.info("No relevant literature found in PubMed for the selected criteria.")
            return

        st.metric("Total articles found", len(all_articles))
        
        sort_recent = st.checkbox("Sort by newest first", value=True, key=f"sort_{tab_id}")
        if sort_recent:
            all_articles.sort(key=lambda x: str(x.get('year', '0')), reverse=True)
            
        st.markdown("<br>", unsafe_allow_html=True)
        
        for article in all_articles:
            with st.container(border=True): 
                pmid = article.get('pmid', '')
                title = article.get('title', 'Untitled')
                
                if pmid and pmid != 'N/A':
                    st.markdown(f"#### [{title}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
                else:
                    st.markdown(f"#### {title}")
                
                st.markdown(f"**Authors:** {article.get('authors')}  \n**Journal:** *{article.get('journal')}* ({article.get('year')})")
                
                with st.expander("Read abstract"):
                    st.write(article.get('abstract'))
                    if pmid and pmid != 'N/A':
                        st.link_button("View on PubMed", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")

def render_literature_view(pubmed_service, df_enzymes: pd.DataFrame, df_transporters: pd.DataFrame, taxonomy_service):
    """Renders the literature view."""

    st.markdown("Explore scientific literature from PubMed.")

    df_e = df_enzymes.copy()
    if not df_e.empty:
        df_e['Protein_Name'] = df_e['Enzyme']
        df_e['Protein_ID'] = df_e['EC number']
        df_e['Type'] = 'Enzyme'

    df_t = df_transporters.copy()
    if not df_t.empty:
        df_t['Protein_Name'] = df_t['Transporter']
        df_t['Protein_ID'] = df_t['TC number']
        df_t['Type'] = 'Transporter'

    dfs_to_concat = [df for df in [df_e, df_t] if not df.empty]
    
    if not dfs_to_concat:
        st.warning("No data available.")
        return
        
    df_unified = pd.concat(dfs_to_concat, ignore_index=True)

    render_literature_tab(
        df=df_unified, 
        item_col="Protein_Name", 
        id_col="Protein_ID", 
        tab_id="all", 
        pubmed_service=pubmed_service, 
        taxonomy_service=taxonomy_service
    )