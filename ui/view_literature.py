
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

def render_literature_view(pubmed_service, df: pd.DataFrame, taxonomy_service):
    """Renders the Articles tab with a clean Class-to-Species multiselect funnel."""

    st.subheader("Scientific literature")
    st.markdown("Explore relevant literature.")

    if df.empty:
        st.warning("No data available.")
        return

    if "Class" not in df.columns:
        with st.spinner("Classifying species for literature..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)

    filtered_df = df.copy()

    col_title, col_btn = st.columns([4, 1])
    with col_title:
        st.markdown("### Search filters")
    with col_btn:
        if st.button("Clear filters", use_container_width=True):
            for key in ["lit_f_class", "lit_f_specie"]:
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
            key="lit_f_class"
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
            key="lit_f_specie"
        )
        if selected_species_list:
            filtered_df = filtered_df[filtered_df["Specie"].isin(selected_species_list)]


    if selected_species_list:
        st.divider()
        st.markdown("### Search results")
        
        all_articles = []
        seen_pmids = set()
        
        for current_species in selected_species_list:
            species_df = df[df['Specie'] == current_species]
            unique_enzymes = species_df[['Enzyme', 'EC number']].drop_duplicates()
            enzymes_to_search = [(row['Enzyme'], row['EC number']) for _, row in unique_enzymes.iterrows()]

            with st.spinner(f"Fetching literature for {current_species}..."):
                for enz_name, enz_ec in enzymes_to_search:
                    arts = pubmed_service.search_articles(current_species, enz_name, enz_ec)
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
        
        sort_recent = st.checkbox("Sort by newest first", value=True)
        if sort_recent:
            all_articles.sort(key=lambda x: str(x.get('year', '0')), reverse=True)
            
        st.markdown("<br>", unsafe_allow_html=True)
        
        for i, article in enumerate(all_articles):
            with st.container(border=True): 
                pmid = article.get('pmid', '')
                title = article.get('title', 'Untitled')
                
                if pmid:
                    st.markdown(f"#### [{title}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
                else:
                    st.markdown(f"#### {title}")
                
                st.markdown(f"**Authors:** {article.get('authors')}  \n**Journal:** *{article.get('journal')}* ({article.get('year')})")
                
                with st.expander("Read Abstract"):
                    st.write(article.get('abstract'))
                    if pmid:
                        st.link_button("View on PubMed", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")