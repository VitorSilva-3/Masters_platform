
import streamlit as st
import pandas as pd

def render_articles_tab(pubmed_service, df: pd.DataFrame):
    """
    Renders the Articles tab in the UI, displaying relevant literature from PubMed.
    """

    st.subheader("Scientific literature (PubMed)")
    st.markdown("Discover the most relevant scientific articles linking specific species to enzymes.")

    if df.empty:
        st.warning("No data available to search for articles.")
        return

    species_list = sorted(df["Specie"].dropna().unique().tolist())

    col1, col2 = st.columns(2)
    
    with col1:
        selected_species = st.selectbox(
            "1. Select a species:",
            options=species_list,
            index=0
        )

    with col2:
        if selected_species:
            species_df = df[df['Specie'] == selected_species]
            
            unique_enzymes = species_df[['Enzyme', 'EC number']].drop_duplicates()
            
            enzyme_dict = {
                f"{row['Enzyme']} (EC {row['EC number']})": (row['Enzyme'], row['EC number']) 
                for _, row in unique_enzymes.iterrows()
            }
            
            selected_enzyme_label = st.selectbox(
                "2. Select an enzyme to search:",
                options=list(enzyme_dict.keys()),
                index=0
            )

    if selected_species and selected_enzyme_label:
        st.divider()
        
        selected_enzyme, selected_ec = enzyme_dict[selected_enzyme_label]
        
        st.markdown(f"### Relevant articles for ***{selected_species}*** and **{selected_enzyme}**")
        
        with st.spinner("Consulting PubMed database..."):
            articles = pubmed_service.search_articles(selected_species, selected_enzyme, selected_ec)
            
        if not articles:
            st.info("No relevant articles were found in PubMed for this specific combination.")
        else:
            st.caption(f"Showing top {len(articles)} most relevant results.")
            
            for i, article in enumerate(articles):
                with st.container():
                    pmid = article.get('pmid', '')
                    title = article.get('title', 'Untitled Article')
                    
                    if pmid:
                        st.markdown(f"#### [{title}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
                    else:
                        st.markdown(f"#### {title}")
                    
                    authors = article.get('authors', 'Unknown Authors')
                    journal = article.get('journal', 'Unknown Journal')
                    year = article.get('year', 'N/A')
                    
                    st.markdown(f"**Authors:** {authors}")
                    st.markdown(f"**Journal:** *{journal}* ({year})")
                    
                    abstract = article.get('abstract', 'No abstract available for this article.')
                    with st.expander("Read Abstract"):
                        st.write(abstract)
                        
                        if pmid:
                            st.link_button(f"View on PubMed (PMID: {pmid})", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
                    
                    if i < len(articles) - 1:
                        st.markdown("---")