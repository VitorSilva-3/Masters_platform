
import streamlit as st
import pandas as pd

def render_data_tab(df: pd.DataFrame):
    """
    Renders the data table and filters in the UI.
    Receives the main DataFrame loaded in app.py.
    """
    st.subheader("Main dataset overview")
    
    if df.empty:
        st.warning("The DataFrame is empty. Please run the data builder script.")
        return

    col1, col2, col3 = st.columns(3)
    
    all_organisms = sorted(df["Specie"].unique().tolist())
    all_enzymes = sorted(df["Enzyme"].unique().tolist())

    with col1:
        selected_org = st.multiselect("Filter by species", options=all_organisms)
    
    with col2:
        selected_enz = st.multiselect("Filter by enzyme", options=all_enzymes)

    with col3:
        selected_status = st.multiselect(
            "Filter by status", 
            options=["🟢 Confirmed", "🟡 Probable/Hypothetical", "🔴 Uncharacterized"]
        )
        
    filtered_df = df.copy()
    
    if selected_org:
        filtered_df = filtered_df[filtered_df["Specie"].isin(selected_org)]
    
    if selected_enz:
        filtered_df = filtered_df[filtered_df["Enzyme"].isin(selected_enz)]

    if selected_status:
        filtered_df = filtered_df[filtered_df["Status"].isin(selected_status)]
        
    st.markdown(f"**Results:** {len(filtered_df)} records found.")
    st.markdown(f"**Unique species:** {filtered_df['Specie'].nunique()}")
    
    display_df = filtered_df.copy()

    if "Source ID (NCBI)" in display_df.columns:
        display_df["Source ID (NCBI)"] = "https://www.ncbi.nlm.nih.gov/protein/" + display_df["Source ID (NCBI)"].astype(str)

    column_cfg = {
        "Source ID (NCBI)": st.column_config.LinkColumn(
            "NCBI accession",
            help="Direct link to the protein entry on NCBI.",
            display_text="https://www.ncbi.nlm.nih.gov/protein/(.*)", 
            width="small"
        ),
    }

    st.dataframe(
        display_df, 
        column_config=column_cfg,
        use_container_width=True,
        hide_index=True
    )

    st.divider()

    st.markdown("""
        ### Definitions of terms frequently found in NCBI annotations:

        > **'Unclassified' / 'sp.'** (e.g., "*Tetraselmis sp.*", "unclassified Leptolyngbya"):
            Refers to organisms that have been isolated and sequenced but lack a formal species description.
        
        > **'Multispecies'** (e.g., "MULTISPECIES: alpha-amylase family glycosyl hydrolase [Spirulina]"):
            Indicates that the protein sequence is 100% identical across multiple species within the same genus. This signifies a highly conserved enzyme structure.
      
        > **'Uncultured'** (e.g., "uncultured *Leptolyngbya sp.*"):
            Sequences derived directly from environmental DNA (metagenomics) without laboratory cultivation of the organism.
            While these sequences demonstrate the genetic potential for enzymatic activity in nature, the specific organism is not currently available in culture collections for industrial application.

        > **'Candidatus'** (e.g., "Candidatus *Synechococcus spongiarum*"):
            Refers to well-characterized taxonomic lineages that remain uncultured or cannot be maintained in pure culture, limiting their immediate biotechnological viability.
        """)