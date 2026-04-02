
import os
import streamlit as st
import pandas as pd
import plotly.express as px
from config import AppConfig

def get_group_for_species(species: str, taxonomy_service) -> str:
    """Helper function to map species to their taxonomic class."""

    lineage = taxonomy_service.fetch_taxonomy_lineage(species)
    if isinstance(lineage, list):
        for node in lineage:
            if node.get("name") in AppConfig.TARGET_TAXA:
                return node.get("name")
        for node in lineage:
            if node.get("rank") == "class":
                return node.get("name")
    return "Other / Unknown"

@st.cache_data(show_spinner=False)
def load_and_merge_predictions(filepath: str, main_df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """Loads prediction data, cleans it, standardizes column names, and merges with the main dataframe."""

    if not os.path.exists(filepath):
        return pd.DataFrame()
    
    try:
        df_pred = pd.read_csv(filepath)
        df_pred = df_pred.loc[:, ~df_pred.columns.str.contains('^Unnamed')]
        
        if 'Localizations' in df_pred.columns:
            df_pred = df_pred.rename(columns={'Localizations': 'Localization'})
            
        if not main_df.empty and "Source ID (NCBI)" in main_df.columns:
            main_subset = main_df[['Source ID (NCBI)', 'Specie', 'Enzyme']].dropna(subset=['Source ID (NCBI)']).drop_duplicates(subset=['Source ID (NCBI)'])
            
            df_pred = pd.merge(
                df_pred, 
                main_subset, 
                left_on=id_col, 
                right_on='Source ID (NCBI)', 
                how='left'
            )
                    
        return df_pred
    except Exception as e:
        st.error(f"Error loading {filepath}: {e}")
        return pd.DataFrame()

def render_dashboard_tab(df: pd.DataFrame, id_col: str, title: str, taxonomy_service):
    """Renders a dashboard for a specific domain with dynamic filtering."""

    if df.empty:
        st.info(f"No prediction data available for {title}.")
        return

    if "Class" not in df.columns:
        with st.spinner("Classifying species for filters..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)

    filtered_df = df.copy()

    col_f1, col_f2, col_f3 = st.columns(3)
    
    with col_f1:
        groups = sorted([g for g in df["Class"].dropna().unique() if g != "Other / Unknown"])
        if "Other / Unknown" in df["Class"].values:
            groups.append("Other / Unknown")
        selected_groups = st.multiselect(f"Select class:", options=groups, key=f"loc_class_{title}")
        if selected_groups:
            filtered_df = filtered_df[filtered_df["Class"].isin(selected_groups)]

    with col_f2:
        available_species = sorted(filtered_df["Specie"].dropna().unique().tolist())
        selected_species = st.multiselect(f"Select species:", options=available_species, key=f"loc_spec_{title}")
        if selected_species:
            filtered_df = filtered_df[filtered_df["Specie"].isin(selected_species)]

    with col_f3:
        available_enzymes = sorted(filtered_df["Enzyme"].dropna().unique().tolist())
        selected_enzymes = st.multiselect(f"Select enzyme:", options=available_enzymes, key=f"loc_enz_{title}")
        if selected_enzymes:
            filtered_df = filtered_df[filtered_df["Enzyme"].isin(selected_enzymes)]

    st.divider()

    if filtered_df.empty:
        st.warning("No data matches the selected filters.")
        return

    total_seqs = len(filtered_df)
    top_loc = filtered_df['Localization'].value_counts().idxmax()
    top_loc_pct = (filtered_df['Localization'].value_counts().max() / total_seqs) * 100

    col1, col2, col3 = st.columns(3)
    col1.metric("Filtered sequences", f"{total_seqs}")
    col2.metric("Predominant localization", top_loc)
    col3.metric("Predominant frequency", f"{top_loc_pct:.1f}%")

    st.divider()

    st.markdown("### Localization distribution")
    
    loc_counts = filtered_df['Localization'].value_counts().reset_index()
    loc_counts.columns = ['Subcellular localization', 'Count']
    
    fig = px.pie(
        loc_counts, 
        values='Count', 
        names='Subcellular localization', 
        hole=0.4, 
        color_discrete_sequence=px.colors.qualitative.Pastel
    )
    fig.update_layout(margin=dict(t=20, b=20, l=0, r=0))
    
    st.plotly_chart(fig, use_container_width=True)

    st.divider()

    st.markdown("### Probability scores")
    st.markdown(f"Showing {len(filtered_df)} records based on your current filters.")
    
    base_exclude = [id_col, 'Localization', 'Signals', 'Source ID (NCBI)', 'Specie', 'Enzyme', 'Class']
    
    prob_cols = []
    extra_text_cols = []
    
    for c in filtered_df.columns:
        if c not in base_exclude:
            if pd.api.types.is_numeric_dtype(filtered_df[c]):
                prob_cols.append(c)
            else:
                extra_text_cols.append(c)
    
    display_cols = ['Specie', 'Enzyme', 'Localization'] + extra_text_cols + prob_cols
    display_df = filtered_df[display_cols]
    
    styled_df = display_df.style.background_gradient(
        subset=prob_cols, 
        cmap="Greens", 
        vmin=0, 
        vmax=1
    ).format({col: "{:.4f}" for col in prob_cols})
    
    col_config = {
        "Specie": st.column_config.TextColumn("Species", width="medium"),
        "Enzyme": st.column_config.TextColumn("Enzyme", width="medium"),
        "Localization": st.column_config.TextColumn("Final prediction", width="medium")
    }
    
    for col in extra_text_cols:
        col_config[col] = st.column_config.TextColumn(str(col).replace("_", " "), width="medium")
        
    st.dataframe(
        styled_df,
        use_container_width=True,
        hide_index=True,
        column_config=col_config
    )

def render_localization_view(df: pd.DataFrame, taxonomy_service):
    """Main function to render the subcellular localization UI."""

    st.markdown("""
    Explore subcellular localization predictions generated by deep learning models developed by the **Technical University of Denmark (DTU)**. 
    
    **DeepLoc 2.0**: Used for predicting the subcellular localization of eukaryotic protein sequences.          
    **DeepLocPro**: Used for predicting the subcellular localization of prokaryotic (archaea and bacteria) protein sequences.           
    
    *Credits: models and base architectures are provided by DTU Health Tech.*
    """)

    euk_path = os.path.join("data", "euk_enzymes_predictions.csv")
    prok_path = os.path.join("data", "prok_enzymes_predictions.csv")

    with st.spinner("Loading predictions..."):
        df_euk = load_and_merge_predictions(euk_path, df, id_col="Protein_ID")
        df_prok = load_and_merge_predictions(prok_path, df, id_col="ACC")

    tab1, tab2 = st.tabs(["Eukaryotes (microalgae)", "Prokaryotes (cyanobacteria)"])

    with tab1:
        render_dashboard_tab(df_euk, id_col="Protein_ID", title="Eukaryotes", taxonomy_service=taxonomy_service)

    with tab2:
        render_dashboard_tab(df_prok, id_col="ACC", title="Prokaryotes", taxonomy_service=taxonomy_service)