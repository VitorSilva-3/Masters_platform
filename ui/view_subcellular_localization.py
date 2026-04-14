
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
                how='inner'
            )
                    
        return df_pred
    except Exception as e:
        st.error(f"Error loading {filepath}: {e}")
        return pd.DataFrame()

def render_dashboard_tab(df: pd.DataFrame, id_col: str, title: str, taxonomy_service):
    """Renders a dashboard for a specific domain."""

    if df.empty:
        st.info(f"No prediction data available for {title}.")
        return

    if "Class" not in df.columns:
        with st.spinner("Classifying species for filters..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)


    key_class = f"loc_class_{title}"
    key_spec = f"loc_spec_{title}"
    key_enz = f"loc_enz_{title}"
    
    filter_keys = [key_class, key_spec, key_enz]
    
    for k in filter_keys:
        if k not in st.session_state:
            st.session_state[k] = []

    def clear_loc_filters():
        for k in filter_keys:
            st.session_state[k] = []

    def get_filtered_df(exclude_col=None):
        temp_df = df.copy()
        if exclude_col != "Class" and st.session_state[key_class]:
            temp_df = temp_df[temp_df["Class"].isin(st.session_state[key_class])]
        if exclude_col != "Specie" and st.session_state[key_spec]:
            temp_df = temp_df[temp_df["Specie"].isin(st.session_state[key_spec])]
        if exclude_col != "Enzyme" and st.session_state[key_enz]:
            temp_df = temp_df[temp_df["Enzyme"].isin(st.session_state[key_enz])]
        return temp_df

    def get_safe_options(col_name, state_key):
        opts = get_filtered_df(col_name)[col_name].dropna().unique().tolist()
        for sel in st.session_state[state_key]:
            if sel not in opts:
                opts.append(sel)
        return opts

    col_title, col_btn = st.columns([4, 1])
    with col_title:
        st.write("") 
    with col_btn:
        st.button("Clear filters", key=f"clear_btn_{title}", on_click=clear_loc_filters, use_container_width=True)

    col_f1, col_f2, col_f3 = st.columns(3)
    
    with col_f1:
        opts_class = get_safe_options("Class", key_class)
        groups = sorted([g for g in opts_class if g != "Other / Unknown"])
        if "Other / Unknown" in opts_class:
            groups.append("Other / Unknown")
        st.multiselect("Select class:", options=groups, key=key_class)

    with col_f2:
        opts_specie = sorted(get_safe_options("Specie", key_spec))
        st.multiselect("Select species:", options=opts_specie, key=key_spec)

    with col_f3:
        opts_enzyme = sorted(get_safe_options("Enzyme", key_enz))
        st.multiselect("Select enzyme:", options=opts_enzyme, key=key_enz)

    filtered_df = get_filtered_df(exclude_col=None)

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