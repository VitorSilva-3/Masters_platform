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
def load_and_merge_all_predictions(dataset_name: str, main_df: pd.DataFrame) -> pd.DataFrame:
    """Loads Euk and Prok predictions, standardizes them, and merges with the main dataframe."""
    
    euk_path = os.path.join("data", f"euk_{dataset_name}_predictions.csv")
    prok_path = os.path.join("data", f"prok_{dataset_name}_predictions.csv")
    
    dfs = []
    
    if os.path.exists(euk_path):
        df_euk = pd.read_csv(euk_path)
        df_euk = df_euk.loc[:, ~df_euk.columns.str.contains('^Unnamed')]
        if 'Localizations' in df_euk.columns: 
            df_euk = df_euk.rename(columns={'Localizations': 'Localization'})
        if 'Protein_ID' in df_euk.columns: 
            df_euk = df_euk.rename(columns={'Protein_ID': 'Model_ID'})
        df_euk['Domain'] = 'Eukaryotes (Microalgae)'
        dfs.append(df_euk)
        
    if os.path.exists(prok_path):
        df_prok = pd.read_csv(prok_path)
        df_prok = df_prok.loc[:, ~df_prok.columns.str.contains('^Unnamed')]
        if 'Localizations' in df_prok.columns: 
            df_prok = df_prok.rename(columns={'Localizations': 'Localization'})
        if 'ACC' in df_prok.columns: 
            df_prok = df_prok.rename(columns={'ACC': 'Model_ID'})
        df_prok['Domain'] = 'Prokaryotes (Cyanobacteria)'
        dfs.append(df_prok)

    if not dfs:
        return pd.DataFrame()
        
    all_preds = pd.concat(dfs, ignore_index=True)
    all_preds['Model_ID'] = all_preds['Model_ID'].astype(str).str.strip()
    
    item_col = "Enzyme" if "Enzyme" in main_df.columns else "Transporter"
    
    if not main_df.empty and "Source ID (NCBI)" in main_df.columns:
        main_subset = main_df[['Source ID (NCBI)', 'Specie', item_col]].dropna(subset=['Source ID (NCBI)']).drop_duplicates(subset=['Source ID (NCBI)'])
        
        merged_df = pd.merge(
            all_preds, 
            main_subset, 
            left_on='Model_ID', 
            right_on='Source ID (NCBI)', 
            how='inner'
        )
        return merged_df
        
    return pd.DataFrame()

def render_dashboard_tab(df: pd.DataFrame, item_col: str, tab_id: str, taxonomy_service):
    """Renders a dashboard dynamically for Enzymes or Transporters."""

    if df.empty:
        st.info(f"No prediction data available.")
        return

    domain_choice = st.radio(
        "Select biological domain to analyze:",
        ["Eukaryotes (Microalgae)", "Prokaryotes (Cyanobacteria)"],
        horizontal=True,
        key=f"radio_domain_{tab_id}"
    )

    working_df = df[df['Domain'] == domain_choice]

    if working_df.empty:
        st.warning(f"No data found for {domain_choice} in this category.")
        return

    if "Class" not in working_df.columns:
        with st.spinner("Classifying species for filters..."):
            unique_species = working_df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            working_df["Class"] = working_df["Specie"].map(species_to_group)

    key_class = f"loc_class_{tab_id}"
    key_spec = f"loc_spec_{tab_id}"
    key_item = f"loc_item_{tab_id}"
    
    filter_keys = [key_class, key_spec, key_item]
    for k in filter_keys:
        if k not in st.session_state:
            st.session_state[k] = []

    def clear_loc_filters():
        for k in filter_keys:
            st.session_state[k] = []

    def get_filtered_df(exclude_col=None):
        temp_df = working_df.copy()
        if exclude_col != "Class" and st.session_state[key_class]:
            temp_df = temp_df[temp_df["Class"].isin(st.session_state[key_class])]
        if exclude_col != "Specie" and st.session_state[key_spec]:
            temp_df = temp_df[temp_df["Specie"].isin(st.session_state[key_spec])]
        if exclude_col != item_col and st.session_state[key_item]:
            temp_df = temp_df[temp_df[item_col].isin(st.session_state[key_item])]
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
        st.button("Clear filters", key=f"clear_btn_{tab_id}", on_click=clear_loc_filters, use_container_width=True)

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
        opts_item = sorted(get_safe_options(item_col, key_item))
        st.multiselect(f"Select {item_col.lower()}:", options=opts_item, key=key_item)

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
    st.success(f"Showing {len(filtered_df)} records based on your current filters.")
    
    base_exclude = ['Model_ID', 'Localization', 'Signals', 'Source ID (NCBI)', 'Specie', item_col, 'Class', 'Domain']
    
    prob_cols = []
    extra_text_cols = []
    
    for c in filtered_df.columns:
        if c not in base_exclude:
            if pd.api.types.is_numeric_dtype(filtered_df[c]) and filtered_df[c].notna().any():
                prob_cols.append(c)
            elif not pd.api.types.is_numeric_dtype(filtered_df[c]) and filtered_df[c].notna().any():
                extra_text_cols.append(c)
    
    display_cols = ['Specie', item_col, 'Localization'] + extra_text_cols + prob_cols
    display_df = filtered_df[display_cols]
    
    styled_df = display_df.style.background_gradient(
        subset=prob_cols, 
        cmap="Greens", 
        vmin=0, 
        vmax=1
    ).format({col: "{:.4f}" for col in prob_cols})
    
    col_config = {
        "Specie": st.column_config.TextColumn("Species", width="medium"),
        item_col: st.column_config.TextColumn(item_col, width="medium"),
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

def render_localization_view(df_enzymes: pd.DataFrame, df_transporters: pd.DataFrame, taxonomy_service):
    """Main function to render the subcellular localization UI with 2 main Tabs."""

    st.markdown("""
    Explore subcellular localization predictions generated by deep learning models developed by the **Technical University of Denmark (DTU)**. 
    
    **DeepLoc 2.0**: Used for predicting the subcellular localization of eukaryotic protein sequences.          
    **DeepLocPro**: Used for predicting the subcellular localization of prokaryotic (archaea and bacteria) protein sequences.           
    
    *Credits: models and base architectures are provided by DTU Health Tech.*
    """)

    with st.spinner("Loading AI predictions..."):
        preds_enzymes = load_and_merge_all_predictions(dataset_name="enzymes", main_df=df_enzymes)
        preds_transporters = load_and_merge_all_predictions(dataset_name="transporters", main_df=df_transporters)

    tab_enz, tab_trans = st.tabs(["Enzymes", "Transporters"])

    with tab_enz:
        render_dashboard_tab(preds_enzymes, item_col="Enzyme", tab_id="enz", taxonomy_service=taxonomy_service)

    with tab_trans:
        render_dashboard_tab(preds_transporters, item_col="Transporter", tab_id="tra", taxonomy_service=taxonomy_service)