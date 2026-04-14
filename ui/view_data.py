
import streamlit as st
import pandas as pd
from config import AppConfig

def get_group_for_species(species: str, taxonomy_service) -> str:
    """Helper function to find which Target Taxa a species belongs to."""
    
    lineage = taxonomy_service.fetch_taxonomy_lineage(species)
    if isinstance(lineage, list):
        for node in lineage:
            if node.get("name") in AppConfig.TARGET_TAXA:
                return node.get("name")
        for node in lineage:
            if node.get("rank") == "class":
                return node.get("name")
    return "Other / Unknown"

def render_data_view(df: pd.DataFrame, taxonomy_service):
    """
    Renders the Data tab, displaying the main dataset.
    """
    st.subheader("Global dataset")

    if df.empty:
        st.warning("No data available.")
        return

    col_m1, col_m2, col_m3, col_m4 = st.columns(4)
    col_m1.metric("Species", df["Specie"].nunique())
    col_m2.metric("Enzymes", df["Enzyme"].nunique())
    
    if "Target sugar" in df.columns:
        col_m3.metric("Target sugars", df["Target sugar"].nunique())
    else:
        col_m3.metric("Target sugars", "N/A")
        
    col_m4.metric("Total records", len(df))
    st.divider()

    if "Class" not in df.columns:
        with st.spinner("Classifying species..."):
            unique_species = df["Specie"].dropna().unique()
            species_to_group = {sp: get_group_for_species(sp, taxonomy_service) for sp in unique_species}
            df["Class"] = df["Specie"].map(species_to_group)

    cols = df.columns.tolist()
    if "Class" in cols:
        cols.insert(cols.index("Specie") + 1, cols.pop(cols.index("Class")))
        df = df[cols]


    filter_keys = ["f_class", "f_specie", "f_enzyme", "f_sugar", "f_status"]
    for k in filter_keys:
        if k not in st.session_state:
            st.session_state[k] = []

    def clear_data_filters():
        for key in filter_keys:
            st.session_state[key] = []

    def get_filtered_df(exclude_col=None):
        temp_df = df.copy()
        if exclude_col != "Class" and st.session_state["f_class"]:
            temp_df = temp_df[temp_df["Class"].isin(st.session_state["f_class"])]
        if exclude_col != "Specie" and st.session_state["f_specie"]:
            temp_df = temp_df[temp_df["Specie"].isin(st.session_state["f_specie"])]
        if exclude_col != "Enzyme" and st.session_state["f_enzyme"]:
            temp_df = temp_df[temp_df["Enzyme"].isin(st.session_state["f_enzyme"])]
        if exclude_col != "Target sugar" and "Target sugar" in df.columns and st.session_state["f_sugar"]:
            temp_df = temp_df[temp_df["Target sugar"].isin(st.session_state["f_sugar"])]
        if exclude_col != "Status" and "Status" in df.columns and st.session_state["f_status"]:
            temp_df = temp_df[temp_df["Status"].isin(st.session_state["f_status"])]
        return temp_df

    def get_safe_options(col_name, state_key):
        opts = get_filtered_df(col_name)[col_name].dropna().unique().tolist()
        for sel in st.session_state[state_key]:
            if sel not in opts:
                opts.append(sel)
        return opts
    
    col_title, col_btn = st.columns([4, 1])
    with col_title:
        st.markdown("### Filters")
    with col_btn:
        st.button("Clear filters", use_container_width=True, on_click=clear_data_filters)

    col1, col2, col3 = st.columns(3, gap="large")

    with col1:
        opts_class = get_safe_options("Class", "f_class")
        groups = sorted([g for g in opts_class if g != "Other / Unknown"])
        if "Other / Unknown" in opts_class:
            groups.append("Other / Unknown")
        st.multiselect("Select class:", options=groups, placeholder="Select...", key="f_class")

    with col2:
        opts_specie = sorted(get_safe_options("Specie", "f_specie"))
        st.multiselect("Select species:", options=opts_specie, placeholder="Select...", key="f_specie")

    with col3:
        opts_enzyme = sorted(get_safe_options("Enzyme", "f_enzyme"))
        st.multiselect("Select enzyme:", options=opts_enzyme, placeholder="Select...", key="f_enzyme")

    col4, col5 = st.columns(2, gap="large")

    with col4:
        if "Target sugar" in df.columns:
            opts_sugar = sorted(get_safe_options("Target sugar", "f_sugar"))
            st.multiselect("Select target sugar:", options=opts_sugar, placeholder="Select...", key="f_sugar")
        else:
            st.multiselect("Select target sugar:", options=[], placeholder="N/A", disabled=True, key="f_sugar_disabled")

    with col5:
        if "Status" in df.columns:
            opts_status = sorted(get_safe_options("Status", "f_status"))
            st.multiselect("Select status:", options=opts_status, placeholder="Select...", key="f_status")
        else:
            st.multiselect("Select status:", options=[], placeholder="N/A", disabled=True, key="f_status_disabled")

    filtered_df = get_filtered_df(exclude_col=None) 

    st.divider()
    st.markdown(f"Showing **{len(filtered_df)}** records based on your current filters.")

    display_df = filtered_df.copy()
    if "Source ID (NCBI)" in display_df.columns:
        display_df["NCBI_URL"] = "https://www.ncbi.nlm.nih.gov/protein/" + display_df["Source ID (NCBI)"].astype(str)

    st.markdown("<br>", unsafe_allow_html=True)
    
    col_config = {
        "Specie": st.column_config.TextColumn("Species", width="medium"),
        "Source ID (NCBI)": None, 
        "NCBI_URL": st.column_config.LinkColumn(
            "Source ID (NCBI)", 
            display_text=r"https://www\.ncbi\.nlm\.nih\.gov/protein/(.*)" 
        )
    }

    st.dataframe(
        display_df,
        use_container_width=True,
        hide_index=True,
        column_config=col_config
    )