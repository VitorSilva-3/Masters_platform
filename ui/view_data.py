
import streamlit as st
import pandas as pd
from config import AppConfig
from st_aggrid import AgGrid, GridOptionsBuilder, ColumnsAutoSizeMode, JsCode

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
    Renders the Data tab, displaying the main dataset with advanced filtering.
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

    filtered_df = df.copy()

    def clear_data_filters():
        for key in ["f_class", "f_specie", "f_enzyme", "f_sugar", "f_status"]:
            st.session_state[key] = []

    col_title, col_btn = st.columns([4, 1])
    with col_title:
        st.markdown("### Filters")
    with col_btn:
        st.button("Clear filters", use_container_width=True, on_click=clear_data_filters)

    col1, col2, col3, col4, col5 = st.columns(5)

    with col1:
        groups = sorted([g for g in df["Class"].unique() if g != "Other / Unknown"])
        if "Other / Unknown" in df["Class"].values:
            groups.append("Other / Unknown")
        selected_groups = st.multiselect("Select class:", options=groups, default=[], placeholder="Select...", key="f_class")
        if selected_groups:
            filtered_df = filtered_df[filtered_df["Class"].isin(selected_groups)]

    with col2:
        available_species = sorted(filtered_df["Specie"].dropna().unique().tolist())
        selected_species = st.multiselect("Select species:", options=available_species, default=[], placeholder="Select...", key="f_specie")
        if selected_species:
            filtered_df = filtered_df[filtered_df["Specie"].isin(selected_species)]

    with col3:
        available_enzymes = sorted(filtered_df["Enzyme"].dropna().unique().tolist())
        selected_enzymes = st.multiselect("Select enzyme:", options=available_enzymes, default=[], placeholder="Select...", key="f_enzyme")
        if selected_enzymes:
            filtered_df = filtered_df[filtered_df["Enzyme"].isin(selected_enzymes)]

    with col4:
        if "Target sugar" in df.columns:
            available_sugars = sorted(filtered_df["Target sugar"].dropna().unique().tolist())
            selected_sugars = st.multiselect("Select target sugar:", options=available_sugars, default=[], placeholder="Select...", key="f_sugar")
            if selected_sugars:
                filtered_df = filtered_df[filtered_df["Target sugar"].isin(selected_sugars)]
        else:
            st.multiselect("Select target sugar:", options=[], default=[], placeholder="N/A", disabled=True, key="f_sugar_disabled")

    with col5:
        if "Status" in df.columns:
            available_status = sorted(filtered_df["Status"].dropna().unique().tolist())
            selected_statuses = st.multiselect("Select status:", options=available_status, default=[], placeholder="Select...", key="f_status")
            if selected_statuses:
                filtered_df = filtered_df[filtered_df["Status"].isin(selected_statuses)]
        else:
            st.multiselect("Select status:", options=[], default=[], placeholder="N/A", disabled=True, key="f_status_disabled")

    st.divider()
    st.markdown(f"**Showing {len(filtered_df)} records** based on your current filters.")

    gb = GridOptionsBuilder.from_dataframe(filtered_df)
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=15) 
    gb.configure_side_bar() 
    gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, sortable=True, filter=True)
    
    column_ncbi = "Source ID (NCBI)"
    if column_ncbi in filtered_df.columns:
        link_jscode = JsCode("""
        class UrlCellRenderer {
            init(params) {
                this.eGui = document.createElement('a');
                this.eGui.innerText = params.value;
                if (params.value && params.value !== 'nan' && params.value.trim() !== '') {
                    this.eGui.setAttribute('href', 'https://www.ncbi.nlm.nih.gov/protein/' + params.value);
                    this.eGui.setAttribute('target', '_blank');
                    this.eGui.setAttribute('style', 'text-decoration: underline; color: #4DA6FF; font-weight: 500;');
                }
            }
            getGui() {
                return this.eGui;
            }
        }
        """)
        
        gb.configure_column(
            column_ncbi, 
            headerName="Source ID (NCBI)", 
            cellRenderer=link_jscode
        )
         
    gridOptions = gb.build()

    st.markdown("<br>", unsafe_allow_html=True)
    AgGrid(
        filtered_df,
        gridOptions=gridOptions,
        enable_enterprise_modules=True,
        allow_unsafe_jscode=True, 
        columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS,
        theme="streamlit" 
    )