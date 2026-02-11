
import streamlit as st
import pandas as pd

def render_data_tab(df: pd.DataFrame):
    """
    Renders the data table and filters in the UI.
    Receives the main DataFrame loaded in app.py.
    """
    st.header("Main dataset overview")
    
    if df.empty:
        st.warning("The DataFrame is empty. Please run the data builder script.")
        return

    col1, col2, col3 = st.columns(3)
    
    all_organisms = sorted(df["Specie"].unique().tolist())
    all_enzymes = sorted(df["Enzyme"].unique().tolist())

    with col1:
        selected_org = st.multiselect("Filter by specie", options=all_organisms)
    
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
    st.markdown(f"**Different species:** {filtered_df['Specie'].nunique()}")
    
    display_df = filtered_df.copy()

    if "Source ID (NCBI)" in display_df.columns:
        display_df = display_df.drop(columns=["Source ID (NCBI)"])
    
    st.dataframe(
        display_df, 
        use_container_width=True,
        hide_index=True
    )