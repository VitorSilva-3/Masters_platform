import streamlit as st
import json
import os
import pandas as pd

def render_feedipedia_view():
    st.title("Agro-industrial waste composition")
    st.markdown("Explore chemical and mineral profiles extracted from the **Feedipedia** database.")

    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(root_dir, "data", "feedipedia_raw_data.json")

    if not os.path.exists(data_path):
        st.error(f"Data file not found at local path: {data_path}. Please execute the acquisition service first.")
        return

    with open(data_path, "r", encoding="utf-8") as f:
        master_data = json.load(f)

    records_map = [
        {
            "main_category": item["Main_Category"],
            "residue_name": item["Residue_Name"],
            "index": idx
        }
        for idx, item in enumerate(master_data)
    ]
    df_structure = pd.DataFrame(records_map)

    st.markdown("### Filters")
    col1, col2= st.columns(2)
    
    with col1:
        main_categories = sorted(df_structure["main_category"].unique())
        selected_main = st.selectbox("Select category", main_categories)
        df_filtered_main = df_structure[df_structure["main_category"] == selected_main]

    with col2:
        residue_names = sorted(df_filtered_main["residue_name"].unique())
        selected_residue = st.selectbox("Select agro-industrial waste", residue_names)

    st.divider()

    target_match = df_filtered_main[df_filtered_main["residue_name"] == selected_residue]
    
    if target_match.empty:
        st.warning("No data found matching the selected criteria.")
        return

    target_index = target_match["index"].values[0]
    selected_record = master_data[target_index]
    
    st.info(f"**Feedipedia registry identifier:** Node ID {selected_record['Node_ID']}")

    
    st.markdown("### Main analysis")
    main_analysis_data = selected_record["Data"].get("Main analysis", {})
    if main_analysis_data:
        df_main = pd.DataFrame.from_dict(main_analysis_data, orient="index")
        st.dataframe(df_main, use_container_width=True)
    else:
        st.info("No nutritional main analysis matrix compiled for this substrate reference.")

    st.write("") 

    st.markdown("### Minerals")
    minerals_data = selected_record["Data"].get("Minerals", {})
    if minerals_data:
        df_minerals = pd.DataFrame.from_dict(minerals_data, orient="index")
        st.dataframe(df_minerals, use_container_width=True)
    else:
        st.info("No specific micro or macro mineral concentrations documented for this substrate.")