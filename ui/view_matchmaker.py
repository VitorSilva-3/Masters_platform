
import streamlit as st
import pandas as pd
import json
import urllib.parse
import re
from ml.matchmaker import BiotransformationMatchmaker

@st.cache_resource
def load_ai_engine():
    with open("data/feedipedia_raw_data.json", "r", encoding="utf-8") as f:
        raw_json = json.load(f)
        
    if isinstance(raw_json, dict):
        for key, val in raw_json.items():
            if isinstance(val, list):
                raw_json = val
                break
                
    flattened_data = []
    for residue in raw_json:
        nome_residuo = residue.get("Residue_Name", "Unknown")
        flat_row = {"Residue_Name": nome_residuo}
        
        data_section = residue.get("Data", {})
        for category, nutrients in data_section.items():
            if isinstance(nutrients, dict):
                for nut_name, nut_info in nutrients.items():
                    if isinstance(nut_info, dict) and "Avg" in nut_info:
                        try:
                            clean_name = str(nut_name).strip()
                            flat_row[clean_name] = float(nut_info["Avg"])
                        except ValueError:
                            flat_row[clean_name] = 0.0
                            
        flattened_data.append(flat_row)

    df_feedipedia = pd.DataFrame(flattened_data)
    df_feedipedia.set_index("Residue_Name", inplace=True)
    df_feedipedia.fillna(0.0, inplace=True)

    # 3. Min-Max Scaling (Crucial for accurate Cosine Similarity)
    for col in df_feedipedia.columns:
        if pd.api.types.is_numeric_dtype(df_feedipedia[col]):
            max_val = df_feedipedia[col].max()
            min_val = df_feedipedia[col].min()
            if max_val > min_val:
                df_feedipedia[col] = (df_feedipedia[col] - min_val) / (max_val - min_val)
            else:
                df_feedipedia[col] = 0.0

    # 4. Loading Main CSVs and Predictions
    try:
        # Main files (with Biological metadata: Specie, EC, Targets)
        df_enz_main = pd.read_csv("data/enzymes_data.csv")
        df_trans_main = pd.read_csv("data/transporters_data.csv")

        # Prediction files (with Spatial Localization)
        df_euk_enz = pd.read_csv("data/euk_enzymes_predictions.csv")
        df_euk_trans = pd.read_csv("data/euk_transporters_predictions.csv")
        df_prok_enz = pd.read_csv("data/prok_enzymes_predictions.csv")
        df_prok_trans = pd.read_csv("data/prok_transporters_predictions.csv")
    except FileNotFoundError as e:
        st.error(f"Missing file in the data/ folder: {e}")
        st.stop()

    # 5. Standardize prediction columns before concatenation
    # Eukaryotes use 'Protein_ID' and 'Localizations', Prokaryotes use 'ACC' and 'Localization'
    df_euk_enz = df_euk_enz.rename(columns={'Protein_ID': 'Source ID (NCBI)', 'Localizations': 'Localization'})
    df_prok_enz = df_prok_enz.rename(columns={'ACC': 'Source ID (NCBI)'})
    
    df_euk_trans = df_euk_trans.rename(columns={'Protein_ID': 'Source ID (NCBI)', 'Localizations': 'Localization'})
    df_prok_trans = df_prok_trans.rename(columns={'ACC': 'Source ID (NCBI)'})

    df_all_enz_preds = pd.concat([df_euk_enz, df_prok_enz], ignore_index=True)
    df_all_trans_preds = pd.concat([df_euk_trans, df_prok_trans], ignore_index=True)

    # 6. Merge Main Metadata with Localization Predictions
    # Left merge to attach the 'Localization' column to our main biological data
    df_all_enzymes = pd.merge(df_enz_main, df_all_enz_preds[['Source ID (NCBI)', 'Localization']], on='Source ID (NCBI)', how='left')
    df_all_transporters = pd.merge(df_trans_main, df_all_trans_preds[['Source ID (NCBI)', 'Localization']], on='Source ID (NCBI)', how='left')
    
    # Rename 'Specie' to 'Strain' to guarantee the Matchmaker finds the correct column
    df_all_enzymes.rename(columns={'Specie': 'Strain'}, inplace=True)
    df_all_transporters.rename(columns={'Specie': 'Strain'}, inplace=True)

    return BiotransformationMatchmaker(df_feedipedia, df_all_enzymes, df_all_transporters)


def render_matchmaker_page():
    st.title("Matchmaker")
    st.markdown("""
    This page utilizes a knowledge-based recommendation system to cross-reference the 
    nutritional profile of inputs with the genotypic functional potential 
    of microalgae and cyanobacteria strains.
    """)
    st.divider()

    try:
        ai_engine = load_ai_engine()
    except Exception as e:
        st.error(f"Error loading the Machine Learning engine or datasets: {e}")
        return

    # Helper function to format and display the results cleanly
    def format_results(df, key_suffix):
        if "Affinity score percent" in df.columns:
            df["Metabolic affinity index"] = df["Affinity score percent"] / 10.0
            
        def categorize_match(score):
            if score >= 8.0: return "🟢 High affinity (complete pathway)"
            elif score >= 4.0: return "🟡 Partial affinity (incomplete pathway)"
            elif score > 0: return "🔴 Low affinity"
            else: return "⚫ No mapped enzymes or transporters"
            
        df["Match level"] = df["Metabolic affinity index"].apply(categorize_match)
        
        # Add a column for GCM search links
        def clean_for_search(name):
            clean_name = re.sub(r"[\[\]\(\)\']", "", str(name))
            parts = clean_name.split()
            if len(parts) >= 2:
                clean_name = f"{parts[0]} {parts[1]}"
            elif len(parts) == 1:
                clean_name = parts[0]
            return clean_name

        search_queries = df["Strain"].apply(clean_for_search)
        sp_query = search_queries.apply(lambda x: urllib.parse.quote(str(x)))
        df["GCM search"] = "https://gcm.wdcm.org/search?search=" + sp_query + "&list=strain"
        
        df_display = df[["Strain", "Match level", "Metabolic affinity index", "Mapped enzymes", "Mapped transporters", "GCM search"]]
        
        available_levels = df_display["Match level"].unique().tolist()
        
        available_levels.sort(reverse=True) 

        selected_levels = st.multiselect(
            "Filter results by Match level:",
            options=available_levels,
            default=available_levels, 
            key=f"filter_match_{key_suffix}"
        )

        if selected_levels:
            df_display = df_display[df_display["Match level"].isin(selected_levels)]
        else:
            df_display = df_display.iloc[0:0]

        st.dataframe(
            df_display,
            column_config={
                "Strain": st.column_config.TextColumn("Microalgae / cyanobacteria strain", width="medium"),
                "Match level": st.column_config.TextColumn("Match level", width="medium"),
                "Metabolic affinity index": st.column_config.NumberColumn(
                    "Affinity index",
                    help="Scale from 0 to 10.",
                    format="%.2f"
                ),
                "Mapped enzymes": st.column_config.TextColumn("Enzymes found", width="medium"),
                "Mapped transporters": st.column_config.TextColumn("Transporters found", width="medium"),
                "GCM search": st.column_config.LinkColumn(
                    "Availability", 
                    display_text="Search GCM",
                    width="small"
                )
            },
            use_container_width=True,
            hide_index=True
        )

    tab_waste, tab_sugar = st.tabs(["Agro-industrial byproducts", "Sugars"])

    with tab_waste:
        st.subheader("1. Agro-industrial byproduct selection")
        residue_options = ai_engine.df_feedipedia.index.tolist()
        
        selected_residue = st.selectbox(
            "Select an agro-industrial byproduct:",
            options=residue_options
        )

        if st.button("Run recommendation for byproduct", type="primary"):
            with st.spinner("Calculating targeted pathway affinities..."):
                st.session_state["results_waste"] = ai_engine.recommend_strains(selected_residue)
                st.session_state["last_residue"] = selected_residue
                
        if "results_waste" in st.session_state and st.session_state.get("last_residue") == selected_residue:
            results_df = st.session_state["results_waste"]
            
            if results_df.empty:
                st.warning(f"No mapped strains found for **{selected_residue}**.")
                st.dataframe(ai_engine.df_feedipedia.loc[selected_residue].replace(0.0, pd.NA).dropna())
            else:
                st.success(f"Analysis complete for **{selected_residue}**!")
                format_results(results_df, "waste")

    with tab_sugar:
        st.subheader("1. Sugar selection")
        sugar_options = ["fructose", "glucose", "galactose", "sucrose", "lactose", "maltose", "starch", "cellulose"]
        
        selected_sugar = st.selectbox(
            "Select a sugar:",
            options=[s.capitalize() for s in sugar_options]
        )

        if st.button("Run recommendation for sugar", type="primary"):
            with st.spinner("Calculating targeted pathway affinities..."):
                st.session_state["results_sugar"] = ai_engine.recommend_strains_for_pure_sugar(selected_sugar.lower())
                st.session_state["last_sugar"] = selected_sugar
                
        if "results_sugar" in st.session_state and st.session_state.get("last_sugar") == selected_sugar:
            results_sugar_df = st.session_state["results_sugar"]
            
            if results_sugar_df.empty:
                st.warning(f"No mapped strains found for {selected_sugar}.")
            else:
                st.success(f"Analysis complete for {selected_sugar}!")
                format_results(results_sugar_df, "sugar")