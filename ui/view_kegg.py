
import streamlit as st
import pandas as pd
from services.kegg_service import KeggService

@st.cache_data(show_spinner=False)
def get_kegg_info(enzyme_name: str, ec_number: str):
    service = KeggService()
    return service.fetch_enzyme_details(enzyme_name, ec_number)

def render_kegg_view(df: pd.DataFrame):
    """
    Renders the KEGG data tab with an analytical dashboard layout.
    """

    st.subheader("KEGG database")
    st.markdown("Explore fundamental biochemical properties and metabolic pathways from the KEGG database.")
    
    if df.empty:
        st.warning("No data available. Please build the main dataset first.")
        return

    unique_enzymes = df[['Enzyme', 'EC number']].dropna().drop_duplicates()
    options = {f"{row['Enzyme']} (EC {row['EC number']})": row for _, row in unique_enzymes.iterrows()}
    
    selected_option = st.selectbox(
        "Select an enzyme to analyze:", 
        options=list(options.keys()),
        index=None,
        placeholder="Search for an enzyme...",
        key="kegg_enzyme_selector"
    )
    
    if selected_option:
        selected_enzyme = options[selected_option]['Enzyme']
        selected_ec = options[selected_option]['EC number']
        
        st.divider()
        
        with st.spinner(f"Fetching KEGG data for EC {selected_ec}..."):
            info = get_kegg_info(selected_enzyme, selected_ec)
            
        if not info or not info.get('name'):
            st.error(f"Could not retrieve KEGG data for EC {selected_ec}. It may not be mapped in the KEGG database yet.")
            return

        ec_parts = selected_ec.split(".")
        class_num = ec_parts[0] if len(ec_parts) > 0 else "-"
        subclass_num = ec_parts[1] if len(ec_parts) > 1 else "-"
        sub_subclass_num = ec_parts[2] if len(ec_parts) > 2 else "-"
        serial_num = ec_parts[3] if len(ec_parts) > 3 else "-"

        official_name = info['name'].split(';')[0] if info['name'] else "Unknown"
        reaction_text = info['reaction'] if info['reaction'] else "Reaction not available"
        kegg_url = f"https://www.genome.jp/entry/{selected_ec}"

        col_title, col_btn = st.columns([4, 1])
        with col_title:
            st.markdown(f"### {official_name}")
        with col_btn:
            st.link_button("View on KEGG", kegg_url, use_container_width=True)

        st.markdown("<br>", unsafe_allow_html=True)

        st.markdown("#### Catalytic reaction")
        st.info(reaction_text)
        
        st.markdown("<br>", unsafe_allow_html=True)

        st.markdown("#### EC number breakdown")
        with st.container(border=True):
            cc1, cc2, cc3, cc4 = st.columns(4)
            cc1.metric("Class", class_num, help="Main enzyme class")
            cc2.metric("Subclass", subclass_num)
            cc3.metric("Sub-subclass", sub_subclass_num)
            cc4.metric("Serial number", serial_num)
            
        st.markdown("<br>", unsafe_allow_html=True)
        
        st.markdown("#### Associated metabolic pathways")
        
        if info['pathways']:
            for pw in info['pathways']:
                pw_code = pw['code']
                pw_desc = pw['description']
                st.markdown(f"- **{pw_code}**: {pw_desc}")
        else:
            st.write("No pathways mapped for this enzyme.")