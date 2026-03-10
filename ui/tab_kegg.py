
import streamlit as st
import pandas as pd
from services.kegg_service import KeggService

@st.cache_data(show_spinner=False)
def get_kegg_info(enzyme_name: str, ec_number: str):
    service = KeggService()
    return service.fetch_enzyme_details(enzyme_name, ec_number)

def render_kegg_tab(df: pd.DataFrame):
    """
    Renders the KEGG data tab.
    """

    st.subheader("KEGG Enzyme integration")
    
    if df.empty:
        st.warning("No data available. Please build the main dataset first.")
        return

    enzyme_data = df[['Enzyme', 'EC number']].drop_duplicates().to_dict('records')
    options = {f"{row['Enzyme']} (EC: {row['EC number']})": row for row in enzyme_data}
    
    st.markdown("Select an enzyme from your dataset to view its fundamental biochemical properties:")
    selected_option = st.selectbox(
        "Choose Enzyme:", 
        options=list(options.keys()),
        label_visibility="collapsed",
        key = "kegg_enzyme_selector"
    )
    
    if selected_option:
        selected_enzyme = options[selected_option]['Enzyme']
        selected_ec = options[selected_option]['EC number']
        
        with st.spinner(f"Fetching KEGG data for {selected_ec}..."):
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
        reaction_text = info['reaction'] if info['reaction'] else "Not available"
        kegg_url = f"https://www.genome.jp/entry/{selected_ec}"

        st.markdown("<br>", unsafe_allow_html=True)
        
        with st.expander(f"{selected_enzyme.capitalize()} (EC {selected_ec})", expanded=True):
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.write(f"**Official name:** {official_name}")
                st.markdown("<br>", unsafe_allow_html=True) 
                st.write(f"**Reaction:** {reaction_text}")
                
            with col2:
                st.write(f"**EC Number:** {selected_ec}")
                st.markdown("<br>", unsafe_allow_html=True) 
                st.write(f"**Class:** {class_num} (Main enzyme class)")
                st.markdown("<br>", unsafe_allow_html=True)
                st.write(f"**Subclass:** {subclass_num}")
                st.markdown("<br>", unsafe_allow_html=True)
                st.write(f"**Sub-subclass:** {sub_subclass_num}")
                st.markdown("<br>", unsafe_allow_html=True)
                st.write(f"**Serial number:** {serial_num}")
                
            st.divider()
            
            st.write("**Associated metabolic pathways:**")
            st.write("") 
            
            if info['pathways']:
                for pw in info['pathways']:
                    st.write(f"• {pw['code']}: {pw['description']}")
            else:
                st.write("• No pathways mapped.")
                
            st.markdown("<br>", unsafe_allow_html=True)
            
            st.markdown(f"[View detailed information on KEGG database]({kegg_url})")