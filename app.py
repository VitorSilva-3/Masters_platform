
import streamlit as st
import time
import streamlit_authenticator as stauth
import bcrypt
from streamlit_lottie import st_lottie
from utils import load_lottiefile, configure_page
from services.mongo_service import MongoService

configure_page("Home")

@st.cache_resource
def get_mongo_service():
    try:
        return MongoService(uri=st.secrets["MONGO_STRING"])
    except Exception as e:
        st.error(f"Error connecting to the database: {e}")
        return None

mongo = get_mongo_service()
if mongo is None:
    st.stop()

credentials = mongo.get_all_users()

authenticator = stauth.Authenticate(
    credentials,
    st.secrets["cookie"]["name"],
    st.secrets["cookie"]["key"],
    st.secrets["cookie"]["expiry_days"]
)

if not st.session_state.get("authentication_status"):
    tab1, tab2 = st.tabs(["Login", "Register new profile"])
    
    with tab1:
        authenticator.login(location='main', fields={'Form name': 'Login', 'Username': 'Email', 'Password': 'Password', 'Login': 'Login'})
        
    with tab2:
        st.markdown("Register new profile")
        with st.form("register_custom_form", clear_on_submit=True):
            new_name = st.text_input("Name").strip()
            new_email = st.text_input("Email").strip()
            new_password = st.text_input("Password", type="password")
            
            submit_register = st.form_submit_button("Register", type="primary")

            if submit_register:
                if new_name and new_email and new_password:
                    salt = bcrypt.gensalt()
                    hashed_password = bcrypt.hashpw(new_password.encode('utf-8'), salt).decode('utf-8')
                    
                    try:
                        mongo.add_user(
                            email=new_email,
                            name=new_name,
                            hashed_password=hashed_password
                        )
                        st.success("Profile registered successfully! You can now log in via the Login tab.")
                    except Exception as e:
                        st.error(f"An error occurred: {e}")
                else:
                    st.warning("Please fill in all fields to register.")

if st.session_state["authentication_status"] is False:
    st.error('Email or password is incorrect')

elif st.session_state["authentication_status"] is None:
    st.warning('Please enter your email and password')

elif st.session_state["authentication_status"]:
    current_email = st.session_state["username"]
    current_role = credentials['usernames'].get(current_email, {}).get('role', 'user')

    with st.sidebar:
        st.markdown(f"**Logged in as:** {st.session_state['name']}")
        st.caption(f"Role: {current_role.capitalize()}")
        authenticator.logout('Logout', 'sidebar')
        
        if current_role == "admin":
            st.divider()
            st.markdown("### Admin panel")
            with st.expander("Manage users"):
                users_list = mongo.get_users_list()
                if len(users_list) == 1:
                    st.write("No other users registered yet.")
                else:
                    for u in users_list:
                        if u['email'] != current_email: 
                            col1, col2 = st.columns([3, 1])
                            
                            col1.markdown(f"**{u['name']}**")
                            
                            if u.get('role') == 'admin':
                                col2.caption("Admin")
                            else:
                                if col2.button("Delete", type="primary", key=f"del_{u['email']}", help="Delete profile and chats"):
                                    mongo.delete_user(u['email'])
                                    st.rerun()

    LOTTIE_DNA = "images/DNA strand.json" 
    LOTTIE_FACTORY = "images/Factory pollution city air and water.json"

    if 'app_loaded' not in st.session_state:
        st.session_state.app_loaded = False

    if not st.session_state.app_loaded:
        hide_sidebar_css = """
            <style>
                [data-testid="stSidebar"] {display: none !important;}
                [data-testid="collapsedControl"] {display: none !important;}
            </style>
        """
        st.markdown(hide_sidebar_css, unsafe_allow_html=True)
        
        st.markdown("<br><br>", unsafe_allow_html=True)
        
        anim_col1, anim_col2 = st.columns(2)
        
        with anim_col1:
            anim1 = load_lottiefile(LOTTIE_DNA)
            if anim1:
                st_lottie(anim1, height=350, key="loading_dna")
                
        with anim_col2:
            anim2 = load_lottiefile(LOTTIE_FACTORY)
            if anim2:
                st_lottie(anim2, height=350, key="loading_factory")
                
        st.markdown("<h2 style='text-align: center; color: #4CAF50;'>Loading MicroValue...</h2>", unsafe_allow_html=True)
        
        bar_col1, bar_col2, bar_col3 = st.columns([1, 2, 1])
        with bar_col2:
            progress_bar = st.progress(0)
            status_text = st.empty()

            loading_steps = [
                "Initializing platform environment...",
                "Loading central databases...",
                "Configuring modules...",
                "Preparing user interface...",
                "Ready."
            ]

            for i, step in enumerate(loading_steps):
                status_text.markdown(f"<p style='text-align: center; color: #888; font-size: 1.1em;'><i>{step}</i></p>", unsafe_allow_html=True)
                progress_bar.progress((i + 1) / len(loading_steps))
                time.sleep(1)

        time.sleep(0.5)
        st.session_state.app_loaded = True 
        st.rerun()

    else:
        st.title(f"Welcome to MicroValue, {st.session_state['name']}!")
        st.markdown("#### *Biotechnological valorization of agro-industrial byproducts through microalgae and cyanobacteria*")
        st.divider()

        st.markdown(
            """
            This platform integrates data on agro-industrial byproducts as well as the enzymes and transporters produced by various microalgae 
            and cyanobacteria organised by species. It identifies which byproducts can serve as substrates for each microalga/cyanobacterium,
            indicating those with the greatest growth potential.
            """
        )

        st.info("👈 **Use the left menu to navigate through the platform's specialized modules.**")

        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### Platform modules")

        with st.container(border=True):
            st.page_link("pages/1_Data.py", label = ":blue[**Data**]", use_container_width = True)
            st.write("Explore the comprehensive datasets cross-referencing species, target sugars, degrading enzymes, and transport mechanisms.")

        with st.container(border=True):
            st.page_link("pages/2_Scientific_literature.py", label = ":blue[**Scientific literature**]", use_container_width = True)
            st.write("Access an automated curation of PubMed publications relevant to the catalogued species and their biochemical repertoire.")
        
        with st.container(border=True):
            st.page_link("pages/3_Taxonomy.py", label = ":blue[**Taxonomy**]", use_container_width = True)
            st.write("Navigate the detailed evolutionary lineage and taxonomic classification of the catalogued microalgae and cyanobacteria.")
        
        with st.container(border=True):
            st.page_link("pages/4_KEGG.py", label = ":blue[**KEGG**]", use_container_width = True)
            st.write("Map metabolic networks and analyze the specific catalytic reactions.")

        with st.container(border=True):
            st.page_link("pages/5_UniProt.py", label = ":blue[**UniProt**]", use_container_width = True)
            st.write("Investigate the general biochemical properties, functional annotations, and structural features of the targeted proteins.")

        with st.container(border=True):
            st.page_link("pages/6_Subcellular_localization.py", label = ":blue[**Subcellular localization**]", use_container_width = True)
            st.write("Examine AI-driven predictions of subcellular compartments utilizing DeepLoc 2.0 and DeepLocPro models.")

        with st.container(border = True):
            st.page_link("pages/7_Agro-industrial_waste_composition.py", label = ":blue[**Agro-industrial waste composition**]", use_container_width = True)
            st.write("Analyze the compositional data of agro-industrial byproducts according to the Feedipedia database.")

        with st.container(border = True):
            st.page_link("pages/8_Matchmaker.py", label = ":blue[**Matchmaker**]", use_container_width = True)
            st.write("Identify the most suitable microalgae and cyanobacteria to be used for the valorization of some agro-industrial byproducts and sugars.")

        with st.container(border = True):
            st.page_link("pages/9_Chatbot.py", label = ":blue[**Chatbot**]", use_container_width = True)
            st.write("Interact with an AI-powered chatbot to ask questions and get insights about the platform's data and features.")