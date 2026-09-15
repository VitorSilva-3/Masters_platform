
import uuid
import streamlit as st
from services.llm_service import LLMService
from services.mongo_service import MongoService

@st.cache_resource
def get_mongo_service():
    try:
        uri = st.secrets["MONGO_STRING"]
        return MongoService(uri=uri)
    except Exception as e:
        st.error(f"Error initializing database: {e}")
        return None

def render_chat_view():
    st.title("Chatbot assistant")

    try:
        api_key = st.secrets["GEMINI_API_KEY"]
    except KeyError:
        st.error("API key for Gemini not found. Please set it in Streamlit secrets.")
        return

    # Inicialização dos Serviços
    if "llm_service" not in st.session_state:
        with st.spinner("Loading local databases and models..."):
            st.session_state.llm_service = LLMService(api_key=api_key)
            
    mongo_service = get_mongo_service()
    if mongo_service is None:
        st.stop() # Para a execução se a BD falhar

    with st.sidebar:
        st.header("User profile")
        
        profiles = ["Filipe", "Pedro", "Vítor"]
        active_user = st.selectbox("Select profile:", profiles)
        
        if "active_user" not in st.session_state or st.session_state.active_user != active_user:
            st.session_state.active_user = active_user
            st.session_state.current_chat_id = str(uuid.uuid4())[:8]
            st.session_state.messages = [{"role": "assistant", "content": f"Hello, {active_user}! How can I help you today?"}]
            st.rerun()

        st.divider()
        st.header("Management of chats")
        
        user_chats = mongo_service.get_user_chats(active_user)
        
        if st.button("New chat", use_container_width=True):
            st.session_state.current_chat_id = str(uuid.uuid4())[:8]
            st.session_state.messages = [{"role": "assistant", "content": f"Hello, {active_user}! How can I help you today?"}]
            mongo_service.save_chat(active_user, st.session_state.current_chat_id, st.session_state.messages)
            st.rerun()

        if user_chats:
            selected_chat = st.selectbox("Load previous chat:", ["Current Session"] + user_chats)
            if selected_chat != "Current Session":
                if st.button("Load selected", use_container_width=True):
                    st.session_state.current_chat_id = selected_chat
                    st.session_state.messages = mongo_service.load_chat(active_user, selected_chat)
                    st.rerun()

        st.divider()
        if st.button("Clear current chat", type="primary", use_container_width=True):
            st.session_state.messages = [{"role": "assistant", "content": f"Hello, {active_user}! How can I help you today?"}]
            mongo_service.save_chat(active_user, st.session_state.current_chat_id, st.session_state.messages)
            st.rerun()

    # --- ÁREA PRINCIPAL DO CHAT ---
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("Ask the assistant..."):
        with st.chat_message("user"):
            st.markdown(prompt)
            
        with st.chat_message("assistant"):
            chat_history = st.session_state.messages
            response_stream = st.session_state.llm_service.get_chat_response_stream(
                user_prompt=prompt, 
                chat_history=chat_history
            )
            full_response = st.write_stream(response_stream)
                
        st.session_state.messages.append({"role": "user", "content": prompt})
        st.session_state.messages.append({"role": "assistant", "content": full_response})
        mongo_service.save_chat(active_user, st.session_state.current_chat_id, st.session_state.messages)