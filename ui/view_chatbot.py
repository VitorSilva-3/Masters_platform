
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
        st.header("Chat history")
        
        user_chats_data = mongo_service.get_user_chats(active_user)
        
        # Botão: Novo Chat
        if st.button("New chat", use_container_width=True):
            st.session_state.current_chat_id = str(uuid.uuid4())[:8]
            st.session_state.messages = [{"role": "assistant", "content": f"Hello, {active_user}! How can I help you today?"}]
            # A linha de save_chat foi removida daqui! 
            # Só guardamos quando o utilizador escrever algo na caixa de texto.
            st.rerun()

        # Dropdown: Carregar chat anterior (com títulos legíveis)
        if user_chats_data:
            chat_options = {"Current session": "Current session"}
            for chat in user_chats_data:
                label = f"{chat['title']} ({chat['id']})"
                chat_options[label] = chat['id']
                
            selected_label = st.selectbox("Chat history:", list(chat_options.keys()))
            selected_chat_id = chat_options[selected_label]
            
            if selected_chat_id != "Current session":
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("Load", use_container_width=True):
                        st.session_state.current_chat_id = selected_chat_id
                        st.session_state.messages = mongo_service.load_chat(active_user, selected_chat_id)
                        st.rerun()
                with col2:
                    if st.button("Delete", type="primary", use_container_width=True):
                        mongo_service.delete_chat(active_user, selected_chat_id)
                        
                        # Se eliminarmos o chat aberto, o ecrã reinicia para um Novo Chat vazio na memória (não na DB)
                        if st.session_state.current_chat_id == selected_chat_id:
                            st.session_state.current_chat_id = str(uuid.uuid4())[:8]
                            st.session_state.messages = [{"role": "assistant", "content": f"Hello, {active_user}! How can I help you today?"}]
                            # A linha de save_chat também foi removida daqui!
                            
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
        
        st.rerun()