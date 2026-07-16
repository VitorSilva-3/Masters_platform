
import streamlit as st
from services.llm_service import LLMService

def render_chat_view():
    st.title("Chatbot assistant")

    try:
        api_key = st.secrets["GEMINI_API_KEY"]
    except KeyError:
        st.error("API key for Gemini not found. Please set it in Streamlit secrets.")
        return

    with st.sidebar:
        st.header("Chatbot assistant")
        if st.button("Clear chat history"):
            st.session_state.messages = []
            st.rerun()

    if "llm_service" not in st.session_state:
        try:
            with st.spinner("Loading local databases and models... This may take a few seconds."):
                st.session_state.llm_service = LLMService(api_key=api_key)
        except Exception as e:
            st.error(f"Error initializing the AI service: {e}")
            return

    if "messages" not in st.session_state or not st.session_state.messages:
        st.session_state.messages = [
            {"role": "assistant", "content": "Hello! Welcome to MicroValue. How can I help you today?"}
        ]

    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("Ask chatbot assistant..."):
        
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