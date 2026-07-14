
import streamlit as st
from services.llm_service import LLMService

def render_chat_view():
    st.title("Chatbot assistant")

    # --- Sidebar: API Configuration ---
    with st.sidebar:
        st.header("AI Engine")
        st.markdown("To activate the assistant, enter your free Google Gemini API key.")
        api_key = st.text_input("Gemini API Key", type="password", help="The key is not saved, it is only used for this active session.")
        st.markdown("[Get Free API Key (Google AI Studio)](https://aistudio.google.com/)")
        
        st.divider()
        if st.button("Clear chat history"):
            st.session_state.messages = []
            st.rerun()

    # Stop rendering if no API key is provided
    if not api_key:
        st.info("Please enter your Gemini API Key in the sidebar to start chatting with the assistant.")
        return

    # --- LLM Service Initialization ---
    if "llm_service" not in st.session_state or st.session_state.get("current_api_key") != api_key:
        try:
            with st.spinner("Loading knowledge base and platform models..."):
                st.session_state.llm_service = LLMService(api_key=api_key)
                st.session_state.current_api_key = api_key
        except Exception as e:
            st.error(f"Error initializing AI service: {e}")
            return

    # --- Chat History Management ---
    if "messages" not in st.session_state or not st.session_state.messages:
        st.session_state.messages = [
            {"role": "assistant", "content": "Hello! Welcome to the Agro-industrial Waste Valorization Platform. You can ask me how the integrated Machine Learning models work, clarify scientific doubts about biotechnology, or request the chemical composition of specific residues in our database (e.g., *'What is the amount of lactose in whey?'*). How can I help you today?"}
        ]

    # Display all past messages
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    # --- User Input Box ---
    if prompt := st.chat_input("Ex: Explain bioinformatics to me or tell me which residue has the highest lignin content."):
        
        # 1. Display user prompt
        with st.chat_message("user"):
            st.markdown(prompt)
            
        # 2. Display thinking indicator
        with st.chat_message("assistant"):
            with st.spinner("Analyzing query..."):
                chat_history = st.session_state.messages
                
                # Call the LLM service
                response = st.session_state.llm_service.get_chat_response(
                    user_prompt=prompt, 
                    chat_history=chat_history
                )
                
                # Display the final response
                st.markdown(response)
                
        # 3. Save messages to session state
        st.session_state.messages.append({"role": "user", "content": prompt})
        st.session_state.messages.append({"role": "assistant", "content": response})