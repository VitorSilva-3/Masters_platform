
import os
import json
from google import genai
from google.genai import types

class LLMService:
    def __init__(self, api_key: str):
        """
        Initializes the LLM Service with the Gemini API key and loads the local database context.
        """
        self.client = genai.Client(api_key=api_key)
        self.data_context = self._load_feedipedia_context()

    def _load_feedipedia_context(self) -> str:
        """
        Reads the local JSON file and compresses it into a string to be injected into the LLM prompt.
        """
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        data_path = os.path.join(project_root, "data", "feedipedia_raw_data.json")
        
        if os.path.exists(data_path):
            try:
                with open(data_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    return json.dumps(data, ensure_ascii=False)
            except Exception as e:
                print(f"[-] Error loading Feedipedia context: {e}")
                
        return "No local database context available."

    def get_chat_response(self, user_prompt: str, chat_history: list) -> str:
        """
        Sends the user's prompt, the conversation history, and the system instructions to the Gemini model.
        """
        # 1. Prepare the Hybrid System Instruction (Platform Guide + Database Expert) in English
        system_instruction = (
            "You are the core intelligent assistant of an advanced bioinformatics and machine learning platform "
            "dedicated to agro-industrial waste valorization.\n"
            "Your main mission is:\n"
            "1. Help users navigate the platform and understand its purpose (predicting and optimizing the use of "
            "residues using artificial intelligence).\n"
            "2. Explain scientific concepts related to biotechnology, bioinformatics, microalgae metabolism, "
            "and predictive models in a clear and accessible way.\n"
            "3. Act as an expert in agro-industrial residues.\n\n"
            "In addition to your general scientific knowledge, you have direct access to the platform's internal database "
            "(extracted from Feedipedia). When the user asks about exact values, metrics, or chemical composition "
            "(e.g., dry matter, crude protein, lactose, lignin) of specific residues, you MUST consult "
            "the provided data below and base your answer on those exact numbers.\n"
            "Always reply in English.\n\n"
            f"--- PLATFORM INTERNAL DATABASE (FEEDIPEDIA) ---\n{self.data_context}\n---------------------------------------"
        )

        # 2. Format the Streamlit chat history into the Gemini API format
        formatted_contents = []
        for msg in chat_history:
            role = "user" if msg["role"] == "user" else "model"
            formatted_contents.append(
                types.Content(role=role, parts=[types.Part.from_text(text=msg["content"])])
            )
        
        # 3. Add the current user prompt to the very end of the contents list
        formatted_contents.append(
            types.Content(role="user", parts=[types.Part.from_text(text=user_prompt)])
        )

        # 4. Make the API Call
        try:
            response = self.client.models.generate_content(
                model='gemini-1.5-flash',
                contents=formatted_contents,
                config=types.GenerateContentConfig(
                    system_instruction=system_instruction,
                    temperature=0.3, # Slightly higher for fluid scientific explanations
                )
            )
            return response.text
            
        except Exception as e:
            return f"Sorry, an error occurred while connecting to the AI server: {str(e)}"