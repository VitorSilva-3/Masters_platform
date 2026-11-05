
import os
import time
import pandas as pd
import streamlit as st
from google import genai
from google.genai import types

class LLMService:
    def __init__(self, api_key: str):
        """Initializes the LLMService with the provided API key and sets up context."""

        self.client = genai.Client(api_key=api_key)
        
        self.static_csv_context = self._load_csv_previews()
        
        self.uploaded_feedipedia = self._get_or_upload_file()

    def _get_or_upload_file(self):
        """Checks if the Feedipedia file is already uploaded. If not, uploads it."""

        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        file_path = os.path.join(project_root, "data", "feedipedia_raw_data.json")
        display_name = "feedipedia_master_data"

        try:
            for f in self.client.files.list():
                if f.display_name == display_name:
                    return f
        except Exception as e:
            print(f"[-] Error checking existing files: {e}")

        if os.path.exists(file_path):
            try:
                print("Uploading Feedipedia to Gemini API...")
                uploaded_file = self.client.files.upload(
                    file=file_path,
                    config={'display_name': display_name}
                )
                return uploaded_file
            except Exception as e:
                print(f"[-] Error uploading file: {e}")
                return None
        return None

    def _load_csv_previews(self) -> str:
        """Loads a subset of CSVs to provide structural context to the LLM."""

        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        data_path = os.path.join(project_root, "data")
        
        context_parts = []
        for filename in ["enzymes_data.csv", "transporters_data.csv"]:
            file_path = os.path.join(data_path, filename)
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path)
                    csv_str = df.head(20).to_csv(index=False)
                    context_parts.append(f"--- FILE: {filename} (Preview) ---\n{csv_str}\n... [TRUNCATED]\n")
                except Exception:
                    pass
        return "\n".join(context_parts)

    def get_chat_response_stream(self, user_prompt: str, chat_history: list):
        """Sends the prompt to Gemini and yields the response in chunks.
           Includes exponential backoff for API overload and smart file attachment."""
        
        system_instruction = (
            "You are the core intelligent assistant of platform "
            "dedicated to agro-industrial waste valorization.\n"
            "Your main mission is:\n"
            "1. Help users navigate the platform and understand its purpose.\n"
            "2. Explain scientific concepts related to biotechnology and biological engineering.\n"
            "3. Act as an expert in agro-industrial residues.\n\n"
            "You have direct access to the Feedipedia JSON file attached to this prompt. "
            "Always consult it carefully to extract exact nutritional data when requested.\n"
            "Always reply in English.\n\n"
            f"--- PLATFORM CSV PREVIEWS ---\n{self.static_csv_context}\n"
            f"---------------------------------------"
        )

        formatted_contents = []
        is_first_message = True
        
        for msg in chat_history:
            role = "user" if msg["role"] == "user" else "model"
            parts = [types.Part.from_text(text=msg["content"])]
            
            if is_first_message and role == "user" and self.uploaded_feedipedia:
                parts.insert(0, types.Part.from_uri(
                    file_uri=self.uploaded_feedipedia.uri,
                    mime_type=self.uploaded_feedipedia.mime_type
                ))
                is_first_message = False
                
            formatted_contents.append(types.Content(role=role, parts=parts))
        
        current_parts = []
        
        if is_first_message and self.uploaded_feedipedia:
            current_parts.append(types.Part.from_uri(
                file_uri=self.uploaded_feedipedia.uri,
                mime_type=self.uploaded_feedipedia.mime_type
            ))
            
        current_parts.append(types.Part.from_text(text=user_prompt))
        formatted_contents.append(types.Content(role="user", parts=current_parts))

        model_name = st.secrets.get("GEMINI_MODEL", "gemini-1.5-flash")
        max_tentativas = 3
        
        for tentativa in range(max_tentativas):
            try:
                response_stream = self.client.models.generate_content_stream(
                    model=model_name,
                    contents=formatted_contents,
                    config=types.GenerateContentConfig(
                        system_instruction=system_instruction,
                        temperature=0.3,
                    )
                )
                
                for chunk in response_stream:
                    yield chunk.text
                return 
                
            except Exception as e:
                erro_str = str(e)
                if ("503" in erro_str or "429" in erro_str or "UNAVAILABLE" in erro_str.upper()):
                    if tentativa < max_tentativas - 1:
                        tempo_espera = 2 ** tentativa  
                        time.sleep(tempo_espera)
                        continue  
                
                yield f"\n\n Could not connect to the AI server. Reason: {erro_str}"
                return