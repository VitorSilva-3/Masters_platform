
import os
import pandas as pd
import streamlit as st
from google import genai
from google.genai import types

class LLMService:
    def __init__(self, api_key: str):
        self.client = genai.Client(api_key=api_key)
        
        # 1. Carrega os CSVs pequenos como texto estático
        self.static_csv_context = self._load_csv_previews()
        
        # 2. Gere o Upload do ficheiro pesado (Feedipedia) para a Google
        self.uploaded_feedipedia = self._get_or_upload_file()

    def _get_or_upload_file(self):
        """
        Verifica se o ficheiro Feedipedia já está na nuvem da Google.
        Se não estiver, faz o upload automaticamente.
        """
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        file_path = os.path.join(project_root, "data", "feedipedia_raw_data.json")
        display_name = "feedipedia_master_data"

        # Tenta encontrar o ficheiro já carregado para não duplicar uploads
        try:
            for f in self.client.files.list():
                if f.display_name == display_name:
                    return f
        except Exception as e:
            print(f"[-] Error checking existing files: {e}")

        # Se não encontrou, faz o upload
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
        """Loads previews of the main CSV files."""
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
        """Sends the prompt to Gemini and YIELDS the response in chunks for real-time streaming."""
        
        system_instruction = (
            "You are the core intelligent assistant of platform "
            "dedicated to agro-industrial waste valorization.\n"
            "Your main mission is:\n"
            "1. Help users navigate the platform and understand its purpose.\n"
            "2. Explain scientific concepts related to biotechnology.\n"
            "3. Act as an expert in agro-industrial residues.\n\n"
            "You have direct access to the Feedipedia JSON file attached to this prompt. "
            "Always consult it carefully to extract exact nutritional data when requested.\n"
            "Always reply in English.\n\n"
            f"--- PLATFORM CSV PREVIEWS ---\n{self.static_csv_context}\n"
            f"---------------------------------------"
        )

        formatted_contents = []
        
        # 1. Reconstrói o histórico do chat
        for msg in chat_history:
            role = "user" if msg["role"] == "user" else "model"
            formatted_contents.append(
                types.Content(role=role, parts=[types.Part.from_text(text=msg["content"])])
            )
        
        # 2. Constrói o novo prompt do utilizador
        current_parts = []
        
        # Anexa o ficheiro diretamente ao prompt atual do utilizador (se o upload funcionou)
        # Anexa o ficheiro usando a referência URI correta da nuvem da Google
        if self.uploaded_feedipedia:
            file_part = types.Part.from_uri(
                file_uri=self.uploaded_feedipedia.uri,
                mime_type=self.uploaded_feedipedia.mime_type
            )
            current_parts.append(file_part)
            
        current_parts.append(types.Part.from_text(text=user_prompt))
        
        formatted_contents.append(
            types.Content(role="user", parts=current_parts)
        )

        try:
            # Usa o modelo mais leve e rápido (fallback para 1.5-flash)
            model_name = st.secrets.get("GEMINI_MODEL", "gemini-1.5-flash")
            
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
                
        except Exception as e:
            yield f"Sorry, an error occurred while connecting to the AI server: {str(e)}"