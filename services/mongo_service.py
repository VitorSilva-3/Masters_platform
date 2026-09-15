
from pymongo import MongoClient

class MongoService:
    def __init__(self, uri: str, db_name="Masters_platform", collection_name="chat_histories"):
        """Initializes the MongoDB connection."""

        try:
            self.client = MongoClient(uri)
            self.db = self.client[db_name]
            self.collection = self.db[collection_name]
            # Testa a ligação (opcional, mas recomendado)
            self.client.admin.command('ping')
        except Exception as e:
            print(f"Error initializing MongoDB connection: {e}")
            raise e

    def save_chat(self, user: str, chat_id: str, messages: list):
        """Saves the chat messages for a specific user and chat ID."""

        self.collection.update_one(
            {"chat_id": chat_id, "user": user},
            {"$set": {"messages": messages}},
            upsert=True
        )

    def load_chat(self, user: str, chat_id: str) -> list:
        """Loads the messages for a specific user and chat ID."""

        document = self.collection.find_one({"chat_id": chat_id, "user": user})
        if document:
            return document.get("messages", [])
        return []

    def get_user_chats(self, user: str) -> list:
        """Returns a list of chat IDs and their corresponding titles for a specific user."""

        chats = self.collection.find({"user": user}, {"chat_id": 1, "messages": 1, "_id": 0})
        chat_list = []
        
        for chat in chats:
            chat_id = chat.get("chat_id")
            messages = chat.get("messages", [])
            
            # Procurar a primeira mensagem do utilizador para usar como título
            title = f"Chat {chat_id}" # Título de segurança
            for msg in messages:
                if msg["role"] == "user":
                    content = msg["content"]
                    # Trunca o texto aos primeiros 35 caracteres e adiciona "..."
                    title = content[:35] + "..." if len(content) > 35 else content
                    break
            
            chat_list.append({"id": chat_id, "title": title})
            
        return chat_list[::-1]

    def delete_chat(self, user: str, chat_id: str):
        """Deletes a specific chat from the database"""

        self.collection.delete_one({"chat_id": chat_id, "user": user})