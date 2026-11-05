
from pymongo import MongoClient

class MongoService:
    def __init__(self, uri: str, db_name="Masters_platform", collection_name="chat_histories"):
        """Initializes the MongoDB connection."""

        try:
            self.client = MongoClient(uri)
            self.db = self.client[db_name]
            self.collection = self.db[collection_name] 
            self.users_collection = self.db["users"]  
            self.client.admin.command('ping')
        except Exception as e:
            print(f"Error initializing MongoDB connection: {e}")
            raise e

    def save_chat(self, user: str, chat_id: str, messages: list):
        """Saves or updates a chat history for a specific user."""

        self.collection.update_one(
            {"chat_id": chat_id, "user": user},
            {"$set": {"messages": messages}},
            upsert=True
        )

    def load_chat(self, user: str, chat_id: str) -> list:
        """Loads a chat history for a specific user."""

        document = self.collection.find_one({"chat_id": chat_id, "user": user})
        if document:
            return document.get("messages", [])
        return []

    def get_user_chats(self, user: str) -> list:
        """Returns a list of all chat histories for a specific user."""

        chats = self.collection.find({"user": user}, {"chat_id": 1, "messages": 1, "_id": 0})
        chat_list = []
        for chat in chats:
            chat_id = chat.get("chat_id")
            messages = chat.get("messages", [])
            title = f"Chat {chat_id}" 
            for msg in messages:
                if msg["role"] == "user":
                    content = msg["content"]
                    title = content[:35] + "..." if len(content) > 35 else content
                    break
            chat_list.append({"id": chat_id, "title": title})
        return chat_list[::-1]

    def delete_chat(self, user: str, chat_id: str):
        """Deletes a specific chat history for a specific user."""

        self.collection.delete_one({"chat_id": chat_id, "user": user})

    def get_all_users(self) -> dict:
        """Returns a dictionary of all users with their credentials for authentication purposes."""

        users = self.users_collection.find({})
        credentials = {"usernames": {}}
        for u in users:
            email = u.get("email")
            credentials["usernames"][email] = {
                "email": email,
                "name": u.get("name"),
                "password": u.get("password"),
                "role": u.get("role", "user")
            }
        return credentials

    def add_user(self, email: str, name: str, hashed_password: str):
        """Adds a new user to the database."""

        is_first_user = self.users_collection.count_documents({}) == 0
        role = "admin" if is_first_user else "user"
        
        self.users_collection.update_one(
            {"email": email},
            {"$set": {
                "email": email,
                "name": name,
                "password": hashed_password,
                "role": role
            }},
            upsert=True
        )

    def delete_user(self, email: str):
        """Deletes a user and all their associated chat histories from the database."""

        self.users_collection.delete_one({"email": email})
        self.collection.delete_many({"user": email})

    def get_users_list(self) -> list:
        """Returns a list of all users with their email, name and role."""
        
        return list(self.users_collection.find({}, {"_id": 0, "email": 1, "name": 1, "role": 1}))