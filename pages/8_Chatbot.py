
import streamlit as st
from ui.view_chatbot import render_chat_view
from utils import configure_page

configure_page("Chatbot assistant")

render_chat_view()