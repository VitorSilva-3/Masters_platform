
import streamlit as st
from ui.view_feedipedia import render_feedipedia_view
from utils import configure_page

configure_page("Agro-industrial waste composition")

render_feedipedia_view()