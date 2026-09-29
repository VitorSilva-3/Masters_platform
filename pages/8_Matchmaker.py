
import streamlit as st
from ui.view_matchmaker import render_matchmaker_page
from utils import configure_page

configure_page("Matchmaker")

render_matchmaker_page()