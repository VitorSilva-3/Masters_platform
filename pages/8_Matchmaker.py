
import streamlit as st
from ui.view_ml import render_ml_page
from utils import configure_page

configure_page("Matchmaker")

render_ml_page()