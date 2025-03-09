import streamlit as st
from streamlit_option_menu import option_menu

# --- Display Sidebar and Pages ---
def display_sidebar_and_pages():
    st.sidebar.title("Welcome!")
    st.sidebar.markdown("#### Select a page to navigate:")

    with st.sidebar:
        page = option_menu(
            menu_title=None,
            options=["Home Page", "Run Tasks"],
            icons=["house", "play-circle"],
            menu_icon="cast",
            default_index=0,
            orientation="vertical",
            styles={
                "nav-link": {
                    "font-size": "16px",
                    "text-align": "left",
                    "margin": "5px",
                    "--hover-color": "#f0f0f0",
                },
                "icon": {
                    "font-size": "18px",
                    "color": "#008cba",
                },
                "container": {
                    "padding": "5px",
                    "background-color": "transparent"
                },
                "nav-link-selected": {
                    "background-color": "rgba(0, 123, 255, 0.15)",
                    "font-weight": "bold",
                    "color": "#008cba",
                },
            },
        )

    # --- Load Pages Based on Sidebar Selection ---
    if page == "Home Page":
        from Home_page import main as show_home_page
        show_home_page()

    elif page == "Run Tasks":
        from Components.Run_mc_parser import main as show_mc_parser_page
        show_mc_parser_page()

# --- Main Function ---
def main():
    st.set_page_config(
        page_title="Missed Cleavage Parser",
        #page_icon="🧬",
        layout="wide",
        initial_sidebar_state="collapsed",
    )

    display_sidebar_and_pages()

if __name__ == '__main__':
    main()
