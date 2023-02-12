### Toxicity filter 2.0 - 2023

# import packages
import streamlit as st
import pandas as pd
import numpy as np
import json
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, DataReturnMode, JsCode
import streamlit.components.v1 as components
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")

with open('dict_mp.json', 'r') as fp:
    dict_mp = json.load(fp)

sig_phenotypes = pd.read_csv('sig_phenotypes.csv')
all_phenotypes = pd.read_feather('all_phenotypes.feather')

# define functions
@st.cache
def get_association_data():
    df = pd.read_csv('results_0202023_confirmednames4streamlit.tsv', sep="\t")
    return df.set_index("gene_symbol")

def aggrid_interactive_table(df: pd.DataFrame):
    options = GridOptionsBuilder.from_dataframe(  df,  enableRowGroup=True,
                                                  enableValue=True,
                                                  enableRangeSelection=True,
                                                  use_checkbox=True,
                                                  enableCellTextSelection=True,
                                                  ensureDomOrder=True  )
    # options.configure_side_bar()
    options.configure_auto_height(True)

    options.configure_selection("single")
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        enableRangeSelection=True,
        enableCellTextSelection=True,
        ensureDomOrder=True,
        gridOptions=options.build(),
        custom_css={
            ".ag-row-hover": {
                "background-color": "‘#ff0000"}},
        theme="alpine",
       # update_mode=GridUpdateMode.MODEL_CHANGED,
        allow_unsafe_jscode=True
    )
    return selection

def aggrid_interactive_table2(df: pd.DataFrame):
    options = GridOptionsBuilder.from_dataframe(  df,  enableRowGroup=True,
                                                  enableValue=True,
                                                  enablePivot=True, use_checkbox=True )
    options.configure_side_bar()
    options.configure_auto_height(False)
    # options.configure_selection("single", pre_selected_rows = [0])
    options.configure_selection("single")
    cellsytle_jscode = JsCode("""
    function(params) {
        if (params.value < 0.5) {
            return {
                'color': 'darkgreen',
                'backgroundColor': 'brighter'
            }
        } else {
            return {
                'color': 'red',
                'backgroundColor': 'white'
            }
        }
    };
    """)
    options.configure_column("Score", cellStyle=cellsytle_jscode)
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        enableCellTextSelection=True,
        fit_columns_on_grid_load=True,
        sizeColumnsToFit=True,
        ensureDomOrder=True,
        gridOptions=options.build(),
        theme="alpine",
        # theme="light",
       # update_mode=GridUpdateMode.MODEL_CHANGED,
        allow_unsafe_jscode=True
    )
    return selection

def aggrid_interactive_table_mp(df: pd.DataFrame):
    options = GridOptionsBuilder.from_dataframe(  df,  enableRowGroup=True,
                                                  enableValue=True,
                                                  enablePivot=True, use_checkbox=True )
    options.configure_side_bar()
    options.configure_auto_height(False)
    # options.configure_selection("single", pre_selected_rows = [0])
    options.configure_selection("single")
    cellsytle_jscode = JsCode("""
    function(params) {
        if (params.value == "significant") {
            return {
                'color': 'red',
                'backgroundColor': 'brighter'
            }
        } 
        if (params.value == "not tested") {
            return {
                'color': 'grey',
                'backgroundColor': 'brighter'
            }
        } 
        
        else {
            return {
                'color': 'darkgreen',
                'backgroundColor': 'white'
            }
        }
    };
    """)
    options.configure_column("Status", cellStyle=cellsytle_jscode)
    selection = AgGrid(
        df,
        enable_enterprise_modules=True,
        enableCellTextSelection=True,
        ensureDomOrder=True,
        sizeColumnsToFit=True,
        fit_columns_on_grid_load=True,
        gridOptions=options.build(),
        theme="alpine",
        # theme="light",
       # update_mode=GridUpdateMode.MODEL_CHANGED,
        allow_unsafe_jscode=True
    )
    return selection

############################### The body of the app
# Protein selection

with st.sidebar:
    df = get_association_data()
    list_proteins = list(set(list(df.index)))
    # list_proteins = list(set(list_proteins))
    list_proteins.sort()
    protein = st.selectbox(
        "Choose the protein", list_proteins
    )
    if not protein:
        st.error("Please select at least one protein.")
    st.write("")
    st.write("")
    st.write("")
    st.write("")
    st.write("")
    st.write("")
        
    st.caption(
        "The detailed info on the established prediction approach could be found on Confluence [page](https://insilicoteam.atlassian.net/l/cp/Jwy5uFHG)")

    # st.caption("##### Toxicity data on 7851 adverse events are available for 911 proteins through the association with 2005 drugs")




st.title(f'Adverse effects predictor for {protein}')
df_selected_protein = df.loc[[protein]]

col1, col2 = st.columns([8, 3])

if protein:
    with col1:
        st.write("##### ADRs predictions")
        df_adrs = df_selected_protein.T.reset_index()
        df_adrs.columns = ["Organ system","Score"]
        new_col = []
        for i in df_adrs.Score.tolist():
            if i.endswith("confirmed"):
                new_col.append(i.split("_")[0] + " (confirmed)")
            else:
                new_col.append(i)
        df_adrs["Score"] = new_col
        selection = aggrid_interactive_table2(df=df_adrs)

sig_phenotypes_protein = sig_phenotypes[sig_phenotypes["gene_symbol"] == protein]
all_phenotypes_protein = all_phenotypes[all_phenotypes["gene_symbol"] == protein]


mp_list = []
mp_status = []
for key , values in zip(dict_mp.keys(),dict_mp.values()) :
        mp_list.append(key)
        for val in values:
            if val in sig_phenotypes_protein.top_level_mp_term_name.unique().tolist():
                mp_status.append("significant")
                break
            if val in list(set(all_phenotypes_protein.top_level_mp_term_name.unique().tolist()) - set(sig_phenotypes_protein.top_level_mp_term_name.unique())):
                mp_status.append("not significant")
                break
            else:
                mp_status.append("not tested")
                break

mp_df = pd.DataFrame(zip(mp_list, mp_status), columns = ["Phenotype","Status"])

if protein:
    with col1:
        st.write("##### Phenotypes data")
        if list(set(mp_df.Status)) == ["not tested"]:
            st.write("Phenotypes are not available")
        else:
            selection2 = aggrid_interactive_table_mp(df=mp_df)
show_sig_mo = st.checkbox(f"Show full data of the significant phenotypes")

if show_sig_mo:
    # sig_phenotypes_protein.drop(columns = ["gene_symbol"])
    st.table(sig_phenotypes_protein)


    @st.experimental_memo
    def convert_df(df):
        return df.to_csv(index=False).encode('utf-8')


    csv = convert_df(sig_phenotypes_protein)

    st.download_button(
        "Press to download the table",
        csv,
        f"{protein}_significant_mouse_phenotypes.csv",
        "text/csv",
        key='download-csv'
    )