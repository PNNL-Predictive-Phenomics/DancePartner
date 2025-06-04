import pandas as pd

def _read_csv_or_txt(file_path: str):
    '''
    A simple function to read either a csv or txt file
    
    Parameters
    ---------
    file_path
        A string indicating the path to the file
    
    Returns
    -------
        A pandas DataFrame
    '''
    file_path = str(file_path)
    if ".csv" in file_path:
        return(pd.read_csv(file_path))
    elif ".txt" in file_path:
        return(pd.read_table(file_path))
    
def _text_clean(string: str):
    '''
    This function cleans up the text of the input string, removing nonalphanumerics.
    This is essential to proper collapsing papers based on title. 

    Parameters
    ---------
    string
        The input string to clean
    
    Return
    ------
        The string with non-alphanumerics removed (with the exception of spaces)

    '''
    return(''.join(char for char in string.lower() if char.isalnum() or char == " "))


def _table_merge(pubmed_info: pd.DataFrame, scopus_info: pd.DataFrame, osti_info: pd.DataFrame):
    '''
    Merge DOI or Title tables from doi_pull and title_pull. This function allows for proper merging, for any 
    combination of PubMed, Scopus, or OSTI information.

    Parameters
    ----------
    pubmed_info
        The output from doi_pull or title_pull from pubmed

    scopus_info
        The output from doi_pull or title_pull from scopus

    osti_info
        The output from doi_pull and title_pull from osti

    Returns
    -------
        A pandas DataFrame        
    '''

    # Merge as long as the dataframe is not none
    summary_tables = [pubmed_info, scopus_info, osti_info]
    exists = [pubmed_info is not None, scopus_info is not None, osti_info is not None]
    position = [pos for pos in range(len(exists)) if exists[pos] == True]

    # If length 1, no merging is required. If length 2, merge the two together. If length
    # 3, merge all three togeter. 
    if sum(exists) == 1: 
        merged = summary_tables[position[0]]
    elif sum(exists) == 2:
        merged = pd.merge(summary_tables[position[0]], summary_tables[position[1]], how = "outer")
    else: 
        merged = pubmed_info.merge(scopus_info, how = "outer").merge(osti_info, how = "outer")

    return(merged)

def deduplicate_papers(pubmed_path: str = None, scopus_path: str = None, osti_path: str = None):
    '''
    Deduplicate papers across databases. 
    
    Parameters
    ----------
    pubmed_path
        The path to the PubMed export of paper information. To obtain, enter the query, hit search, hit save, select "All results" and "csv".
    
    scopus_path
        The path to the Scopus export of paper information. To obtain, enter the query, hit search, hit export, select "CSV" and keep all defaults checked.
    
    osti_path
        The path to the OSTI export of paper information. To obtain, enter the query, hit search, and save results as a "CSV"
    
    Returns
    -------
        A table with deduplicated papers
    '''

    # If provided, pull essential information from each text file and clean paper titles.
    pubmed, scopus, osti = None, None, None
    if pubmed_path is not None:
        pubmed = _read_csv_or_txt(pubmed_path)
        pubmed = pubmed[["PMID", "Title", "DOI"]].rename({"PMID":"pubmed"}, axis = 1).sort_values(by = "DOI")
        pubmed["Title"] = [_text_clean(x) for x in pubmed["Title"]]
    if scopus_path is not None:
        scopus = _read_csv_or_txt(scopus_path)
        scopus = scopus[["Title", "DOI", "EID"]].rename({"EID":"scopus"}, axis = 1).sort_values(by = "DOI")
        scopus["Title"] = [_text_clean(x) for x in scopus["Title"]]
    if osti_path is not None:
        osti = _read_csv_or_txt(osti_path)
        osti = osti[["TITLE", "DOI", "OSTI_IDENTIFIER"]].rename({"TITLE":"Title", "OSTI_IDENTIFIER":"osti"}, axis = 1).sort_values(by = "DOI")
        osti["Title"] = [_text_clean(x) for x in osti["Title"]]

    # Merge tables 
    return _table_merge(pubmed, scopus, osti).reset_index(drop = True)[["Title", "DOI", "pubmed", "scopus", "osti"]]

def litportal_deduplicate_papers(pubmed_path: str = None, scopus_path: str = None, osti_path: str = None):
    '''
    We recommend using the "generic" depulicate papers, as this special function is for an internal tool
    called "LitPortal." This function takes up to 3 tables from LitPortal, removes duplicate publications, and gives
    the deduplicated table to the pull_papers function to properly pull the data. These files are
    all exports from LitPortal.

    Parameters
    ----------
    pubmed_path
        The path to the LitPortal file exported from a PubMed search. Optional.
    
    scopus_path
        The path to the LitPortal file exported from a Scopus search. Optional.
    
    osti_path
        The path to the LitPortal file exported from an OSTI search. Optional.  
    
    Returns
    -------
        A table with deduplicated papers
    '''
    
    ################
    ## READ FILES ##
    ################

    # Read each file. If OSTI is provided, fix the DOI annotation, remove patents, and blank abstracts. 
    pubmed, scopus, osti = None, None, None
    if pubmed_path is not None:
        pubmed = _read_csv_or_txt(pubmed_path)
    if scopus_path is not None:
        scopus = _read_csv_or_txt(scopus_path)
    if osti_path is not None:
        osti = _read_csv_or_txt(osti_path)
        osti["DOI"] = osti["DOI"].str.replace("https://doi.org/", "")
        osti = osti[(osti["Publication Type"] != "Patent") & (osti["Abstract"] != "") & (osti["Abstract"] != "Not provided.")]

    ##########################
    ## PULL DOIs AND TITLES ##
    ##########################

    def doi_pull(summary_table: pd.DataFrame, table_name: str):
        '''
        This function pulls DOIs for papers from the LitPortal table 

        Parameters
        ----------
        summary_table
            The table exported from LitPortal
        
        table_name
            The desired name of the table in the column

        Returns
        -------
            A pandas DataFrame with DOIs 
        '''
        if summary_table is None:
            return None
        else:
            return(
                pd.DataFrame({
                    table_name : [str(x) for x in summary_table["OriginId"][summary_table["DOI"].isna() == False]],
                    "DOI": summary_table["DOI"][summary_table["DOI"].isna() == False]
                })
            )

    def title_pull(summary_table: pd.DataFrame, table_name: str):
        '''
        This function pulls titles for papers from the LitPortal table for whenever a DOI is not available. 

        Parameters
        ----------
        summary_table
            The table exported from LitPortal
        
        table_name
            The desired name of the table in the column

        Returns
        -------
            A pandas DaataFrame with paper titles
        '''
        if summary_table is None:
            return None
        else:
            return(
                pd.DataFrame({
                    table_name: [str(x) for x in summary_table["OriginId"][summary_table["DOI"].isna()]],
                    "Title": [_text_clean(x) for x in summary_table["Title"][summary_table["DOI"].isna()]]
                })
            )

    # Create the table of papers 
    paper_table = pd.concat([
        _table_merge(doi_pull(pubmed, "pubmed"), doi_pull(scopus, "scopus"), doi_pull(osti, "osti")),
        _table_merge(title_pull(pubmed, "pubmed"), title_pull(scopus, "scopus"), title_pull(osti, "osti"))
    ]).reset_index(drop = True)[["Title", "DOI", "pubmed", "scopus", "osti"]]

    # Return the paper table
    return(paper_table)
    





