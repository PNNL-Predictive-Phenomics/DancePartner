import pandas as pd

def deduplicate_papers(pubmed_path = None, scopus_path = None, osti_path = None):
    '''
    This function takes up to 3 tables from LitPortal, removes duplicate publications, and gives
    the deduplicated table to the pull_papers function to properly pull the data. These files are
    all exports from LitPortal.

    pubmed_path: (str) The path to the LitPortal file exported from a PubMed search. Optional.
    scopus_path: (str) The path to the LitPortal file exported from a Scopus search. Optional.
    osti_path: (str) The path to the LitPortal file exported from an OSTI search. Optional.  
    '''
    
    ################
    ## READ FILES ##
    ################

    def read_csv_or_txt(file_path):
        '''A simple function to read either a csv or txt file'''
        file_path = str(file_path)
        if ".csv" in file_path:
            return(pd.read_csv(file_path))
        elif ".txt" in file_path:
            return(pd.read_table(file_path))

    # Read each file. If OSTI is provided, fix the DOI annotation, remove patents, and blank abstracts. 
    pubmed, scopus, osti = None, None, None
    if pubmed_path is not None:
        pubmed = read_csv_or_txt(pubmed_path)
    if scopus_path is not None:
        scopus = read_csv_or_txt(scopus_path)
    if osti_path is not None:
        osti = read_csv_or_txt(osti_path)
        osti["DOI"] = osti["DOI"].str.replace("https://doi.org/", "")
        osti = osti[(osti["Publication Type"] != "Patent") & (osti["Abstract"] != "") & (osti["Abstract"] != "Not provided.")]

    ##########################
    ## PULL DOIs AND TITLES ##
    ##########################

    def doi_pull(summary_table, table_name):
        '''
        This function pulls DOIs for papers from the LitPortal table 

        summary_table: (pandas DataFrame) The table exported from LitPortal
        table_name: (str) The desired name of the table in the column
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
    
    def text_clean(string):
        '''
        This function cleans up the text of the input string, removing nonalphanumerics.
        This is essential to proper collapsing papers based on title. 

        string: (str) The input string to clean
        
        '''
        return(''.join(char for char in string.lower() if char.isalnum() or char == " "))
    

    def title_pull(summary_table, table_name):
        '''
        This function pulls titles for papers from the LitPortal table for whenever
        a DOI is not available. 

        summary_table: (pandas DataFrame) The table exported from LitPortal
        table_name: (str) The desired name of the table in the column
        '''
        if summary_table is None:
            return None
        else:
            return(
                pd.DataFrame({
                    table_name: [str(x) for x in summary_table["OriginId"][summary_table["DOI"].isna()]],
                    "Title": [text_clean(x) for x in summary_table["Title"][summary_table["DOI"].isna()]]
                })
            )
    
    def table_merge(pubmed_info, scopus_info, osti_info):
        '''
        Merge DOI or Title tables from doi_pull and title_pull. This function allows
        for proper merging, for any combination of PubMed, Scopus, or OSTI LitPortal tables.

        pubmed_info: (DataFrame) The output from doi_pull or title_pull from pubmed
        scopus_info: (DataFrame) The output from doi_pull or title_pull from scopus
        osti_info: (DataFrame) The output from doi_pull and title_pull from osti        
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

    # Create the table of papers 
    paper_table = pd.concat([
        table_merge(doi_pull(pubmed, "pubmed"), doi_pull(scopus, "scopus"), doi_pull(osti, "osti")),
        table_merge(title_pull(pubmed, "pubmed"), title_pull(scopus, "scopus"), title_pull(osti, "osti"))
    ]).reset_index(drop = True)

    # Return the paper table
    return(paper_table)
    





