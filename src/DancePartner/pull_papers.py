import os
import datetime
import pandas as pd
from .pull_pubmed import __pull_pubmed
from .pull_scopus import __pull_scopus
from .pull_osti import __pull_osti

def pull_papers(output_directory: str, 
                pubmed_ids: list[str] = None, 
                scopus_ids: list[str] = None, 
                osti_ids: list[str] = None, 
                deduped_table: pd.DataFrame = None,
                type: str = "both", 
                include_summary_file: bool = True, 
                tarball_path: str = None, 
                scopus_api_key: str = None):
    """ 
    Given a list of IDs referencing a literature database, pull available text, prioritizing full text whenever available, then titles and abstracts. 
    A summary file of what was pulled is also generated. 
    
    Parameters
    ----------
    output_directory
        A string indicating the directory path for where to write results.

    pubmed_ids
        A list of PubMed IDs as strings. Only one of pubmed_ids, scopus_ids, or osti_ids can be provided. If wanting to use multiple databases, upload a dedeuped_table.

    scopus_ids
        A list of DOIs. Only one of pubmed_ids, scopus_ids, or osti_ids can be provided. If wanting to use multiple databases, upload a dedeuped_table.

    osti_ids
        A list of OSTI IDs. A list of PubMed IDs. Only one of pubmed_ids, scopus_ids, or osti_ids can be provided. If wanting to use multiple databases, upload a dedeuped_table.
    
    deduped_table
        A pandas DataFrame o deduplicated table of papers from deduplicate_papers. pubmed_ids, scopus_ids, and osti_ids should all be None to use a deduped table. 
    
    type 
        Either "full text" to pull only full text, "abstract" to pull only abstracts, or "both" to first prioritize full text, and then prioritize abstracts. 
        
    include_summary_file
        A boolean where True will write a summary .txt file desbring number of papers found from each pull_ranking method.

    tarball_path 
        An optional string for the path where to write the (large) tarball files to. Can also be used to specify a tarball path where a previous    
        function run may have saved articles to, which can reduce run time.

    scopus_api_key 
        A string API key for Scopus-Elselvier. Only needed when pulling papers from Scopus. See https://dev.elsevier.com/.

    
    Returns
    -------
        Extracted papers as text files in folders, with a summary file

    """

    # First detect IDs to use
    ID_List = [pubmed_ids is not None, scopus_ids is not None, osti_ids is not None]

    # If a database is specified, use the specified database
    if sum(ID_List) > 1:
        raise Exception("Detected multiple sets of IDs. If using multiple databases, please use deduplicate_papers() first.")
    
    # Pull papers by the specification
    if pubmed_ids is not None:
        pubmed_ids = [str(x) for x in pubmed_ids]
        total_papers = len(pubmed_ids)
        counts_dictionary = __pull_pubmed(ids = pubmed_ids, output_directory = output_directory, type = type, tarball_path = tarball_path)
    if scopus_ids is not None:
        total_papers = len(scopus_ids)
        counts_dictionary = __pull_scopus(ids = scopus_ids, output_directory = output_directory, type = type, scopus_api_key = scopus_api_key)
    if osti_ids is not None:
        osti_ids = [str(x) for x in osti_ids]
        total_papers = len(osti_ids)
        counts_dictionary = __pull_osti(ids = osti_ids, output_directory = output_directory, type = type)

    # Otherwise, this is the deduplicated example
    if deduped_table is not None:

        counts_dictionary = {"full": [], "abstract": []}
        total_papers = len(deduped_table)

        # Now to iterate through the deduped table
        for row in range(len(deduped_table)):

            # Pull each value
            pubmed = str(deduped_table.at[row, "pubmed"])
            scopus = str(deduped_table.at[row, "scopus"])
            osti = str(deduped_table.at[row, "osti"])

            # Determine if full text pubmed is an option. If so, try to pull full text. If not, move on.
            if pubmed != "nan":
                single_counts = __pull_pubmed(ids = [pubmed], output_directory = output_directory, type = "full", tarball_path = tarball_path)
                if single_counts is not None:
                    counts_dictionary["full"].extend(single_counts["full"])
                    continue

            # Determine if full text from scopus is an option. If so, pull full text. If not, move on.
            if scopus != "nan":
                single_counts = __pull_scopus(ids = [scopus], output_directory = output_directory, type = "full", scopus_api_key = scopus_api_key)
                if single_counts is not None:
                    counts_dictionary["full"].extend(single_counts["full"])
                    continue
            
            # Determine if full text from OSTI is an option. If so, pull full text. If not, move on.
            if osti != "nan":
                single_counts = __pull_osti(ids = [osti], output_directory = output_directory, type = "full")
                if single_counts is not None:
                    counts_dictionary["full"].extend(single_counts["full"])
                    continue

            # Determine if abstract from pubmed is an option. If not, move on.
            if pubmed != "nan":
                single_counts = __pull_pubmed(ids = [pubmed], output_directory = output_directory, type = "abstract", tarball_path = tarball_path)
                if single_counts is not None:
                    counts_dictionary["abstract"].extend(single_counts["abstract"])
                    continue

            # Determine if abstract from scopus is an option. If not, move on.
            if scopus != "nan":
                single_counts = __pull_scopus(ids = [scopus], output_directory = output_directory, type = "abstract", scopus_api_key = scopus_api_key)
                if single_counts is not None:
                    counts_dictionary["abstract"].extend(single_counts["abstract"])
                    continue
            
            # Determine if abstract from OSTI is an option. If not, move on.
            if osti != "nan":
                single_counts = __pull_osti(ids = [osti], output_directory = output_directory, type = "abstract")
                if single_counts is not None:
                    counts_dictionary["abstract"].extend(single_counts["abstract"])
                    continue
    
    # Set minimum number 
    if counts_dictionary is None:
        counts_dictionary = {"full": [], "abstract": []}
    
    # Calculate number found
    found = len(counts_dictionary["full"]) + len(counts_dictionary["abstract"])

    # Write summary file in the output_directory
    if include_summary_file:
        with open(os.path.join(output_directory, "output_summary.txt"), "w") as f:
            f.write("Output Summary for Pulling Papers\n")
            f.write("Created: " + str(datetime.datetime.now()) + "\n")
            f.write("Total Num. Articles: " + str(total_papers) + "\n")
            f.write("Total Num. Articles Found: " + str(found) + "\n")
            f.write("Number of Full Text: " + str(len(counts_dictionary["full"])) + "\n")
            f.write("Number of Title & Abstracts: " + str(len(counts_dictionary["abstract"])) + "\n")
            f.write("Number Missing: " + str(total_papers - found) + "\n")
    return None