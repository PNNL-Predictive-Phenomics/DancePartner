import os
import pandas as pd
import re
import networkx as nx

def __get_ome_df(ome_path: str, 
                 delim: str = ",", 
                 ome_type: str = None,
                 lower: bool = True):
    '''
    A support function to pull an ome files, and parse the file to be a pandas dataframe 

    Parameters
    ----------
    omes_path
        The path to the omes file formatted with the first column as the ID and the second as the Synonyms
    
    delim
        The delimiter for the file. Default is a ',' for a csv. 

    ome_type
        An optional identifier for the ome type

    lower
        If True, convert all text to lowercase before processing.

    Return
    ------
    A three column file with the ID, Synonym, and Type (lipid, metabolite, or gene product)
    '''

    # Read the ome csv
    ome = pd.read_csv(ome_path, sep = delim)

    # Instantiate the ome dictionary
    ome_dict = {}

    for row in range(len(ome)):
        
        terms = str(ome["Synonyms"][row]).split("; ")
        if lower:
            terms = [re.sub(r'[^a-zA-Z0-9]', '', term.strip().lower()) for term in terms]
        else:
            terms = [re.sub(r'[^a-zA-Z0-9]', '', term.strip()) for term in terms]

        if not isinstance(terms, list):
            terms = list(terms)

        for do_not_use in ["", "nan"]:
            if do_not_use in terms:
                terms.remove(do_not_use)

        # Only return if there is one match
        if len(terms) > 0:
            ome_dict[ome.iloc[row, 0]] = list(set(terms))

    # Make dictionary
    df = pd.DataFrame(ome_dict.items()).explode(1).rename({0: "ID", 1: "Synonym"}, axis = 1)
    if ome_type is not None:
        df["Type"] = ome_type

    return df

def __get_label_priority(ome_list: list):
    '''
    A support function to collapse a list omes to a single entry, rank the order by lipid, metabolite, 
    and gene product
    '''
    
    if "lipid" in ome_list:
        return "lipid"
    elif "metabolite" in ome_list:
        return "metabolite"
    else:
        return "gene product"

def make_synonym_table(omes_folder: str, 
                       proteome_filename: str = None, 
                       genome_filename: str = None, 
                       min_length: int = 3,
                       drop_numerics: bool = True,
                       output_directory: str = None,
                       lower: bool = True):
    '''
    Generate a complete synonym table with the DancePartner Supergroup ID, a synonym, an ID, and the type.
    The omes folder must have LipidMaps_Lipdome.csv, CHEBI_Metabolome.txt, and stop_words_english.txt
    
    Parameters
    ----------
    omes_folder
        Path to the omes folder which must hold the "LipidMaps_Lipidome.csv" and "CHEBI_Metabolome.txt" files. Required. 
    
    proteome_filename
        Name of the proteome file within the omes folder. Use the full file name. Optional. Use "pull_proteome" to make the
        proper proteome file format for synonym detection.

    genome_filename
        Name fo the genome file within the omes folder. Use the full file name. Optional. Use "pull_genome_from_GTF" to make
        the proper genome file format for synonym detection.
    
    min_length 
        Minimum number of characters in a term. Default is 3.  

    drop_numerics
        If True, drop synonym terms that are numeric. Default is True.

    output_directory
        A path to a directory where results will be written. Optional.

    lower
        If True, convert synonym text to lowercase before cleaning and matching.
        If False, preserve the original letter case. Default is True.
    
    Returns
    -------
    A four column file with the DancePartner Supergroup ID, a synonym, an ID, and the type. Supergroups are made to group 
    synonyms and IDs that refer to the same biomolecule.
    '''

    ## Read and gather files---------------------------------------------------------------------------------------------

    # Parse lipidome
    lipidome = __get_ome_df(os.path.join(omes_folder, "LipidMaps_Lipidome.csv"), ome_type = "lipid", lower = lower)

    # Parse metabolome
    metabolome = __get_ome_df(os.path.join(omes_folder, "CHEBI_Metabolome.txt"), "\t", ome_type = "metabolite", lower = lower)

    # Start the universal ome
    ome = pd.concat([lipidome, metabolome])

    # Add the proteome 
    if proteome_filename is not None:
        proteome = __get_ome_df(os.path.join(omes_folder, proteome_filename), "\t", "gene product", lower = lower)
        ome = pd.concat([ome, proteome])

    # Add the genome, if applicable
    if genome_filename is not None:
        genome = __get_ome_df(os.path.join(omes_folder, genome_filename), "\t", "gene product", lower = lower)
        ome = pd.concat([ome, genome])

    # Read the stop words file
    stopwords = pd.read_csv(os.path.join(omes_folder, "stop_words_english.txt"))["stopwords"].tolist()

    ## Clean terms------------------------------------------------------------------------------------------

    # Clean terms that are too short or are a stopword
    unique_terms = list(set(ome["Synonym"].to_list()))
    cleaned_terms = [term for term in unique_terms if len(term) >= min_length and term not in stopwords]
    ome = ome[ome["Synonym"].isin(cleaned_terms)].reset_index(drop = True)

    # Remove numerics if requested
    def is_not_numeric(value):
        """Checks if a value is numeric or can be converted to a numeric type."""
        try:
            float(value) 
            return False
        except (ValueError, TypeError):
            return True
        
    if drop_numerics:
        numeric_terms = [term for term in unique_terms if is_not_numeric(term)]
        ome = ome[ome["Synonym"].isin(numeric_terms)].reset_index(drop = True)

    ## Define supergroups-----------------------------------------------------------------------------------

    # Build unique Synonym Group Names
    Synonym_Groups = ome.groupby("Synonym").agg({"ID": list, "Type": list}).reset_index()
    Synonym_Groups["SynGroup"] = ["SynGroup" + str(x) for x in range(len(Synonym_Groups))]

    # Build unique ID Group Names
    ID_Groups = ome.groupby("ID").agg({"Synonym": list, "Type": list}).reset_index()
    ID_Groups["IDGroup"] = ["IDGroup" + str(x) for x in range(len(ID_Groups))]

    # Define and add supergroups using unique components in a network
    ome = pd.merge(pd.merge(ome, Synonym_Groups[["Synonym", "SynGroup"]]), ID_Groups[["ID", "IDGroup"]])

    # Build the network of connected groups
    Supergroups = nx.Graph()
    Supergroups.add_edges_from(ome[["SynGroup", "IDGroup"]].itertuples(index=False, name=None))

    # Find connected nodes
    connected_components = list(nx.connected_components(Supergroups))

    # Create a new supergroup mapping
    supergroup_mapping = {}
    for i, component in enumerate(connected_components):
        for node in component:
            supergroup_mapping[node] = i

    # Apply supergroup mapping
    ome["DancePartnerID"] = ome["SynGroup"].map(supergroup_mapping).combine_first(ome["IDGroup"].map(supergroup_mapping))

    # Build out the synonym table
    synonym_table = ome[["DancePartnerID", "ID", "Synonym", "Type"]].groupby("DancePartnerID").agg({"ID": list, "Synonym": list, "Type": list})

    # Set a typing priority
    synonym_table["Type"] = [__get_label_priority(ome_list) for ome_list in synonym_table["Type"]]
    
    # Unnest the lists and make the table in long format
    synonym_table = synonym_table.explode(["ID", "Synonym"]).reset_index()

    # Output file or return it
    if output_directory is not None:
        synonym_table.to_csv(os.path.join(output_directory, "synonym_table.txt"), index=False, sep = "\t")
    else:
        return(synonym_table)


def map_synonyms(term_list: list[str], 
                 synonym_table: pd.DataFrame,
                 add_missing: bool = False, 
                 lower: bool = True,
                 output_directory: str = None):
    '''
    Map synonyms to IDs

    Parameters
    ----------
    term_list 
        List of terms to map to lipidome, metabolome, and proteome. 

    synonym_table
        A synonym table as created in make_synonym_table
        
    add_missing
        If True, add terms that were not mapped to synonyms. Default is False.

    lower
        If True, convert input terms to lowercase before cleaning and matching.
        If False, preserve the original letter case. Default is True.
    
    output_directory
        A path to a directory where results will be written. Optional.
    
    Returns
    -------
    A table with the DancePartner ID group, a synonym, its ID, and the type (gene product, lipid, metabolite)
    '''

    # Format terms
    if lower:
        term_list = [re.sub(r'[^a-zA-Z0-9]', '', term.strip().lower()) for term in term_list]
    else:
        term_list = [re.sub(r'[^a-zA-Z0-9]', '', term.strip()) for term in term_list]

    # Remove the ID column
    synonym_table = synonym_table[["DancePartnerID", "Synonym", "Type"]]

    # Subset the large table. Drop duplicates. 
    SynonymTable = synonym_table[synonym_table["Synonym"].isin(term_list)].drop_duplicates().reset_index(drop = True)

    # Add missing if applicable
    if add_missing:
        missing = [term for term in term_list if term not in SynonymTable["Synonym"].tolist()]
        if len(missing) > 0:
            SynonymTable = pd.concat([
                SynonymTable,
                pd.DataFrame({"DancePartnerID": ["" for x in range(len(missing))],
                              "Synonym": missing, 
                              "Type": ["" for x in range(len(missing))]})
            ]).reset_index(drop = True)
    
    if output_directory is not None:
        SynonymTable.to_csv(os.path.join(output_directory, "synonym_table.txt"), index=False, sep = "\t")
    else:
        return(SynonymTable)
