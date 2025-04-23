import os
import pandas as pd
import re

def __get_ome_df(ome_path, delim = ","):
    '''
    Pull the lipidome csv or metabolome txt, and parse the file to be a pandas dataframe 
    '''

    # Read the ome csv
    ome = pd.read_csv(ome_path, sep = delim)

    # Instantiate the ome dictionary
    ome_dict = {}

    for row in range(len(ome)):
        
        terms = str(ome["Synonyms"][row]).split("; ")
        terms = [re.sub(r'[^a-zA-Z0-9]', '', term.strip().lower()) for term in terms]

        if not isinstance(terms, list):
            terms = list(terms)

        for do_not_use in ["", "nan"]:
            if do_not_use in terms:
                terms.remove(do_not_use)

        # Is iloc ok here?
        if len(terms) > 0:
            ome_dict[ome.iloc[row, 0]] = list(set(terms))

    return pd.DataFrame(ome_dict.items()).explode(1).rename({0: "ID", 1: "Synonym"}, axis = 1)

def list_synonyms(omes_folder, proteome_filename, min_length = 3):
    '''
    List all possible synonyms to match.
    
    Args:
        omes_folder (String): Path to the omes folder. Required. 
        proteome_filename (String): Name of the proteome file within the omes folder. Use the full file name. Required.
        min_length (Integer): Minimum number of characters in a term. Default is 3.  
    '''

    # Parse lipidome
    lipidome = __get_ome_df(os.path.join(omes_folder, "LipidMaps_Lipidome.csv"))

    # Parse metabolome
    metabolome = __get_ome_df(os.path.join(omes_folder, "CHEBI_Metabolome.txt"), "\t")

    # Parse proteome
    proteome = __get_ome_df(os.path.join(omes_folder, proteome_filename), "\t")

    # List terms of interest
    toi = lipidome["Synonym"].tolist()
    toi.extend(metabolome["Synonym"].tolist())
    toi.extend(proteome["Synonym"].tolist())

    # Remove stop words 
    stopwords = pd.read_csv(os.path.join(omes_folder, "stop_words_english.txt"))["stopwords"].tolist()
    toi = [term for term in toi if term not in stopwords ]
    toi = [term for term in toi if len(term) >= min_length]

    # Return terms of interest
    return toi


def map_synonyms(term_list, omes_folder, proteome_filename, add_missing = False, output_directory = None):
    '''
    Map synonyms to IDs in the order of lipids, metabolites, and finally gene products. 

    Args:
        term_list (List[String]): List of terms to map to lipidome, metabolome, and proteome. 
        organism_id (Integer): The UniProt organism ID.
        omes_folder (String): Path to the omes folder. Required. 
        proteome_filename (String): Name of the proteome file within the omes folder. Use the full file name. Required.
        add_missing (Optional, Logical): If True, add terms that weren't mapped to synonyms. Optional.
        output_directory (Optional, String): A path to a directory for where to write results to.
    '''

    # Format terms
    term_list = [re.sub(r'[^a-zA-Z0-9]', '', term.strip().lower()) for term in term_list]

    # Start data.frame to hold all terms
    terms = pd.DataFrame({"Synonym": term_list})

    #################
    ## FIND LIPIDS ##
    #################

    # Parse lipidome
    lipidome = __get_ome_df(os.path.join(omes_folder, "LipidMaps_Lipidome.csv"))

    # Map lipids
    found_lipids = terms.merge(lipidome)
    found_lipids["Type"] = "lipid"

    ######################
    ## FIND METABOLITES ##
    ######################

    # Parse metabolome
    metabolome = __get_ome_df(os.path.join(omes_folder, "CHEBI_Metabolome.txt"), "\t")

    # Map metabolites
    found_metabolites = terms.merge(metabolome)
    found_metabolites["Type"] = "metabolite"

    ########################
    ## FIND GENE PRODUCTS ##
    ########################

    # Parse proteome
    proteome = __get_ome_df(os.path.join(omes_folder, proteome_filename), "\t")

    # Map proteins
    found_proteins = terms.merge(proteome)
    found_proteins["Type"] = "gene product"

    ################
    ## PAIR TERMS ##
    ################

    # Merge and take first entry per synonym
    SynonymTable = pd.concat([
        found_lipids,
        found_metabolites,
        found_proteins
    ]).groupby('Synonym').first().reset_index().dropna()

    # Add missing if applicable
    if add_missing:
        missing = [term for term in term_list if term not in SynonymTable["Synonym"].tolist()]
        if len(missing) > 0:
            SynonymTable = pd.concat([
                SynonymTable,
                pd.DataFrame({"Synonym":missing, "ID": ["" for x in range(len(missing))], "Type": ["" for x in range(len(missing))]})
            ]).reset_index(drop = True)
    
    if output_directory is not None:
        SynonymTable.to_csv(os.path.join(output_directory, "synonym_table.txt"), index=False, sep = "\t")
    else:
        return(SynonymTable)
