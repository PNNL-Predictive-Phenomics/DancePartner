import json
import re
import os
import pandas as pd
import ppi_pipeline as ppi

def parse_lipid_relationships(json_file, output_directory):
    '''
    Pull relationships from the LipidMaps relationships json

    Args:
        json_file (String): Path to the json file. Download from here: https://www.lipidmaps.org/resources/tools/reactions
        output_directory (String): List the folder you would like to write the "LipidMaps_Relationships.txt" file to
    Returns:
        A dataframe denoting relationships in 7 columns: Synonym1, ID1, Type1, Synonym1, ID2, Type2, Source
    '''

    # Parse the json file
    with open(json_file, "r") as file:
        data = json.load(file)

    # Hold the source and target information
    source = []
    target = []

    # Pull all elements
    for el in data["elements"]["edges"]:
        source.append(el["data"]["source"])
        target.append(el["data"]["target"])

    # Format source nad target
    source = [re.sub(r'[^a-zA-Z0-9]', '', term.strip().lower()) for term in source]
    target = [re.sub(r'[^a-zA-Z0-9]', '', term.strip().lower()) for term in target]

    # Map all synonyms to IDs
    all = source
    all.extend(target)
    all = list(set(all))
    syns = ppi.map_synonyms(all, "../omes", "UP000001570_proteome.txt")

    # Merge IDs
    lipid_rel = pd.DataFrame({"Source": source, "Target": target})
    lipids = pd.merge(lipid_rel.rename(columns = {"Source":"Synonym"}), syns).rename(columns = {"Synonym":"Synonym1", "ID":"ID1", "Type":"Type1"})
    lipids = pd.merge(lipids.rename(columns = {"Target":"Synonym"}), syns).rename(columns = {"Synonym":"Synonym2", "ID":"ID2", "Type":"Type2"})
    lipids["Source"] = "database"

    # Format file   
    lipids = lipids[["Synonym1", "ID1", "Type1", "Synonym2", "ID2", "Type2", "Source"]]

    # Export file
    outpath = os.path.join(output_directory,  "LipidMaps_Relationships.txt") 
    lipids.to_csv(outpath, index = False, sep = "\t")
    return(None)

