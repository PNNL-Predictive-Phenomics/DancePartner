import requests
from bs4 import BeautifulSoup
from io import StringIO

import pywikipathways as pwpw

import os
import io
import re

import json
import pandas as pd
import numpy as np

from .create_synonym_table import map_synonyms

pd.options.mode.chained_assignment = None 

##############################
## DATA PREPARING FUNCTIONS ##
##############################

def __parse_url(url):
    '''
    Simple parsing function to extract data from an API call

    Args:
        url (String): The url to parse
    '''
    req = requests.get(url,  timeout = 500)
    return BeautifulSoup(req.content, "html.parser")

def __upper_triangle_meshgrid(x):
    """
    Generates all combinations of x and x, and returns only the upper triangle,
    meaning the not self-self relationships and the unique relationships.

    Args:
        x (Numpy Array): A list of terms to make the meshgrid out of
    Returns:
        A pandas DataFrame of only the unique combinations
    """

    # Make a meshgrid of x and y 
    X, Y = np.meshgrid(x, x)

    # Take the upper triangle
    mask = np.triu(np.ones_like(X, dtype=bool), k=0)
    X[~mask] = np.nan

    # Drop NA values
    combos = pd.DataFrame([X.flatten(), Y.flatten()]).T.dropna()

    # Remove self mapping
    combos = combos[combos[0] != combos[1]]
    return(combos.reset_index(drop = True))

def remove_relationship_duplicates(network_table, remove_self_relationships = True): 
    '''
    Remove all duplicates from a network table.
    
    Args:
        network_table (Pandas DataFrame): Output of build_network_table, pull_protein_protein_interactions, 
           etc. Use pd.concat to concatenate multiple tables together. 
        remove_self_relationships (Logical): True to remove any relationships to self, and False to maintain them. Default is True.
    Returns:
        A file with unique interacting biomolecules
    '''

    # Remove self-relationships if possible
    if remove_self_relationships:
        network_table = network_table[network_table["ID1"] != network_table["ID2"]]

    # Remove duplicates
    network_table["Groups"] = network_table.apply(lambda row: " ".join(sorted([str(row["ID1"]), str(row["ID2"])])), axis = 1)
    network_table = network_table.drop_duplicates(subset = ["Groups"])
    return network_table[["Synonym1", "ID1", "Type1", "Synonym2", "ID2", "Type2", "Source"]].reset_index(drop = True)

##################
## INTERACTOMES ##
##################

def pull_uniprot(species_id, output_directory = None, remove_self_relationships = True, verbose = True):
    """
    Function that pulls protein-protein and protein-metabolite interactions for a species. 

    Args:
        species_id (String or Numeric): The taxon ID for the organism of interest.
        output_directory (String): Path specifying where to write the result.
        remove_self_relationships (Logical): True to remove any relationships to self, and False to maintain them. Default is True.
        verbose (Logical): Whether progress messages should be written or not. Default is False.
    Returns:
        A dataframe denoting relationships in 7 columns: Synonym1, ID1, Type1, Synonym1, ID2, Type2, Source
    """

    # Write progress message if requested
    if verbose:
        print("...pulling organism metadata")

    # Build the url
    url = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_id%2Ccc_interaction%2Ccc_cofactor&format=tsv&query=%28" + str(species_id) + "%29"

    # Read the data.frame
    interactions = pd.read_csv(StringIO(str(__parse_url(url))), sep = "\t")

    ## Protein-Protein Interactions---------------------------------------------------------------------------------

    # Write progress message if requested
    if verbose:
        print("...parsing protein-protein relationships")

    # First, build protein-protein interactions
    pre_prot_prot = interactions.dropna(subset = ["Interacts with"])[["Entry", "Interacts with"]].reset_index(drop = True).rename({"Entry":"ID1", "Interacts with":"ID2"}, axis = 1)

    # Define a list to hold protein-protein interactions 
    ppi_list = []

    for row in range(len(pre_prot_prot)):
        terms = pre_prot_prot.loc[row, "ID2"].split("; ")

        if not isinstance(terms, list):
            terms = list(terms)

        for term in terms:
            ppi_list.append({"ID1": pre_prot_prot.loc[row, "ID1"], "ID2": term})

    # Save interactions as a table
    prot_prot = pd.DataFrame(ppi_list)

    # Add missing columns
    prot_prot["Synonym1"], prot_prot["Synonym2"] = "Not needed", "Not needed"
    prot_prot["Type1"], prot_prot["Type2"] = "gene product", "gene product"
    prot_prot["Source"] = "database"
    prot_prot = prot_prot[["Synonym1", "ID1", "Type1", "Synonym2", "ID2", "Type2", "Source"]]

    ## Protein-Metabolite Interactions------------------------------------------------------------------------------

    # Write progress message if requested
    if verbose:
        print("...parsing protein-metabolite relationships")

    # First, extract unformatted protein-metabolite interactions
    pre_prot_metab = interactions.dropna(subset = ["Cofactor"])[["Entry", "Cofactor"]].reset_index(drop = True).rename({"Entry":"ID1", "Cofactor":"ID2"}, axis = 1)

    # Pull CHEBI Idnetifiers
    CHEBI_list = []

    # Extract out all CHEBI IDs
    for x in pre_prot_metab["ID2"]:
        sub_CHEBI = []
        for y in x.split("ChEBI:")[1:]:
            sub_CHEBI.append(y.split(";")[0])
        CHEBI_list.append("; ".join(sub_CHEBI))
    pre_prot_metab["ID2"] = CHEBI_list

    # Format protein-metabolite interactions
    pmi_list = []

    for row in range(len(pre_prot_metab)):
        terms = pre_prot_metab.loc[row, "ID2"].split("; ")

        if not isinstance(terms, list):
            terms = list(terms)

        for term in terms:
            pmi_list.append({"ID1": pre_prot_metab.loc[row, "ID1"], "ID2": term})
    prot_metab = pd.DataFrame(pmi_list)

    # Add missing columns 
    prot_metab["Synonym1"], prot_metab["Synonym2"] = "Not needed", "Not needed"
    prot_metab["Type1"], prot_metab["Type2"] = "gene product", "metabolite"
    prot_metab["Source"] = "database"
    prot_metab = prot_metab[["Synonym1", "ID1", "Type1", "Synonym2", "ID2", "Type2", "Source"]]

    ## Combine Datasets------------------------------------------------------------------------------

    # Concatenate tables and remove duplicates
    relationships = pd.concat([prot_prot, prot_metab])
    final_relationships = remove_relationship_duplicates(relationships, remove_self_relationships)
    final_relationships = final_relationships.dropna().reset_index(drop = True)

    # Write or return output
    if output_directory is not None:
        final_relationships.to_csv(os.path.join(output_directory, str(species_id) + "_uniprot.txt"), index=False, sep = "\t")
        return None
    else:
        return final_relationships

########################
## METABOLIC NETWORKS ##
########################

def pull_wikipathways(species_name, species_id, omes_folder, proteome_filename, output_directory = None, remove_self_relationships = True, verbose = False):
    '''
    Extract relationships from metabolic networks stored in WikiPathways

    Args:
        species_name (String): The name for the species. Select species from here: https://www.wikipathways.org/browse/organisms.html.
            Use proper Genus species format.
        species_id (String or Numeric): The taxon ID for the organism of interest. Required.
        omes_folder (String): Path to the omes folder. Required. 
        proteome_filename (String): Path to the proteom
        output_directory (String): Path specifying where to write the result.
        remove_self_relationships (Logical): True to remove any relationships to self, and False to maintain them. Default is True.
        verbose (Logical): Whether progress messages should be written or not. Default is False. 
    Returns:
        A dataframe denoting relationships in 7 columns: Synonym1, ID1, Type1, Synonym1, ID2, Type2, Source
    '''

    # Extract all pathways
    pathways = pwpw.list_pathways(species_name)

    # Make a list to hold all relationships
    relationships = []

    # Iterate through pathways 
    for p_id in pathways["id"].tolist():

        if verbose:
            print("Extracting entities for: " + p_id)

        # Build URL and read the JSON file through the API  
        try:
            url = "https://www.wikipathways.org/wikipathways-assets/pathways/" + p_id + "/" + p_id + ".json"
            data = json.loads(str(__parse_url(url)))

            # Hold entities 
            wp_data = []

            # Extract all entities 
            for entry in data["entitiesById"]:
                content = data["entitiesById"][entry]["type"]
                if "DataNode" in content:
                    wp_data.append(content)

            # Synonyms need to be parsed and collapsed 
            pre_nodes = pd.DataFrame(wp_data).drop_duplicates()

            # If there is more than 6 columns, collapse the outside columns 
            if len(pre_nodes.columns) > 6: 
                col5 = []
                for row in range(len(pre_nodes)):
                    col5.append(" & ".join(pre_nodes.iloc[row, 5:].dropna().tolist()))
                pre_nodes[5] = col5 

            # Remove any cases of column 3 having a missing value
            pre_nodes = pre_nodes.dropna(subset = [3])

            # Map Terms
            if verbose:
                print("...Mapping terms to standardized IDs")

            # Map IDs to our list of standardized terms 
            IDs = []
            Types = []

            # Reset the index of pre_nodes
            pre_nodes = pre_nodes.reset_index(drop = True)

            for row in range(len(pre_nodes)):

                terms = []
                terms.append(pre_nodes.loc[row, 2])
                terms.extend(pre_nodes.loc[row, 5].split(" & "))
                terms = [x.split(":")[-1] for x in terms if x not in ["CHEBI:"]]

                syns = map_synonyms(terms, omes_folder, proteome_filename)
                try:
                    IDs.append(syns[syns["ID"] != ""]["ID"][0])
                    Types.append(syns[syns["ID"] != ""]["Type"][0])
                except:
                    IDs.append("")
                    Types.append("")

            # Add official ID and Type
            pre_nodes[2] = IDs
            pre_nodes[3] = Types
            pre_nodes = pre_nodes[pre_nodes[2] != ""]

            if len(pre_nodes) == 0:
                continue

            # Pull node information. An edge will be drawn between every node. 
            nodes = pre_nodes.loc[:,[2,3,5]].rename({2: "ID1", 3: "Type1", 5: "Synonym1"}, axis = 1)

            # Extract upper triangle
            ut = __upper_triangle_meshgrid(nodes["ID1"]).rename({0:"ID1", 1:"ID2"}, axis = 1)

            # Left join both columns 
            relationship_table = pd.merge(ut[["ID1"]], nodes).reset_index(drop = True).join(
                pd.merge(ut[["ID2"]], nodes.rename({"ID1":"ID2", "Type1":"Type2", "Synonym1":"Synonym2"}, axis = 1))
            )
            relationship_table["Source"] = "database"
            relationships.append(relationship_table)

        except:
            print(p_id + " not found")
            continue

    # Remove duplicates
    final_relationships = pd.concat(relationships).reset_index(drop = True)
    final_relationships = remove_relationship_duplicates(final_relationships, remove_self_relationships)

    if output_directory is not None:
        final_relationships.to_csv(os.path.join(output_directory, str(species_id) + "_wikipathways.txt"), index=False, sep = "\t")
        return None
    else:
        return final_relationships
    
def pull_kegg(kegg_species_id, omes_folder, proteome_filename, output_directory = None, flatten_module = False, remove_self_relationships = True, verbose = False):
    '''
    Extract relationships from metabolic networks (modules) stored in KEGG

    Args:
        kegg_species_id (String): The name for the species. Select species from here: https://rest.kegg.jp/list/organism
        species_id (String or Numeric): The taxon ID for the organism of interest. Required.
        omes_folder (String): Path to the omes folder. Required. 
        proteome_filename (String): Path to the proteome
        output_directory (String): Path specifying where to write the result.
        flatten_module (Logical): If True, everything in a module will be considered related to everything else in a module. If False, 
            metabolic relationships as defined by KEGG will be preserved. Default: False. 
        remove_self_relationships (Logical): True to remove any relationships to self, and False to maintain them. Default is True.
        verbose (Logical): Whether progress messages should be written or not. Default is False. 
    Returns:
        A dataframe denoting relationships in 7 columns: Synonym1, ID1, Type1, Synonym1, ID2, Type2, Source
    '''

    ## Pull Organism--------------------------------------------------------------------------------------

    # If applicable, write message
    if verbose:
        print("...extracting organism information")
    
    # Pull all pathways
    species_url = "https://rest.kegg.jp/list/pathway/" + kegg_species_id
    pathways = list(set(pd.read_csv(io.StringIO(str(__parse_url(species_url))), sep = "\t", header = None)[0].to_list()))

    ## Pull Pathways---------------------------------------------------------------------------------------

    # If applicable, write message
    if verbose:
        print("...pulling pathways. ETA: " + str(len(pathways) * 1.5) + " seconds")

    # Store modules as they're pulled
    modules = []

    for pathway in pathways:

        # Construct url
        pathway_url = "https://rest.kegg.jp/get/" + pathway

        # Extract data
        text = str(__parse_url(pathway_url))

        # Set module flag to false
        module_flag = False

        # Extract modules if possible
        if "\nMODULE" in text:

            for line in text.split("\n"):

                if line.startswith("MODULE"):
                    module_flag = True

                if module_flag:
                    if "MODULE" in line:
                        modules.append(line.replace("MODULE", "").strip().split()[0].split("_")[-1])
                    elif line.startswith("  ") == False:
                        module_flag = False
                    else:
                        modules.append(line.strip().split()[0].split("_")[-1])

    modules = list(set(modules))

    ## Pull Modules---------------------------------------------------------------------------------------

    # If applicable, write message
    if verbose:
        print("...pulling modules. ETA: " + str(len(modules) * 4) + " seconds")

    # Store all relationships
    relations = []

    # Determine whether to attempt the flatten module route or not
    if flatten_module:

        if verbose:
            print("...flatten module mode selected.")

        # Define flags and list to hold biomolecules and relationships
        orthology_flag = False
        compound_flag = False
        biomolecules = []

        # Iterate through modules 
        for module in modules: 

            try:

                # Construct the module url
                module_url = "https://rest.kegg.jp/get/" + module
                text = str(__parse_url(module_url))

                if "\nORTHOLOGY" in text or "\nCompound" in text:

                    for line in text.split("\n"):

                        if line.startswith("ORTHOLOGY"):
                            orthology_flag = True
                        elif line.startswith("COMPOUND"):
                            compound_flag = True

                        if orthology_flag:
                            if "ORTHOLOGY" in line:
                                biomolecules.append(line.replace("ORTHOLOGY", "").strip().split()[1].replace(";", ""))
                            elif line.startswith("  ") == False:
                                orthology_flag = False
                            else:
                                biomolecules.append(line.strip().split()[1].replace(";", ""))

                        if compound_flag:
                            if "COMPOUND" in line:
                                biomolecules.append(line.replace("COMPOUND", "").strip().split()[1])
                            elif line.startswith("  ") == False:
                                compound_flag = False
                            else:
                                biomolecules.append(line.strip().split()[1])

                    # Collapse any potential duplicates
                    biomolecules = list(set(biomolecules))

                    if len(biomolecules) > 0:

                        # Pull node information. An edge will be drawn between every node
                        nodes = map_synonyms(biomolecules, omes_folder, proteome_filename).rename({"ID": "ID1", "Type": "Type1", "Synonym": "Synonym1"}, axis = 1)

                        # Extract upper triangle
                        ut = __upper_triangle_meshgrid(nodes["ID1"]).rename({0:"ID1", 1:"ID2"}, axis = 1)

                        # Left join both columns 
                        relationship_table = pd.merge(ut[["ID1"]], nodes).reset_index(drop = True).join(
                            pd.merge(ut[["ID2"]], nodes.rename({"ID1":"ID2", "Type1":"Type2", "Synonym1":"Synonym2"}, axis = 1))
                        )
                        relationship_table["Source"] = "database"

                        # Save result
                        relations.append(relationship_table.dropna())

            except:
                print(module + " module not formatted correctly, or no synonyms were mapped")
                continue

    else: 

        for module in modules:

            try: 

                # Construct the module url
                module_url = "https://rest.kegg.jp/get/" + module
                text = str(__parse_url(module_url))

                orthology = []
                orthology_flag = False
                reaction = []
                reaction_flag = False
                compound = []
                compound_flag = False

                # Pull relevant text---------------------------------------
                for line in text.split("\n"):

                    if line.startswith("ORTHOLOGY"):
                        orthology_flag = True
                    elif line.startswith("REACTION"):
                        reaction_flag = True
                    elif line.startswith("COMPOUND"):
                        compound_flag = True

                    # Pull orthology 
                    if orthology_flag:
                        if "ORTHOLOGY" in line or line.startswith(" "):
                            orthology.append(line.replace("ORTHOLOGY", " ").strip())
                        else:
                            orthology_flag = False

                    # Pull reactions
                    if reaction_flag:
                        if "REACTION" in line or line.startswith(" "):
                            reaction.append(line.replace("REACTION", " ").strip())
                        else:
                            reaction_flag = False

                    # Pull compounds
                    if compound_flag:
                        if "COMPOUND" in line or line.startswith(" "):
                            compound.append(line.replace("COMPOUND", " ").strip())
                        else:
                            compound_flag = False
                            
                # Pull protein-protein relationships----------------------------------

                # Extract proteins
                proteins = [protein.split("  ")[1].split(" [")[0] for protein in orthology]

                # Extract protein-protein relationships
                protein_IDs = []

                for el in range(len(proteins)-1):
                    protein_IDs.append(proteins[el] + " & " + proteins[el+1])

                protein_IDs = list(set(protein_IDs))

                # Pull metabolite-metabolite relationships-----------------------------

                # Extract metabolites
                metabolite_rellys = []

                # Extract common names
                KIDs = pd.DataFrame([x.split("  ") for x in compound]).rename({0:"KEGG", 1:"Name"}, axis = 1)

                # Extract metaboilte relationships
                for react in reaction:
                    metabs = react.split("  ")[1].replace("-&gt;", "+").split("+")
                    metabs = [KIDs[KIDs["KEGG"] == metab.strip()]["Name"].tolist()[0] for metab in metabs]
                    metabolite_rellys.append(__upper_triangle_meshgrid(metabs))

                # Clean duplicates
                metabolite_rellys = pd.concat(metabolite_rellys).reset_index(drop = True).replace("nan", np.nan).dropna().reset_index(drop = True)
                metabolite_IDs = []
                for x in range(len(metabolite_rellys)):
                    terms = metabolite_rellys.loc[x, :].tolist()
                    terms = [term.strip() for term in terms]
                    metabolite_IDs.append(" & ".join(terms))
                metabolite_IDs = list(set(metabolite_IDs))

                # Pull metabolite-protein relationships---------------------------------

                metab_prot_IDs = []

                # Pull associated reaction per protein
                for ortho in orthology:
                    protein = ortho.split("  ")[1].split(" [")[0]
                    rns = ortho.split("[RN:")[1].split(" ")
                    rns = [rn.replace("]", "").strip() for rn in rns]

                    # Pull metabolites in a reaction
                    for rn in rns:
                        selected = [x for x in reaction if rn in x][0]
                        metabs = selected.split("  ")[1].replace("-&gt;", "+").split("+")
                        metabs = [KIDs[KIDs["KEGG"] == metab.strip()]["Name"].tolist()[0] for metab in metabs]
                        metab_prot_IDs.extend([metab + " & " + protein for metab in metabs])

                # Map synonyms and return relationships---------------------------------

                # Extract all IDs
                pre_all_IDs = protein_IDs
                pre_all_IDs.extend(metabolite_IDs)
                all_IDs = []
                for id in pre_all_IDs:
                    all_IDs.extend(id.split(" & "))
                all_IDs = list(set(all_IDs))

                # MAP IDs to synonyms
                syns = map_synonyms(all_IDs, omes_folder, proteome_filename)

                # Construct table
                module_table = []

                # Pull all module relationships
                all_rels = protein_IDs
                all_rels.extend(metabolite_IDs)
                all_rels.extend(metab_prot_IDs)

                for rel in all_rels:

                    # Extract terms
                    term1 = re.sub(r'[^a-zA-Z0-9]', '', rel.split(" & ")[0].strip().lower())
                    term2 = re.sub(r'[^a-zA-Z0-9]', '', rel.split(" & ")[1].strip().lower())

                    # Make table
                    sub_table = pd.concat([syns[syns["Synonym"] == term1].rename({"Synonym":"Synonym1", "ID":"ID1", "Type":"Type1"}, axis = 1).reset_index(drop = True), 
                                           syns[syns["Synonym"] == term2].rename({"Synonym":"Synonym2", "ID":"ID2", "Type":"Type2"}, axis = 1).reset_index(drop = True)], axis = 1)
                    
                    module_table.append(sub_table)

                module_table = pd.concat(module_table).dropna().reset_index(drop = True)
                module_table["Source"] = "database"
                relations.append(module_table)

            except:
                print(module + " module not formatted correctly, or no synonyms were mapped")
                continue

    ## Remove Duplicates-------------------------------------------------------------------------------------

    # Remove duplicates
    relationships = pd.concat(relations)
    final_relationships = remove_relationship_duplicates(relationships, remove_self_relationships)

    if output_directory is not None:
        final_relationships.to_csv(os.path.join(output_directory, str(kegg_species_id) + "_kegg.txt"), index=False, sep = "\t")
        return None
    else:
        return final_relationships