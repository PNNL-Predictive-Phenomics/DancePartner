import os
import requests
import re
import pandas as pd
from bs4 import BeautifulSoup
from requests.auth import HTTPBasicAuth
import zipfile
import io
import json

def pull_proteome(proteome_id: str, output_directory: str):
    """
    Function that pulls a proteome and its synonyms for a species. 

    Parameters
    ----------
    proteome_id
        Search for a proteome ID here: https://www.uniprot.org/proteomes/. It starts with "UP"
        
    output_directory
        Path specifying where to write the result within the current directory. 
    
    Returns
    -------
    Protein IDs and their synonyms in a text file
    """

    # Construct the URL
    try:
        url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=%28%28proteome%3A" + str(proteome_id) + "%29%29"

        # Read the filev
        req = requests.get(url, stream = True)
        soup = BeautifulSoup(req.content, 'html.parser')
        myjson = json.loads(str(soup))

        # Iterate through file 
        ID = []
        Synonyms = []

        for entry in range(len(myjson["results"])):

            # Pull the subdata
            subdata = myjson["results"][entry]

            # Extract accession
            try:
                ID.append(subdata["primaryAccession"])
            except:
                ID.append("")

            # Extract the synonyms
            try:

                # Pull values
                split = str(subdata["proteinDescription"]).split(",")

                # Clean splits
                clean_splits = []

                # Extract out synonym
                for x in split:
                    if "'value':" in x:
                        x = x.replace("'value':", "")
                        x = re.sub(r'[^a-zA-Z0-9 .]', '', x.replace("'value':", ""))
                        x = x.replace('recommendedName', '')
                        x = x.replace('fullName', '')
                        x = x.replace('alternativeNames', '')
                        x = x.strip()
                        if len(x) >= 3:
                            clean_splits.append(x)

                # Add to Synonym table
                Synonyms.append("; ".join(clean_splits))

            except:
                Synonyms.append("")

        # Extract out simplfied gene name
        for el in range(len(Synonyms)):
            Synonyms[el] = "; ".join([Synonyms[el], Synonyms[el].split("; ")[0].split(" ")[-1]])

        # Build proteome
        proteome = pd.DataFrame([ID, Synonyms]).T.rename({0:"UniProtID", 1:"Synonyms"}, axis = 1)

        if output_directory is not None:
            proteome.to_csv(os.path.join(output_directory, proteome_id + "_proteome.txt"), index=False, sep = "\t")
        else:
            return(proteome)
    except:
        print(proteome_id + " is not recognized as a proper proteome_id")

def pull_genome_from_GTF(GTF_file: str, output_directory: str):
    '''
    Extract the Gene IDs and Gene Names from a GTF file

    Parameters
    ----------
    GFT_file
        Path to the GTF genome file downloaded from ENSEMBL: https://www.ensembl.org/index.html
    
    output_directory
        Path specifying where to write the result within the current directory. 
    
    Returns
    -------
    Gene IDs and their synonyms in a text file
    '''

    # Hold all the gene names
    gene_id = []
    gene_names = []

    # Extract all the gene names 

    # Open the file
    with open(GTF_file, "r") as file:
        
        # Iterate through the file, searching for all gene names
        for line in file:

            # If the item is a gene name, extract the value
            if "gene_id" in line and "gene_name" in line:

                # Add identifier
                gene_id.append(line.split("gene_id")[1].strip().split('"')[1])

                # Find all genes
                genes = line.split("gene_name")[1].strip().split('"')[1].split("{}")

                # Remove numerics
                gene_list = []
                for gene in genes:
                    try:
                        float(gene)
                    except:
                        gene_list.append(gene.split(":")[-1])

                # Add all genes that are not numberics       
                gene_names.append("; ".join(gene_list))

    # Start dataframe
    genome = pd.DataFrame({
        "Gene": gene_id,
        "Synonyms": gene_names
    }).drop_duplicates().reset_index(drop = True)

    # Drop any cases where there is not a synonym or ID
    genome = genome[(genome["Gene"] != "") & (genome["Synonyms"] != "")].reset_index(drop = True)

    # Write genome
    if output_directory is not None:
        genome.to_csv(os.path.join(output_directory, "genome.txt"), index=False, sep = "\t")
    else:
        return(genome)