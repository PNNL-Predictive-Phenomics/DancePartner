import os
import requests
import re
import pandas as pd
from bs4 import BeautifulSoup
from requests.auth import HTTPBasicAuth
import zipfile
import io
import json

def pull_proteome(proteome_id, output_directory):
    """
    Function that pulls a proteome and its synonyms for a species. 

    Args:
        proteome_id (String): Search for a proteome ID here: https://www.uniprot.org/proteomes/. It starts with "UP"
        output_directory (String): Path specifying where to write the result within the current directory. 
    Returns:
        Protein IDs and their synonyms in a text file
    """

    # Construct the URL
    try:
        url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=%28%28proteome%3A" + str(proteome_id) + "%29%29"

        # Read the file
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
    
    
def pull_genome(species_id, ncbi_api_key, output_directory):
    """
    Function that pulls a genome for a species.

    Args:
        species_id (String or Numeric): The taxon ID for the organism of interest.
        ncbi_api_key (String): A String API key to the NCBI database. See: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ 
        output_directory (String): Path specifying where to write the result.
    Returns:
        Text file with all genes from genome
    Example Code:
        ncbi_key = open("example_data/ncbi_key.txt")
        pull_genome(2097, "test_omes", ncbi_key.read())
    """
    
    if ncbi_api_key is None:
        raise ValueError("Must specify `ncbi_api_key` for proper use of function. See documentation.")
    auth = HTTPBasicAuth('api-key', ncbi_api_key)

    # Use species ID to find accession IDs
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/" + str(species_id) + "/dataset_report?filters.assembly_level=chromosome&filters.assembly_level=complete_genome&table_fields=assminfo-accession&table_fields=assminfo-name"
    req = requests.get(url, headers={'Accept':'application/json'}, auth=auth)
    response = req.json()
    if len(response['reports']) > 0:
        accession_id = response['reports'][0]['accession']
    else:
        print("No accessions found for species ID.")
        return(None)
    
    # Use accession id to download fasta file (cds)
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/" + str(accession_id) + \
        "/download?chromosomes=1&chromosomes=2&chromosomes=3&chromosomes=X&chromosomes=Y&chromosomes=MT&include_annotation_type=GENOME_GFF&include_annotation_type=GENOME_GBFF&include_annotation_type=GENOME_GTF&include_annotation_type=CDS_FASTA"
    data = requests.get(url, headers={'Accept':'application/zip'}, auth=auth)
    z = zipfile.ZipFile(io.BytesIO(data.content))

    # Look for correct file in zip folder
    try:
        fasta_file = next(x for x in z.namelist() if ".fna" in x)
    except StopIteration:
        print("No Fasta file found in request.")
        return(None)
    
    # Parse fasta file  (keep lines that aren't nucleotide sequences, i.e line begins with '>') amd write results to .txt file
    with open(os.path.join(output_directory, str(species_id) + "_ncbi_genes.txt"), "w") as f:
        for line in z.open(fasta_file, "r").readlines():
            line_string = line.decode() #downloaded zip file reads in bytes --> convert to strings
            if line_string[0] == ">":

                # Finally extract the gene name if possible. If not, proceed.
                try:
                    f.write(re.search("\[gene=(.*?)\]", line_string).group(1))
                except AttributeError:
                    continue

    return(None)
    