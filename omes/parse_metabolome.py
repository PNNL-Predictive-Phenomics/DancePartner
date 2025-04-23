import pandas as pd
import os

def parse_metabolome(sdf_file, output_directory):
    '''
    Function to parse a metabolome from CHEBI. Takes about 3 seconds to run.

    Args:
        sdf_file (String): Path to unzipped SDF file. Get it here https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/. Pull the complete one.
        output_directory (String): List the folder you would like to write the "CHEBI_Metabolome.txt" to.
    Returns:
        A file with RefMet_ID, Synonyms. 
    '''

    # Hold all mined inputs
    CHEBI_ID = []
    Synonyms = []

    # Store important flags
    CHEBI_Flag = False
    First_Pass = False
    Synonym_Flag = False

    # Parse out IDs and Synonyms
    with open(sdf_file, "r") as file:
        for line in file:        

            # If the CHEBI flag is true
            if CHEBI_Flag:

                # Save flag on first pass
                if First_Pass:
                    CHEBI_ID.append(line.strip("\n"))
                    First_Pass = False

                # Not every entry has a synonym. If it doesn't, then return a blank
                if "> <ChEBI ID>" in line:
                    Synonyms.append("")
                
                # If there is a list of synonyms, let's capture them  
                elif "> <Synonyms>" in line:
                    Synonym_Flag = True
                    Synonym_List = []
            
                # We will continue to capture synonyms until there isn't anymore. Then we will append. 
                elif Synonym_Flag:

                    if line != "\n":
                        Synonym_List.append(line.strip("\n"))
                    else:
                        Synonyms.append("; ".join(Synonym_List))
                        Synonym_Flag = False
                        CHEBI_Flag = False
        
            # Set CHEBI flag to true
            if "> <ChEBI ID>" in line:
                CHEBI_Flag = True
                First_Pass = True
    
    # Write txt file
    pd.DataFrame({
        "CHEBI": CHEBI_ID,
        "Synonyms": Synonyms
    }).to_csv(os.path.join(output_directory, "CHEBI_Metabolome.txt"), sep = "\t", index = False)
    return None

