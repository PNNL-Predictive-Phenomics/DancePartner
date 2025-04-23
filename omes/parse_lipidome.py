import os
import pandas as pd

def parse_lipidome(sdf_file, csv_file, output_directory):
    '''
    Function to parse a lipidome from output files in lipidmaps

    Args:
        sdf_file (String): Path to the SDF file. Download and unzip the SDF file of the entire database from here: https://www.lipidmaps.org/databases/lmsd/download
        csv_file (String): Path to the CSV file. Go to https://www.lipidmaps.org/data/structure/LMSDSearch.php and enter a blank search. Download the entire database as a csv.
        output_directory (String): List the folder you would like to write the "LipidMaps_Lipidome.csv" file to
    Returns:
        A file with the LMID, Synonyms, Category, Class, and Abbreviation from LipidMaps
    '''

    ##############
    ## SDF FILE ##
    ##############

    # Let's build a dictionary with a list of IDs and their synonyms, as well as a separate lists to track class names
    Lipids = {}
    Categories = []
    Classes = []

    # Add trackers for each item 
    id_tracker = False
    common_tracker = False
    systematic_tracker = False
    synonym_tracker = False
    category_tracker = False
    class_tracker = False

    # Read file 
    sdf = open(sdf_file, "r")

    # Parsing through the file line by line
    for line in sdf:

        # If any trackers are true, add the information to the dictionary
        if id_tracker:
            theID = line.strip("\n")
            Lipids[theID] = ""
            id_tracker = False
        
        # These entries are all synonyms and appended with semicolons
        if common_tracker or systematic_tracker or synonym_tracker:
            theName = line.strip("\n").strip('"') + "; "
            Lipids[theID] = Lipids[theID] + theName
            common_tracker = False
            systematic_tracker = False
            synonym_tracker = False

        # These entries all list lipid categories, classes, and subclasses
        if category_tracker:
            Categories.append(line.strip("\n").strip('"'))
            category_tracker = False
        if class_tracker:
            Classes.append(line.strip("\n").strip('"'))
            class_tracker = False

        # Set trackers to true when the headers are detected
        if "> <LM_ID>" in line:
            id_tracker = True
        if "> <COMMON_NAME>" in line:
            common_tracker = True
        if "> <SYSTEMATIC_NAME>" in line: 
            systematic_tracker = True
        if "> <SYNONYMS>" in line:
            synonym_tracker = True
        if "> <CATEGORY>" in line:
            category_tracker = True
        if "> <MAIN_CLASS>" in line:
            class_tracker = True

    # Make a data.frame 
    LipExt = pd.DataFrame.from_dict(Lipids, orient = "index").rename(columns = {0:"Synonyms"}).rename_axis("LMID").reset_index()
    LipExt["Category"] = Categories
    LipExt["Class"] = Classes

    ##############
    ## CSV FILE ##
    ##############

    # Read in CSV data
    LipAcr = pd.read_csv(csv_file)

    # Extract the two needed columns
    LipAcr = LipAcr[["LM_ID", "ABBREV"]].rename(columns = {"LM_ID":"LMID", "ABBREV":"Abbreviation"})

    # Filter out blank abbreviations
    LipAcr = LipAcr[LipAcr["Abbreviation"] != "-"].reset_index(drop = True)

    # Add abbreviation
    LipExt = LipExt.merge(LipAcr, how = "left")

    ##################
    ## WRITE OUTPUT ##
    ##################

    outpath = os.path.join(output_directory,  "LipidMaps_Lipidome.csv") 
    LipExt.to_csv(outpath, index = False)
    return(None)
