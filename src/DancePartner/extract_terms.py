import os
import pandas as pd
import re

# Build a support function to get directories 
def __get_all_files(directory):
    '''
    Walk through a direcotry and pull all txt files. Remove the summary file and gzipped files. 
    '''

    # Make a list to hold all filepaths 
    file_paths = []

    # Iterate through directories 
    for dirpath, dirnames, filenames in os.walk(directory):

        # Iterate throughfiles
        for filename in filenames:

            # Get the file path
            file_path = os.path.join(dirpath, filename)

            # Append a file path only if it is not the output_summary file or a gz file
            if ("output_summary.txt" in file_path) is False and (".gz" in file_path) is False:
                file_paths.append(file_path)


    return file_paths

# Extract unique terms from papers 
def extract_terms_scispacy(paper_directory, 
                           omes_folder, 
                           tags = ["GENE_OR_GENE_PRODUCT", "SIMPLE_CHEMICAL", "AMINO_ACID"],
                           additional_stop_words = None,
                           min_length = 3, 
                           max_length = 100,
                           verbose = False):
    '''
    Extract terms from papers

    Args:
        paper_directory (String): Directory to papers in txt format. Subdirectories are searched, and .gz and 
            output_summary.txt files are ignored.
        omes_folder (String): Path to the omes folder where "stop_words_english.txt" is stored. Required. 
        tags (List[String]): A list of tags from the `en_ner_bionlp13cg_md` model. See https://allenai.github.io/scispacy/ 
        additional_stop_words (List[String]): Add more words to be removed from consideration. Default is None. 
        min_length (Integer): The minimum number of non-whitespace characters required. Default is 3.
        max_length (Integer): The maximum number of characters allowed. Default is 100. 
        verbose (Logical): Indicate whether a message should be printed as each file is processed. Default is "FALSE"
    Returns:
        A list of unique terms found in papers 
    '''

    # Pull papers
    paper_list = __get_all_files(paper_directory)

    try:
        import spacy
        import scispacy
        import en_ner_bionlp13cg_md
    except:
        raise Exception("To use this function, install spacy and scispacy. You must also pull the en_ner_bionlp13cg_md model. See the README for more details.")


    # Load the model
    nlp_b13 = en_ner_bionlp13cg_md.load()

    # Load stop words
    stop_words = pd.read_csv(os.path.join(omes_folder, "stop_words_english.txt"))["stopwords"].tolist()

    if additional_stop_words is not None:
        stop_words.extend(additional_stop_words)

    # Hold a list of terms 
    identified_terms = []

    for file in paper_list:

        if verbose:
            print("On paper...", file)
        
        with open(file, "r") as info:

            # Read entire paper
            content = info.read()

            # Search for terms throughout paper
            b13_apply = nlp_b13(text = content)
            
            # Extract terms and labels
            b13_terms = []
            b13_labels = []

            # Unwrap bionlp13cg labels 
            for ent in b13_apply.ents:
                b13_terms.append(ent.text)
                b13_labels.append(ent.label_)

            # Create data.frame
            res = pd.DataFrame({"term":b13_terms, "label":b13_labels})
            terms = res[res["label"].isin(tags)]["term"].unique().tolist()

            # Clean terms 
            for term in terms:
                term = re.sub("[^\s\d\w]|\n", '', term.lower().strip())
                if len(term) >= min_length and len(term) <= max_length and term not in stop_words:
                    identified_terms.append(term)

    # Return unique list
    return list(set(identified_terms))