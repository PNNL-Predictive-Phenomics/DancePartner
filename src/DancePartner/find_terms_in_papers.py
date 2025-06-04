import os
import pandas as pd
import nltk
import re
from pathlib import Path

def find_terms_in_papers(paper_directory: str, terms: list[str], output_directory: str = None, 
                         n_gram_max: int = 3, max_char_length: int = 250, padding: int = 10, 
                         verbose: bool = False):   
    """
    This function searches through sentences of papers to extract biomolecule pairs present in each sentence.
    It utilizes a set-intersection method on the n-grams of the sentences with the already-found biomolecule synonyms.

    Parameters
    ----------
    paper_directory
        A directory path pointing to the list of papers to be parsed through.
    
    terms
        List of terms to find in papers
    
    output_directory
        An optional path to a directory for where to write results to. Otherwise, the function will return the table.

    n_gram_max
        The number of n_grams to consider when combing the papers. (e.g. n_grams=2 will catch "protein A" but n_grams=1 will not). If unsure, use the default.
    
    max_char_length 
        The number of maximum characters that can be in a segment containing the pair of biomolecules
        
    padding
        The amount of padding (in characters) to surround the terms in a segment by at minimum.
    
    verbose
        If True, print status messages
    
    Returns
    -------
        A Pandas DataFrame of the resulting data.
    """
    
    # Modify term matches
    terms = [re.sub("[^\s\d\w]|\n", '', term.lower()) for term in terms]

    matches = []
    for root, _, files in os.walk(paper_directory):
        for file in files:
            if file.endswith(".txt") is False:
                continue
            if verbose:
                print("On file " + file)
            file_path = os.path.join(root, file)
            file_id = Path(file_path).stem
            with open(file_path, "r") as f:
                sentences = [re.sub("[^\s\d\w]|\n", '', x.lower()) for x in nltk.sent_tokenize(f.read())]
                for sentence_ind, sentence in enumerate(sentences):
                    words = nltk.word_tokenize(sentence)
                    joined_word_ngrams = [' '.join(x) for x in list(nltk.everygrams(words, 1, n_gram_max))]
                    found_terms = list(set(joined_word_ngrams).intersection(terms))
                    if len(found_terms) > 1:
                        pairs = [(a,b) if a < b else (b,a) for idx, a in enumerate(found_terms) for b in found_terms[idx + 1:] if a != b]
                        for pair in set(pairs): 
                            try:
                                index1 = re.search(rf"\b{pair[0]}\b", sentence).span()[0]
                                index2 = re.search(rf"\b{pair[1]}\b", sentence).span()[0]
                            except AttributeError:
                                continue
                            # Assign term based off of which occurs first in sentence
                            if index1 < index2:
                                term1 = pair[0]
                                term2 = pair[1]
                            else:
                                term1 = pair[1]
                                term2 = pair[0]
                            # Create segment of sentence with terms if len(sentence) is too long
                            if len(sentence) < max_char_length:
                                segment = sentence
                            else:
                                first_index, second_index = min(index1, index2), max(index1, index2)
                                segment = sentence[max(0, first_index-padding):
                                                    min(len(sentence)-1, second_index + max(len(term1), len(term2)) +padding)]
                                if len(segment) > max_char_length:
                                    continue

                            matches.append([file_id, term1, term2, '_'.join([file_id]), sentence_ind, len(segment), segment])
                
    # Wrap up matches and write to a CSV file       
    column_names = ['paper_id','term_1','term_2','id','sentence_index', 'segment_length','segment']
    matches_df = pd.DataFrame(matches, columns = column_names).sort_values('paper_id')

    # Clean up any nonalphanumerics
    matches_df["segment"] = matches_df["segment"].str.replace(r'[^a-zA-Z0-9 ]', '', regex=True)

    if output_directory is not None:
        matches_df.to_csv(os.path.join(output_directory, "sentence_biomolecule_pairs.csv"), index=False)
    else:
        return(matches_df)