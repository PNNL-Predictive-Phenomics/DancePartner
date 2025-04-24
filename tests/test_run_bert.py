import shutil
import os
import pandas as pd
import DancePartner as dance

## How to calculate coverage (from within main package directory): 
# coverage run --source=DancePartner -m pytest -x tests/* -W ignore
# coverage report
# coverage html

# This function tests both find names in papers, running BERT, building networks, and constructing metrics
def test_run_bert():

    # Find terms to test in BERT in a paper
    output_directory = os.path.join(os.getcwd(), "testfind")
    os.mkdir(output_directory)
    terms = dance.list_synonyms("omes", "UP000001940_proteome.txt")

    # Format for BERT and run BERT
    dance.find_terms_in_papers("example_data/papers", terms, output_directory = output_directory)
    dance.run_bert(os.path.join(output_directory, "sentence_biomolecule_pairs.csv"), model_path = "biobert", output_directory = output_directory, segment_col_name = "segment")

    # Subset BERT
    BERT = pd.read_table(os.path.join(output_directory, "bert_results.txt"))
    TrueBERT = BERT[BERT["True Positive"] >= 0.5]

    # Pull terms 
    all_found_terms = TrueBERT["term_1"].tolist()
    all_found_terms.extend(TrueBERT["term_2"].to_list())
    all_found_terms = list(set(all_found_terms))

    # Map synonyms
    dance.map_synonyms(
        term_list = all_found_terms,
        omes_folder = "omes",
        proteome_filename = "UP000001940_proteome.txt",
        add_missing = True,
        output_directory = output_directory
    )
    synonyms = pd.read_csv(os.path.join(output_directory, "synonym_table.txt"), sep = "\t")

    # Build network
    network_table = dance.build_network_table(TrueBERT, synonyms)
    network = dance.visualize_network(network_table)
    metrics = dance.calculate_network_metrics(network)

    shutil.rmtree(output_directory)