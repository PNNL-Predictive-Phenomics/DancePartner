import os
import shutil
import pytest
import DancePartner as dance

## How to calculate coverage (from within main package directory): 
# coverage run --source=DancePartner -m pytest -x tests/* -W ignore
# coverage report
# coverage html

# This test is for pulling papers when only a list of PMIDs is given
def test_pull_papers_pmids():

    # Create output directories 
    output_directory = os.path.join(os.getcwd(), "pubmed_test")
    os.mkdir(output_directory, mode = 0o777)

    # Pull three papers (one abstract, one pdf, on clean text) and one non-paper
    paper_ids = [9851916, 9531567, 16803962, 16803962333333]

    # Pull papers
    dance.pull_papers(pubmed_ids = paper_ids, output_directory = output_directory)

    # Remove directory
    shutil.rmtree(output_directory)

# Test scopus 
def test_pull_scopus():

    # Create output directories 
    output_directory_one = os.path.join(os.getcwd(), "scopus_test_1")
    os.mkdir(output_directory_one)

    # Read key
    api_key_file = open("example_data/scopus_key.txt")
    api_key = api_key_file.read()
    api_key_file.close()

    # Pull a clean scopus file
    dance.pull_papers(output_directory = output_directory_one, scopus_ids = ["10.1016/j.bbamem.2020.183488"], scopus_api_key = api_key)

    # Remove files
    shutil.rmtree(output_directory_one)

# Test OSTI
def test_pull_osti():

    # Create output diectories
    output_directory_uno = os.path.join(os.getcwd(), "osti_test_1")
    os.mkdir(output_directory_uno)

    # Pull a clean scopus file
    dance.pull_papers(output_directory = output_directory_uno, osti_ids = [1629838])

    # Remove files
    shutil.rmtree(output_directory_uno)

# This test is for pulling papers when only three small files are given
def test_pull_papers_deduplicate():

    # Create an output directory 
    output_directory = os.path.join(os.getcwd(), "deduplication_test")
    os.mkdir(output_directory)

    # Ensure a scopus key has been made
    if os.path.isfile("example_data/scopus_key.txt") == False:
        raise Exception("Request a scopus key at https://dev.elsevier.com/ and save the key in a file called 'example_data/scopus_key.txt' to run this function.")

    # Read key
    api_key_file = open("example_data/scopus_key.txt")
    api_key = api_key_file.read()
    api_key_file.close()

    # Pull two papers and one non-paper
    pubmed_path = "example_data/pubmed_ecoli_short.txt"
    scopus_path = "example_data/scopus_ecoli_short.txt"
    osti_path = "example_data/osti_ecoli_short.txt"

    # Deduplicate papers
    deduped = dance.deduplicate_papers(pubmed_path = pubmed_path, scopus_path = scopus_path, osti_path = osti_path)

    # Download papers
    dance.pull_papers(output_directory, deduped_table = deduped, scopus_api_key=api_key)
    shutil.rmtree(output_directory)