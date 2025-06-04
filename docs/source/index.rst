########################
Welcome to DancePartner!
########################

*DancePartner* is a python package for mining multi-omics relationship networks from literature and databases.
Though *DancePartner* may be organized into an intuitive pipeline, it should be thought of as a toolbox of 
functions for building multi-omics networks for various needs, whether those be networks derived only from literature,
databases, or a mix of the two. We welcome additions to our package and are happy to collaborate with anyone willing 
to add code to this framework. For simplicity, we will present the code in three different pipelines: 

`A. Mining from Literature`_

`B. Mining from Databases`_

`C. Building & Combining Networks`_

#########################
A. Mining from Literature
#########################

The literature mining pipeline may be summarized in ? key steps: `1. Pulling Publications`_, `2. Identifying Entities`_,
`3. Extracting Relationships`_, `4. Collapsing Synonyms`_, and `5. Building Networks`_.

***********************
1. Pulling Publications
***********************

Publications may be pulled from any of three databases: `PubMed <https://pubmed.ncbi.nlm.nih.gov/>`_, `Scopus <https://www.scopus.com/search/form.uri?display=basic#basic>`_,
and `OSTI <https://www.osti.gov/>`_. There is also a function to deduplicate papers across databases (which returns a table) that may be passed to the paper pulling function.

Documentation
=============

.. autoclass:: DancePartner.pull_papers.pull_papers

PubMed Example
==============

PubMed requires a list of PubMed IDs called PMIDs. To obtain PMIDs, simply enter a query into the search bar of PubMed, click “Save”, select “All results”, 
and output the format as “PMID”. 

.. code-block:: python

    pull_papers(database = "pubmed", ids = [9851916], output_directory = "papers")

Scopus Example
==============

Scopus uses DOIs to identify papers. To obtain these DOIs, enter a query into the search bar, click “Export”, select the desired format, select all documents, 
and then export at least the DOI column. Scopus also requires string API key. See https://dev.elsevier.com/.

.. code-block:: python

    # Save the scopus key as a variable
    pull_papers(
        scopus_ids = ["10.1186/s40168-021-01035-8", "10.1002/bit.26296", "10.1002/pmic.200300397", "10.1074/mcp.M115.057117"],
        output_directory = "scopus_papers", 
        scopus_api_key = scopus_api_key
    )

OSTI Example
============

To obtain OSTI IDs, enter a query and click “Save Results”, and the resulting file will contain the OSTI IDs.

.. code-block:: python

    pull_papers(
        osti_ids = ["2229172", "1629838", "1766618", "1379914"], 
        output_directory = "osti_papers"
    )


Deduplication
=============

To deduplicate papers across databases, the following steps must be followed:

**PubMed:** Enter the query, hit search, hit save, select "All results" and "csv"

**Scopus:** Enter the query, hit search, hit export, select "CSV" and keep all defaults checked.

**OSTI:** Enter the query, hit search, and save results as a "CSV"

.. autoclass:: DancePartner.deduplicate_papers.deduplicate_papers

.. code-block:: python

    # Save the results as a deduplicated table
    deduped_table = deduplicate_papers(pubmed_path, scopus_path, osti_path)

    # Then pull publications using the deduped table. Use the saved scopus_api_key
    pull_papers(deduped_table = deduped_table, output_directory = "deduped_example", scopus_api_key = scopus_api_key)

***********************
2. Identifying Entities
***********************

First, biological terms (entities) need to be defined and found in articles. An option is to use ScispaCy, or 
simply use synonym files, which is the recommended approach. To install ScispaCy, see the "Optional ScispaCy Model"
section of the `README <https://github.com/pnnl-predictive-phenomics/DancePartner/tree/main>`_.

ScispaCy
========

.. autoclass:: DancePartner.extract_terms.extract_terms_scispacy

.. code-block:: python

    # Extracting terms requires a path to the papers and a path to the omes folder
    extract_terms_scispacy(paper_directory = paper_directory, omes_folder = "../omes")

Synonym Files
=============

First, pull a proteome.

.. autoclass:: DancePartner.pull_ome.pull_proteome

.. code-block:: python

    # Pull a proteome (protein and its synonyms) and place it in the omes folder
    pull_proteome(proteome_id = "UP000001940", output_directory = "../omes")

List Synonyms
=============

Then, list all synonyms.

.. autoclass:: DancePartner.create_synonym_table.list_synonyms

.. code-block:: python

    # List the omes folder and the proteome to use in the omes folder
    list_synonyms("../omes", "UP000001940_proteome.txt")

***************************
3. Extracting Relationships
***************************

First, sentences with terms need to be extracted and formatted for the downstream BERT model.

.. autoclass:: DancePartner.find_terms_in_papers.find_terms_in_papers

.. code-block:: python

    # Supply this function with the directory with the papers, a list of terms, and the output directory
    find_terms_in_papers(
        paper_directory = paper_directory,
        terms = my_terms,
        n_gram_max = 5,
        padding = 50,
        output_directory = output_directory
    )

Next, BERT can be run. Extract the BERT model from `here <https://huggingface.co/david-degnan/BioBERT-RE/tree/main>`_. 
Place in the top level directory of this repo in a folder called "biobert". Pull the config.json, the pytorch_model.bin, and the training_args.bin files.

.. autoclass:: DancePartner.bert_functions.run_bert

.. code-block:: python

    # Create a variable to the output directory. Make sure the biobert model is pulled and in a location that can be referenced by model_path.
    run_bert(
        input_path = "/path/to/sentence_biomolecule_pairs.csv",
        model_path = "../biobert", # Update to the path of your BERT model if necessary. We recommend placing biobert in the top directory.
        output_directory = output_directory,
        segment_col_name = "segment",
        use_cpu = True
    )

**********************
4. Collapsing Synonyms
**********************

