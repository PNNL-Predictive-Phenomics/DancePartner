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
