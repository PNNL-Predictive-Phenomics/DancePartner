import requests
from xml.etree import ElementTree as ET
import time
import pandas as pd

def query_pubmed(query: str, pubmed_api_key: str = None, max_results: int = 10000):
    """
    Pulls PubMed IDs, Titles, and DOIs for articles matching a given PubMed query. This function uses the E-utilities API to first search for PubMed IDs matching the query, and then fetches the details of those articles in batches. The output is a DataFrame containing the PubMed IDs, Titles, and DOIs, which can be passed to downstream deduplication functions in the DancePartner package.
    
    Parameters
    -----------
    query
        the PubMed query to search for. See the README on the github repository for suggestions on constructing specific PubMed queries. 
    pubmed_api_key
        the PubMed API key to use (optional). Providing an API key allows for faster and more reliable querying, and also allows for a higher rate limit. 
    max_results
        the maximum number of results to return (default is 10000). More than this number is not recommended, as it may lead to long query times and potential timeouts. If you need more than this number of results, consider breaking up your query into smaller subqueries.

    Returns
    -------
    pd.DataFrame
        a DataFrame containing the PubMed IDs, Titles, and DOIs. This output table can be passed to downstream deduplication functions in the DancePartner package
    """

    ## STEP 1: EQuery PubMed IDs-------------------------------------

    # Start a list of parameters
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "usehistory": "y"
    }

    # Add the API key if provided, for faster and more reliable querying
    if pubmed_api_key is not None:
        params["api_key"] = pubmed_api_key

    # Make the request and return the list of PubMed IDs
    search_response = requests.get(url, params = params)
    search_response.raise_for_status() # Automatically raise an error if the request was unsuccessful
    search_data = search_response.json()["esearchresult"]

    # If no PMIDs returned, return a blank table
    pmids = search_data.get("idlist", [])
    if not pmids:
        return pd.DataFrame(columns=["PMID", "Title", "DOI"])

    # Now, extract the webenv and query_key from the search data, which will be used in the next step to fetch the details of the articles using EFetch
    webenv = search_data.get("webenv")
    query_key = search_data.get("querykey")

    ## Step 2: Efetch article details using the webenv and query_key-------------------------------------

    # Extrac the records in batches, while following guidance from PubMed
    records = []
    batch_size = 200
    sleep_time = 0.11 if pubmed_api_key else 0.34  # Respect 2026 rate limits

    # Iterate through the list of PMIDs in batches
    for start in range(0, len(pmids), batch_size):
        fetch_params = {
            "db": "pubmed",
            "retmode": "xml",
            "retstart": start,
            "retmax": batch_size,
            "WebEnv": webenv,
            "query_key": query_key,
        }
        if pubmed_api_key is not None:
            fetch_params["api_key"] = pubmed_api_key

        fetch_response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params = fetch_params)
        fetch_response.raise_for_status() # Once again, automatically raise an error if the request was unsuccessful

        root = ET.fromstring(fetch_response.content)
        for article in root.findall(".//PubmedArticle"):

            # Extract the PMID
            pmid_el = article.find(".//MedlineCitation/PMID")
            pmid = pmid_el.text if pmid_el is not None else ""

            # Extract the Title (use itertext to capture italics, etc.)
            title_el = article.find(".//ArticleTitle")
            title = "".join(title_el.itertext()) if title_el is not None else ""

            # Extract the DOI
            doi = ""
            for id_el in article.findall(".//ArticleIdList/ArticleId"):
                if id_el.attrib.get("IdType") == "doi":
                    doi = id_el.text or ""
                    break

            records.append({"PMID": pmid, "Title": title, "DOI": doi})

        time.sleep(sleep_time)

    papers = pd.DataFrame(records)
    if len(papers) > max_results:
        papers = papers.iloc[:max_results]
    return papers

def query_scopus(query: str, scopus_api_key: str, max_results: int = 5000):
    """
    Pulls Scopus records matching a query. Returns a DataFrame of Titles, DOIs, and EIDs.

    Parameters
    -----------
    query
        the Scopus query to search for. Scopus uses its own query syntax
        (e.g. 'TITLE-ABS-KEY("e coli proteomics") AND PUBYEAR > 1999').
        See the README on the github repository for suggestions on constructing
        specific Scopus queries.
    scopus_api_key
        a string API key for Scopus-Elsevier. See the Elsevier Developer Portal
        for assistance acquiring one. Required (unlike PubMed, Scopus has no
        unauthenticated access).
    max_results
        the maximum number of results to return (default is 5000). Note that
        Scopus limits standard search results to 5000 per query — for larger
        result sets, break the query up by year, topic, etc. Per Elsevier's
        terms, weekly quotas also apply.

    Returns
    -------
    pd.DataFrame
        a DataFrame containing the Titles, DOIs, and EIDs (Scopus's unique
        identifier). This output table can be passed to downstream
        deduplication functions in the DancePartner package.
    """

    if scopus_api_key is None:
        raise ValueError("scopus_api_key must be provided. See Elsevier Dev Portal for assistance")

    url = "https://api.elsevier.com/content/search/scopus"
    headers = {
        "X-ELS-APIKey": scopus_api_key,
        "Accept": "application/json",
    }

    # Scopus returns results in pages — the API max is 25 per page on the
    # standard view (200 on the COMPLETE view, but that requires additional
    # institutional access). We'll use 25 as a safe default.
    records = []
    page_size = 25
    start = 0
    sleep_time = 0.15  # Respect Scopus rate limits (~9 req/sec is the typical cap)

    while start < max_results:
        params = {
            "query": query,
            "count": min(page_size, max_results - start),
            "start": start,
        }
        response = requests.get(url, headers=headers, params=params)
        response.raise_for_status()  # Automatically raise an error if the request was unsuccessful
        data = response.json().get("search-results", {})

        # Total number of hits available in Scopus for this query
        total_hits = int(data.get("opensearch:totalResults", 0))
        if total_hits == 0:
            return pd.DataFrame(columns=["Title", "DOI", "EID"])

        # Iterate through entries on this page
        entries = data.get("entry", [])
        if not entries:
            break  # No more results

        for entry in entries:
            # Scopus returns "error" entries when no result is found on a page
            if "error" in entry:
                continue
            records.append({
                "Title": entry.get("dc:title", ""),
                "DOI": entry.get("prism:doi", ""),
                "EID": entry.get("eid", ""),
            })

        # Stop if we've pulled all available results
        start += page_size
        if start >= total_hits:
            break
        time.sleep(sleep_time)

    return pd.DataFrame(records)


def query_osti(query: str, max_results: int = 10000):
    """
    Pulls OSTI records matching a query. Returns a DataFrame of TITLE, DOI, and OSTI_IDENTIFIER for each record.

    Parameters
    -----------
    query
        the OSTI query to search for. OSTI supports full-text search across titles,
        abstracts, and other fields. See the README on the github repository for
        suggestions on constructing specific OSTI queries.
    max_results
        the maximum number of results to return (default is 10000). More than this
        number is not recommended, as it may lead to long query times and potential
        timeouts. If you need more than this number of results, consider breaking
        up your query into smaller subqueries.

    Returns
    -------
    pd.DataFrame
        a DataFrame containing the OSTI_IDENTIFIERs, Titles, and DOIs. This output table can
        be passed to downstream deduplication functions in the DancePartner package.
    """
    url = "https://www.osti.gov/api/v1/records"

    # OSTI returns results in pages of up to 100 records, controlled by the
    # `rows` and `page` parameters. We paginate until we hit max_results or
    # run out of records.
    records = []
    page_size = 100
    page = 0
    sleep_time = 0.2  # Be polite — OSTI doesn't publish a rate limit, but throttling avoids issues

    while len(records) < max_results:
        params = {
            "q": query,
            "rows": min(page_size, max_results - len(records)),
            "page": page,
        }
        response = requests.get(url, params=params)
        response.raise_for_status()  # Automatically raise an error if the request was unsuccessful
        data = response.json()

        # OSTI returns a list of records directly (not wrapped in an envelope)
        if not isinstance(data, list) or not data:
            break  # No more results

        for entry in data:
            records.append({
                "TITLE": entry.get("title", ""),
                 "DOI": entry.get("doi", ""),
                "OSTI_IDENTIFIER": str(entry.get("osti_id", ""))
            })

        # If OSTI returned fewer than a full page, we've hit the end
        if len(data) < page_size:
            break

        page += 1
        time.sleep(sleep_time)

    # If we somehow over-shot max_results, trim
    records = records[:max_results]
    if not records:
        return pd.DataFrame(columns=["TITLE", "DOI", "OSTI_IDENTIFIER"])
    return pd.DataFrame(records)