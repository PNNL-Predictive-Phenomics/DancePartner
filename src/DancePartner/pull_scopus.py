import os
import requests
from bs4 import BeautifulSoup
import re

def __pull_scopus_clean(ids, output_directory, scopus_api_key):
    """
    Function that pulls paper abstracts from Scopus. Writes them to a directory.
    Args:
        ids (List): A list of DOIs - scopus pulls based off the DOI. 
        output_directory (String): Path specifying where to write the papers to.
        scopus_api_key (String): A string API key for Scopus-Elselvier. See documentation on how to acquire.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """

    if scopus_api_key is None:
        raise ValueError("scopus_api_key must be provided. See Elsevier Dev Portal for assistance")
    found_ids = []
    write_path = os.path.join(output_directory, "scopus_clean")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path)
    for paper_id in ids:
        req = requests.get("https://api.elsevier.com/content/article/doi/" + str(paper_id) + "?apiKey=" + scopus_api_key)
        if req.status_code == 200:
            soup = BeautifulSoup(req.content, features="lxml")
            try:
                # Remove non-paper text
                for x in soup.find_all("ce:bibliography"):
                    x.decompose()
                for x in soup.find_all("xocs:references"):
                    x.decompose()
                for x in soup.find_all("ce:author"):
                    x.decompose()
                for x in soup.find_all("ce:affiliation"):
                    x.decompose()
                for x in soup.find_all("object"):
                    x.decompose()
                for x in soup.find_all("xocs:attachments"):
                    x.decompose()
                # Grab text and remove empty lines
                parsed_text = soup.get_text()
                parsed_text_lines = [x for x in parsed_text.split("\n") if x.split()]
                with open(os.path.join(write_path, re.sub("[./-]", "_", paper_id) + ".txt"), "w") as f:
                    f.write('\n'.join(parsed_text_lines))
                found_ids.append(paper_id)
            except AttributeError:
                continue
    return(found_ids)

def __pull_scopus_abstracts(ids, output_directory, scopus_api_key, abstract_include_title=True):
    """
    Function that pulls paper abstracts from Scopus. Writes them to a directory.
    Args:
        ids (List): A list of Scopus IDs.
        output_directory (String): Path specifying where to write the papers to.
        scopus_api_key (String): A string API key for Scopus-Elselvier. See documentation on how to acquire.
        abstract_include_title (Boolean, default=True): Whether to include the paper's title as the first sentence of the text.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """

    if scopus_api_key is None:
        raise ValueError("scopus_api_key must be provided. See Elsevier Dev Portal for assistance")
    found_ids = []
    write_path = os.path.join(output_directory, "scopus_abstracts")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path)
    for paper_id in ids:
        req = requests.get("https://api.elsevier.com/content/abstract/doi/" + str(paper_id) + "?apiKey=" + scopus_api_key)
        soup = BeautifulSoup(req.content, 'html.parser')
        try:
            abstract = soup.find("abstract").find("ce:para").get_text()
            if abstract is None:
                continue
            with open(os.path.join(write_path, re.sub("[./-]", "_", str(paper_id)) + ".txt"), "w") as f:
                if abstract_include_title:
                    f.write(soup.find("dc:title").get_text() + ". ")
                f.write(abstract)
            found_ids.append(paper_id)
        except AttributeError:
            continue
    return(found_ids)

def __pull_scopus(ids, output_directory, type, scopus_api_key):
    """
    Function that pulls paper abstracts from Scopus. Writes them to a directory.
    Args:
        ids (List): A list of Scopus IDs.
        output_directory (String): Path specifying where to write the papers to.
        scopus_api_key (String): A string API key for Scopus-Elselvier. See documentation on how to acquire.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """

    # If pulling full text, first pass through clean text and then pdfs 
    if type == "full text":
        return({"full": __pull_scopus_clean(ids, output_directory, scopus_api_key), "abstract": []})
    elif type == "abstract":
        return({"full": [], "abstract": __pull_scopus_abstracts(ids, output_directory, scopus_api_key)})
    elif type == "both":
        found_ids_clean = __pull_scopus_clean(ids, output_directory, scopus_api_key)
        remaining_ids = [the_id for the_id in ids if the_id not in found_ids_clean]
        found_ids_abstract = __pull_scopus_abstracts(remaining_ids, output_directory, scopus_api_key)
        return({"full": found_ids_clean, "abstract": found_ids_abstract})