import os
import requests
from bs4 import BeautifulSoup
import json

def __pull_osti_clean(ids, output_directory):
    """
    Function that pulls clean text papers from OSTI. Writes papers to a directory.
    Args:
        ids (List): A list of OSTI IDs.
        output_directory (String): Path specifying where to write the papers to.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """
    found_ids = []
    write_path = os.path.join(output_directory, "osti_clean")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path)
    for id in ids:
        try:
            paper_id = str(int(float(id)))
            req = requests.get("https://www.osti.gov/api/v1/records/" + paper_id )
            data = json.loads(req.content)
            # Data should be list, otherwise cannot find
            if not isinstance(data, list):
                continue
            try:
                if "links" in data[0].keys():
                    for link in data[0]['links']:
                        if "fulltext" in link.values():
                            text_url = link['href']
                            text_req = requests.get(text_url)
                            text_soup = BeautifulSoup(text_req.content, features="lxml")
                            fulltext_url = text_soup.find("a", {"title":"Document DOI URL", "data-product-type":"Journal Article"}).get_text()
                            fulltext_req = requests.get(fulltext_url)
                            #print(fulltext_req.status_code)
                            fulltextsoup = BeautifulSoup(fulltext_req.content)
                            # Attempt to Remove Non Document Text
                            for x in fulltextsoup.find_all("div", {"class": "References"}):
                                x.decompose()
                            for x in fulltextsoup.find_all("form"):
                                x.decompose()
                            for x in fulltextsoup.find_all("select"):
                                x.decompose()
                            for x in fulltextsoup.find_all("section", {"data-title": "References"}):
                                x.decompose()
                            parsed_text = fulltextsoup.get_text()
                            parsed_text_lines = [x for x in parsed_text.split("\n") if x.split()]
                            # Bad request (empty or redirected)
                            if len(parsed_text_lines) < 5:
                                continue
                            filename = os.path.join(write_path, paper_id + ".txt")
                            with open(filename, 'w') as f:
                                f.write("\n".join(parsed_text_lines))
                            found_ids.append(paper_id)
                            #break from for loop over links once fulltext found
                            break
            except AttributeError:
                continue
        except:
            continue

    return(found_ids)

def __pull_osti_abstracts(ids, output_directory, abstract_include_title=True):
    """
    Function that pulls paper abstracts from OSTI. Writes them to a directory.
    Args:
        ids (List): A list of OSTI IDs.
        output_directory (String): Path specifying where to write the papers to.
        abstract_include_title (Boolean, default=True): Whether to include the paper's title as the first sentence of the text.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """
    found_ids = []
    write_path = os.path.join(output_directory, "osti_abstracts")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path)
    for id in ids:
        try:
            paper_id = str(int(float(id)))
            req = requests.get("https://www.osti.gov/api/v1/records/" + paper_id)
            try:
                data = json.loads(req.content)[0]
                abstract = BeautifulSoup(data['description'], features="lxml").find("p").get_text()
                if abstract is None:
                    continue
                with open(os.path.join(write_path, paper_id + ".txt"), "w") as f:
                    if abstract_include_title:
                        f.write(BeautifulSoup(data['title'], features="lxml").get_text() + ". ")
                    f.write(abstract)
                found_ids.append(paper_id)
            except AttributeError:
                continue
            except IndexError:
                continue
            except KeyError:
                continue
        except:
            continue
    return(found_ids)

def __pull_osti(ids, output_directory, type):
    """
    Function that pulls text from OSTI.
    Args:
        ids (List): A list of OSTI IDs.
        output_directory (String): Path specifying where to write the papers to.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """

    # If pulling full text, first pass through clean text and then pdfs 
    if type == "full text":
        return({"full": __pull_osti_clean(ids, output_directory), "abstract": []})
    elif type == "abstract":
        return({"full": [], "abstract": __pull_osti_abstracts(ids, output_directory)})
    elif type == "both":
        found_ids_clean = __pull_osti_clean(ids, output_directory)
        remaining_ids = [the_id for the_id in ids if the_id not in found_ids_clean]
        found_ids_abstract = __pull_osti_abstracts(remaining_ids, output_directory)
        return({"full": found_ids_clean, "abstract": found_ids_abstract})