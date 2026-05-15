import os
import requests
from bs4 import BeautifulSoup
import tarfile
import pypdf
import io
import urllib3
import importlib


def __get_ncbi_params(pubmed_api_key: str = None, **params):
    if pubmed_api_key is not None:
        params["api_key"] = pubmed_api_key
    return params


def __get_metapub_module():
    return importlib.import_module("metapub")


def __fetch_pubmed_record(pmid: str, pubmed_api_key: str = None):
    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params = __get_ncbi_params(pubmed_api_key, db = "pubmed", id = pmid, retmode = "xml"),
        timeout = 30,
    )
    response.raise_for_status()
    return BeautifulSoup(response.content, "xml")


def __extract_article_identifier(soup: BeautifulSoup, identifier_type: str):
    for tag in soup.find_all(lambda entry: entry.name is not None and entry.name.lower() in {"articleid", "article-id"}):
        attrs = {str(key).lower(): str(value) for key, value in tag.attrs.items()}
        if attrs.get("idtype") == identifier_type or attrs.get("pub-id-type") == identifier_type:
            return tag.get_text().strip()
    return None


def __extract_article_title(soup: BeautifulSoup):
    for tag_name in ["ArticleTitle", "article-title"]:
        tag = soup.find(tag_name)
        if tag is not None:
            return tag.get_text(" ", strip = True)
    return None


def __extract_article_abstract(soup: BeautifulSoup):
    abstract_nodes = []
    for tag_name in ["AbstractText", "abstracttext"]:
        abstract_nodes.extend(soup.find_all(tag_name))

    abstract_parts = [node.get_text(" ", strip = True) for node in abstract_nodes if node.get_text(strip = True)]
    if len(abstract_parts) == 0:
        return None
    return " ".join(abstract_parts)


def __pull_pubmed_clean(ids: list[str], output_directory: str, tarball_path: str, pubmed_api_key: str = None):
    """
    Function that pulls paper abstracts from PubMed. Writes them to a directory.

    Parameters
    ----------
    ids
        A list of PubMed IDs.
    
    output_string
        Path specifying where to write the papers to.
    
    tarball_path
        An optional path of where to write the (large) tarball files to. Can also be used to specify a tarball path where a previous function 
        run may have saved articles to, which can reduce run time.

    pubmed_api_key
        An optional NCBI API key for PubMed E-utilities requests.

    Returns
    -------
        List of IDs that were found. A subset of the `ids` argument.
    """
    
    notfound_count = 0
    found_ids = []

    pmc_list = []
    if tarball_path is None:
        tarball_path = os.path.join(output_directory, "pubmed_tarballs")
        if os.path.exists(tarball_path) == False:
            os.mkdir(tarball_path, mode = 0o777)
    else:
        # Create a list of pre-written PMC names if tarball_path has been pre-specified
        for _, _, files in os.walk(tarball_path):
            # If files are present, we will make a list of what's in there already and then append new files to this directory instead of 
            #   a subfolder of output_dire
            # Else, we will write new files to this directory instead of writing into our output_dir
            if len(files) > 0:
                for file in files:
                    if ".tar.gz" in file and "PMC" in file:
                        pmc_list.append(file.split(".")[0])

    write_path = os.path.join(output_directory, "pubmed_clean")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path, mode = 0o777)

    # First find tarballs and download them into an internal directory # Delete them?
    for id in ids:
        pmid = str(int(float(id)))
        try:
            soup = __fetch_pubmed_record(pmid, pubmed_api_key)
            pmcid = __extract_article_identifier(soup, "pmc")
            if pmcid is None:
                continue
            if pmcid in pmc_list:
                continue

            response = requests.get(
                "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi",
                params = __get_ncbi_params(pubmed_api_key, id = pmcid),
                timeout = 30,
            )
            tgz_link = BeautifulSoup(response.content, 'html.parser').find("link", attrs={"format":"tgz"})
            if tgz_link is None:
                notfound_count += 1
                continue

            tgz_url = "https://" + tgz_link.get("href")[6:]
            response = requests.get(tgz_url, stream=True, timeout = 60)
            if response.status_code == 200:
                filename = os.path.join(tarball_path, pmcid + ".tar.gz")
                with open(filename, 'wb') as f:
                    f.write(response.raw.read())
            else:
                notfound_count += 1
        except AttributeError:
            notfound_count += 1
        except TimeoutError:
            notfound_count += 1
        except requests.exceptions.RequestException:
            notfound_count += 1
        except urllib3.exceptions.ProtocolError:
            notfound_count += 1

    # Now grab the text from the tarballs
    for _, _, files in os.walk(tarball_path):
        for file in files:
            # Grab .nxml file in each tarball
            if ".tar.gz" in file:
                try:
                    tar = tarfile.open(os.path.join(tarball_path, file), encoding = "utf8")
                    for member in tar.getmembers():
                        # Each tarball should have one .nxml file that contains the full article
                        if ".nxml" in member.name:
                            f = tar.extractfile(member)
                            content = f.read()  
                            # Create text file from xml (html parsed with Beautiful Soup)
                            soup = BeautifulSoup(content, "html.parser")
                            # Remove tables and certain math objects from xml
                            for x in soup.find_all('table-wrap'):
                                x.decompose()
                            for x in soup.find_all('mml:annotation'):
                                x.decompose()
                            pmid = soup.find("article-id", attrs={"pub-id-type":"pmid"}).get_text()
                            file_name = os.path.join(write_path, str(pmid) + ".txt")
                            with open(file_name, "w", encoding = "utf8") as f:
                                for p in soup.find_all("p", recurisve=False):
                                    f.write(p.get_text())
                            # success --> append id (Integer type) to found-list
                            found_ids.append(pmid.strip())
                            #.nxml found and txt written, go to next tarball
                            tar.close()
                            break                         
                except tarfile.ReadError:
                    # Some cases where a tarball downloaded, but it's empty ??
                    notfound_count += 1
    
    return(found_ids)

def __pull_pubmed_pdfs(ids: list[str], output_directory: str, pubmed_api_key: str = None):
    """
    Function that pulls PDF papers from PubMed. Writes them to a directory.

    Parameters
    ----------
    ids
        A list of PubMed IDs as strings.
    
    output_directory
        Path specifying where to write the papers to.

    pubmed_api_key
        An optional NCBI API key. When provided, it is forwarded to metapub via
        the NCBI_API_KEY environment variable.
    
    Returns
    -------
        List of IDs that were found. A subset of the `ids` argument.
    """

    notfound_count = 0
    write_path = os.path.join(output_directory, "pubmed_pdfs")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path, mode = 0o777)
    found_ids = []

    original_api_key = os.environ.get("NCBI_API_KEY")
    if pubmed_api_key is not None:
        os.environ["NCBI_API_KEY"] = pubmed_api_key

    try:
        metapub = __get_metapub_module()

        # Iterate through list and try to scan the pdf and save to a folder
        for id in ids:
            pmid = str(int(float(id)))
            try:
                src = metapub.FindIt(str(pmid))
                req = requests.get(src.url, timeout = 60)
                pdf = io.BytesIO(req.content)
                reader = pypdf.PdfReader(pdf)
                filename = os.path.join(write_path, str(pmid) + ".txt")
                with open(filename, 'w', encoding = "utf8") as f:
                    for i in range(len(reader.pages)):
                        f.write(" ".join(reader.pages[i].extract_text().split("\n"))) 
                # success --> append id (Integer type) to found-list
                found_ids.append(pmid)

            except requests.exceptions.MissingSchema:
                notfound_count += 1
            except pypdf._utils.PdfStreamError:
                notfound_count += 1
            except pypdf.generic._data_structures.PdfReadError:
                notfound_count += 1
            except metapub.exceptions.InvalidPMID:
                notfound_count += 1
            except metapub.exceptions.MetaPubError:
                notfound_count += 1
            except AttributeError:
                notfound_count += 1
            except TypeError:
                notfound_count += 1
            except UnicodeEncodeError:
                notfound_count += 1
            except requests.exceptions.RequestException:
                notfound_count += 1
    finally:
        if pubmed_api_key is not None:
            if original_api_key is None:
                del os.environ["NCBI_API_KEY"]
            else:
                os.environ["NCBI_API_KEY"] = original_api_key
    return(found_ids)

def __pull_pubmed_abstracts(ids: str, output_directory: str, abstract_include_title: bool = True, pubmed_api_key: str = None):
    """
    Function that pulls paper abstracts from PubMed. Writes them to a directory.
    

    Parameters
    ----------
    ids
        A list of PubMed IDs.
        
    output_directory
        Path specifying where to write the papers to.
    
    abstract_include_title
        Whether to include the paper's title as the first sentence of the text.

    pubmed_api_key
        An optional NCBI API key for PubMed E-utilities requests.
    
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """ 

    write_path = os.path.join(output_directory, "pubmed_abstracts")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path, mode = 0o777)
    found_ids = []
    notfound_count = 0

    for id in ids:
        pmid = str(int(float(id)))
        try:
            soup = __fetch_pubmed_record(pmid, pubmed_api_key)
            abstract = __extract_article_abstract(soup)
            title = __extract_article_title(soup)
            if abstract is None:
                notfound_count += 1
                continue
            with open(os.path.join(write_path, str(pmid) + ".txt"), "w", encoding = "utf8") as f:
                if abstract_include_title and title is not None:
                    f.write(title)
                    f.write(". ")
                f.write(abstract)
            # success --> append id (Integer type) to found-list
            found_ids.append(pmid)
        except AttributeError:
            notfound_count += 1
        except requests.exceptions.RequestException:
            notfound_count += 1
        
    return(found_ids)

def __pull_pubmed(ids: list[str], output_directory: str, type: str, tarball_path: str, pubmed_api_key: str = None):
    """
    Function to pull papers from PubMed.

    Parameters
    ----------
    ids
        A list of PubMed IDs.
    
    output_directory
        Path specifying where to write the papers to.
        
    type
        Either "full text" to pull only full text, "abstract" to pull only abstracts, or "both" to first prioritize full text, and then prioritize abstracts. 
    
    tarball_path
        An optional path of where to write the (large) tarball files to. Can also be used to specify a tarball path where a previous function run may have saved 
        articles to, which can reduce run time.

    pubmed_api_key
        An optional NCBI API key for PubMed E-utilities requests.

    Returns
    -------
        List of IDs that were found. A subset of the `ids` argument.
    """

    # If pulling full text, first pass through clean text and then pdfs 
    if type == "full text":
        found_ids_clean = __pull_pubmed_clean(ids, output_directory, tarball_path, pubmed_api_key = pubmed_api_key)
        remaining_ids = [the_id for the_id in ids if the_id not in found_ids_clean]
        found_ids_pdf = __pull_pubmed_pdfs(remaining_ids, output_directory, pubmed_api_key = pubmed_api_key)
        found_ids_clean.extend(found_ids_pdf)
        return({"full": found_ids_clean, "abstract": []})
    elif type == "abstract":
        return({"full": [], "abstract": __pull_pubmed_abstracts(ids, output_directory, pubmed_api_key = pubmed_api_key)})
    elif type == "both":
        found_ids_clean = __pull_pubmed_clean(ids, output_directory, tarball_path, pubmed_api_key = pubmed_api_key)
        remaining_ids = [the_id for the_id in ids if the_id not in found_ids_clean]
        found_ids_pdf = __pull_pubmed_pdfs(remaining_ids, output_directory, pubmed_api_key = pubmed_api_key)
        remaining_ids = [the_id for the_id in remaining_ids if the_id not in found_ids_pdf]
        found_ids_abstract = __pull_pubmed_abstracts(remaining_ids, output_directory, pubmed_api_key = pubmed_api_key)
        found_ids_clean.extend(found_ids_pdf)
        return({"full": found_ids_clean, "abstract": found_ids_abstract})