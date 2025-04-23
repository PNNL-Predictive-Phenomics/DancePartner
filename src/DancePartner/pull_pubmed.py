import os
import requests
from bs4 import BeautifulSoup
import tarfile
import metapub
import pypdf
import io
import urllib3


def __pull_pubmed_clean(ids, output_directory, tarball_path):
    """
    Function that pulls paper abstracts from PubMed. Writes them to a directory.
    Args:
        ids (List): A list of PubMed IDs.
        output_directory (String): Path specifying where to write the papers to.
        tarball_path (Optional, String): An optional path of where to write the (large) tarball files to. Can also be used to specify a tarball path where a previous function run may have saved articles to, which can reduce run time.
    Returns:
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
        # Find if the PubMed Article has a corresponding PubMed Central ID and page
        req = requests.get("https://pubmed.ncbi.nlm.nih.gov/" + str(pmid) + "/")
        soup = BeautifulSoup(req.content, 'html.parser')
        pmc_url = soup.find_all("a", class_="id-link", attrs={"data-ga-action":"PMCID"})
        if len(pmc_url) > 0:
            try:
                # Use that PubMedCentral ID to find where the article is stored on FTP
                pmcid = pmc_url[0].get_text().strip()
                if pmcid in pmc_list:
                    # tarball has already been downloaded to the `tarball_path`. Don't re-download. Break from loop to next pmid.
                    continue
                link = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=" + pmcid
                tgz_url = "https://" + BeautifulSoup(requests.get(link).content, 'html.parser').find("link", attrs={"format":"tgz"}).get("href")[6:]
                response = requests.get(tgz_url, stream=True)
                # Download the tarball from the FTP location
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
            except urllib3.exceptions.ProtocolError:
                notfound_count += 1

    # Now grab the text from the tarballs
    for _, _, files in os.walk(tarball_path):
        for file in files:
            # Grab .nxml file in each tarball
            if ".tar.gz" in file:
                try:
                    tar = tarfile.open(os.path.join(tarball_path, file))
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
                            with open(file_name, "w") as f:
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

def __pull_pubmed_pdfs(ids, output_directory):
    """
    Function that pulls PDF papers from PubMed. Writes them to a directory.
    Args:
        ids (List): A list of PubMed IDs.
        output_directory (String): Path specifying where to write the papers to.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """

    notfound_count = 0
    write_path = os.path.join(output_directory, "pubmed_pdfs")
    if os.path.exists(write_path) == False:
        os.mkdir(write_path, mode = 0o777)
    found_ids = []

    # Iterate through list and try to scan the pdf and save to a folder
    for id in ids:
        pmid = str(int(float(id)))
        try:
            src = metapub.FindIt(str(pmid))
            req = requests.get(src.url)
            pdf = io.BytesIO(req.content)
            reader = pypdf.PdfReader(pdf)
            filename = os.path.join(write_path, str(pmid) + ".txt")
            with open(filename, 'w') as f:
                for i in range(len(reader.pages)):
                    f.write(" ".join(reader.pages[i].extract_text().split("\n"))) 
            # success --> append id (Integer type) to found-list
            found_ids.append(pmid)

        except requests.exceptions.MissingSchema:
            #print("Invalid URL for article {}".format(pmid))
            notfound_count += 1
        except pypdf._utils.PdfStreamError:
            #print("PDF Stream Error with article {}".format(pmid))
            notfound_count += 1
        except pypdf.generic._data_structures.PdfReadError:
            #print("PDF Read Error with article {}".format(pmid))
            notfound_count += 1
        except metapub.exceptions.InvalidPMID:
            #print("PubMed invalid article error for article {}".format(pmid))
            notfound_count += 1
        except AttributeError:
            #print("Attribute Error with article {}".format(pmid))
            notfound_count += 1
        except TypeError:
            #print("Type Error with article {}".format(pmid))
            notfound_count += 1
        except UnicodeEncodeError:
            #print("Encoding Error with article {}".format(pmid))
            notfound_count += 1
    return(found_ids)

def __pull_pubmed_abstracts(ids, output_directory, abstract_include_title=True):
    """
    Function that pulls paper abstracts from PubMed. Writes them to a directory.
    Args:
        ids (List): A list of PubMed IDs.
        output_directory (String): Path specifying where to write the papers to.
        abstract_include_title (Boolean, default=True): Whether to include the paper's title as the first sentence of the text.
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
        url = "https://pubmed.ncbi.nlm.nih.gov/" + str(pmid) + "/"
        req = requests.get(url)
        soup = BeautifulSoup(req.content, "html.parser")
        try:
            abstract = soup.find(id="eng-abstract").get_text().strip()
            with open(os.path.join(write_path, str(pmid) + ".txt"), "w")as f:
                if abstract_include_title:
                    f.write(soup.find("meta", {"name":"citation_title"})['content'])
                    f.write(". ")
                f.write(abstract)
            # success --> append id (Integer type) to found-list
            found_ids.append(pmid)
        except AttributeError:
            try:
                #print("{} for article {}".format(soup.find(class_="empty-abstract").get_text(), pmid))
                notfound_count += 1
            except AttributeError:
                #print("Error with article {}".format(pmid))
                notfound_count += 1
        
    return(found_ids)

def __pull_pubmed(ids, output_directory, type, tarball_path):
    """
    Function to pull papers from PubMed.
    Args:
        ids (List): A list of PubMed IDs.
        output_directory (String): Path specifying where to write the papers to.
        type (String): Either "full text" to pull only full text, "abstract" to pull only abstracts, or "both" to first prioritize full text,
        and then prioritize abstracts. 
        tarball_path (Optional, String): An optional path of where to write the (large) tarball files to. Can also be used to specify a tarball path where a previous function run may have saved articles to, which can reduce run time.
    Returns:
        List of IDs that were found. A subset of the `ids` argument.
    """

    # If pulling full text, first pass through clean text and then pdfs 
    if type == "full text":
        found_ids_clean = __pull_pubmed_clean(ids, output_directory, tarball_path)
        remaining_ids = [the_id for the_id in ids if the_id not in found_ids_clean]
        found_ids_pdf = __pull_pubmed_pdfs(remaining_ids, output_directory)
        found_ids_clean.extend(found_ids_pdf)
        return({"full": found_ids_clean, "abstract": []})
    elif type == "abstract":
        return({"full": [], "abstract": __pull_pubmed_abstracts(ids, output_directory)})
    elif type == "both":
        found_ids_clean = __pull_pubmed_clean(ids, output_directory, tarball_path)
        remaining_ids = [the_id for the_id in ids if the_id not in found_ids_clean]
        found_ids_pdf = __pull_pubmed_pdfs(remaining_ids, output_directory)
        remaining_ids = [the_id for the_id in remaining_ids if the_id not in found_ids_pdf]
        found_ids_abstract = __pull_pubmed_abstracts(remaining_ids, output_directory)
        found_ids_clean.extend(found_ids_pdf)
        return({"full": found_ids_clean, "abstract": found_ids_abstract})