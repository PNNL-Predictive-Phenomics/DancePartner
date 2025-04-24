import os
import shutil
import DancePartner as dance

## How to calculate coverage (from within main package directory): 
# coverage run --source=DancePartner -m pytest -x tests/* -W ignore
# coverage report
# coverage html

def test_uniprot():
    output_directory = os.path.join(os.getcwd(), "test_uniprot")
    os.mkdir(output_directory)
    dance.pull_uniprot(6239, output_directory, verbose = False)
    shutil.rmtree(output_directory)