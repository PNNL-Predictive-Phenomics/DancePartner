import os
import shutil
import DancePartner as dance

# Create output directory for pulling omes
output_directory = os.path.join(os.getcwd(), "test_rellypull")
os.mkdir(output_directory)

def test_uniprot():

    # Pull a small interactome
    dance.pull_uniprot(6239, output_directory, verbose = False)
    shutil.rmtree(output_directory)