import os
import shutil
import DancePartner as dance

# Create output directory for pulling omes
output_directory = os.path.join(os.getcwd(), "test_omes")
os.mkdir(output_directory)

## How to calculate coverage (from within main package directory): 
# coverage run --source=DancePartner -m pytest -x tests/* -W ignore
# coverage report
# coverage html

def test_pull_proteome():

    # Pull a small proteome - like the one for a butterfly
    dance.pull_proteome("UP000464024", output_directory)

    # Try a nonsensical pull
    dance.pull_proteome("ugly_carrot", output_directory)
    shutil.rmtree(output_directory)


