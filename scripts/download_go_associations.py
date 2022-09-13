import os
import gzip
import shutil
import requests

def download_file(url):
    """
    Download file from given url.

    Returns
    -------
    str :
        Path and name of the downloaded file.
    """
    base = os.path.basename(url)

    with open(base, "wb") as f:
        r = requests.get(url)
        f.write(r.content)

    return os.path.join(os.path.abspath(os.getcwd()), base)

def unpack_file(file_in, file_out):
    """ Unpack given file to output file. """
    with gzip.open(file_in, 'rb') as f_in:
        with open(file_out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

if __name__ == "__main__":
    # import pdb; pdb.set_trace() # uncomment for interactive debugging

    for key in snakemake.params.keys():
        url = snakemake.params[key]

        file = download_file(url)

        if file.endswith(".gz"):
            unpack_file(file, snakemake.output[key])
        else:
            shutil.copyfile(file, snakemake.output[key])
