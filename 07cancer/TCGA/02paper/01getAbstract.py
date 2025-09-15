#  python 01getAbstract.py pmid.txt abstract.txt
from Bio import Entrez, Medline
import sys
input = sys.argv[1] # input = "pmid.txt"
output = sys.argv[2] # output = "abstract.txt"

# Set your Entrez email address
Entrez.email = "sky.alex.jww@gmail.com"

# allAbstract={}
with open(input, "r") as file, open(output, "w") as w:
    for line in file:
        pubmed_id = line.strip() #  pubmed_id = "15118073"
    # Retrieve records for a PubMed ID via efetch
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
        records = Medline.parse(handle)
    # Extract abstract from the record
        abstract = None
        for record in records:
            abstract = record.get("AB", None)
            # print(f"PubMed ID: {pubmed_id}\nAbstract: {abstract}")
        if abstract is not None:
            w.write(pubmed_id+"\t"+abstract+"\n")

