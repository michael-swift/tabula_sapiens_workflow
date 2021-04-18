import sys
import requests
from lxml import html
from datetime import datetime


def main(gene, species='Homo+sapiens', outfile=None):
    """ 
    Script for extracting fasta from IMGT website text
    Args
        gene: e.g. IGHV, IGHJ
    """

    if not outfile:
        today = datetime.now().date()
        outfile = '%s_%s%.02d%02d.fa' % (gene, today.year, today.month, today.day)

    imgt_url = 'http://www.imgt.org/genedb/GENElect?query=7.14+%s&species=%s' % (gene, species)

    response = requests.get(imgt_url)
    tree = html.fromstring(response.content)
    fasta_text = tree.xpath('//*[@id="GENE-DB"]/pre[2]/text()')[0]

    # strip handles leading and trailing whitespace (e.g. newlines)
    with open(outfile, 'w') as out:
        out.write(fasta_text.strip() + '\n')

if __name__ == "__main__":
    main(sys.argv[1])
    
