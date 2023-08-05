import Bio.Blast
from Bio.Blast import NCBIWWW
result_handle= NCBIWWW.qblast('blastn', 'nt', 'TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG')

from Bio.Blast import NCBIXML
blast_record= NCBIXML.read(result_handle)


E_value_thresh= 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
            if hsp.expect < E_value_thresh:
                print('***Alignment***')
                print('Sequence:', alignment.title)
                print('Length:', alignment.length)
                print('E Value:', hsp.expect)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)
