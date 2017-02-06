"""Searches PubMed for orchid related articles and prints out the
first 5 such Medline records.
"""
import sys
sys.path.insert(0,'/biol/programs/python/lib/python/')

from Bio import PubMed
from Bio import Medline
import string

rec_parser = Medline.RecordParser()
medline_dict = PubMed.Dictionary(parser = rec_parser)



search_term = input('Please eneter a search term:   ')
orchid_ids = PubMed.search_for(search_term)

no_records = int(input('How many Records would you like returned?   '))

for id in orchid_ids[0:no_records]:
    cur_record = medline_dict[id]
    print'title:', string.rstrip(cur_record.title)
    print'authors:', cur_record.authors
    print'source:', string.strip(cur_record.source)
    print

