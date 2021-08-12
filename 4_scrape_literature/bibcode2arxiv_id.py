import os
import re
import sys
import numpy as np
import argparse

# extract arxiv IDs from bibtex file

parser = argparse.ArgumentParser(description='Arxiv')
parser.add_argument('--planet', '-p')

args = parser.parse_args()

planet = args.planet

def main(planet):
	with open(f'/Users/kate/Desktop/arxiv/article_database/{planet}/{planet}.bib') as bibtex_file:
	    bibtex_str = bibtex_file.read()


	regex_search_term = '{'
	regex_replacement = ''
	text_after = re.sub(regex_search_term, regex_replacement, bibtex_str)
	regex_search_term = '}'
	regex_replacement = ''
	text_after = re.sub(regex_search_term, regex_replacement, text_after)
	regex_search_term = ','
	regex_replacement = ''
	text_after = re.sub(regex_search_term, regex_replacement, text_after)

	splitted_bibtex = text_after.split()


	IDs = []
	not_on_arxiv = []
	for i in range(2,len(splitted_bibtex)):
	  if splitted_bibtex[i-1] == '=' and splitted_bibtex[i-2] == 'eprint':
	    #print('ID: ', splitted_bibtex[i])
	    IDs.append(splitted_bibtex[i])
	 # else:
	 # 	not_on_arxiv.append(splitted_bibtex[i])


	np.savetxt(f'/Users/kate/Desktop/arxiv/article_database/{planet}/arxiv_ids.txt', IDs, fmt='%s')
	#np.savetxt(f'/Users/kate/Desktop/arxiv/article_database/{planet}/not_on_arxiv.txt', not_on_arxiv, fmt='%s')