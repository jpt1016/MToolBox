#!/usr/bin/env python
# coding: utf-8

###Description: Removes blank annotations from Annovar 
#Version 1.0

#import modules
import csv
import sys
import re
import argparse
import string
import types

parser = argparse.ArgumentParser(description='Manipulate Annovar VCF')
#parser.add_argument('infile', '-i', type=argparse.FileType('r'), help='Path to the input file to transform')
parser.add_argument('-f', help='path to file')
parser.add_argument('-o', help='output file prefix')
args=parser.parse_args()

with open(args.f, 'rb') as inFile:
	
	reader = csv.reader(inFile, delimiter='\t')
	outFile = open(str(args.o) +".vcf", 'wb')
	writer = csv.writer(outFile, quoting=csv.QUOTE_NONE, delimiter='\t', escapechar='\n', lineterminator="\n") 
	
	for ln in reader:
#		print ln
		if ln[0].startswith("##"):
			row = ''.join(ln)
#			print row
			outFile.write(row + "\n")
			continue

		elif ln[0].startswith("#"):
			row = '\t'.join(ln)
			outFile.write(row + "\n")
#			print row 
			continue

		else:
			annotation=[]
			annotation = ln[7].split(";")
			#print "og"
			#print annotation
			
			#for word in annotation:
				
			annotation = [x for x in annotation if "=." not in x]
				#if '=.' in word:
				#	print word
				#	annotation.remove(word)
				#	print "does not contain " + word 

				#	continue
					#print "this word does not contain =. " + word
			#print annotation
			#for x in annotation:
			#	if ' ' in x:
			#		print x
			#		x = x.replace(" ", "")
					#print x
			#print annotation
			#final_anno = [annotation.strip(' ') for annotation in final_anno]
			ln[7] = ';'.join(annotation)
			writer.writerow(ln)

inFile.close()
outFile.close()
