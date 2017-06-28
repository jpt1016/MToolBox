#!/usr/bin/env python

import subprocess
import argparse
import os
import csv
import re
import copy

#defines CL args
parser = argparse.ArgumentParser(description='Script that runs the MToolbox pipeline. Can annotate, manipulate and combine the MToolbox VCFs')
parser.add_argument('-v', help='if this argument is passed, filtering will be done to remove variants not meeting this threshold')
parser.add_argument('-p', help='path of input files')
parser.add_argument('-i', help='the MToolbox config file')
parser.add_argument('-a', help='path to annovar databases')
parser.add_argument('-o', help='output file prefix')
parser.add_argument('-c', help='include path to existing VCF if merging new VCF with previously generated VCF')

#stores CL arguments
args=parser.parse_args()
config_file = args.i
variant_threshold = args.v
input_path = args.p
db_path = args.a

#runs the MToolbox script
if os.path.isfile(config_file) and os.path.exists(input_path):
	os.chdir(input_path)
	os.system("bash /RQusagers/jtsaip/soft/packages/MToolBox/MToolBox/MToolBox.sh -i " + str(config_file) + " >> MToolbox.out 2>> MToolbox.out")

	#changes chrMT to MT
	if os.path.isfile(VCF_file.vcf):
		os.system("sed -i 's/^chrMT/MT/g' VCF_file.vcf")

	else:
		print "Error: MToolbox VCFOutput failed. Exiting."
		quit()	
#At least one path/file did not exist
else:
	print "Error: config file and/or input files do not exist. Exiting"
	quit()

#GATK command to combine two VCFs
def combine(existing, new):
	
	if os.path.isfile(existing):
		os.system("java -jar /RQexec/dionnela/soft/packages/GenomeAnalysisTK-3.7/dist/GenomeAnalysisTK.jar -T CombineVariants -R /RQexec/dionnela/data/pipeline.svn/data/reference/human_g1k_v37_MT.fasta --variant " + str(existing) + "--variant " + str(new) + " -o final_combined.vcf")
	else:
		print "Error: existing VCF file not found. Exiting."
		quit()

#filter the VCF file
def dofilter(threshold):
	
	#first opens file just for writing headers
	with open(str(input_path) + 'VCF_file.vcf', 'rb') as inFile:
		plainFile = open("MToolbox.vcf", 'wb')

		reader = csv.reader(inFile, delimiter='\t')
	
		headers = []
		y=0
		for ln in reader:

			#checks for additional header information
			if ln[0].startswith("##"):
				line = ""

				for element in ln:
					line = line + str(element)
					plainFile.write(line + "\n")
				
					y+=1

					#places the correct format entry in the correct place
					if y == 7:
						plainFile.write('##FORMAT=<ID=VC,Number=1,Type=String,Description="Variant copy number">\n')		

			#checks for the header with column names 
			elif (ln[0].startswith("#CHROM")):
				plainFile.close()
				outFile = open("MToolbox.vcf", 'a')
				writer = csv.writer(outFile, quoting=csv.QUOTE_NONE, delimiter='\t', escapechar='\n', lineterminator="\n")
				writer.writerow(ln)

			#starts writting body of file
			else:
				i = 0
				for element in ln:
					#changes '0' into './.'
					if element == str(0):
						
						ln[i] = "./."
						#writer.writerow(ln)
					
					#adds 'VC' into format field 
					elif "GT:" in element:
						if threshold:
							newList = element.split(":")
							newList.append("VC")
							ln[i]=":".join(newList)

					#parses and removes certain calls
					elif re.search(r':\d', element):
						if threshold:
							newList = element.split(":")

							#checks for calls with multiple alleles (denoted with a comma)
							if "," in newList[2]:
                                                        	print "Row corresponding with: " + str(ln[i]) + " contains multiple variant alleles."
                                                        	het_freqs = copy.copy(newList[2])
                                                        	VCList = []
                                                        	for freq in het_freqs.split(','):
                                                               		VC = float(freq) * int(newList[1])
                                                                	VCList.append(int(round(VC)))
                                                        	count = 0
                                                        	for vc in VCList:
                                                                	count += vc
                                                        	if int(count) < int(threshold):
                                                                	print "Row corresponding with: " + str(ln[i]) + " did not meet variant threshold \"" + str(threshold) + "\" count= " + str(count) + ", genotype will be dropped at this position."
                                                                	ln[i]="./."
                                                        	else:
                                                                	VCList_str = ','.join(str(e) for e in VCList)
                                                                	newList.append(VCList_str)
                                                                	ln[i]=":".join(newList)
							else:

								#drops the calls that do not meet the threshold
								variant_reads = int(round(float(newList[1]) * float(newList[2])))

								if int(variant_reads) < int(threshold):
									print "Row corresponding with: " + str(ln[i]) + " did not meet variant threshold \"" + str(threshold) + "\", \"" + str(variant_reads) + "\" genotype will be dropped at this position."
									ln[i]="./."

								#the call meets all thresholds and will not be removed
								else:
									newList.append(str(variant_reads))
									ln[i]=":".join(newList)

					i+=1
				#writes the correctly formated line
				writer.writerow(ln)
	inFile.close()
	outFile.close()



#Outputs annotated VCF file 
def annotate(path, prefix):
	output_name = str(prefix)
	os.system("perl ~/runs/annovar/table_annovar.pl new.vcf " + str(path) + " -vcfinput -buildver hg19 -protocol CodingControl.csv,ControlRegion.csv,PolymorphismsCoding.csv,RNA.csv,somatic.csv,mitimpact27,HmtDB_healthy,HmtDB_pa -operation f,f,f,f,f,f,f,f -nastring .")
#	os.system("python ~/runs/annovar/scripts/AnnovarEmptyAnoRemover.py -f new.vcf.hg19_multianno.vcf -o " + output_name)
	if os.path.isfile('new.vcf.hg19_multianno.vcf'):
		with open(new.vcf.hg19_multianno.vcf, 'rb') as inFile:
			reader = csv.reader(inFile, delimiter='\t')
			outFile = open(output_name +".vcf", 'wb')
			writer = csv.writer(outFile, quoting=csv.QUOTE_NONE, delimiter='\t', escapechar='\n', lineterminator="\n")

			for ln in reader:
				if ln[0].startswith("##"):
					row = ''.join(ln)
					outFile.write(row + "\n")
					continue

				elif ln[0].startswith("#"):
					row = '\t'.join(ln)
					outFile.write(row + "\n")
					continue

				else:
					annotation=[]
					annotation = ln[7].split(";")
					annotation = [x for x in annotation if "=." not in x]
					ln[7] = ';'.join(annotation)
					writer.writerow(ln)

		inFile.close()
		outFile.close()
	else:
		print "Error: Annotation step failed, no annotated VCF file found. Exiting"
		quit()

#calls the function if CL argument provided
if (variant_threshold):
	dofilter(variant_threshold)
	os.system("java -jar /RQexec/dionnela/soft/packages/GenomeAnalysisTK-3.7/dist/GenomeAnalysisTK.jar -T SelectVariants -R /RQexec/dionnela/data/pipeline.svn/data/reference/human_g1k_v37_MT.fasta -o new.vcf --excludeNonVariants -V MToolbox.vcf")


#calls the anotate function, removes intermediate vcf
if args.a and os.path.exists(args.a):
	annotate(db_path, args.o)
	os.system('rm new.vcf* MToolbox.vcf*')

#combines vcf with additional samples with vcf of existing samplnees 
if args.c and os.path.isfile(args.c):
	combine(output_name + ".vcf", str(args.c))

#Either filtering, annotating or combining not flagged, or the annovar database directory was not found, exits the script
else:
	print "Annovar database directory not found, exiting." 
	quit()
