import sys
import os
import argparse
import pysam
import logging

logging.basicConfig(filename='FP_SNPs.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
failed_num = 0

parser = argparse.ArgumentParser(description='This script define which one of two allele variants is reference or alternate using human genome reference file (GRCh38). It raise warning when allele from the reference file is not found in the input file.')

parser.add_argument('-r', '--reference', type=str,
                    help='Name of genome reference file (fasta format)')
parser.add_argument('-i', '--input', type=str,
                    help='Name of file with alleles variants (tsv format)') 
parser.add_argument('-o', '--output', type=str,
                    help='Name of output file (if not stated output file will be named as input file with "REF_ALT" postfix in tsv format)')

args = parser.parse_args()
logging.info('Script execution started with settings: ' + ' '.join(sys.argv[1:]))

try:

    log_file = f"{os.path.splitext(args.input)[0]}_log.txt"
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%D-%M-%Y %H:%M:%S')

    reference_genome = pysam.FastaFile(args.reference) 
    input_file = args.input
    if not args.output:
        name = input_file.split('.')
        output_file = name[0] + '_REF_ALT.tsv'
    else:
        output_file = args.output
    
    logging.info(f'Reference file: {args.reference}')
    logging.info(f'Input file: {args.input}')
    logging.info(f'Output file: {output_file}')

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write("#CROM\tPOS\tID\tREF\tALT\n")
        header = infile.readline()
        line_num = 1
        if header[:28] != "#CROM\tPOS\tID\tallele1\tallele2":
            sys.stderr.write('ERROR. Check the input file header format.\n')
            logging.error(f'Header is not in #CROM\tPOS\tID\tallele1\tallele2 format\n')
            sys.exit(1)
        for line in infile:
            line_num += 1
            if not line.startswith('chr'):
                logging.warning(f'Wrong format in line {line_num}, failed to process.')
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            variant_id = fields[2]
            allele1 = fields[3]
            allele2 = fields[4]
            
            ref_allele = reference_genome.fetch(chrom, pos - 1, pos).upper() 
                
            if allele1.upper() == ref_allele:
                alt_allele = allele2
            elif allele2.upper() == ref_allele:
                alt_allele = allele1
            else: 
                sys.stderr.write(f'Failed to define reference allele for variant {variant_id}, line {line_num}.\n')
                logging.warning(f'Failed to define reference allele for variant {variant_id}, line {line_num}. {ref_allele} is not {allele1}/{allele2}.')
                ref_allele = alt_allele = "NA" 
                failed_num += 1
            
            outfile.write(f"{chrom}\t{pos}\t{variant_id}\t{ref_allele}\t{alt_allele}\n")

    reference_genome.close()
    infile.close()
    outfile.close()
    print(f'The resulting file {output_file} is ready. Failed to process {failed_num} cases.')
    logging.info(f'Script execution completed successfully. Failed to process {failed_num} cases.\n')
except Exception as e:
        sys.stderr.write(f'ERROR. Try -h for script usage.\n')
        logging.error(f'An error occurred during script execution: {e}\n')
        sys.exit(1)


