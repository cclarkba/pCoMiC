####### MAIN FILE TO RUN #######

from find_mutations import find_muts
from induce_mutations import induce_muts
import argparse, sys

def main():
    parser = argparse.ArgumentParser(description = '''"Find" module identifies mutations in patient CFTR CDS sequence when comparing to reference sequence. Compares mutations to input database and calcualtes risk factor.\n
                                     "Induce" module produces new fasta file by inducing mutation on original reference CDS fasta. Does mutate original file.''')
    parser.add_argument("--find", help = "For find module (T if yes) Ex. --find T")
    parser.add_argument("--induce", help = "For induce module (T if yes) Ex. --induce T")
    parser.add_argument("-p", "--patient", help = "For find, file path of patient CFTR CDS sequence in fasta file format.")
    parser.add_argument("-r","--reference", help = "For find, file path of reference CFTR CDS sequence in fasta file format.")
    parser.add_argument("-d", "--database", help = "Database containing mutations as .xlsx or .csv file. See readme file for correct formatting of file.")
    parser.add_argument("-m", "--mutation", help = "Mutation in variant cDNA format. Only required for induce module. Exs. c.2052_2054del, c.981_982insAA, c.3042C>A")
    parser.add_argument("-o", "--output", help = "Optional: Name of folder to contain output.", default = "CFTR_mutants")
    
    argument = parser.parse_args()

    if argument.find and argument.induce:
        print("Please indicate only one module to run at a time.")
        sys.exit()
    if not argument.find and not argument.induce:
        print("Please indicate a module to run. Ex. --find T")
        sys.exit()
    if argument.find:
        if not argument.patient or not argument.reference or not argument.database:
            print("Please enter input for --patient, --reference and --database.")
            sys.exit()
        if not argument.patient.lower().endswith(('.fa','.fasta')) or not argument.reference.lower().endswith(('.fa','.fasta')) or not argument.database.lower().endswith(('.xlsx','.csv')):
            print("Please ensure patient and reference files are fasta files and database is .xlsx or .csv.")
            sys.exit()
        find_muts(argument.patient,argument.reference,argument.databse,argument.output)
    if argument.induce:
        if not argument.mutation or not argument.reference:
            print("Please enter input for both --reference and --mutation")
            sys.exit()
        if not argument.reference.lower().endswith('.fa','.fasta'):
            print("Please ensure reference file is a fasta file.")
            sys.exit()
        induce_muts(argument.mutation,argument.reference,argument.output)

if __name__ == '__main__':
    main()

