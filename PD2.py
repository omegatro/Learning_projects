import sys, argparse, re, logging, urllib, pandas as pd, igraph
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq


class Interface:
    def __init__(self, alphabet, valid_input_type, accession_pattern):
        self.alphabet = alphabet
        self.valid_input_type = valid_input_type
        self.accession_pattern = accession_pattern
        
    def parse_arguments(self):
        '''Method is used to parse command-line arguments, allowing to interact with the database from script.'''
        parser = argparse.ArgumentParser(description='This script can be used to work with local database of bacterial reference genomes in fasta format.')
        req_arg_grp = parser.add_argument_group('Required arguments')
        arg_dict = {
            "--add-fasta": f"Add new record(s) from locally stored fasta file & metadata file of {self.valid_input_type} format. Expected: path_to_fasta path_to_{self.valid_input_type}",
            "--add-ncbi": f"Add new record from NCBI RefSeq database by correct accession number. Expected: valid_accession_number",
            "--add-ncbi-list": f"Add new record(s) from NCBI RefSeq database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--exp-fasta": f"Create a fasta file containing sequences from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--exp-meta": f"Create a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--exp-records": f"Create a fasta file containing sequences & a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--rm-record": f"Remove records from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--ch-header": f"Replace an accession number that exists in local database with user-provided accession number if provided accession number was found in local database. Expected: valid_accession_from_local_db",
            "--ch-tax": f"Replace a taxonomy string that exists in local database with user-provided taxonomy string if provided accession number was found in local database. Expected: valid_accession_from_local_db",
            "--view-data": "Visualize the contents of the local database as a tree chart, showing the number and percentage of records that belong to each taxonomic group."
        }

        for arg in arg_dict.keys():
            if arg != "--view-data":
                req_arg_grp.add_argument(arg, metavar='\b', help = arg_dict[arg],required=False)  
            else:
                req_arg_grp.add_argument(arg, help = arg_dict[arg], action='store_true', required=False)

        if len(sys.argv)==1:
            parser.print_help(sys.stderr)
            sys.exit(1)
        self.args = parser.parse_args()

    def check_format(self, raw_accession):
        '''Method checks if provided accession number matches RefSeq accession number conventions.'''
        if re.match(self.accession_pattern, raw_accession): return True
        else: return False

    def check_alphabet(self, raw_sequence):
        '''Method checks if sequence that corresponds to the id contains only allowed characters.'''
        for letter in str(raw_sequence):
            if letter not in self.alphabet:
                return False
        return True

    def check_file_type(self, file_name):
        '''Method checks if provided file type matches allowed file type.'''
        return file_name.endswith(f'.{self.valid_input_type}')

    def read_local(self, valid_accession):
        #Database class required to test.
        local_fasta = ""
        taxonomy = ""
        return local_fasta, taxonomy

    def print_tax(self, taxonomy):
        #Database class required to test.
        pass


class Database:
    def __init__(self, db_name):
        self.sequence_file = Path(f'{db_name}.fasta')
        self.taxomony_file = Path(f'{db_name}.csv')
        
    def create_db_files(self):
        '''Method is used to create database file in the directory where the script is.'''
        self.sequence_file.touch(exist_ok=True)
        self.taxomony_file.touch(exist_ok=True)
        tax_file_columns = pd.DataFrame(columns=["accession_number","taxonomy_string"])
        tax_file_columns.to_csv(self.taxomony_file, index=False, header=True)

    def add_record(self, header, sequence, taxonomy, description):
        '''Method is used to add new record to database files (taxonomy file & fasta file).'''
        new_record_seq = SeqIO.SeqRecord(Seq(sequence),id=header, description=description)
        new_record_tax = pd.DataFrame.from_dict({"accession_number":[header],"taxonomy_string":[taxonomy]})
        with open(self.sequence_file, "a+") as seq_db:
            SeqIO.write(new_record_seq, seq_db, "fasta")
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        tax_db = tax_db.append(new_record_tax)
        tax_db.to_csv(self.taxomony_file, index=False)
        

    def calculate_content(self, count_by_group_path=None):
        '''Method is used to calculate current database content by taxonomic group and store in a file.'''
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        tax_db['taxonomy_string'].str.split(expand=True)
        pass

    def find_id(self, valid_accession):
        '''Method is used to check if given id exists in the database.'''
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        if valid_accession in list(tax_db["accession_number"]): return True
        else: return False

    def find_tax(self, valid_accession):
        '''Method is used to get taxonomy information of a record in the database.'''
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        accession_ind = tax_db.index[tax_db["accession_number"] == valid_accession].to_list()[0]
        return tax_db.iloc[accession_ind].loc['taxonomy_string']

    def rm_record(self, valid_accession):
        '''Method is used to remove a record from local database based on accession number.'''
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        seq_db = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        tax_db = tax_db[tax_db["accession_number"] != valid_accession]
        del seq_db[valid_accession]
        with open(self.sequence_file, "w+") as seq_db_file:
            SeqIO.write(seq_db.values(), seq_db_file, "fasta")
        tax_db.to_csv(self.taxomony_file, index=False)

    def write_tax(self, valid_accession, new_taxonomy):
        '''Method is used to replace a taxonomy information of a record in the database.'''
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        accession_ind = tax_db.index[tax_db["accession_number"] == valid_accession].to_list()[0]
        tax_db.at[accession_ind,'taxonomy_string']=new_taxonomy
        tax_db.to_csv(self.taxomony_file, index=False)

    def write_id(self, old_accession, new_accession):
        '''Method is used to replace an accession of a record in the database.'''
        tax_db = pd.read_csv(self.taxomony_file, header=[0])
        seq_db = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        accession_ind = tax_db.index[tax_db["accession_number"] == old_accession].to_list()[0]
        tax_db.at[accession_ind,'accession_number']=new_accession
        seq_db[new_accession] = seq_db[old_accession]
        del seq_db[old_accession]
        seq_db[new_accession].id = new_accession
        with open(self.sequence_file, "w+") as seq_db_file:
            SeqIO.write(seq_db.values(), seq_db_file, "fasta")
        tax_db.to_csv(self.taxomony_file, index=False)


    def export_fasta(self, valid_accession):
        #CREATES "datestamp-export-sequences.fasta"
        pass           
    def export_meta(self, valid_accession):
        #CREATES "datestamp-export-metadata.csv"
        pass
    def export_record(self, valid_accession):
        #CREATES "datestamp-export-metadata.csv"
        #CREATES "datestamp-export-sequences.fasta"
        pass


class Query_ncbi:
    #Finished (unit)
    def __init__(self, database, email, alphabet):
        self.database = database
        self.email = email
        self.alphabet = alphabet

    def check_output(self, valid_accession):
        '''Method checks if record exists in NCBI database given accession number.'''
        Entrez.email = self.email
        try:
            handle = Entrez.efetch(db=self.database, id = valid_accession, rettype="acc") #Attempts to open connection to NCBI database
            data = handle.read() #Reads up-to-date accession from NCBI
            handle.close() #Closes connection
            if valid_accession in data: #Checks if found accession contains query
                return True
            else: 
                return False
        except urllib.error.HTTPError as e: #If query is invalid, the exception is thrown
            if str(e) == "HTTP Error 400: Bad Request":
                return False

    def get_fasta(self, valid_accession):
        '''Method attempts to get fasta sequence (header + sequence separately) given accession number.'''
        #generate request
        Entrez.email = self.email
        handle = Entrez.efetch(db=self.database, id = valid_accession, rettype="fasta") #Attempts to open connection to NCBI database
        record = SeqIO.read(handle, "fasta") #Reads sequence in fasta format from ncbi
        handle.close() #Closes connection
        return record.seq, record.id

    def get_taxonomy(self, valid_accession):
        '''Method attempts to get taxonomy information(list) from NCBI database given accession number.'''
        Entrez.email = self.email
        handle = Entrez.efetch(db=self.database, id = valid_accession, rettype="gb") #Attempts to open connection to NCBI database
        tax_data = SeqIO.read(handle, format='genbank').annotations['taxonomy']
        handle.close() #Closes connection
        return tax_data[0:8]


class Logger:
    #Finished (unit)
    def __init__(self, operation_type, valid_accession=None,  taxonomy=None,  old_accession=None,  old_taxonomy=None,  new_accession=None,  new_taxonomy=None, log_file='PD2.log'):
        self.operation_type = operation_type
        self.log_file = log_file
        self.log_dict = {
        "add_fasta":f"New record was added from fasta file (id: {valid_accession}).",
        "add_ncbi":f"New record was added from ncbi database (id: {valid_accession})",
        "rm_record":f"Record was removed (id: {valid_accession})",
        "ch_header":f"ID of the record was changed from {old_accession} to {new_accession}",
        "ch_tax":f"Taxonomy information for {valid_accession} was changed from {old_taxonomy} to {new_taxonomy}"
    }

    def log_change(self):
        #Create log file if it does not exist
        log_file = Path(self.log_file)
        log_file.touch(exist_ok=True)
        
        #Log changes - full documentation link - https://docs.python.org/3/howto/logging.html
        pd2_logger = logging.getLogger(self.operation_type)
        pd2_logger.setLevel(logging.DEBUG)
        log_file_handler = logging.FileHandler(filename=self.log_file)
        log_file_formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S')
        log_file_handler.setFormatter(log_file_formatter)
        pd2_logger.addHandler(log_file_handler)
        pd2_logger.info(self.log_dict[self.operation_type])
        log_file_handler.close() #Close file after logging is complete


class Plotter:
    def __init__(self, name, age):
        self.name = name
        self.age = age
    def display(self, path_to_count_file):
        pass


if __name__ == "__main__":
    pass