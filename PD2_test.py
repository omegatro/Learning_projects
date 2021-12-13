import unittest, os, pandas as pd
from datetime import datetime
from Bio import SeqIO
from PD2 import Interface, Logger, Query_ncbi, Database

class TestInterfaceMethods(unittest.TestCase):
    global my_interface 
    my_interface = Interface("ACGT","csv", r"^[A-z]{2}_[A-Z0-9]*$")
    def test_Interface_check_format(self):
        valid_accession_list = [
        "AC_1234567", "NC_1234567", "NG_1234567", "NT_1234567",
        "NW_1234567", "NZ_1234567", "NM_1234567", "NR_1234567",
        "XM_1234567", "XR_1234567", "AP_1234567", "NP_1234567",
        "YP_1234567", "XP_1234567", "WP_1234567", "NZ_OADK00000000"
        ]
        for ac_number in valid_accession_list:
            self.assertTrue(my_interface.check_format(ac_number))
        invalid_accession = "XYX_124#5"
        self.assertFalse(my_interface.check_format(invalid_accession))

    def test_Interface_check_file_type(self):
        self.assertTrue(my_interface.check_file_type("file.csv"))
        self.assertFalse(my_interface.check_file_type("file.tsv"))

    def test_Interface_check_alphabet(self):
        my_interface = Interface("ACGT","csv", r"^[A-z]{2}_[0-9]*$")
        self.assertTrue(my_interface.check_alphabet("ACGTTTTGCCATGGTAC"))
        self.assertFalse(my_interface.check_alphabet("RGHATGCNTKLFR"))


class TestLoggerMethods(unittest.TestCase):
    def test_Logger_log_change_convention(self):
        op_types = ['add_fasta',"add_ncbi","rm_record","ch_header","ch_tax"]
        for op_type in op_types:    
            my_logger = Logger(op_type, log_file="PD2_test.log", valid_accession="NC_157326", old_accession="NL_1456", new_accession="NL_14561", old_taxonomy="old", new_taxonomy="new")
            my_logger.log_change()
            with open(f"PD2_test.log", "r+") as log_file:
                self.assertEqual(log_file.readlines()[-1].strip(), f"{datetime.now().strftime('%Y-%m-%d %I:%M:%S')} {my_logger.log_dict[op_type]}")
            os.remove("PD2_test.log")
        

class TestQueryMethods(unittest.TestCase):
    global my_query
    my_query = Query_ncbi('nucleotide', "jb17048@edu.lu.lv","ACGT")
    def test_Query_check_output(self):
        self.assertTrue(my_query.check_output("NC_045512"))
        self.assertFalse(my_query.check_output("ABCDEFGH"))

    def test_Query_get_fasta(self):
        valid_output = my_query.get_fasta("NC_045512")
        self.assertRegex(str(valid_output[0]), f'[{my_query.alphabet}]*')
        id_check = "NC_045512" in str(valid_output[1])
        self.assertTrue(id_check)

    def test_Query_get_taxonomy(self):
        bact_tax = ['Bacteria', 'Proteobacteria', 'Betaproteobacteria', 'Neisseriales', 'Neisseriaceae', 'Neisseria']
        self.assertEqual(my_query.get_taxonomy("NZ_OADK00000000"), bact_tax)


class TestDatabaseMethods(unittest.TestCase):
    global my_database
    my_database = Database("test_database")       

    def test_Database_create_db_files(self):
        my_database.create_db_files()
        seq_file_test = os.path.exists(my_database.sequence_file) and os.path.isfile(my_database.sequence_file)
        tax_file_test = os.path.exists(my_database.taxomony_file) and os.path.isfile(my_database.taxomony_file)
        self.assertTrue(seq_file_test and tax_file_test)
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxomony_file)


    def test_Database_add_record(self):
        my_database.create_db_files()
        my_database.add_record("test_header","ACGT", "test_taxonomy_string","test_description")
        test_id = "test_header"
        test_description = "test_description"
        test_sequence = "ACGT"
        test_tax_record = ["test_header","test_taxonomy_string"]
        seq_record = SeqIO.read(my_database.sequence_file, format="fasta")
        tax_record = pd.read_csv(my_database.taxomony_file, header=[0])
        self.assertEqual(seq_record.seq, test_sequence)
        self.assertEqual(seq_record.id, test_id)
        self.assertEqual(seq_record.description.split(" ")[-1], test_description)
        self.assertEqual(tax_record["accession_number"][0], test_tax_record[0])
        self.assertEqual(tax_record["taxonomy_string"][0], test_tax_record[1])
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxomony_file)
        
    
    def test_Database_find_id(self):
        my_database.create_db_files()
        my_database.add_record("test_header","ACGT", "test_taxonomy_string","test_description")
        self.assertTrue(my_database.find_id("test_header"))
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxomony_file)

    def test_Database_rm_record(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        my_database.add_record("test_header_1","ACGT_1", "test_taxonomy_string_1","test_description_1")
        my_database.rm_record("test_header")
        test_tax_record = ["test_header_1","test_taxonomy_string_1"]
        with open(my_database.sequence_file, "r+") as seq_db:
            seq_record = SeqIO.read(my_database.sequence_file, format="fasta")
            self.assertEqual(seq_record.id, "test_header_1")
            self.assertEqual(seq_record.seq, "ACGT_1")
            self.assertEqual(seq_record.description.split(" ")[-1], "test_description_1")
        tax_db = pd.read_csv(my_database.taxomony_file, header=[0])
        self.assertEqual(tax_db["accession_number"][0], test_tax_record[0])
        self.assertEqual(tax_db["taxonomy_string"][0], test_tax_record[1])
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxomony_file)

if __name__ == '__main__':
    unittest.main()
    