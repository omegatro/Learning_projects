import unittest, os, pandas as pd
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
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
        Path('file.csv').touch(exist_ok=True)
        Path('file.fasta').touch(exist_ok=True)
        self.assertTrue(my_interface.check_file_type("file"))
        os.remove('file.csv')
        os.remove('file.fasta')

    def test_Interface_check_alphabet(self):
        my_interface = Interface("ACGT","csv", r"^[A-z]{2}_[0-9]*$")
        self.assertTrue(my_interface.check_alphabet("ACGTTTTGCCATGGTAC"))
        self.assertFalse(my_interface.check_alphabet("RGHATGCNTKLFR"))
    
    def test_Interface_read_local(self):
        test_output_df = pd.DataFrame.from_dict({'accession_number':['NZ_CP021325'],'taxonomy_string':['Bacteria|Firmicutes|Bacilli|Bacillales|Listeriaceae|Listeria']})
        test_output_fasta = SeqIO.SeqRecord(Seq('GTGTTTGGATAACCTTATCCATAGCTTTTTCTATCTGTGGATAACTTTATAGCATCCATTTACATTACATAAAAAGGGGGGGTACTAGTGCAATCAATTGAAGACATCTGGCAGGAAACA'),id='NZ_CP021325', description='test_description')
        test_output_df.to_csv('test_file.csv', index=False, header=True)
        with open('test_file.fasta', "a+") as seq_db:
            SeqIO.write(test_output_fasta, seq_db, "fasta")
        output = my_interface.read_local("test_file")
        self.assertTrue(output[0]['NZ_CP021325'].id == test_output_fasta.id)
        self.assertTrue(output[0]['NZ_CP021325'].seq == test_output_fasta.seq)
        self.assertTrue(output[0]['NZ_CP021325'].description == 'NZ_CP021325 test_description')
        self.assertTrue(output[1].equals(test_output_df))
        os.remove('test_file.csv')
        os.remove('test_file.fasta')

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
        tax_file_test = os.path.exists(my_database.taxonomy_file) and os.path.isfile(my_database.taxonomy_file)
        self.assertTrue(seq_file_test and tax_file_test)
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_add_record(self):
        my_database.create_db_files()
        my_database.add_record("test_header","ACGT", "test_taxonomy_string","test_description")
        test_id = "test_header"
        test_description = "test_description"
        test_sequence = "ACGT"
        test_tax_record = ["test_header","test_taxonomy_string"]
        seq_record = SeqIO.read(my_database.sequence_file, format="fasta")
        tax_record = pd.read_csv(my_database.taxonomy_file, header=[0])
        self.assertEqual(seq_record.seq, test_sequence)
        self.assertEqual(seq_record.id, test_id)
        self.assertEqual(seq_record.description.split(" ")[-1], test_description)
        self.assertEqual(tax_record["accession_number"][0], test_tax_record[0])
        self.assertEqual(tax_record["taxonomy_string"][0], test_tax_record[1])
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)
        
    def test_Database_find_id(self):
        my_database.create_db_files()
        my_database.add_record("test_header","ACGT", "test_taxonomy_string","test_description")
        self.assertTrue(my_database.find_id("test_header"))
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_find_tax(self):
        my_database.create_db_files()
        my_database.add_record("test_header","ACGT","test_taxonomy_string","test_description")
        self.assertEqual(my_database.find_tax("test_header"), "test_taxonomy_string")
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_rm_record(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        my_database.add_record("test_header_1","ACGT_1", "test_taxonomy_string_1","test_description_1")
        my_database.rm_record("test_header")
        test_tax_record = ["test_header_1","test_taxonomy_string_1"]
        with open(my_database.sequence_file, "r+") as seq_db:
            seq_record = SeqIO.read(seq_db, format="fasta")
            self.assertEqual(seq_record.id, "test_header_1")
            self.assertEqual(seq_record.seq, "ACGT_1")
            self.assertEqual(seq_record.description.split(" ")[-1], "test_description_1")
        tax_db = pd.read_csv(my_database.taxonomy_file, header=[0])
        self.assertEqual(tax_db["accession_number"][0], test_tax_record[0])
        self.assertEqual(tax_db["taxonomy_string"][0], test_tax_record[1])
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_write_tax(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        my_database.write_tax("test_header", "replaced_test_taxonomy_string")
        test_tax_record = ["test_header","replaced_test_taxonomy_string"]
        tax_db = pd.read_csv(my_database.taxonomy_file, header=[0])
        self.assertEqual(tax_db["accession_number"][0], test_tax_record[0])
        self.assertEqual(tax_db["taxonomy_string"][0], test_tax_record[1])
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_write_id(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        my_database.write_id("test_header", "replaced_test_header")
        test_tax_record = "replaced_test_header"
        with open(my_database.sequence_file, "r+") as seq_db:
            seq_record = SeqIO.read(seq_db, format="fasta")
            self.assertEqual(seq_record.id, "replaced_test_header")
        tax_db = pd.read_csv(my_database.taxonomy_file, header=[0])
        self.assertEqual(tax_db["accession_number"][0], test_tax_record)
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_export_fasta(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        exported_seq = my_database.export_fasta("test_header")
        self.assertEqual(exported_seq.id, "test_header")
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_export_meta(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        exported_tax = my_database.export_tax("test_header")
        self.assertEqual(exported_tax["taxonomy_string"], "test_taxonomy_string")
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

    def test_Database_export_record(self):
        my_database.create_db_files()
        my_database.add_record("test_header", "ACGT", "test_taxonomy_string","test_description")
        export_record = my_database.export_record("test_header")
        self.assertEqual(export_record[1].id, "test_header")
        self.assertEqual(export_record[0]["taxonomy_string"], "test_taxonomy_string")
        os.remove(my_database.sequence_file)
        os.remove(my_database.taxonomy_file)

if __name__ == '__main__':
    unittest.main()   