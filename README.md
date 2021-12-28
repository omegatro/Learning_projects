# Reference sequence database management tool
- A script is designed to create and manage a local database of reference sequences in fasta format.
- Taxonomic information for each sequence is stored in the metadata table.
- Visualization function allows to plot the contents of the database using the taxonomic information stored in the metadata table.
- The script will generate:
    - 1 fasta file to store genome sequences 
      - can take up much space, depending on the size of stored sequences
      - not more than 300Mb for test cases
    - 1 csv file to store taxonomic information for each sequence 
    - 4 json files, storing information required to visualize the current content of the database
    - PD2.log file to save information about any changes to the database files (fasta & csv)
- The files will be created in the same directory where the script is stored.
- Change the email in line 432 of the script to valid email if using with ncbi options.
## Special Python packages (tested on Python 3.7.7)
```
pip install biopython plotly igraph pandas
```
## Options
- --add-fasta: Add new record(s) from locally stored fasta file & metadata file of csv format 
  - Expected: file_name (same file name for csv & fasta)
- --add-ncbi: Add new record from NCBI RefSeq database by correct accession number 
  - Expected: valid_accession_number
- --add-ncbi-list: Add new record(s) from NCBI RefSeq database based on a list of correct accession numbers provided as a file of csv format (single column with one header row)
  - Expected: path_to_csv
- --exp-fasta: Create a fasta file containing sequences from local database based on a list of correct accession numbers provided as a file of csv format 
  - Expected: path_to_csv
- --exp-meta: Create a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of csv format
  - Expected: path_to_csv
- --exp-records: Create a fasta file containing sequences & a csv file containing taxonomic information from local database
     based on a list of correct accession numbers provided as a file of csv format 
  - Expected: path_to_csv
- --rm-record: Remove records from local database if given correct accession number. 
  - Expected: valid_accession_number_from_local_db
- --ch-header: Replace an accession number that exists in local database with user-provided accession number if provided accession number was found in local database 
  - Expected: valid_accession_from_local_db,new_accession
- --ch-tax: Replace a taxonomy string that exists in local database with user-provided taxonomy string if provided accession number was found in local database 
  - Expected: valid_accession_from_local_dbnew_taxonomy_string
- --view-data: Visualize the contents of the local database as a tree chart, showing the number of records that belong to each taxonomic group
## Visualization example
![The image will be here shortly](https://github.com/omegatro/UNPG/blob/datz5032_final/newplot.jpg?raw=true)

## Examples
- Before you start, move all files you intend to work with to the directory where PD2.py script is locates, e.g.
```
cp test_files/* ./
```
- Creating a test database by downloading sequences 
from https://www.ncbi.nlm.nih.gov/refseq/ and viewing the result
```
python PD2.py --add-ncbi-list ncbi_headers.csv
python PD2.py --view-data
```
- Adding sequence from local file
```
python PD2.py --add-fasta sequence
python PD2.py --view-data
```
- Adding sequence from ncbi by [accession number](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#:~:text=An%20accession%20number%20applies%20to,the%20type%20of%20sequence%20record.)
```
python PD2.py --add-ncbi NZ_CP033718
python PD2.py --view-data
```
- Changing accession number for a record
```
python PD2.py --ch-header NZ_CP033718,NZ_CP033719
python PD2.py --view-data
```
- Changing taxonomy information for a record
```
python PD2.py --ch-tax NZ_CP033719,Bacteria|Actinobacteria|Corynebacteriales|Mycobacteriaceae|Mycobacteroides|Mycobacteroides abscessus
python PD2.py --view-data
```
- Removing a record from the database
```
python PD2.py --rm-record NZ_CP033719
python PD2.py --view-data
```
- Exporting: sequences, taxonomy information or both (output files will be generated in the folder where PD2.py is located)
```
python PD2.py --exp-fasta sequence.csv
python PD2.py --exp-meta sequence.csv
python PD2.py --exp-records sequence.csv
```
- Big test - 110 sequences will be downloaded from ncbi database - around 275Mb - may take some time to download (~ 20 min depending on the connection speed).
```
python PD2.py --add-ncbi-l big_test.csv
```
