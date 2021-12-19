# Reference sequence database management tool
- A script is designed to create and manage a local database of reference sequences in fasta format, along with annotation table of taxonomic data for each sequence and visualization tool.
- The script will generate 1 fasta file to store genome sequences (can get large), 1 csv file to store taxonomic information for each sequence and 4 json files, storing information required to visualize the current content of the database.
- Additionally the PD2.log file will be created, where changes to the database files (fasta & csv) will be logged.
- The files will be created in the same directory where the script is stored.
## Options
- --add-fasta: Add new record(s) from locally stored fasta file & metadata file of csv format Expected: file_name (same file name for csv & fasta)
- --add-ncbi: Add new record from NCBI RefSeq database by correct accession number Expected: valid_accession_number
- --add-ncbi-list: Add new record(s) from NCBI RefSeq database based on a list of correct accession numbers provided as a file of csv format (single column with one header row) Expected: path_to_csv
- --exp-fasta: Create a fasta file containing sequences from local database based on a list of correct accession numbers provided as a file of csv format Expected: path_to_csv
- --exp-meta: Create a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of csv format Expected: path_to_csv
- --exp-records: Create a fasta file containing sequences & a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of csv format Expected: path_to_csv
- --rm-record: Remove records from local database based on a list of correct accession numbers provided as a file of csv format Expected: path_to_csv
- --ch-header: Replace an accession number that exists in local database with user-provided accession number if provided accession number was found in local database Expected: valid_accession_from_local_dbnew_accession
- --ch-tax: Replace a taxonomy string that exists in local database with user-provided taxonomy string if provided accession number was found in local database Expected: valid_accession_from_local_dbnew_taxonomy_string
- --view-data: Visualize the contents of the local database as a tree chart, showing the number of records that belong to each taxonomic group
## Visualization example
![The image will be here shortly](https://github.com/omegatro/UNPG/blob/datz5032_final/newplot.jpg?raw=true)
