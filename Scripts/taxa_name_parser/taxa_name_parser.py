# taxa name parser script
# Run the script with "$ python taxa_name_parser.py [name_of_archive]"
from zipfile import ZipFile # To unzip archive
from sys import argv # To get user input
import csv # To handle tab separated value files
import subprocess # To access shell scripting
import timeit # To time processes

###Functions:
def get_vernacular_dict(file): # A function that takes a tab-separated file and returns it as a list of dictionaries
    """
        A function that takes a csv-file ('VernacularName.csv' specifically) and returns it as a list of dictionaries.
    """
    with open(file, encoding='utf-8-sig') as file_tsv:
        file_dictreader = csv.DictReader(file_tsv, delimiter='\t', quoting=csv.QUOTE_NONE) # Convert tab separated file into a DictReader-object
        list_of_dicts = list([]) # Initiate list that will contain the dictionaries
        for row in file_dictreader: # For every dictionary in the DictReader object
            if row['isPreferredName'] == 'true' and row['language'] == 'sv': # If the dictionary is relevant
                list_of_dicts.append(row) # Add the dictionary to the list of dictionaries
        return list_of_dicts # Return list of dictionaries

def get_taxa_dict(file):
    """
        A function that takes a csv-file ('Taxon.csv' specifically) and returns it as a list of dictionaries.
    """
    with open(file, encoding='utf-8-sig') as file_tsv:
        file_dictreader = csv.DictReader(file_tsv, delimiter='\t', quoting=csv.QUOTE_NONE) # Convert tab separated file into a DictReader-object
        list_of_dicts = list([]) # Initiate list that will contain the dictionaries
        for row in file_dictreader: # For every dictionary in the DictReader object
            if row['taxonRank'] == 'species' and row['taxonomicStatus'] == 'accepted' and row['nomenclaturalStatus'] == 'valid': # If the dictionary is relevant
                list_of_dicts.append(row) # Add the dictionary to the list of dictionaries
        return list_of_dicts # Return list of dictionaries

###Main:
input = argv[1] if len(argv) > 1 else None # Check if the user supplied a file to parse
if input == None: # If user did not supply a file
    print("Please supply a zip-file to parse.") # Inform user that no file was given
    quit() # Terminate program

print("Unzipping...")
zipped_file = ZipFile(input, 'r') # Unzip user supplied file
zipped_file.extract('VernacularName.csv',path='./tmp_to_be_removed') # Extract file with vernacular names to a temporary location
zipped_file.extract('Taxon.csv',path='./tmp_to_be_removed') # Extract file with latin names to a temporary location

print("Extracting...")
verna_dict_list = get_vernacular_dict('./tmp_to_be_removed/VernacularName.csv') # Get list of dictionaries from tsv-file
taxa_dict_list = get_taxa_dict('./tmp_to_be_removed/Taxon.csv') # Get list of dictionaries from tsv-file

print("Matching latin names with vernacular names...")
starttime = timeit.default_timer() # Start timer to time how long it takes to match all names
with open("latin_vernacular.csv", "w", encoding='utf-8-sig') as output: # Initiate output-file
    output.write("Latin,Swedish\n")
    for item in taxa_dict_list: # For every dictionary with latin names
        for thing in verna_dict_list: # For every dictionary with vernacular names
            if item['acceptedNameUsageID'] == thing['taxonId']: # If the ID in both dictionaries match
                latin = item["scientificName"] # Get latin name
                swedish = thing["vernacularName"] # Get vernacular name
                output.write(latin + "," + swedish + "\n") # Write the names in the output file csv-style
                print("Matching: " + latin + (" "*(80-len(latin))) , end='\r') # Print which latin name is currently being matched
print("Matching: Done." + (" "*80))
stoptime = timeit.default_timer() # Stop timer
print("Matching names took a total of: " + str(round(stoptime-starttime, 2)) + " seconds.") # Print how many seconds it took to match names

print("Cleaning up...")
subprocess.run(["rm -r tmp_to_be_removed/"], shell=True) # Remove temporary directory containing unzipped files
print("Done.")
