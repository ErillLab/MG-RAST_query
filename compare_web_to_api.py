# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 11:08:03 2013

@author: Talmo


MG-RAST has two interfaces for searching:
    1. The main website, which has a huge table that can be filtered.
        URL: http://metagenomics.anl.gov/?page=MetagenomeSelect
    2. A REST API that returns JSON-encoded data.
        URL: http://api.metagenomics.anl.gov/1/api.html

They may not be outputting data from the same backend, since their total counts
do not match. This script compares the output from the two.

To get the web table, go to the main website (URL above), click the
"clear table filters" at the top of the table and display all of the items.
At the time the table was last retrieved, there were 14837 total items.

NOTE: The web table actually has ALL of the data in the HTML mark-up so the
page is actually ~9.1 MB. The "filtering" and TSV generation is done via
JavaScript.
"""
import csv
import requests
from pprint import pprint
from time import time

# The base API URL. Use the api() function to build a quick query URL, ex:
#   >>> print api("metagenome?limit=20&order=name")
#   http://api.metagenomics.anl.gov/1/metagenome?limit=20&order=name
api_url = "http://api.metagenomics.anl.gov/"
api = lambda action: api_url + action

# Initialize stuff for requests
api_ids = []
offset = 0
max_items_per_query = 10000

# Python's ghetto do-while loop
while True:
    # Retrieves all the entries with minimal verbosity, no filtering
    query = "metagenome?verbosity=minimal&limit=%d&order=id&direction=asc&match=all&status=both&offset=%d" % (max_items_per_query, offset)
    
    # Send the request and download the response (may take a while)
    s = time()
    print "Requesting items %d-%d..." % (offset + 1, offset + 1 + max_items_per_query),
    r = requests.get(api(query))
    print "Done. [%.2fs]" % (time() - s)
    
    # Save IDs
    for data in r.json()["data"]:
        api_ids.append(str(data["id"][3:]))
    
    # Stop once we've fetched all the items
    if offset + max_items_per_query >= r.json()["total_count"]:
        break
    
    # Iterate to next chunk of items
    offset += max_items_per_query

print "IDs found from API:", len(api_ids)


# Load the IDs from the downloaded table
web_ids = []
with open('web_table.tsv','rb') as tsv:
    web_table = csv.reader(tsv, delimiter='\t')
    web_header = web_table.next() # header row
    for row in web_table:
        web_ids.append(row[0])
print "IDs loaded from table:", len(web_ids)

# Cross-check the two:
not_in_api = [id for id in web_ids if id not in api_ids]
not_in_web = [id for id in api_ids if id not in web_ids]

print "Entries not in API:", len(not_in_api)
print "Entries not in web:", len(not_in_web)

# Double-check by looking up missing IDs in the API:
not_in_api_check = []
for id in not_in_api[:10]: # first 10 anyway (this takes a while)
    not_in_api_check.append("ERROR" in requests.get(api("metagenome/" + id)).json())
    print id, "-", "Not found" if not_in_api_check[-1] else "Found"

