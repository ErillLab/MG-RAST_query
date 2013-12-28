# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 18:07:46 2013

@author: Talmo

This script will query the MG-RAST database and obtain metagenomes that should
be assembled.

TODO:
- Query the database using the API to confirm metagenomes are assemblies.

"""

import csv
import numpy as np
import requests
from prettytable import PrettyTable

######################
##### Parameters #####
######################

data_table = "web_table.tsv"

min_avg_seq_length = 1000
max_avg_seq_length = 20000
min_bps = 1e6

######################
print "Parameters:"
params = PrettyTable()
params.add_row(["min_avg_seq_length", min_avg_seq_length])
params.add_row(["max_avg_seq_length", max_avg_seq_length])
params.add_row(["min_bps", min_bps])
print params.get_string(header=False)
print

# CLean up data table before we read it (get rid of double tabs)
with open(data_table, "r") as f:
    _table = f.read()
    _table = _table.replace("\t\t", "\t")
    _table = _table.replace("&nbsp;", "_")
with open(data_table, "w") as f:
    f.write(_table)


# Load data from table
all_mgs = []
with open(data_table, "r") as tsv:
    tsv_table = csv.reader(tsv, delimiter="\t")
    
    # Read header row
    table_header = tsv_table.next()
    
    # Parse every row into dictionary
    for row in tsv_table:
        all_mgs.append({key:value for key, value in zip(table_header, row)})

# Convert numerical columns to appropriate data types
for mg in all_mgs:
    mg["bps"] = int(mg["bps"])
    mg["sequences"] = int(mg["sequences"])
    mg["avg_seq_length"] = float(mg["avg_seq_length"])

print "Loaded %d metagenomes from %s." % (len(all_mgs), data_table)
print


# Filter for metagenomes meeting criteria
assemblies = []
for mg in all_mgs:
    if mg["avg_seq_length"] >= min_avg_seq_length and mg["avg_seq_length"] <= max_avg_seq_length and mg["bps"] >= min_bps:
        assemblies.append(mg)

print "Detected %d assemblies based on parameters." % len(assemblies)

# Make a table from the assemblies list
assemblies_table = PrettyTable(table_header)
for mg in assemblies:
    assemblies_table.add_row([mg[field] for field in table_header])
print assemblies_table.get_string(fields=["id", "avg_seq_length", "bps", "sequencing_method", "sequencing_type"], sortby="avg_seq_length", reversesort=True)
print

# Calculate statistics on metagenomes
print "Statistics:"
assemblies_stats = PrettyTable(["Field", "Min", "Mean", "Max"])
for field in ["bps", "sequences", "avg_seq_length"]:
    col = np.array([mg[field] for mg in assemblies])
    assemblies_stats.add_row([field, min(col), np.mean(col), max(col)])
print assemblies_stats
print

# Display count per method
print "Metagenomes per sequencing_method:"
assemblies_methods = [mg["sequencing_method"] for mg in assemblies]
all_mgs_methods = [mg["sequencing_method"] for mg in all_mgs]
methods = set(all_mgs_methods)
methods_table = PrettyTable(["Method", "Assemblies", "All"])
for method in methods:
    methods_table.add_row([method, assemblies_methods.count(method), all_mgs_methods.count(method)])
methods_table.add_row(["Total", len(assemblies), len(all_mgs)])
print methods_table

# To download each MG:
# http://api.metagenomics.anl.gov/download/mgm4544224.3?file=050.2
# 050.2 = original sequence, gzipped
# >>> mg = requests.get("http://api.metagenomics.anl.gov/download/mgm4544224.3?file=050.2")
# >>> mg.headers
# CaseInsensitiveDict({'content-length': '5509921', 'content-disposition': 'attachment;filename=111857.fna.gz', 'expires': 'Sat, 28 Dec 2013 01:16:15 GMT', 'keep-alive': 'timeout=90', 'server': 'nginx', 'connection': 'keep-alive', 'cache-control': 'max-age=0', 'date': 'Sat, 28 Dec 2013 01:16:15 GMT', 'access-control-allow-origin': '*', 'content-type': 'application/x-download'})
# >>> len(mg.content)
# 5509921
