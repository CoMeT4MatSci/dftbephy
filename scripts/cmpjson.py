import sys
import numpy as np
import json

def compare_approximate(first, second):
    """Return whether two dicts of arrays are roughly equal"""
    assert(first.keys() == second.keys())
    
    for key in first:
        try:
            np.testing.assert_almost_equal(first[key], second[key], decimal=5, err_msg='Issue in field %s: ' % key)
        except TypeError:
            assert(first[key] == second[key])


# check if we got two command line arguments
if len(sys.argv) != 3:
    print('call with two json files for comparison: python cmpjson.py file1.json file2.json')
    exit()

# open and read json files
with open(sys.argv[1]) as jfile:
    data = json.load(jfile)
data = json.loads(data)

with open(sys.argv[2]) as jfile:
    data2 = json.load(jfile)
data2 = json.loads(data2)

# do the comparison and throw an exception
# if it fails
compare_approximate(data, data2)
