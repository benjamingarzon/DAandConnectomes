#!/usr/bin/env python
"""
usage: simplify_samples.py input_file output_file param1 param2 ...
Removes unnecessary parameters from stan sample file for easier handling.

"""
# Author: Benjamin Garzon <benjamin.garzon@gmail.com>
# License: BSD 3 clause

import string
import sys
import argparse
import re
import numpy as np

def simplify_samples(args):

    n_pars = len(args.params)

    f_in = open(args.input_file, 'r')
    f_out = open(args.output_file, 'w')

    block = 0	# 0 : first comments, 1 : warmup values, 2 : non-warmup values
    l = 0
    for line in f_in:
        l = l + 1 
#        print("Reading line %d of %s."%(l, args.input_file))

	if block == 0: 	

	    if line[0] == "#": 	
#                print("Commented line: %s"%line[:-1])
                f_out.write(line)
	
	    else:
                block += 1
		header = line
                params = string.split(header, sep = ',')
                indices = np.array([i for i, param in enumerate(params) if all([re.search(r"\b%s"%myparam, param)==None for myparam in args.params]) ])
                new_header = ",".join( params[i] for i in indices )
                f_out.write(new_header)
#                f_out.write("\n")

        elif ( line[0] != "#" ) and ( len(line) > 1 ): 	

            values = np.array(string.split(line, sep=","))
            f_out.write(",".join(values[indices]))
#            f_out.write("\n")

	else:
#            print("Commented line: %s"%line[:-1])
            f_out.write(line)

    print("Total number of parameters is %d."%len(params))
    print("Reduced number of parameters is %d."%len(indices))
    
    f_in.close()
    f_out.close()

def main():

    parser = argparse.ArgumentParser(description='simplify_samples.py: \
    remove parameters from stan samples file.')
    
    parser.add_argument('input_file',
                   help='Input File.')
                   
    parser.add_argument('output_file',
                   help='Output File.')

    parser.add_argument('params', type=str, nargs='+',
                    help='Parameter names')           
        
    simplify_samples(parser.parse_args())

if __name__ == "__main__":
    main()	  
