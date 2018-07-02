#!/usr/bin/env python
"""
usage: reduce_samples.py input_file output_file param1 param2 ...
Removes unnecessary parameters from stan sample file for easier handling.

"""
# Author: Benjamin Garzon <benjamin.garzon@gmail.com>
# License: BSD 3 clause

import string
import sys
import argparse
import re
import numpy as np

def reduce_samples(args):

    n_pars = len(args.params)

    f_in = open(args.input_file, 'r')
    f_out = open(args.output_file, 'w')
    mcmc = args.mcmc

    block = 0	# 0 : first comments, 1 : warmup values, 2 : non-warmup values
    l = 0
    for line in f_in:
        l = l + 1 
        print l, block
        print("Reading line %d of %s."%(l, args.input_file))

	if block == 0: 	

	    if line[0] == "#": 	
	         print("Skipping line: %s"%line[:-1])
		
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
            if mcmc and block == 2:
    	        f_out.write(",".join(values[indices]))
#                f_out.write("\n")
            elif not mcmc:
    	        f_out.write(",".join(values[indices]))
#                f_out.write("\n")

	else:
             print("Skipping line: %s"%line[:-1])	

             if re.search("terminated", line)!=None:
                 block += 1		
                   

    print("Total number of parameters is %d."%len(params))
    print("Reduced number of parameters is %d."%len(indices))
    
    f_in.close()
    f_out.close()

def main():

    parser = argparse.ArgumentParser(description='reduce_samples.py: \
    remove parameters from stan samples file.')
    
    parser.add_argument('input_file',
                   help='Input File.')
                   
    parser.add_argument('output_file',
                   help='Output File.')

    parser.add_argument('-mcmc',
                   help='The input file containt MCMC samples (default=FALSE)', action='store_true')

    parser.add_argument('params', type=str, nargs='+',
                    help='Parameter names')           
        
    reduce_samples(parser.parse_args())

if __name__ == "__main__":
    main()	  
