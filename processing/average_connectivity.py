#!/usr/bin/env python
import argparse
import numpy as np

def average_fc(file_list, output_filename):
    """
    Average connectivity matrices.
    """
    all_fc = list()
    filenames = file_list.split(" ")
    for filename in filenames:
    	fc = np.loadtxt(filename, delimiter=" ")
 	all_fc.append(fc)

    fc = np.dstack(all_fc)
    fc_mean = np.mean(fc, axis = 2)
    np.savetxt(output_filename, fc_mean, delimiter=" ")

def main():

    parser = argparse.ArgumentParser(description='FC average.')
    
    parser.add_argument('output_file',
                   help='Output file.')

    parser.add_argument('file_list',
                   help='FC files.')

    args = parser.parse_args()
    average_fc(args.file_list.strip(), args.output_file)

if __name__ == "__main__":
    main()


