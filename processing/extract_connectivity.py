#!/usr/bin/env python
from __future__ import division
import argparse
from sklearn.covariance import ledoit_wolf
from nilearn.input_data import NiftiLabelsMasker, NiftiMasker
from scipy.spatial.distance import squareform, pdist
import numpy as np
import pandas as pd

def compute_fc(X, TYPE):
    """
    Compute different types of FC estimators. 
    """
    if TYPE == 'lw':
        fc_mat, _ = ledoit_wolf(X)
        np.fill_diagonal(fc_mat, 0)
        #fc = squareform(fc_mat)

    if TYPE == 'full':
        fc = 1 - pdist(X.T, metric = 'correlation')
#        fc = np.arctanh(fc)
        fc_mat = squareform(fc)

    return(fc_mat) 

def get_fc(func_mni_filename, atlas_filename, confounds_filename, output_filename, 
    TR, TYPE, FWHM):
    """
    Extract connectivity matrix given atlas and processed fMRI data.
    """
    if FWHM == 0:
        FWHM = None

    confounds = pd.read_csv(confounds_filename, sep='\t')
    global_signal = confounds['GlobalSignal'].values
    FD = confounds['FramewiseDisplacement'].values
#    masker = NiftiLabelsMasker(labels_img = atlas_filename, 
#                     smoothing_fwhm = FWHM,
#                     t_r = TR, 
#                     memory_level = 1, 
#                     memory = 'nilearn_cache')

#    X = masker.fit_transform(func_mni_filename, confounds = global_signal)
#    nvols = X.shape[0]
#    fc = compute_fc(X, TYPE)
#    fc = np.arctanh(fc)

#    np.savetxt(output_filename, fc, delimiter=" ")

    # print average and max FD
    print("{};{};{};{}".format(np.mean(FD[1:]), np.max(FD[1:]), np.sum(FD[1:] > 0.2)/FD.size, np.sum(FD[1:] > 0.3)/FD.size ))

def main():

    parser = argparse.ArgumentParser(description='FC extractor.')
    
    parser.add_argument('func_file',
                   help='Functional file.')

    parser.add_argument('atlas_file',
                   help='Atlas file.')

    parser.add_argument('conf_file',
                   help='Confounds file.')

    parser.add_argument('output_file',
                   help='Output file.')

    parser.add_argument('TR',
                   help='TR')

    parser.add_argument('TYPE',
                   help='TYPE')
                   
    parser.add_argument('FWHM',
                   help='FWHM') 
                    
    args = parser.parse_args()
    get_fc(args.func_file, args.atlas_file, 
    args.conf_file, args.output_file, float(args.TR), args.TYPE, float(args.FWHM))

if __name__ == "__main__":
    main()


