#!/home/benjamin.garzon/Software/anaconda2/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import imsave
from PIL import Image, ImageFont, ImageDraw
import os

#widthBW = 3430 # 8.7cm x 600 dpi
width = 3425 # 8.7cm x 1000 dpi
width_h = 7007 # 17.8cm x 1000 dpi

font = ImageFont.truetype("/usr/share/fonts/liberation/LiberationSerif-Bold.ttf", size = 150)

letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]

def combine_figs(fig_list, nrows, ncols, outputname, width = 0, letters = letters, with_letter = True):

    rows = []
    for row in range(nrows):
        cols = []
        for col in range(ncols):
            newim = Image.open(fig_list[row*ncols + col]).convert('RGB') 
            draw = ImageDraw.Draw(newim) 
            if with_letter:
                 print letters[row*ncols + col], fig_list[row*ncols + col]
                 draw.text((20, 10), letters[row*ncols + col], (0,0,0), font=font)
            del draw
#            print(np.array(newim).shape)
            cols.append(np.array(newim)[:, :, :3])      

        rows.append(np.concatenate(tuple(cols), axis=1))

    im = np.concatenate(tuple(rows), axis=0)
    im = Image.fromarray(im)

    if width !=0:
        wpercent = (width/float(im.size[0]))
        height = int((float(im.size[1])*float(wpercent)))
        im = im.resize((width, height), Image.ANTIALIAS)
        print("Resizing image")
        print(im.size)
#    fig = plt.imshow(im)
#    fig.axes.get_xaxis().set_visible(False)
#    fig.axes.get_yaxis().set_visible(False)
#    plt.figure(facecolor='white')
    #plt.show()
    imsave(outputname, im)

# definitions here

ext='.tif'

# -------------------- new figure
fig_list_1 = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/Figure1a_Cartoon_HD.tif'
]

outputname_1 = '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/Figure1' + ext
combine_figs(fig_list_1, 1, 1, outputname_1, width = width_h, with_letter = False)

# -------------------- new figure
fig_list_2 = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/FigureSimilarity_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/FigureMDS_GNG.png'
#  '/home/benjamin.garzon/Data/DAD/analysis_variability/Figure1d_Graph_Model_HD.tif'
]

outputname_2 = '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/Figure2' + ext
combine_figs(fig_list_2, 1, 2, outputname_2, width = width_h)


# -------------------- new figure

fig_list_3a = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_WM_scores.png'
]

fig_list_3b = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_WM_BetaMu_GNG_modules20.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_WM_BetaMu_GNG_modules70.png'
]

outputname_3 = '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/Figure3' + ext
combine_figs(fig_list_3a + fig_list_3b, 1, 3, outputname_3, width = width_h, letters = ["A", "B", "C"])

# -------------------- new figure
fig_list_4 = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_FxBP_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_FxVBM_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_FxBP_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_FxVBM_GNG.png'
]

outputname_4 = '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/Figure4' + ext
combine_figs(fig_list_4, 2, 2, outputname_4, width = width_h)


# -------------------- new figure
fig_list_s1 = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_MNI_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_MNI_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_MNI_RS.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_MNI_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_MNI_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_MNI_RS.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_BM20_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_BM20_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaMu_BM20_RS.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_BM20_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_BM20_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_BetaSigma_BM20_RS.png',
]

outputname_s1 = '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/FigureS1' + ext
combine_figs(fig_list_s1, 4, 3, outputname_s1, width = width_h)

# -------------------- new figure
fig_list_s2 = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaMu_MNI_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaMu_MNI_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaMu_MNI_RS.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaSigma_MNI_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaSigma_MNI_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaSigma_MNI_RS.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaMu_BM20_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaMu_BM20_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaMu_BM20_RS.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaSigma_BM20_GNG.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaSigma_BM20_TAB.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsGSR/Figure_BetaSigma_BM20_RS.png',
]

outputname_s2 = '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/FigureS2' + ext
combine_figs(fig_list_s2, 4, 3, outputname_s2, width = width_h)



# -------------------- new figure
fig_list_s3 = [
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_Removal_GNG-a.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_Removal_GNG-b.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_Removal_GNG-c.png',
  '/home/benjamin.garzon/Data/DAD/analysis_variability/figsNOGSR/Figure_Removal_GNG-d.png'
]

outputname_s3= '/home/benjamin.garzon/Data/DAD/analysis_variability/finalfigs/FigureS3' + ext
combine_figs(fig_list_s3, 2, 2, outputname_s3, width = width_h)
