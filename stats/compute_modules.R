
# Average the parameters over anatomical and BrainMap modules
WORKDIR="~/Data/DAD/processed/fmriprep"
PARCDIR="~/Data/DAD/parcellations"

FIGS_DIR = file.path(WORKDIR, "connectomes/modules/")
unlink(paste0(FIGS_DIR, '*'))
fig_count = 0
HCTHRESH = 0.42 # hierachical clustering threshold
ONLY_BEHAV = TRUE # use ońly behavioural domains

# ------------------------------------------------------------------------
# paths

MODULES_FILE="modules.RData"
MODULES_FILE_MNI=file.path(PARCDIR, "shen/modules_MNI.RData")
MODULES_FILE_MNI_RED=file.path(PARCDIR, "shen/modules_MNI_reduced.RData")

MODULES_FILE_20="~/Data/DAD/ICNs/Laird/modules.RData"
MODULES_FILE_70="~/Data/DAD/ICNs/Ray/modules.RData"

INPUT_FILE.TAB=file.path(WORKDIR, 'connectomes/TAB/zFC_all_150.mat')
INPUT_FILE.GNG=file.path(WORKDIR, 'connectomes/GNG/zFC_all_150.mat')
INPUT_FILE.RS=file.path(WORKDIR, 'connectomes/RS/zFC_all_150.mat')

INPUT_FILE.TAB.valid=file.path(WORKDIR, 'connectomes/TAB/zFC_all_150_valid.mat')
INPUT_FILE.GNG.valid=file.path(WORKDIR, 'connectomes/GNG/zFC_all_150_valid.mat')
INPUT_FILE.RS.valid=file.path(WORKDIR, 'connectomes/RS/zFC_all_150_valid.mat')

INPUT_FILE.TAB.mod.20=file.path(WORKDIR, 'connectomes/TAB/modules/zFC_all_150_TAB_20.mat')
INPUT_FILE.GNG.mod.20=file.path(WORKDIR, 'connectomes/GNG/modules/zFC_all_150_GNG_20.mat')
INPUT_FILE.RS.mod.20=file.path(WORKDIR, 'connectomes/RS/modules/zFC_all_150_RS_20.mat')

INPUT_FILE.TAB.mod.70=file.path(WORKDIR, 'connectomes/TAB/modules/zFC_all_150_TAB_70.mat')
INPUT_FILE.GNG.mod.70=file.path(WORKDIR, 'connectomes/GNG/modules/zFC_all_150_GNG_70.mat')
INPUT_FILE.RS.mod.70=file.path(WORKDIR, 'connectomes/RS/modules/zFC_all_150_RS_70.mat')

ROI_VOL_FILE=file.path(PARCDIR, "shen/parc_shen_150.volume.csv")
LABELS_FILE=file.path(PARCDIR, "shen/fconn_150_labels.txt")
labels = read.csv(LABELS_FILE, header=TRUE, sep='\t')

if( T ) {
  
# ------------------------------------------------------------------------
# sync valid indices
info.TAB = readMat(INPUT_FILE.TAB)
info.GNG = readMat(INPUT_FILE.GNG)
info.RS = readMat(INPUT_FILE.RS)

merged_matrices = abind(
  info.TAB$merged.matrices, 
  info.GNG$merged.matrices, 
  info.RS$merged.matrices, 
  along=3)

# ------------------------------------------------------------------------

# select those that are valid for all tasks
valid = !apply(merged_matrices, c(1,2), function(x) any(is.na(x)))
diag(valid) = F
valid = apply(valid, 2, any)
valid.indices = which(valid)
n.rois = length(valid)

# write out 
writeMat(INPUT_FILE.TAB.valid,
         subjects=unlist(info.TAB$subjects), 
         merged.matrices = info.TAB$merged.matrices, 
         merged.matrices.mat=flatten(info.TAB$merged.matrices, valid.indices), 
         labels=unlist(info.TAB$labels))
writeMat(INPUT_FILE.GNG.valid,
         subjects=unlist(info.GNG$subjects), 
         merged.matrices = info.GNG$merged.matrices, 
         merged.matrices.mat=flatten(info.GNG$merged.matrices, valid.indices), 
         labels=unlist(info.GNG$labels))
writeMat(INPUT_FILE.RS.valid, 
         subjects=unlist(info.RS$subjects), 
         merged.matrices = info.RS$merged.matrices, 
         merged.matrices.mat=flatten(info.RS$merged.matrices, valid.indices), 
         labels=unlist(info.RS$labels))

#info.TAB.valid = readMat(INPUT_FILE.TAB.valid)
#info.GNG.valid = readMat(INPUT_FILE.GNG.valid)
#info.RS.valid = readMat(INPUT_FILE.RS.valid)

# ------------------------------------------------------------------------
# get MNI modules
map_name = list(
  "Brain-Stem"="BStem",
  "Caudate"="Caud",
  "Cerebellum"="Cereb",
  "Frontal Lobe"="Front",
  "Insula"="Ins",
  "Occipital Lobe"="Occip",
  "Parietal Lobe"="Pariet",
  "Putamen"="Put",       
  "Temporal Lobe"="Temp",
  "Thalamus"="Thal"
)

order_partition = list(
  "BStem"=10, 
  "LCaud"=2, "LCereb"=1, "LFront"=3, "LIns"=8, "LOccip"=9, "LPariet"=7,
  "LPut"=6, "LTemp"=5, "LThal"=4 ,
  "RCaud"=12, "RCereb"=11, "RFront"=13, "RIns"=18, "ROccip"=19, "RPariet"=17,
  "RPut"=16, "RTemp"=15, "RThal"=14
  
)

names_partition = paste0(substring(labels$Hemisphere,1,1), unlist(map_name[labels$LabelMNI]))
#partition = as.numeric(as.factor(names_partition))
partition = unlist(order_partition[as.factor(names_partition)])

names(partition) = names_partition
module_names = names(order_partition)[order(unlist(order_partition))]
save(partition, n.rois, valid.indices, module_names, file = MODULES_FILE_MNI )

# MNI modules for small structures
partition[grep("Temp|Occip|Pariet|Front", names_partition)] = 0
partition = as.numeric(as.factor(partition))-1
names(partition) = names_partition
module_names = module_names[grep("Temp|Occip|Pariet|Front", module_names, invert=TRUE)]
save(partition, n.rois, valid.indices, module_names, file = MODULES_FILE_MNI_RED )

# ------------------------------------------------------------------------
# BrainMap modules

# modules for 20 networks
WD="~/Data/DAD/ICNs/Laird/"
artifacts = c(19, 20)
ICNnames = formatC(seq(20), width=2, flag="0")
compute_modules(WD, artifacts, ICNnames)

# modules for 70 networks
WD="~/Data/DAD/ICNs/Ray/"
artifacts = c(66, 69, 70)
ICNnames = formatC(seq(70), width=2, flag="0")
compute_modules(WD, artifacts, ICNnames, figname="FigureSupp2d")

# average FC in modules
average_over_modules(INPUT_FILE.TAB.valid, MODULES_FILE_20, INPUT_FILE.TAB.mod.20)
average_over_modules(INPUT_FILE.GNG.valid, MODULES_FILE_20, INPUT_FILE.GNG.mod.20)
average_over_modules(INPUT_FILE.RS.valid, MODULES_FILE_20, INPUT_FILE.RS.mod.20)

average_over_modules(INPUT_FILE.TAB.valid, MODULES_FILE_70, INPUT_FILE.TAB.mod.70)
average_over_modules(INPUT_FILE.GNG.valid, MODULES_FILE_70, INPUT_FILE.GNG.mod.70)
average_over_modules(INPUT_FILE.RS.valid, MODULES_FILE_70, INPUT_FILE.RS.mod.70)
}
