
# remove some ICN weights??
FIGS_DIR = "~/Data/DAD/processed/analysis_variability/modules/"
unlink(paste0(FIGS_DIR, '*'))
fig_count = 0


HCTHRESH = 0.42 # hierachical clustering threshold
ONLY_BEHAV = TRUE # use ońly behavioural domains


# -------------
MODULES_FILE="modules.R"
MODULES_FILE_MNI="/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/parcellations/shen/modules_MNI.R"
MODULES_FILE_MNI_RED="/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/parcellations/shen/modules_MNI_reduced.R"

MODULES_FILE_20="~/Data/DAD/ICNs/Laird/modules.R"
MODULES_FILE_70="~/Data/DAD/ICNs/Ray/modules.R"

INPUT_FILE.TAB='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/TAB/Connectome0.3/zFC_all_150_0.3_TAB.mat'
INPUT_FILE.GNG='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/GNG/Connectome0.3/zFC_all_150_0.3_GNG.mat'
INPUT_FILE.RS='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/RS/Connectome0.4/zFC_all_150_0.4_RS.mat'
#INPUT_FILE.0.4.TAB='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/TAB/Connectome0.4/zFC_all_150_0.4_TAB.mat'
#INPUT_FILE.0.4.GNG='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/GNG/Connectome0.4/zFC_all_150_0.4_GNG.mat'

INPUT_FILE.TAB.valid='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/TAB/Connectome0.3/zFC_all_150_0.3_TAB_valid.mat'
INPUT_FILE.GNG.valid='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/GNG/Connectome0.3/zFC_all_150_0.3_GNG_valid.mat'
INPUT_FILE.RS.valid='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/RS/Connectome0.4/zFC_all_150_0.4_RS_valid.mat'
#INPUT_FILE.0.4.GNG.valid='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/GNG/Connectome0.4/zFC_all_150_0.4_GNG_valid.mat'
#INPUT_FILE.0.4.TAB.valid='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/TAB/Connectome0.4/zFC_all_150_0.4_TAB_valid.mat'

INPUT_FILE.TAB.mod.20='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/TAB/modules/zFC_all_150_0.3_TAB_20.mat'
INPUT_FILE.GNG.mod.20='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/GNG/modules/zFC_all_150_0.3_GNG_20.mat'
INPUT_FILE.RS.mod.20='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/RS/modules/zFC_all_150_0.4_RS_20.mat'

INPUT_FILE.TAB.mod.70='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/TAB/modules/zFC_all_150_0.3_TAB_70.mat'
INPUT_FILE.GNG.mod.70='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/GNG/modules/zFC_all_150_0.3_GNG_70.mat'
INPUT_FILE.RS.mod.70='/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/processed/RS/modules/zFC_all_150_0.4_RS_70.mat'

ROI_VOL_FILE="/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/parcellations/shen/parc_shen_150.volume.csv"
LABELS_FILE="/home/ALDRECENTRUM/benjamin.garzon/Data/DAD/parcellations/shen/fconn_150_labels.txt"
labels = read.csv(LABELS_FILE, header=TRUE, sep='\t')
# -------------

# sync valid indices
info.TAB = readMat(INPUT_FILE.TAB)
info.GNG = readMat(INPUT_FILE.GNG)
info.RS = readMat(INPUT_FILE.RS)
#info.0.4.TAB = readMat(INPUT_FILE.0.4.TAB)
#info.0.4.GNG = readMat(INPUT_FILE.0.4.GNG)

merged_matrices = abind(
  info.TAB$merged.matrices, 
  info.GNG$merged.matrices, 
  info.RS$merged.matrices, 
  #        info.0.4.TAB$merged.matrices, 
  #        info.0.4.GNG$merged.matrices,
  along=3)

valid = !apply(merged_matrices, c(1,2), function(x) any(is.na(x)))
valid = apply(valid, 2, any)
valid.indices = which(valid)
n.rois = length(valid)

writeMat(INPUT_FILE.TAB.valid,
         subjects=unlist(info.TAB$subjects), merged.matrices = info.TAB$merged.matrices, merged.matrices.mat=flatten(info.TAB$merged.matrices, valid.indices), labels=unlist(info.TAB$labels))
writeMat(INPUT_FILE.GNG.valid,
         subjects=unlist(info.GNG$subjects), merged.matrices = info.GNG$merged.matrices, merged.matrices.mat=flatten(info.GNG$merged.matrices, valid.indices), labels=unlist(info.GNG$labels))
writeMat(INPUT_FILE.RS.valid, 
         subjects=unlist(info.RS$subjects), merged.matrices = info.RS$merged.matrices, merged.matrices.mat=flatten(info.RS$merged.matrices, valid.indices), labels=unlist(info.RS$labels))
#writeMat(INPUT_FILE.0.4.TAB.valid,
#	subjects=unlist(info.0.4.TAB$subjects), merged.matrices.mat=flatten(info.0.4.TAB$merged.matrices, valid.indices), labels=unlist(info.0.4.TAB$labels))
#writeMat(INPUT_FILE.0.4.GNG.valid,
#	subjects=unlist(info.0.4.GNG$subjects), merged.matrices.mat=flatten(info.0.4.GNG$merged.matrices, valid.indices), labels=unlist(info.0.4.GNG$labels))


# -------------
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

# modules for 20 regions
WD="~/Data/DAD/ICNs/Laird/"
artifacts = c(19, 20)
ICNnames = formatC(seq(20), width=2, flag="0")
compute_modules(WD, artifacts, ICNnames)

# modules for 70 regions
WD="~/Data/DAD/ICNs/Ray/"
artifacts = c(66, 69, 70)
ICNnames = formatC(seq(70), width=2, flag="0")
compute_modules(WD, artifacts, ICNnames, figname="FigureSupp2d")

if( F ) {
# average FC in modules
average_over_modules(INPUT_FILE.TAB.valid, MODULES_FILE_20, INPUT_FILE.TAB.mod.20)
average_over_modules(INPUT_FILE.GNG.valid, MODULES_FILE_20, INPUT_FILE.GNG.mod.20)
average_over_modules(INPUT_FILE.RS.valid, MODULES_FILE_20, INPUT_FILE.RS.mod.20)

average_over_modules(INPUT_FILE.TAB.valid, MODULES_FILE_70, INPUT_FILE.TAB.mod.70)
average_over_modules(INPUT_FILE.GNG.valid, MODULES_FILE_70, INPUT_FILE.GNG.mod.70)
average_over_modules(INPUT_FILE.RS.valid, MODULES_FILE_70, INPUT_FILE.RS.mod.70)
}
