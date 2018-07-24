
rm(list=ls())
library(ggplot2)

# GNG: 210 vols x 3 runs x 2 sec
# TAB: 330 vols x 2 runs x 2 sec
# GNG: 170 vols x 1 runs x 2 sec

motion_file = "/home/benjamin.garzon/Data/DAD/processed/fmriprep/connectomes/FramewiseDisplacement.csv"
motion = read.csv(motion_file, sep = ';')

maxplot = ggplot(motion, aes(as.factor(RUN), FD_MAX)) + geom_boxplot() + facet_grid( .~ TASK)
print(maxplot)

meanplot = ggplot(motion, aes(as.factor(RUN), FD_MEAN)) + geom_boxplot() + facet_grid( .~ TASK)
print(meanplot)

propplot = ggplot(motion, aes(as.factor(RUN), FD_2)) + geom_boxplot() + facet_grid( .~ TASK)
print(propplot)


lenient = subset(motion, FD_MEAN < 0.55)
print(table(lenient$TASK))

stringent = subset(motion, FD_MEAN < 0.3 & FD_MAX < 5 & FD_3 < 0.3)
print(table(stringent$TASK))

stringent.GNG = subset(stringent, TASK == 'GNG')
stringent.TAB = subset(stringent, TASK == 'TAB')
stringent.RS = subset(stringent, TASK == 'RS')

print(
  c(length(unique(stringent.GNG$SUBJECT)), 
    length(unique(stringent.RS$SUBJECT)),
    length(unique(stringent.TAB$SUBJECT))
  )
)
