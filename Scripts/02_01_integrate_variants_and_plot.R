##################################################################################################################
### Description: Integrate classifier predictions to obtain final mutation set
###              

library(dplyr)
library(ggplot2)
library(VennDiagram)
workdir = '/rocker-build/gmkf_nbl_somatic/'


#import mutations files
consensusVarsFile = paste0(workdir,'Data/coding_mutations_annotated.tsv')
uniqueVarsFile = paste0(workdir,'Data/unique_coding_mutations_full.tsv')
classifierPredictionsFile = paste0(workdir,'Data/classifier_results_vars_predicted.tsv')

consensusVars = read.csv(consensusVarsFile,sep='\t')
uniqueVars = read.csv(uniqueVarsFile,sep='\t')
classifierPredictions = read.csv(classifierPredictionsFile,sep='\t')
colnames(classifierPredictions)[1] = 'mutid'

#attach classifier predictions
uniqueVars = dplyr::left_join(uniqueVars,classifierPredictions,by='mutid')
consensusVars = dplyr::left_join(consensusVars,classifierPredictions,by='mutid')

#filter consensus variants
'''
the final call set include any variant call made by at least 3 callers
any call made by 2 or fewer calls need a somatic prediction by the classifier
to be included in the final call set.
'''
finalConsVars = consensusVars %>% filter(numCallers > 2 | Prediction != 'fail')
finalUnVars = uniqueVars %>% filter(Prediction == 'somatic')


#--------------------------------------------------------------------------------------------
#plot - venn diagram of overlapping callers in final variant call set
'''
callers:
m - mutect2
s - strelka2
l - lancet
v - vardict2
'''
table(finalUnVars$caller)
m = finalUnVars %>% filter(caller == 'mutect') %>% tally() %>% pull(n)
s = finalUnVars %>% filter(caller == 'strelka') %>% tally() %>% pull(n)
l = finalUnVars %>% filter(caller == 'lancet') %>% tally() %>% pull(n)
v = finalUnVars %>% filter(caller == 'vardict') %>% tally() %>% pull(n)
table(consensusVars$caller)
ms = consensusVars %>% filter(caller == 'mutect,strelka') %>% tally() %>% pull(n)
ml = consensusVars %>% filter(caller == 'mutect,lancet') %>% tally() %>% pull(n)
mv = consensusVars %>% filter(caller == 'mutect,vardict') %>% tally() %>% pull(n)
sl = consensusVars %>% filter(caller == 'lancet,strelka') %>% tally() %>% pull(n)
sv = consensusVars %>% filter(caller == 'vardict,strelka') %>% tally() %>% pull(n)
lv = consensusVars %>% filter(caller == 'lancet,vardict') %>% tally() %>% pull(n)
msl = consensusVars %>% filter(caller == 'mutect,strelka,lancet') %>% tally() %>% pull(n)
msv = consensusVars %>% filter(caller == 'vardict,mutect,strelka') %>% tally() %>% pull(n)
mlv = consensusVars %>% filter(caller == 'mutect,lancet,vardict') %>% tally() %>% pull(n)
slv = consensusVars %>% filter(caller == 'strelka,lancet,vardict') %>% tally() %>% pull(n)
mslv = consensusVars %>% filter(caller == 'mutect,strelka,lancet,vardict') %>% tally() %>% pull(n)

#make venn diagram
callerOverlap = setNames(as.numeric(c(m,s,l,v,ms,ml,mv,sl,sv,lv,msl,msv,mlv,slv,mslv)),
                         nm = c('mutect','strelka','lancet','vardict',
                                'mutect,strelka','mutect,lancet','mutect,vardict',
                                'strelka,lancet','strelka,vardict','lancet,vardict',
                                'mutect,strelka,lancet','mutect,strelka,vardict',
                                'mutect,lancet,vardict','strelka,lancet,vardict',
                                'mutect,strelka,lancet,vardict'))

grid.newpage()
draw.quad.venn(
  area1 = sum(callerOverlap[c(1, 5, 6, 7, 11, 12, 13, 15)]),  # mutect total
  area2 = sum(callerOverlap[c(2, 5, 8, 9, 11, 12, 14, 15)]),  # strelka total
  area3 = sum(callerOverlap[c(3, 6, 8, 10, 11, 13, 14, 15)]), # lancet total
  area4 = sum(callerOverlap[c(4, 7, 9, 10, 12, 13, 14, 15)]), # vardict total
  n12 = callerOverlap["mutect,strelka"] + callerOverlap["mutect,strelka,lancet"] + callerOverlap["mutect,strelka,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n13 = callerOverlap["mutect,lancet"] + callerOverlap["mutect,strelka,lancet"] + callerOverlap["mutect,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n14 = callerOverlap["mutect,vardict"] + callerOverlap["mutect,strelka,vardict"] + callerOverlap["mutect,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n23 = callerOverlap["strelka,lancet"] + callerOverlap["mutect,strelka,lancet"] + callerOverlap["strelka,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n24 = callerOverlap["strelka,vardict"] + callerOverlap["mutect,strelka,vardict"] + callerOverlap["strelka,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n34 = callerOverlap["lancet,vardict"] + callerOverlap["mutect,lancet,vardict"] + callerOverlap["strelka,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n123 = callerOverlap["mutect,strelka,lancet"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n124 = callerOverlap["mutect,strelka,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n134 = callerOverlap["mutect,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n234 = callerOverlap["strelka,lancet,vardict"] + callerOverlap["mutect,strelka,lancet,vardict"],
  n1234 = callerOverlap["mutect,strelka,lancet,vardict"],
  category = c("Mutect2", "Strelka2", "Lancet", "VarDict2"),
  fill = c("red", "blue", "green", "yellow"),
  alpha = c(0.4, 0.4, 0.3,0.5),
  cex = 2,
  cat.cex = 2
)


#--------------------------------------------------------------------------------------------
#plot - Unit circle filled proportionally by variant caller overlap in final variant call set
# Input counts and labels
counts = c(m+s+l+v,ms+ml+mv+sl+sv+lv,msl+msv+mlv+slv,mslv)
labels = c("Unique", "2 Callers", "3 Callers", "Full Consensu")
props = counts / sum(counts)

# Compute left/right x-boundaries starting from -1 to 1
x_bounds = c(-1, -1 + cumsum(props) * 2)  # because width of full circle is 2
x_start = head(x_bounds, -1)
x_end = tail(x_bounds, -1)

# Generate vertical wedges that follow the circle boundary
make_wedge <- function(xmin, xmax, n = 100) {
  x <- seq(xmin, xmax, length.out = n)
  y_top <- sqrt(1 - x^2)
  y_bot <- -y_top
  data.frame(
    x = c(x, rev(x)),
    y = c(y_top, rev(y_bot))
  )
}

# Combine data for all wedges
circle_data <- bind_rows(lapply(1:length(counts), function(i) {
  wedge <- make_wedge(x_start[i], x_end[i])
  wedge$group <- labels[i]
  return(wedge)
}))

circle_outline <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 200)),
  y = sin(seq(0, 2 * pi, length.out = 200))
)

fills = c('#D4D146','#39A865','#6B82D1','#C75B3E')

# Plot with black outline
ggplot(circle_data, aes(x = x, y = y, fill = group)) +
  geom_polygon(aes(group = group)) +
  geom_path(data = circle_outline, aes(x = x, y = y), inherit.aes = FALSE, color = "black", size = 1) +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = fills) +
  ggtitle("Proportions of variant caller overlap in final consensus pool") +
  theme(legend.position = "none")

#-------------------------------------------
#Ensure final variant count:
finalConsCount = nrow(finalConsVars) + nrow(finalUnVars)
print(finalConsCount)