rm(list=ls())
setwd('/home/nzhou/git/CAFApi/post_ddl/assessment/aggregate/')
#nametable = read.csv('./name_table.csv', colClasses = rep("character",5))
#toggle below
#row=5
#toggle above
# metric = nametable[row,'metric']
# taxon = nametable[row,'taxon']
# taxon_name = nametable[row, 'taxon_name']
# go = nametable[row, 'go']
# go_name = nametable[row,'go_name']
metric = 'AUC'
taxon = '208963'
taxon_name = 'Pseudomonas'
go = '0042710'
go_name = 'biofilm formation'

agg<-read.table(paste('./pi_',metric,'_',taxon,'_',go,'_sorted.txt',sep=""), header=F)
colnames(agg)<-c('teamname', 'modelnum', 'metric')
#baseline methods are dealt with separately
bsl = c('expression', 'blast', 'blastcomp')

  
#Use the top n teams (Use the top model from that team)
teams_indices = c()
teams_name = c()
n = 5
#And then put the baseline models at the end
#Just use loops because why not
for (i in 1:nrow(agg)){
  if (!agg$teamname[i] %in% bsl) {
    if (!agg$teamname[i]%in%teams_name){
      teams_indices  = c(teams_indices,i)
      teams_name = c(teams_name, as.character(agg$teamname[i]))
    }
  }
  if (length(teams_indices)==n) break
}

bsl_indices=which(agg$teamname %in% bsl)
#manually change order of bsl_indices to be "blast1, blast2, blastcomp1, blastcomp2"
agg[bsl_indices,]
bsl_indices=c(2,1,40,30,16,17)
indices = c(teams_indices, bsl_indices)
cols = rep("#999999",length(teams_indices))
cols=c(cols, rep(c('#fe9929','#C4302B','#00539F'), each=length(bsl_indices)/3)) #yellow for expression,red for blast, blue for blastcomp, 

#adaptive y_lim
y_lb = floor(max(min(agg[indices,"metric"])-0.2, 0)*10)/10
y_ub = ceiling(min(max(agg[indices,"metric"])+0.2, 1)*10)/10
adpt = F
png(filename = paste0("barplot_",metric,"_", go,"_", taxon, ".png"), 
    width = 2000, height=1500, res=300)
par(mar = c(6,4.1,6.1,2.1))
if (adpt){
  ylim= c(y_lb, y_ub)
} else{
  ylim = c(0.0,0.8)  
}
xx<-barplot(agg[indices, "metric"], ylim = ylim, las=2, xpd=F
            ,col = cols, border = F, xaxt='n', ylab = "AUC")
abline(h=0, col=1)
title(paste0(go_name," in ",taxon_name), line=1)
#axis(side=1,at=xx,labels=agg[indices,1])
textcols = c(rep(1, length(teams_indices)), rep(0, length(bsl_indices)))
text(x = xx, y = agg[indices,"metric"]-0.09, label = round(agg[indices,"metric"],2), 
     pos = 3, cex = 0.9, col = textcols)
labels=c(as.character(agg$teamname[teams_indices]), paste0(as.character(agg$teamname[bsl_indices]), agg$modelnum[bsl_indices]))
axis(1, at=xx, labels = F, tick=FALSE, line=-0.5, cex.axis=1.1, srt = 45, adj = 1)
text(x=xx, y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), labels=labels,
     srt = 45, adj = 1, xpd=T, cex = 0.9)
#manually choose which is the highest bsl
abline(h = agg$metric[bsl_indices[2]], col = cols[6], lwd = 2, lty = 2)
dev.off()

     