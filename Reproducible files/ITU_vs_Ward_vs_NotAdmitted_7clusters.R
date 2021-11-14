library(fst) #load the package
library(ccdf)
path <- paste0("/gpfs/home/mgauth911e/data", ".fst")
ft <- read_fst(path) # load the data as an fst object

metadata <- read.delim("GSE153931_cd8_t24_processed_data_annotations.txt", header = TRUE, sep = "\t", dec = ".")

# ID of task
slar_taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# remove patient void
void <- which(as.factor(metadata$orig.donor)=="void")
ft <- ft[,-as.numeric(void)]
metadata <- metadata[-void,]

# keep covid
illness <- rep(1,length(metadata$orig.donor))
for (k in 1:8){
  illness[which(metadata$orig.donor==unique(metadata$orig.donor)[24+k])] <- rep(0,length(which(
    metadata$orig.donor==unique(metadata$orig.donor)[24+k])))
}
covid <- metadata$orig.unit[which(illness==1)] # X
covid_cluster <- metadata$RNA_snn_res.0.2[which(illness==1)] # Z
ft <- ft[,which(illness==1)] # matrix

# arrange X
covid[which(as.factor(covid)=="void")] <- "Not_admitted"

covid_bin <- rep(1,length(covid))
covid_bin[which(as.factor(covid)=="Not_admitted")] <- rep(0,length(which(as.factor(covid)=="Not_admitted")))

# arrange Z

keep_cluster <- which(as.factor(covid_cluster)!=7)

pval <- rep(NA,14)
task <- rep(NA,14)
stat <- rep(NA,14)
j <- 0
for (i in (14*(slar_taskid-1)+1):(slar_taskid*14)){
  j <- j+1
  if (i>nrow(ft)){
    break
  }
  else{
    print(i)
    res <- test_asymp(as.numeric(ft[i,keep_cluster]),data.frame(X=factor(covid,levels=c("ITU","Ward","Not_admitted"))[keep_cluster]), 
                      #data.frame(Z=factor(covid_cluster,levels=c("0","1","2","3","4","5","6"))[keep_cluster]),
                      space_y = TRUE, number_y = 244, log = FALSE)
    pval[j] <- res$raw_pval
    stat[j] <- res$Stat
    task[j] <- i
  }
}

res <- data.frame(task=task,raw_pval=pval, test_statistic=stat)

### sauvegarde du fichier

filename_r <- paste0("result_r_2/CD8+ Analysis/244seuils/noZ/pval", "_", slar_taskid, ".csv")
write.csv(res, file = filename_r, row.names = F)