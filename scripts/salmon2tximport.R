mainDir <- paste0(system('echo $HOME/', intern = T), 'Documents/Saeed/')

#----Make a coldata matrix----
run <- list.files(paste0(mainDir, 'data/'))
condition <- c(rep('IRE1_KO', 6), rep('IRE1_WT', 6))
batch <- rep((1:6), 2)

samples <- data.frame(run, condition, batch, row.names = run)

#----Locate counts data----
files <- file.path(paste0(mainDir, 'quants'), paste0(samples$run, '_quant'), 'quant.sf')             

names(files) <- samples$run
all(file.exists(files))

library(tximport)
library(ensembldb)
library(readr)

#----Generate transcript database from annotation file----
txdb <- makeTxDbFromGFF(file = './gencode.v42.annotation.gtf.gz')
saveDb(x=txdb, file = "./gencode.v42.annotation.TxDb")

#txdb <- loadDb(file = './gencode.v42.annotation.TxDb')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- select(txdb, k, "GENEID", 'TXNAME')

# Remove unwanted string in transcript database----
tx2gene$TXNAME <- gsub("\\..*", "", tx2gene$TXNAME)
tx2gene$GENEID <- gsub("\\..*", "", tx2gene$GENEID)

#----Read counts for each sample----
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = T)

#----session----
sessionInfo()