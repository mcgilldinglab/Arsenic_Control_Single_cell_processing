library(cisTopic)
library(Matrix)

pathTo10X <- './'
fragments <- paste0(pathTo10X, 'data_ARS/atac_fragments.tsv.gz')
metrics <- paste0(pathTo10X, 'data_ARS/per_barcode_metrics.csv')
bed <- paste0(pathTo10X, 'data_ARS/atac_peaks.bed')
cisTopicObject <- createcisTopicObjectFrom10X( fragments, bed, metrics , project.name='ARS')



cisTopicObject <- runModels(cisTopicObject, topic=c(100), seed=987, nCores=16, burnin = 120, iterations = 150, addModels=FALSE)
logLikelihoodByIter(cisTopicObject, select = c(2,5,10,15,20,30,40))
cisTopicObject <- selectModel(cisTopicObject, select=100)
saveRDS(cisTopicObject,file = "./cistopic_ARS_100.Rds")

A <- t(cisTopicObject@selected.model$document_expects)
write.csv(A,file = "./Cistopic_ARS.csv",row.names = cisTopicObject@cell.names)

getBedFiles(cisTopicObject, path='/group_homes/abc/home/u2600703/cistopic/output/cisBed100')

