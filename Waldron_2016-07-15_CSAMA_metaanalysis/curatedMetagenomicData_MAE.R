## curatedMetagenomicData and ExperimentHub demo
## requires bioc-devel (2.4)

library(ExperimentHub)
eh = ExperimentHub()
myquery = query(eh, "curatedMetagenomicData")

View(mcols(myquery))

candela.eh = display(myquery)

##candela.ids = c("EH2", EH6", "EH7")
candela.list = lapply(names(candela.eh), function(id) candela.eh[[id]])
names(candela.list) = candela.eh$title

library(MultiAssayExperiment)
mae = MultiAssayExperiment(ExperimentList = candela.list, 
                           pData = pData(candela.list[[2]]))

experiments(mae)
rownames(mae)
colnames(mae)
mae[1:6, , ]
mae[, pData(mae)$gender == "male", ]
mae[, , 1:2]

## Subsetting by CharacterList