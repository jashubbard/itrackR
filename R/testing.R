#
#
setwd('~/Dropbox/data_archive/switchrate_eug/')


alledfs <- list.files(pattern=glob2rx('4*.edf'))
edfs <- alledfs[1:2]


test <- itrackr(edfs=edfs,samples=F)
# #
# # plot(test,zoom=T)
# # # test <- load_edfs(test,path='~/Dropbox/data_archive/switchrate_eug/',pattern='4*.edf')
# # #
# # #
# # # allbatch <- edfR::edf.batch()
#
#
pat <- '^BLOCK [0-9]*'


obj <- test

as.numeric(gsub("[^0-9]", "",obj$messages$message[grepl(pat,obj$messages$message)]))


source('~/Dropbox/analysis_lib/r_helpers.R')
beh <- load_files('~/Dropbox/finalized_experiments/switchrate_eug/data','^4*.txt')

test<- index_vars(test,c('Block','Trial'),c('^BLOCK [0-9]*','^TRIAL [0-9]*'),numeric.only=T)
test <- find_messages(test,c('STIMONSET','RESPONSE'),c('STIMONSET','RESPONSE'),timestamp = T)
test <- add_behdata(test,beh)


test2 <- epoch_samples(test,'STIMONSET',epoch=c(-1100,200))

