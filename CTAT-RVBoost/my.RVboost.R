#!/usr/bin/env Rscript

stdin <- commandArgs(TRUE) 

if(length(stdin) < 1 || length(stdin) > 2) {
	stop("ERROR! Incorrect number of arguments. \nUSAGE: RVboost.R input_VCF [attrs=\"DJ,PctExtPos,ReadPosRankSum,QD,FS,ED\"] \n\n")
}

###arguments
inputVCF <- stdin[1]


### use these attributes from VCF file to make a model
sel.attri <- c("DJ","PctExtPos","ReadPosRankSum","QD","FS","ED") 

if (length(stdin) == 2) {
    sel.attri = strsplit(stdin[2], ",")[[1]]
}
message("Using attribute list: ", sel.attri)

##library <- "/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/Rlibs/"
##.libPaths(library)

hapmap <- "/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/resources/hapmap.ids.txt.gz"
model <- "adaboost"
output <- "rvboost_outdir"


if (! file.exists(output)) {
    dir.create(output);
}

### setting R library

library(VariantAnnotation)
library(gbm)
library(mgcv)


parseVCF <- function(VCF.filename,
                     genome.ver = "hg19", 
                     sel.info.attr = c("QD","ReadPosRankSum","DP","FS","MQ"),
                     with.fliter = TRUE)
{

    message("parsing vcf file:", VCF.filename)
    tmp.vcf <- suppressWarnings(readVcf(file=VCF.filename,genome=genome.ver))
    
    info.res <- as.data.frame(info(tmp.vcf))

    message("Existing and available attributes include: ")
    print(colnames(info.res))
    
    check.flag <- is.element(sel.info.attr,colnames(info.res))
    if(!all(check.flag))
    {
        stop(paste("the following attribute(s) not available in provided VCF: \n",
                   sel.info.attr[which(check.flag==FALSE)])
             )
    }
    
    
    N.attri <- length(sel.info.attr)
    N.var <- nrow(info.res)
    attri.mtx <- mat.or.vec(N.var,N.attri)
    colnames(attri.mtx) <- sel.info.attr
    rownames(attri.mtx) <- rownames(info.res)
    for(k in 1:N.attri)
    {
        attri.mtx[,k] <- as.numeric(info.res[[sel.info.attr[k]]])
    }
    
    ret.res <- list()
    ret.res$attri.mtx <- attri.mtx
    
    if(with.fliter) # with VQSR results, usually for comparison purpose
    {
        ret.res$filt.vec <- filt(tmp.vcf)
        ret.res$VQSR.score <- info.res[["VQSLOD"]]
    }
    
    ret.res$GT.vec <- geno(tmp.vcf)$GT
    ##ret.res$pos.vec <- paste(rowData(tmp.vcf)@seqnames, start(rowData(tmp.vcf)@ranges),sep=":")
    ret.res$pos.vec <- paste(seqnames(tmp.vcf), start(ranges(tmp.vcf)),sep=":")  # from vadir
    
    return(ret.res)
}


fitRVmodel <- function( input.mtx, DB.filename, DB.ID="rsid", pos.vec=NULL,
                       model.type="adaboost", ada.n=5e3)
{

    message("Fitting RV model")
    model.type <- tolower(model.type)
    
    provided.model.type.vec <- 
        c("glm.bino","glm.gaus","gam.bino","gam.gaus","adaboost","huberized","bernoulli")
    
    if(!is.element(model.type, provided.model.type.vec))  {
        stop("Input model type is not supported! \n")
    }
    
    
    tmp.read <- unlist(read.delim(gzfile(DB.filename),header=FALSE,colClasses="character"))
    if(DB.ID=="rsid")       {
        train.label <- as.numeric(is.element(rownames(input.mtx),tmp.read))   }
    if(DB.ID=="pos" & !is.null(pos.vec)) {
        train.label <- as.numeric(is.element(pos.vec,tmp.read))  }      
    
    cat("\n distribution of variants in given SNP database (1 = in, 0 = not in) \n")
    print(table(train.label))
    
    
    if (model.type=="glm.bino")  { RV.res <- glm(train.label~input.mtx,family=binomial) }
    if (model.type=="glm.gaus")  { RV.res <- glm(train.label~input.mtx,family=gaussian) }
    if (model.type=="gam.bino")  { RV.res <- gam(train.label~input.mtx,family=binomial) }
    if (model.type=="gam.gaus")  { RV.res <- gam(train.label~input.mtx,family=gaussian) }
    
    if (model.type=="adaboost" | model.type=="huberized" | model.type=="bernoulli"){  
        RV.res <- gbm.fit(x=input.mtx,
                          y=train.label,
                          n.trees=ada.n,
                          interaction.depth=2,
                          distribution=model.type,
                          verbose=FALSE) 
        
        RV.res$fitted.values <-  plogis(RV.res$fit) 
                                        # convert it to the 0-1 scale since the adaboost method gives the predictions on logit scale. 
                                        # http://stats.stackexchange.com/questions/37497/how-to-use-r-gbm-with-distribution-adaboost
    }
    
    
    RV.res$ID <- rownames(input.mtx)
    RV.res$train.label <- train.label
    
    return(RV.res)
    
}


parseRNA.res <- parseVCF(inputVCF,
                         sel.info.attr=sel.attri)

tmp.mtx <- parseRNA.res$attri.mtx
tmp.mtx[,"ReadPosRankSum"] <- abs(tmp.mtx[,"ReadPosRankSum"]) # make ReandPosRank monotonical

                                        #=== imputation of missing values
for(k in 1:ncol(tmp.mtx)){
    sel.na.idx <- which(is.na(tmp.mtx[,k]))
    sel.nan.idx <- which(!is.na(tmp.mtx[,k]))
    if(!is.null(sel.na.idx)){
        tmp.mtx[sel.na.idx,k] <- median(tmp.mtx[sel.nan.idx,k])
    }
}


write.table(tmp.mtx, file=paste(output, "tmp.mtx", sep="/"), quote=F, sep="\t")

                                        #=== fit adaboost model
fit.res <- fitRVmodel(input.mtx=tmp.mtx,
                      DB.filename=hapmap,
                      DB.ID="pos",
                      pos.vec=parseRNA.res$pos.vec,
                      ada.n=2e4)


## === compute SNP.conf score
RVboost.ECDF <- ecdf(fit.res$fitted.values[which(fit.res$train.label==1)])
RVboost.Q.score <- RVboost.ECDF(fit.res$fitted.values)
score <- paste (output,"original_score.txt",sep="/");
write.table(fit.res$fitted.values,file=score,col.names=FALSE,row.names=FALSE,quote=FALSE)
Qscore <- paste (output,"RV.Qscore.txt",sep="/"); 
write.table(RVboost.Q.score,file=Qscore,col.names=FALSE,row.names=FALSE,quote=FALSE)
#=== Final outputs to export 
#=== 1. fit.res$fitted.values : orginial adaboost scores 
#=== 2. RVboost.Q.score : SNP% confidence score
