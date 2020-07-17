

#Get variant summary
#' @export
summaryVariant<-function(maf,vIdCol=c("Hugo_Symbol","Protein_Change"),outCol=colnames(maf@data),includeSilent=FALSE) {
	mafDat=maf@data

	mafSummaryVariantCount=mafDat[, .N, , vIdCol]
	setkeyv(mafSummaryVariantCount,vIdCol)
	setkeyv(mafDat,vIdCol)
	mafSummaryVariantOut=mafDat[mafSummaryVariantCount,mult = "first",nomatch = 0L]

	if (includeSilent) {
		message(paste0("summaryVariant: includeSilent==TRUE. Silent variants were added to summary table"))
		mafDat=maf@maf.silent

		mafSummaryVariantCount=mafDat[, .N, , vIdCol]
		setkeyv(mafSummaryVariantCount,vIdCol)
		setkeyv(mafDat,vIdCol)
		mafSummaryVariantOut=rbind(mafSummaryVariantOut,mafDat[mafSummaryVariantCount,mult = "first",nomatch = 0L])
	}

	outCol=c(vIdCol,"N",setdiff(outCol,vIdCol))
	message(paste0("summaryVariant: Column ", paste(vIdCol,collapse="+")," were used as Variants ID"))
#	if (nrow(mafSummaryVariantOut)>outNum) {
#   message(paste0(nrow(mafSummaryVariantOut)," variants in total and only top ",outNum," variants will be reported"))
#	}
	return(mafSummaryVariantOut[order(-N)][,..outCol])
}

#Prepare maf object for VariantOncoPlot
#' @export
prepareVariantOncoPlot<-function(mafObj,AACol="HGVSp_Short") {
  dataForPlot=mafObj
  temp=paste0(dataForPlot@data[,Hugo_Symbol],"_",dataForPlot@data[,get(AACol)])
  dataForPlot@data$Hugo_Symbol=temp
  dataForPlot@gene.summary =
    maftools:::summarizeMaf(maf = subsetMaf(maf = dataForPlot, fields = 'Hugo_Symbol', includeSyn = FALSE,mafObj = FALSE,genes=unique(temp)), chatty = FALSE)[[c("gene.summary")]]
  return(dataForPlot)
}

#Filter MAF file
#' @export
filterMaf=function(mafFile,
                   mafCol=c("AF","ExAC_AF","gnomAD_AF","gnomAD_exome_AF","gnomAD_genome_AF"),mafMax=0.001,
                   OrDefinedKeepCol=NULL,
                   ExAC_FILTER=TRUE,ExAC_FILTER_Col=c("ExAC_FILTER","gnomAD_exome_FILTER","gnomAD_genome_FILTER"),
                   variantInSampleMaxPercent=0.9) {
  mafTable=fread(cmd=paste0('grep -v "#" ',mafFile),quote="")
  maf=read.maf(mafTable,verbose=TRUE)

  matSubsetObj=maf
  mafRemovedTable=NULL

  #MAF max
  mafCol=intersect(mafCol,colnames(maf@data))
  if (!is.null(mafMax) & length(mafCol)>0) {
    print(paste0(mafCol))
    for (i in 1:length(mafCol)) {
      if (class(matSubsetObj@data[[mafCol[i]]])=="character") {
        #change maf data type
        temp=suppressWarnings(sapply(strsplit(matSubsetObj@data[[mafCol[i]]],"\\|"),function(x) max(as.numeric(x))))
        temp[is.infinite(temp)]=NA
        matSubsetObj@data[[mafCol[i]]]=temp
        #change maf.silent data type
        temp=suppressWarnings(sapply(strsplit(matSubsetObj@maf.silent[[mafCol[i]]],"\\|"),function(x) max(as.numeric(x))))
        temp[is.infinite(temp)]=NA
        matSubsetObj@maf.silent[[mafCol[i]]]=temp
      }
    }
    temp1=paste0(mafCol,"<",mafMax)
    temp2=paste0("is.na(",mafCol,") | (",mafCol,"=='')")
    temp3=paste0("(",temp1," | ",temp2,")")
    queryExp=paste0("(",paste(temp3,collapse=" & "),")")
  } else {
    queryExp=""
  }

  if (!is.null(OrDefinedKeepCol[1])) {
    #Somatic (Cosmic)
    #OrDefinedKeepCol=="SOMATIC"
    ##queryExp=paste0(queryExp," | ",'(SOMATIC!="")')
    #queryExp=paste0(queryExp," | (",OrDefinedKeepCol,'!="")')
  }

  ExAC_FILTER_Col=intersect(ExAC_FILTER_Col,colnames(maf@data))
  if (ExAC_FILTER & length(ExAC_FILTER_Col)>0) {
    #Exac Filter
    #matSubsetObj@data$ExAC_FILTER
    for (ExAC_FILTER_Col_One in ExAC_FILTER_Col) {
      queryExpOne=paste0('(',ExAC_FILTER_Col_One,'=="" | ',ExAC_FILTER_Col_One,'=="PASS" | is.na(',ExAC_FILTER_Col_One,'))')
      queryExp=paste0(queryExp," & ",queryExpOne)
    }
    #    queryExp=paste0(queryExp," & ",'(ExAC_FILTER=="" | ExAC_FILTER=="PASS" | is.na(ExAC_FILTER))')
  }

  print(queryExp)
  #matDataSubset=subsetMaf(maf = dataForReport[["maf"]], query = queryExp,includeSyn=TRUE)

  matSubsetObjKept=subsetMaf(maf = matSubsetObj, query = queryExp,includeSyn=TRUE,mafObj=TRUE)
  if (nrow(matSubsetObj@data)>nrow(matSubsetObjKept@data) | nrow(matSubsetObj@maf.silent)>nrow(matSubsetObjKept@maf.silent)) {
    temp=subsetMaf(maf = matSubsetObj, query =paste0("!(",queryExp,")") ,includeSyn=TRUE,mafObj=FALSE)
    mafRemovedTable<-rbind(mafRemovedTable,
                           cbind(temp,RemovedQuery=rep(queryExp,nrow(temp)),RemovedTag=rep("PublicDatabase",nrow(temp)))
    )
  }
  matSubsetObj=matSubsetObjKept

  if (!is.null(variantInSampleMaxPercent)) {
    #remove variants within more than 90% samples
    dataSampleN=as.integer(maf@summary["Samples",summary,on="ID"])
    vIdCol=c("Hugo_Symbol","HGVSp_Short","Protein_Change","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2")
    vIdCol=intersect(vIdCol,colnames(maf@data))

    #Get variant summary function
    ##Need to load summaryVariant function in Rmd file
    selectedVariantsDT=summaryVariant(matSubsetObj,vIdCol=vIdCol,includeSilent=TRUE)
    selectedVariantsDT=selectedVariantsDT[which(selectedVariantsDT$N>=max(dataSampleN*variantInSampleMaxPercent,2)),..vIdCol]
    #dim(selectedVariantsDT)
    if (nrow(selectedVariantsDT)>0) {
      print(paste0(nrow(selectedVariantsDT)," variants removed because exist in more than ",as.integer(variantInSampleMaxPercent*100),"% samples"))
      for (i in 1:nrow(selectedVariantsDT)) {
        queryExp=paste0("!(",paste(paste0(colnames(selectedVariantsDT),"=='",selectedVariantsDT[i,],"'"),collapse=" & "),")")
        #print(queryExp)

        matSubsetObjKept=subsetMaf(maf = matSubsetObj, query = queryExp,includeSyn=TRUE,mafObj=TRUE)
        if (nrow(matSubsetObj@data)>nrow(matSubsetObjKept@data) | nrow(matSubsetObj@maf.silent)>nrow(matSubsetObjKept@maf.silent)) {
          temp=subsetMaf(maf = matSubsetObj, query =paste0("!(",queryExp,")") ,includeSyn=TRUE,mafObj=FALSE)
          mafRemovedTable<-rbind(mafRemovedTable,
                                 cbind(temp,RemovedQuery=rep(queryExp,nrow(temp)),RemovedTag=rep("HighFrequencyInData",nrow(temp)))
          )
        }
        matSubsetObj=matSubsetObjKept

      }
    } else {
      print(paste0("variantInSampleMaxPercent defined at ",as.integer(variantInSampleMaxPercent*100),"% but no variant was more than it. All variants kept."))
    }
  }

  matDataSubset<-rbind( matSubsetObj@data,matSubsetObj@maf.silent)
  write.table(matDataSubset,paste0(tools::file_path_sans_ext(basename(mafFile)),".subset.maf"),sep="\t",quote=FALSE,row.names=FALSE)

  temp=table(mafRemovedTable$RemovedTag)
  names(temp)=paste0("Removed: ",names(temp))
  filterMafSummaryTable=c(totalRecords=(nrow(maf@data)+nrow(maf@maf.silent)),subsetRecords=nrow(matDataSubset),temp)
  cat("filterMafSummary:\n")
  print(filterMafSummaryTable)
  write.table(as.data.frame(filterMafSummaryTable),paste0(tools::file_path_sans_ext(basename(mafFile)),".filterMafSummary.txt"),sep="\t",quote=FALSE,col.names=FALSE)
  write.table(mafRemovedTable,paste0(tools::file_path_sans_ext(basename(mafFile)),".removed.maf"),sep="\t",quote=FALSE,row.names=FALSE)

  invisible(matDataSubset)
}





#modified based on read.maf to summary silent maf and return data.table
#' @export
summarySilentMaf<-function(maf,vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins",
                                                "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation",
                                                "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
                                                "Missense_Mutation")) {
  maf.silent = maf[!Variant_Classification %in% vc.nonSilent]
  if (nrow(maf.silent) > 0) {
    maf.silent.vc = maf.silent[, .N, .(Tumor_Sample_Barcode,
                                       Variant_Classification)]
    maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc,
                                           formula = Tumor_Sample_Barcode ~ Variant_Classification,
                                           fill = 0, value.var = "N")
    summary.silent = data.table::data.table(ID = c("Samples",
                                                   colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                            N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,
                                                                                                       2:ncol(maf.silent.vc.cast), with = FALSE])))
    #		maf = maf[Variant_Classification %in% vc.nonSilent]

    message(paste0("silent variants: ", nrow(maf.silent)))
    return(summary.silent)
  } else {
    message(paste0("silent variants: ", 0))
    return(NULL)
  }
}

#' @export
makeDataTable<-function(x) {
  DT::datatable(x, extensions = 'Buttons',
                rownames=FALSE,
                options = list(
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  scrollX = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                )
  )

}
