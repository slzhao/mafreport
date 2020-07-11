
#' @export
initialize_maf_report_parameter<-
  function(mafFile,
           reportOut=paste0(basename(mafFile),".report.html"),
           reportOutDir=dirname(mafFile),
           reportTemplate=system.file('templates', "report.Rmd", package = 'mafreport'),
           reportModules=c("Initialize.Rmd","SummaryTables.Rmd","VariantsVisualization.Rmd","HighFreqVariantsCheck.Rmd","VariantsAnalysis.Rmd"),
           genome="hg19",
           vafCol=c("t_vaf","tumor_f"),
           vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                       "Nonsense_Mutation",	"Nonstop_Mutation",
                       #"In_Frame_Del", "In_Frame_Ins",
                       "Missense_Mutation"),
           dndscv.refdb=genome,
           interestedGenes=NULL,
           clinicalData=NULL,
           clinicalFeatures=NULL,
           performGroupSummary=NULL,
           showCode=TRUE
           ) {

  dataForReport=list()
  dataForReport[["mafFile"]]=mafFile
  dataForReport[["reportOut"]]=reportOut
  dataForReport[["reportOutDir"]]=reportOutDir
  dataForReport[["reportTemplate"]]=reportTemplate
  dataForReport[["reportModules"]]=reportModules
  dataForReport[["genome"]]=genome
  dataForReport[["vafCol"]]=vafCol
  dataForReport[["vc_nonSyn"]]=vc_nonSyn
  dataForReport[["dndscv.refdb"]]=dndscv.refdb
  dataForReport[["interestedGenes"]]=interestedGenes
  dataForReport[["clinicalData"]]=clinicalData
  dataForReport[["clinicalFeatures"]]=clinicalFeatures
  dataForReport[["performGroupSummary"]]=performGroupSummary
  dataForReport[["showCode"]]=showCode

  # checkAndAddParameter=function(dataForReport,ParameterName,defaultValue) {
  #   if (is.null(ParameterName)) {
  #     dataForReport[[deparse(substitute(ParameterName))]]=defaultValue
  #   } else {
  #     dataForReport[[deparse(substitute(ParameterName))]]=reportOut
  #   }
  #   return(dataForReport)
  # }
  # dataForReport=checkAndAddParameter(dataForReport,reportOut,paste0(basename(dataForReport[["mafFile"]]),".report.html"))

  return(dataForReport)
}


#' @export
make_maf_report<-function(dataForReport) {
  dataForReport=dataForReport
  render(dataForReport[["reportTemplate"]],output_dir=dataForReport[["reportOutDir"]], output_file =dataForReport[["reportOut"]],
         intermediates_dir=tempdir())
  return(dataForReport)
}
