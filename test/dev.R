##################################################
#install mafreport package from github
##################################################
#library(devtools)
#install_github("slzhao/mafreport")


##################################################
#get OS type to determin if it is on local PC or ACCRE server
##################################################
osInfo=.Platform$OS.type
if (osInfo=="unix") { #ACCRE
  #load the latest version by devtools from ACCRE
  devtools::load_all("~/source/mafreport")
  reportOutDir="/scratch/cqs/zhaos/temp"
} else if (osInfo=="windows") { #local PC
  #load the latest version by devtools from local
  devtools::load_all("d:/source/mafreport")
  reportOutDir="D:\\temp\\"
} else {
  stop(paste0("Unknown OS type: ",osInfo))
}
#Or load the installed version
#library(mafreport)


##################################################
#Prepare. Four examples. Only need one of them
##################################################

#Simple Example. TCGA data from maftools package
mafFile=system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
clinicalData=system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
clinicalFeatures="FAB_classification"
#reportTemplate="d:/source/mafreport/inst/templates/report.Rmd" for test purpose
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
                                               vafCol="i_TumorVAF_WU",
                                               clinicalData=clinicalData,
                                               clinicalFeatures=clinicalFeatures,
                                               reportModules=c("Initialize.Rmd","SummaryTables.Rmd","VariantsVisualization.Rmd"))

#Example for GATK. Tim's Exome sequence project
mafFile="/scratch/cqs/baiy7/Tim_proj/IPF_new_WES/analysis/_IPF_new_WES_result_reseq_exLCR/bwa_refine_gatk4_SNV_09_toMAF/result/IPF_new_WES.freq0.001.filtered.tsv.maf"
clinicalData="/gpfs23/scratch/cqs/baiy7/Tim_proj/IPF_new_WES/analysis/_IPF_new_WES_result_reseq_exLCR/bwa_refine_gatk4_SNV_10_report/result/familyInfo.txt"
clinicalFeatures="Family"
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
                                               vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                                                           "Nonsense_Mutation",	"Nonstop_Mutation", "Missense_Mutation"),
                                               interestedGenes=c("TERT","RTEL1","DKC1","PARN"),
                                               clinicalData=clinicalData,
                                               clinicalFeatures=clinicalFeatures,
                                               reportModules=c("Initialize.Rmd","SummaryTables.Rmd","GroupSummary.Rmd")
)


#Example for GATK. Tiger's SNV project
mafFile="D:\\OneDriveWork\\OneDriveVanderbilt\\work\\software_QC\\test\\linton_exomeseq_2118.maf_filtered.vcf.maf.subset.maf"
clinicalData="D:\\OneDriveWork\\OneDriveVanderbilt\\work\\software_QC\\test\\linton_exomeseq_2118.pass.family.txt"
clinicalFeatures="family"
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
                                               vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                                                           "Nonsense_Mutation",	"Nonstop_Mutation", "Missense_Mutation"),
                                               interestedGenes=c("LDLR","APOB","PCSK9","LDLRAP1","STAP1","LIPA","ABCG5","ABCGB","APOE","LPA","PNPLA5","CH25H","INSIG2","SIRT1"),
                                               clinicalData=clinicalData,
                                               clinicalFeatures=clinicalFeatures
                                               #		reportModules=c("Initialize.Rmd","SummaryTables.Rmd","GroupSummary.Rmd")
)

#Example for GATK. Tiger's SNV project in ACCRE
mafFile="/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/bwa_refine_gatk4_hc_gvcf_vqsr_filterMAF_annovar_filter_toMAF/result/linton_exomeseq_3321.freq0.001.filtered.tsv.maf"
clinicalData="/scratch/cqs/shengq2/macrae_linton/20180913_linton_exomeseq_2118_human_cutadapt/linton_exomeseq_2118.pass.family.txt"
clinicalFeatures="family"
interestedGeneStr=c("LDLR","APOB","PCSK9","LDLRAP1","STAP1","LIPA","ABCG5","ABCGB","APOE","LPA","PNPLA5","CH25H","INSIG2","SIRT1")
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
                                              vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                                                          "Nonsense_Mutation",	"Nonstop_Mutation", "Missense_Mutation"),
                                              interestedGenes=interestedGenes
#                                              clinicalData=clinicalData,
#                                              clinicalFeatures=clinicalFeatures,
#                                              reportModules=c("Initialize.Rmd","SummaryTables.Rmd","GroupSummary.Rmd")
)

#Example for mouse. Jing's WES somatic mutation project
mafFile="/data/CRC_expression/JingYang/result/2221/bwa_refine_muTect_combined_annovar_toMAF/result/2221_mm10_pass.combined.vcf.annovar.final.tsv.maf"
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
                                               genome="mm10",dndscv.refdb="/scratch/cqs/references/dndscv/RefCDS_mouse_GRCm38.p2.rda")

##Example for MuTect2 and hg38. Pierre's WES project
mafFile="D:\\OneDriveWork\\OneDriveVanderbilt\\work\\Pierre\\20180830_WESTestResult\\20181206_LungWesVcfToMaf.all.subset.maf"
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,genome="hg38",
                                               #		reportModules=c("Initialize.Rmd","HighFreqVariantsCheck.Rmd"),
                                               dndscv.refdb="D:\\OneDriveWork\\OneDriveVanderbilt\\work\\software_QC\\test\\RefCDS_human_GRCh38.p12.rda",
                                               clinicalData="D:\\OneDriveWork\\OneDriveVanderbilt\\work\\Pierre\\20180830_WESTestResult\\wes_sample_grp.txt",
                                               clinicalFeatures=c("PatientID","SampleGroup")
                                               #		reportTemplate="d:/source/mafreport/inst/templates/report.Rmd"
)

##Example for JarVarDict and hg38. Pierre's WES project
mafFile="/gpfs23/scratch/cqs/zhaos/Pierre/WES/20190509_VarDict_vcf2maf/20190522_LungWesUmiVcfToMaf.all.subset.noError.maf"
reportOutDir="/gpfs23/scratch/cqs/zhaos/Pierre/WES/20190509_VarDict_vcf2maf/"
dataForReport=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,genome="hg38",
                                               #reportModules=c("Initialize.Rmd","SummaryTables.Rmd","VariantsVisualization.Rmd","HighFreqVariantsCheck.Rmd"),
                                               dndscv.refdb="/scratch/cqs/references/dndscv/RefCDS_mouse_GRCm38.p2.rda"
                                               #                                               clinicalData="/gpfs23/scratch/cqs/zhaos/Pierre/WES/wes_sample_grp.txt",
                                               #                                               clinicalFeatures=c("PatientID","SampleGroup")
                                               #		reportTemplate="d:/source/mafreport/inst/templates/report.Rmd"
)

##################################################
#Make report
##################################################
dataForReport=make_maf_report(dataForReport)


##################################################
#Result
##################################################
#Html Report: *.maf.report.html
#RData to reproduce Report: *.maf.report.html.RData





