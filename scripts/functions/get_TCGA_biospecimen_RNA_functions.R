
## Loading R data objects
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

################################
## FUNCTIONS TO RETRIEVE DATA ##
################################

get_biospecimen_tcga.legacy <- function(selectedcancer){
  print(paste0("Getting biospecimen data for ", selectedcancer ))
  query <- GDCquery(project = selectedcancer, 
                    data.category = "Biospecimen",
                    file.type = "xml", legacy=TRUE
  )
  if (is.null(query)){
    print(paste0("No biospecimen legacy data for ", selectedcancer))
    message(paste0("No biospecimen legacy data for ", selectedcancer))
  }else{
    #df_clin_data <- query %>% getResults
    GDCdownload(query,files.per.chunk = 50)
    biospecimen1 <- distinct(GDCprepare_clinic(query,clinical.info = c("sample")))#"bcr_patient_barcode",sample_type_id,sample_type,project,bcr_sample_barcode   
    biospecimen1_sub <- biospecimen1[,c("bcr_patient_barcode","sample_type_id","sample_type","project","bcr_sample_barcode")]
    biospecimen2 <- distinct(GDCprepare_clinic(query,clinical.info = c("aliquot")))
    biospecimen2_sub <- biospecimen2[,c("bcr_patient_barcode","bcr_aliquot_barcode")]
    biospecimen2_sub$bcr_sample_barcode <-substr(biospecimen2_sub$bcr_aliquot_barcode,1,16)
    #bcr_patient_barcode,bcr_aliquot_barcode
    # biospecimen3 <- GDCprepare_clinic(query,clinical.info = c("bio_patient"))
    biospecimen4 <- distinct(GDCprepare_clinic(query,clinical.info = c("analyte"))) 
    biospecimen4_sub <- biospecimen4[,c("analyte_type_id","bcr_analyte_barcode")]
    biospecimen4_sub$bcr_patient_barcode <- substr(biospecimen4_sub$bcr_analyte_barcode,1,12)
    biospecimen4_sub$bcr_sample_barcode <- substr(biospecimen4_sub$bcr_analyte_barcode,1,16)
    # analyte_type_id,bcr_analyte_barcode
    #Merge 1 and 2 on sample barcode
    biospecimen_12 <- merge(biospecimen1_sub,biospecimen2_sub,by ="bcr_sample_barcode", all=TRUE)
    # Add analyte type
    biospecimen_12$bcr_analyte_barcode <- substr(biospecimen_12$bcr_aliquot_barcode,1,20)
    
    biospecimen_12$analyte_type_id <- substr(biospecimen_12$bcr_aliquot_barcode,20,20)
    
    # Merge 12 with 4 to have also analyte type
    biospecimen_all <- merge(biospecimen_12,biospecimen4_sub,by ="bcr_analyte_barcode", all=TRUE)
    # cleanup duplicate columns
    biospecimen_all_final <- biospecimen_all[,c("bcr_sample_barcode.x","sample_type_id","sample_type","project","bcr_aliquot_barcode","analyte_type_id.x","bcr_analyte_barcode","bcr_patient_barcode.x")]
    # 
    colnames(biospecimen_all_final)[c(1,6,8)] <-c("bcr_sample_barcode","analyte_type_id", "bcr_patient_barcode") 
    
    biospecimen_tcga.legacy  <- rbind(biospecimen_tcga.legacy ,biospecimen_all_final)
    biospecimen_tcga.legacy 
  }
}

get_biospecimen_tcga.nonlegacy <- function(selectedcancer){
  print(paste0("Getting biospecimen data for ", selectedcancer ))
  query <- GDCquery(project = selectedcancer, 
                    data.category = "Biospecimen",
                    file.type = "xml"
  )
  if (is.null(query)){
    print(paste0("No biospecimen data for ", selectedcancer))
  }else{
    #df_clin_data <- query %>% getResults
    GDCdownload(query)
    biospecimen1 <- distinct(GDCprepare_clinic(query,clinical.info = c("sample")))#"bcr_patient_barcode",sample_type_id,sample_type,project,bcr_sample_barcode   
    biospecimen1_sub <- biospecimen1[,c("bcr_patient_barcode","sample_type_id","sample_type","project","bcr_sample_barcode")]
    biospecimen2 <- distinct(GDCprepare_clinic(query,clinical.info = c("aliquot")))
    biospecimen2_sub <- biospecimen2[,c("bcr_patient_barcode","bcr_aliquot_barcode")]
    biospecimen2_sub$bcr_sample_barcode <-substr(biospecimen2_sub$bcr_aliquot_barcode,1,16)
    #bcr_patient_barcode,bcr_aliquot_barcode
    # biospecimen3 <- GDCprepare_clinic(query,clinical.info = c("bio_patient"))
    biospecimen4 <- distinct(GDCprepare_clinic(query,clinical.info = c("analyte"))) 
    biospecimen4_sub <- biospecimen4[,c("analyte_type_id","bcr_analyte_barcode")]
    biospecimen4_sub$bcr_patient_barcode <- substr(biospecimen4_sub$bcr_analyte_barcode,1,12)
    biospecimen4_sub$bcr_sample_barcode <- substr(biospecimen4_sub$bcr_analyte_barcode,1,16)
    # analyte_type_id,bcr_analyte_barcode
    #Merge 1 and 2 on sample barcode
    biospecimen_12 <- merge(biospecimen1_sub,biospecimen2_sub,by ="bcr_sample_barcode", all=TRUE)
    # Add analyte type
    biospecimen_12$bcr_analyte_barcode <- substr(biospecimen_12$bcr_aliquot_barcode,1,20)
    
    biospecimen_12$analyte_type_id <- substr(biospecimen_12$bcr_aliquot_barcode,20,20)
    
    # Merge 12 with 4 to have also analyte type
    biospecimen_all <- merge(biospecimen_12,biospecimen4_sub,by ="bcr_analyte_barcode", all=TRUE)
    # cleanup duplicate columns
    biospecimen_all_final <- biospecimen_all[,c("bcr_sample_barcode.x","sample_type_id","sample_type","project","bcr_aliquot_barcode","analyte_type_id.x","bcr_analyte_barcode","bcr_patient_barcode.x")]
    # 
    colnames(biospecimen_all_final)[c(1,6,8)] <-c("bcr_sample_barcode","analyte_type_id", "bcr_patient_barcode") 
    
    biospecimen_tcga.nonlegacy  <- rbind(biospecimen_tcga.nonlegacy ,biospecimen_all_final)
    biospecimen_tcga.nonlegacy 
  }
}

##################
## PROCESS DATA ##
##################
process_biospecimenData <- function(biospecimenList, msiDF){
  # Merge data
  biospecimen.tcga  <-bind_rows(biospecimenList , .id = "column_label")
  # Extract Primary tumor and Normal Solid Tissue samples: looking to assign as tumor (TP) with 01 code and normal (NT) witth code 11
  biospecimen.tcga.tpnt  <-biospecimen.tcga %>% dplyr::filter(sample_type_id %in% c(1,11)) %>% distinct()
  biospecimen.tcga.tpnt$sample_type <- ifelse(biospecimen.tcga.tpnt$sample_type_id==1,"TP","NT")
  # Add msi Information
  biospecimen.tcga.tpnt$MSI_status <-msiDF$mononucleotide_and_dinucleotide_marker_panel_analysis_status[match(biospecimen.tcga.tpnt$bcr_patient_barcode,msiDF$bcr_patient_barcode)]
  biospecimen.tcga.tpnt
}

process_biospecimenData_TM.NT <- function(biospecimenList, msiDF){
  # biospecimenList <-biospecimen.legacy
  # Merge data
  biospecimen.tcga  <-bind_rows(biospecimenList , .id = "column_label")
  # Extract Primary tumor and Normal Solid Tissue samples: looking to assign as tumor (TP) with 01 code and normal (NT) witth code 11
  biospecimen.tcga.tmnt  <-biospecimen.tcga %>% dplyr::filter(sample_type_id %in% c(6,11)) %>% distinct()
  biospecimen.tcga.tmnt$sample_type <- ifelse(biospecimen.tcga.tmnt$sample_type_id==6,"TM","NT")
  # Add msi Information
  biospecimen.tcga.tmnt$MSI_status <-msiDF$mononucleotide_and_dinucleotide_marker_panel_analysis_status[match(biospecimen.tcga.tmnt$bcr_patient_barcode,msiDF$bcr_patient_barcode)]
  biospecimen.tcga.tmnt
}

add_brcaSubtypes <- function(biospecimenDF,molSubtypesDF){
  # First subset to brca
  tcga_brca <- biospecimenDF[biospecimenDF$type=='BRCA',]
  tcga_brca$PatientID <- substr(tcga_brca$bcr_aliquot_barcode,1,12)
  tcga_brca.types <-merge(tcga_brca,molSubtypesDF[,c("PatientID","Subtype_mRNA")], by = "PatientID")
  
  tcga_brca.types$type <- tcga_brca.types$Subtype_mRNA 
  
  tcga_brca.types$type <- ifelse(grepl("LumA",tcga_brca.types$type),"BRCA.LuminalA",
                                 ifelse(grepl("LumB",tcga_brca.types$type),"BRCA.LuminalB",
                                        ifelse(grepl("Basal",tcga_brca.types$type),"BRCA.TripleNegative",
                                               ifelse(grepl("Her2",tcga_brca.types$type),"BRCA.Her2",
                                                      ifelse(grepl("Normal",tcga_brca.types$type),"BRCA.Normal",NA)))))
  # Add extra brca subtypes table to main table excluding the subtype column
  biospecimenDF.f <- rbind(biospecimenDF,tcga_brca.types[,-c(1,dim(tcga_brca.types)[2])])
  biospecimenDF.f
}

add_coadSubtypes <- function(biospecimenDF){
  # First subset to COAD
  tcga_coad <- biospecimenDF[biospecimenDF$type=='COAD',]
  tcga_coad$type <- ifelse(grepl("MSI-H",tcga_coad$MSI_status),"COAD.MSI_High",
                           ifelse(grepl("MSI-L",tcga_coad$MSI_status),"COAD.MSI_Low",
                                  ifelse(grepl("MSS",tcga_coad$MSI_status),"COAD.MSS",NA)))
  # Remove NAs
  tcga_coad <- tcga_coad[!is.na(tcga_coad$type),]
  #
  # Add extra brca subtypes table to main table excluding the subtype column
  biospecimenDF.f <- rbind(biospecimenDF,tcga_coad)
  biospecimenDF.f
}


add_brcaSubtypesFilt <- function(biospecimenDFpre,biospecimenDFpost,molSubtypesDF){
  # First subset to brca
  tcga_brca <- biospecimenDFpre[biospecimenDFpre$type=='BRCA',]
  tcga_brca$PatientID <- substr(tcga_brca$bcr_aliquot_barcode,1,12)
  tcga_brca.types <-merge(tcga_brca,molSubtypesDF[,c("PatientID","Subtype_mRNA")], by = "PatientID")
  
  tcga_brca.types$type <- tcga_brca.types$Subtype_mRNA 
  
  tcga_brca.types$type <- ifelse(grepl("LumA",tcga_brca.types$type),"BRCA.LuminalA",
                                 ifelse(grepl("LumB",tcga_brca.types$type),"BRCA.LuminalB",
                                        ifelse(grepl("Basal",tcga_brca.types$type),"BRCA.TripleNegative",
                                               ifelse(grepl("Her2",tcga_brca.types$type),"BRCA.Her2",
                                                      ifelse(grepl("Normal",tcga_brca.types$type),"BRCA.Normal",NA)))))
  # Filter barcodes
  tcga_brca.types_barcodes <- tcga_replicateFilter(tcga_brca.types$bcr_aliquot_barcode , analyte_target = "RNA",filter_FFPE=TRUE, full_barcode=TRUE)
  # Subset to filtered
  tcga_brca.types.f <- subset(tcga_brca.types,tcga_brca.types$bcr_aliquot_barcode %in% tcga_brca.types_barcodes)
  # Add extra brca subtypes table to main table excluding the subtype column
  biospecimenDF.f <- rbind(biospecimenDFpost,tcga_brca.types.f[,-c(1,dim(tcga_brca.types)[2])])
  biospecimenDF.f
}
add_coadSubtypesFilt <- function(biospecimenDFpre,biospecimenDFpost){
  # First subset to COAD
  tcga_coad <- biospecimenDFpre[biospecimenDFpre$type=='COAD',]
  tcga_coad$type <- ifelse(grepl("MSI-H",tcga_coad$MSI_status),"COAD.MSI_High",
                           ifelse(grepl("MSI-L",tcga_coad$MSI_status),"COAD.MSI_Low",
                                  ifelse(grepl("MSS",tcga_coad$MSI_status),"COAD.MSS",NA)))
  # Remove NAs
  tcga_coad <- tcga_coad[!is.na(tcga_coad$type),]
  # Filter barcodes
  tcga_coad_barcodes <- tcga_replicateFilter(tcga_coad$bcr_aliquot_barcode , analyte_target = "RNA",filter_FFPE=TRUE, full_barcode=TRUE)
  # Subset to filtered
  tcga_coad.f <- subset(tcga_coad,tcga_coad$bcr_aliquot_barcode %in% tcga_coad_barcodes)
  # Add extra brca subtypes table to main table excluding the subtype column
  biospecimenDF.f <- rbind(biospecimenDFpost,tcga_coad.f)
  biospecimenDF.f
}
