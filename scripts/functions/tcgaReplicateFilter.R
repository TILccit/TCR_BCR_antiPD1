
######################################################
## TCGA Replicate filter function -                 ##
## input is a vector of TCGA barcodes,              ##
## based on analyte target, it selects RNA OR DNA - ##
## FUNCTION CHANGED COMPARED TO ORIGINAL            ##
######################################################
##TCGA Replicate filter function - input is a vector of TCGA barcodes, based on analyte target, it selects RNA OR DNA - FUNCITON CHANGED COMPARED TO ORIGINAL
tcga_replicateFilter = function(tsb, analyte_target=c("DNA","RNA"), decreasing=TRUE, analyte_position=20, plate=c(22,25), portion=c(18,19), filter_FFPE=FALSE, full_barcode=FALSE){
  # basically, user provide tsb and analyte_target is fine. If you
  # want to filter FFPE samples, please set filter_FFPE and full_barcode
  # all to TRUE, and tsb must have nchar of 28
  
  analyte_target = match.arg(analyte_target)
  print(analyte_target)
  # Strings in R are largely lexicographic
  # see ??base::Comparison
  
  # filter FFPE samples
  # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html> 
  if(full_barcode & filter_FFPE){
    ffpe = c("TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01", 
             "TCGA-44-2656-01B-06D-A276-05", "TCGA-44-2656-01B-06D-A27C-26", 
             "TCGA-44-2656-01B-06R-A277-07", "TCGA-44-2662-01B-02D-A271-08", 
             "TCGA-44-2662-01B-02D-A273-01", "TCGA-44-2662-01B-02R-A277-07", 
             "TCGA-44-2665-01B-06D-A271-08", "TCGA-44-2665-01B-06D-A273-01", 
             "TCGA-44-2665-01B-06D-A276-05", "TCGA-44-2665-01B-06R-A277-07", 
             "TCGA-44-2666-01B-02D-A271-08", "TCGA-44-2666-01B-02D-A273-01", 
             "TCGA-44-2666-01B-02D-A276-05", "TCGA-44-2666-01B-02D-A27C-26", 
             "TCGA-44-2666-01B-02R-A277-07", "TCGA-44-2668-01B-02D-A271-08", 
             "TCGA-44-2668-01B-02D-A273-01", "TCGA-44-2668-01B-02D-A276-05", 
             "TCGA-44-2668-01B-02D-A27C-26", "TCGA-44-2668-01B-02R-A277-07", 
             "TCGA-44-3917-01B-02D-A271-08", "TCGA-44-3917-01B-02D-A273-01", 
             "TCGA-44-3917-01B-02D-A276-05", "TCGA-44-3917-01B-02D-A27C-26", 
             "TCGA-44-3917-01B-02R-A277-07", "TCGA-44-3918-01B-02D-A271-08", 
             "TCGA-44-3918-01B-02D-A273-01", "TCGA-44-3918-01B-02D-A276-05", 
             "TCGA-44-3918-01B-02D-A27C-26", "TCGA-44-3918-01B-02R-A277-07", 
             "TCGA-44-4112-01B-06D-A271-08", "TCGA-44-4112-01B-06D-A273-01", 
             "TCGA-44-4112-01B-06D-A276-05", "TCGA-44-4112-01B-06D-A27C-26", 
             "TCGA-44-4112-01B-06R-A277-07", "TCGA-44-5645-01B-04D-A271-08", 
             "TCGA-44-5645-01B-04D-A273-01", "TCGA-44-5645-01B-04D-A276-05", 
             "TCGA-44-5645-01B-04D-A27C-26", "TCGA-44-5645-01B-04R-A277-07", 
             "TCGA-44-6146-01B-04D-A271-08", "TCGA-44-6146-01B-04D-A273-01", 
             "TCGA-44-6146-01B-04D-A276-05", "TCGA-44-6146-01B-04D-A27C-26", 
             "TCGA-44-6146-01B-04R-A277-07", "TCGA-44-6146-01B-04R-A27D-13", 
             "TCGA-44-6147-01B-06D-A271-08", "TCGA-44-6147-01B-06D-A273-01", 
             "TCGA-44-6147-01B-06D-A276-05", "TCGA-44-6147-01B-06D-A27C-26", 
             "TCGA-44-6147-01B-06R-A277-07", "TCGA-44-6147-01B-06R-A27D-13", 
             "TCGA-44-6775-01C-02D-A271-08", "TCGA-44-6775-01C-02D-A273-01", 
             "TCGA-44-6775-01C-02D-A276-05", "TCGA-44-6775-01C-02D-A27C-26", 
             "TCGA-44-6775-01C-02R-A277-07", "TCGA-44-6775-01C-02R-A27D-13", 
             "TCGA-A6-2674-01B-04D-A270-10", "TCGA-A6-2674-01B-04R-A277-07", 
             "TCGA-A6-2677-01B-02D-A270-10", "TCGA-A6-2677-01B-02D-A274-01", 
             "TCGA-A6-2677-01B-02D-A27A-05", "TCGA-A6-2677-01B-02D-A27E-26", 
             "TCGA-A6-2677-01B-02R-A277-07", "TCGA-A6-2684-01C-08D-A270-10", 
             "TCGA-A6-2684-01C-08D-A274-01", "TCGA-A6-2684-01C-08D-A27A-05", 
             "TCGA-A6-2684-01C-08D-A27E-26", "TCGA-A6-2684-01C-08R-A277-07", 
             "TCGA-A6-3809-01B-04D-A270-10", "TCGA-A6-3809-01B-04D-A274-01", 
             "TCGA-A6-3809-01B-04D-A27A-05", "TCGA-A6-3809-01B-04D-A27E-26", 
             "TCGA-A6-3809-01B-04R-A277-07", "TCGA-A6-3810-01B-04D-A270-10", 
             "TCGA-A6-3810-01B-04D-A274-01", "TCGA-A6-3810-01B-04D-A27A-05", 
             "TCGA-A6-3810-01B-04D-A27E-26", "TCGA-A6-3810-01B-04R-A277-07", 
             "TCGA-A6-5656-01B-02D-A270-10", "TCGA-A6-5656-01B-02D-A274-01", 
             "TCGA-A6-5656-01B-02D-A27A-05", "TCGA-A6-5656-01B-02D-A27E-26", 
             "TCGA-A6-5656-01B-02R-A277-07", "TCGA-A6-5656-01B-02R-A27D-13", 
             "TCGA-A6-5659-01B-04D-A270-10", "TCGA-A6-5659-01B-04D-A274-01", 
             "TCGA-A6-5659-01B-04D-A27A-05", "TCGA-A6-5659-01B-04D-A27E-26", 
             "TCGA-A6-5659-01B-04R-A277-07", "TCGA-A6-6650-01B-02D-A270-10", 
             "TCGA-A6-6650-01B-02D-A274-01", "TCGA-A6-6650-01B-02D-A27A-05", 
             "TCGA-A6-6650-01B-02D-A27E-26", "TCGA-A6-6650-01B-02R-A277-07", 
             "TCGA-A6-6650-01B-02R-A27D-13", "TCGA-A6-6780-01B-04D-A270-10", 
             "TCGA-A6-6780-01B-04D-A274-01", "TCGA-A6-6780-01B-04D-A27A-05", 
             "TCGA-A6-6780-01B-04D-A27E-26", "TCGA-A6-6780-01B-04R-A277-07", 
             "TCGA-A6-6780-01B-04R-A27D-13", "TCGA-A6-6781-01B-06D-A270-10", 
             "TCGA-A6-6781-01B-06D-A274-01", "TCGA-A6-6781-01B-06D-A27A-05", 
             "TCGA-A6-6781-01B-06R-A277-07", "TCGA-A6-6781-01B-06R-A27D-13", 
             "TCGA-A7-A0DB-01C-02D-A272-09", "TCGA-A7-A0DB-01C-02R-A277-07", 
             "TCGA-A7-A0DB-01C-02R-A27D-13", "TCGA-A7-A13D-01B-04D-A272-09", 
             "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A13D-01B-04R-A27D-13", 
             "TCGA-A7-A13E-01B-06D-A272-09", "TCGA-A7-A13E-01B-06R-A277-07", 
             "TCGA-A7-A13E-01B-06R-A27D-13", "TCGA-A7-A26E-01B-06D-A272-09", 
             "TCGA-A7-A26E-01B-06D-A275-01", "TCGA-A7-A26E-01B-06D-A27B-05", 
             "TCGA-A7-A26E-01B-06R-A277-07", "TCGA-A7-A26E-01B-06R-A27D-13", 
             "TCGA-A7-A26J-01B-02D-A272-09", "TCGA-A7-A26J-01B-02D-A275-01", 
             "TCGA-A7-A26J-01B-02D-A27B-05", "TCGA-A7-A26J-01B-02D-A27F-26", 
             "TCGA-A7-A26J-01B-02R-A277-07", "TCGA-A7-A26J-01B-02R-A27D-13", 
             "TCGA-B2-3923-01B-10D-A270-10", "TCGA-B2-3923-01B-10R-A277-07", 
             "TCGA-B2-3923-01B-10R-A27D-13", "TCGA-B2-3924-01B-03D-A270-10", 
             "TCGA-B2-3924-01B-03D-A274-01", "TCGA-B2-3924-01B-03D-A27A-05", 
             "TCGA-B2-3924-01B-03D-A27E-26", "TCGA-B2-3924-01B-03R-A277-07", 
             "TCGA-B2-3924-01B-03R-A27D-13", "TCGA-B2-5633-01B-04D-A270-10", 
             "TCGA-B2-5633-01B-04D-A274-01", "TCGA-B2-5633-01B-04D-A27A-05", 
             "TCGA-B2-5633-01B-04D-A27E-26", "TCGA-B2-5633-01B-04R-A277-07", 
             "TCGA-B2-5633-01B-04R-A27D-13", "TCGA-B2-5635-01B-04D-A270-10", 
             "TCGA-B2-5635-01B-04D-A274-01", "TCGA-B2-5635-01B-04D-A27A-05", 
             "TCGA-B2-5635-01B-04D-A27E-26", "TCGA-B2-5635-01B-04R-A277-07", 
             "TCGA-B2-5635-01B-04R-A27D-13", "TCGA-BK-A0CA-01B-02D-A272-09", 
             "TCGA-BK-A0CA-01B-02D-A275-01", "TCGA-BK-A0CA-01B-02D-A27B-05", 
             "TCGA-BK-A0CA-01B-02D-A27F-26", "TCGA-BK-A0CA-01B-02R-A277-07", 
             "TCGA-BK-A0CA-01B-02R-A27D-13", "TCGA-BK-A0CC-01B-04D-A272-09", 
             "TCGA-BK-A0CC-01B-04D-A275-01", "TCGA-BK-A0CC-01B-04D-A27B-05", 
             "TCGA-BK-A0CC-01B-04R-A277-07", "TCGA-BK-A0CC-01B-04R-A27D-13", 
             "TCGA-BK-A139-01C-08D-A272-09", "TCGA-BK-A139-01C-08D-A275-01", 
             "TCGA-BK-A139-01C-08D-A27B-05", "TCGA-BK-A139-01C-08D-A27F-26", 
             "TCGA-BK-A139-01C-08R-A277-07", "TCGA-BK-A139-01C-08R-A27D-13", 
             "TCGA-BK-A26L-01C-04D-A272-09", "TCGA-BK-A26L-01C-04D-A275-01", 
             "TCGA-BK-A26L-01C-04D-A27B-05", "TCGA-BK-A26L-01C-04D-A27F-26", 
             "TCGA-BK-A26L-01C-04R-A277-07", "TCGA-BK-A26L-01C-04R-A27D-13", 
             "TCGA-BL-A0C8-01B-04D-A271-08", "TCGA-BL-A0C8-01B-04D-A273-01", 
             "TCGA-BL-A0C8-01B-04D-A276-05", "TCGA-BL-A0C8-01B-04D-A27C-26", 
             "TCGA-BL-A0C8-01B-04R-A277-07", "TCGA-BL-A0C8-01B-04R-A27D-13", 
             "TCGA-BL-A13I-01B-04D-A271-08", "TCGA-BL-A13I-01B-04D-A276-05", 
             "TCGA-BL-A13I-01B-04R-A277-07", "TCGA-BL-A13I-01B-04R-A27D-13", 
             "TCGA-BL-A13J-01B-04D-A271-08", "TCGA-BL-A13J-01B-04D-A273-01", 
             "TCGA-BL-A13J-01B-04D-A276-05", "TCGA-BL-A13J-01B-04D-A27C-26", 
             "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13")
    
    tsb = setdiff(tsb, tsb[which(tsb %in% ffpe)])
  }
  # Get only samples for analyte in question
  if(analyte_target == "DNA"){
    tsb_dna = c()
    for(x in tsb){
      barcode = x
      analyte = substr(barcode, 
                       start = analyte_position,
                       stop = analyte_position)
      if(analyte %in% c("D","X","G","W")){
        tsb_dna = c(tsb_dna, barcode)
      }
    }
    tsb <-tsb_dna
  }else{
    tsb_rna = c()
    for(x in tsb){
      barcode = x
      analyte = substr(barcode, 
                       start = analyte_position,
                       stop = analyte_position)
      if(analyte %in% c("R","T","H")){
        tsb_rna = c(tsb_rna, barcode)
      }
    }
    tsb <-tsb_rna
  }
  #tsb
  ## find repeated samples
  # Extract sample IDs
  sampleID = substr(tsb, start = 1, stop = 15)
  # Duplicated samples
  dp_samples = unique(sampleID)
  # Not Duplicated samples
  un_samples = setdiff(sampleID,dp_samples)
  ################################
  ## NO DUPLICATES
  ################################
  if(length(dp_samples)==0){
    uniq_tsb = tsb[! sampleID %in% dp_samples]
    
    add_tsb = c()
    # Have not found duplicated values, but check if input is analyte chosen, so filter out whatever is not analyte chosen
    print("ooo Not find any duplicated barcodes, return original input..")
    
  }else{
    ################################
    ## DUPLICATES FOUND
    ################################
    # Unique, not duplicated aliquote barcodes
    uniq_tsb = tsb[! sampleID %in% dp_samples]
    # Duplicated aliquot barcodes-belonging to the same sample & analyte
    dp_tsb = setdiff(tsb, uniq_tsb)
    
    add_tsb = c()
    
    # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
    # if analyte_target = "DNA"
    # analyte:  D > G,W,X
    if(analyte_target == "DNA"){
      uniq_tsb_dna = c()
      print("Grabbing DNA...")
      # Loop through samples having multiple aliquots
      for(x in dp_samples){
        # find the multiple aliquots for given sample
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        # find analytes for multiple aliquotes of given sample
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        # D is preferred over G,W,X
        if(any(analytes == "D") & !(all(analytes == "D"))){
          aliquot = mulaliquots[which(analytes == "D")]
          add_tsb = c(add_tsb, aliquot)
          dp_tsb = setdiff(dp_tsb, mulaliquots)
        }
        
      }
      # Filter duplicate aliquots o be sure they are DNA
      dp_tsb_dna = c()
      # Loop through aliquots
      for(x in dp_tsb){
        # Extract analyte barcode
        barcode = substr(x,1,20)
        analyte = substr(barcode, 
                         start = analyte_position,
                         stop = analyte_position)
        # Check if analyte is DNA
        if(analyte %in% c("D","X","G","W")){
          dp_tsb_dna = c(dp_tsb_dna, barcode)
        }
      }
      # Filter the uniq barcodes to select only for analyte chosen
      # print(un_samples)
      # print(uniq_tsb)
      for(x in uniq_tsb){
        # Extract analyte barcode
        barcode = substr(x,1,20)
        analyte = substr(barcode, 
                         start = analyte_position,
                         stop = analyte_position)
        if(analyte %in% c("D","X","G","W")){
          uniq_tsb_dna = c(uniq_tsb_dna, barcode)
        }
      }
      dp_tsb <- dp_tsb_dna
    }else{
      print("Grabbing RNA...")
      # Filter duplicate and unique aliquots to be sure they are RNA
      dp_tsb_rna = c()
      uniq_tsb_rna = c()
      # Duplicate aliquots
      for(x in dp_tsb){
        barcode = x
        analyte = substr(barcode, 
                         start = analyte_position,
                         stop = analyte_position)
        if(analyte %in% c("R","T","H")){
          dp_tsb_rna = c(dp_tsb_rna, barcode)
        }
      }
      # Filter unique for RNA
      #print(un_samples)
      for(x in uniq_tsb){
        barcode = x
        analyte = substr(barcode, 
                         start = analyte_position,
                         stop = analyte_position)
        if(analyte %in% c("R","T","H")){
          uniq_tsb_rna = c(uniq_tsb_rna, barcode)
        }
      }
      dp_tsb <- dp_tsb_rna
      # Loop through duplicate samples
      for(x in dp_samples){
        # Get multiple aliquots for given sample
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        # Get analytes for the multiple aliquots of given sample
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        # H is preferred to R and T
        if(any(analytes == "H") & !(all(analytes == "H"))){
          aliquot = mulaliquots[which(analytes == "H")]
          if (length(which(analytes == "H"))>1){# If more than one analyte ==H
            # Remove the analyte(s) not R from duplicates
            dp_tsb = setdiff(dp_tsb,mulaliquots[which(analytes != "H")])
          }else{ #If only one analyte H
            # Add analyte to unique analytes
            add_tsb = c(add_tsb, aliquot)
            # Remove all mulaliquots from duplicates
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }
        }else if(any(analytes == "R") & !(all(analytes == "R"))){
          aliquot = mulaliquots[which(analytes == "R")]
          if (length(which(analytes == "R"))>1){# If more than one analyte ==R
            # Remove the analyte(s) not R from duplicates
            dp_tsb = setdiff(dp_tsb,mulaliquots[which(analytes != "R")])
          }else{ #If only one analyte R
            # Add analyte to unique analytes
            add_tsb = c(add_tsb, aliquot)
            # Remove all mulaliquots from duplicates
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }
          
        }else if(any(analytes == "T") & !(all(analytes == "T"))){
          aliquot = mulaliquots[which(analytes == "T")]
          if (length(which(analytes == "T"))>1){# If more than one analyte ==R
            # Remove the analyte(s) not R from duplicates
            dp_tsb = setdiff(dp_tsb,mulaliquots[which(analytes != "T")])
          }else{ #If only one analyte R
            # Add analyte to unique analytes
            add_tsb = c(add_tsb, aliquot)
            # Remove all mulaliquots from duplicates
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }
        }
      }
      
    }
    # if analyte_target = "RNA"
    # analyte: H > R > T 
    # else{
    #     
    # }
    # if analyte_target = "RNA"
    # analyte: H > R > T 
    # else{
    #     
    # }
    
    if(length(dp_tsb) == 0){
      print("ooo Filter barcodes successfully!")
      if(analyte_target == "DNA"){
        #print("Show DNA barcodes")
        uniq_tsb <- uniq_tsb_dna
      }else{
        #print("Show RNA barcodes")
        uniq_tsb <- uniq_tsb_rna
        
      }
      c(uniq_tsb, add_tsb)
    }else{
      print("Filter according to portion number")
      # filter according to portion number
      sampleID_res = substr(dp_tsb, start=1, stop=15)
      dp_samples_res = unique(sampleID_res)
      
      for(x in dp_samples_res){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        portion_codes = substr(mulaliquots,
                               start = portion[1],
                               stop = portion[2])
        portion_keep = sort(portion_codes, decreasing = decreasing)[1]
        # If not all aliquots have the max portion code
        if(!all(portion_codes == portion_keep)){
          # If there is one aliquote matching the max portion code
          if(length(which(portion_codes == portion_keep)) == 1){
            # Add this aliquote to the unique aliquots
            add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
            # Unique best aliquot selected so Remove all multiple aliquots of this sample from duplicates
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }else{# If there are more than one aliquot with the max portion code
            # Remove aliquots that do not have the max portion code from duplicates
            dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)])
          }
          
        }
      }
      
      if(length(dp_tsb)==0){
        print("ooo Filter barcodes successfully!")
        c(uniq_tsb, add_tsb)
      }else{
        print("Filter according to plate number")
        # filter according to plate number
        sampleID_res = substr(dp_tsb, start=1, stop=15)
        dp_samples_res = unique(sampleID_res)
        for(x in dp_samples_res){
          mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
          # if ("TCGA-AA-3663-11A-01R-1723-07" %in% mulaliquots){
          #   print(mulaliquots)
          # }
          plate_codes = substr(mulaliquots,
                               start = plate[1],
                               stop = plate[2])
          plate_keep = sort(plate_codes, decreasing = decreasing)[1]
          add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)])
          dp_tsb = setdiff(dp_tsb, mulaliquots)
          
          
        }
        
        if(length(dp_tsb)==0){
          print("ooo No more duplicates,filter barcodes successfully!")
          c(uniq_tsb, add_tsb)
        }else{
          print(paste0("ooo There are duplicates left.Barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned."))
          c(uniq_tsb, add_tsb)
        }
      }
    }
  }
}