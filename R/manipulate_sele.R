#' Compile allele ages file
#'
#' This function aggregates mut, freq, and sele files inferred by Relate into an allele_ages file.
#'
#' @param mut data.table. mut file
#' @param freq data.frame. freq file
#' @param sele data.frame. sele file
#' @return Returns a data table.
#' @examples
#' get.allele_ages(mut, freq, sele)
#'
get.allele_ages <- function(mut, freq, sele){

  if( any(colnames(mut) == "CHR") & any(colnames(freq) == "CHR") & any(colnames(sele) == "CHR") ){

    if(any(colnames(mut) == "upstream")){
      mut <- mut[,c("CHR","pos_of_snp", "rs-id", "age_begin", "age_end", "ancestral_allele/alternative_allele", "upstream", "downstream")]
      colnames(mut) <- c("CHR","BP", "ID", "lower_age", "upper_age", "ancestral/derived", "upstream", "downstream")
    }else{
      mut <- mut[,c("CHR","pos_of_snp", "rs-id", "age_begin", "age_end", "ancestral_allele/alternative_allele")]
      colnames(mut) <- c("CHR","BP", "ID", "lower_age", "upper_age", "ancestral/derived")
    }

    sele <- sele[,c("CHR","pos", "when_mutation_has_freq2")]
    colnames(sele) <- c("CHR","BP", "pvalue")

    freq <- freq[,c("CHR","pos", "TreeFreq")]
    colnames(freq) <- c("CHR","BP", "DAF")

    allele_ages <- merge(mut, freq, by = c("CHR","BP"), all.x = T)
    allele_ages <- merge(allele_ages, sele, by = c("CHR","BP"), all.x = T)

  }else{

    warning("No column named CHR. Mergin data frame using BP only.")

    if(any(colnames(mut) == "upstream")){
      mut <- mut[,c("pos_of_snp", "rs-id", "age_begin", "age_end", "ancestral_allele/alternative_allele", "upstream", "downstream")]
      colnames(mut) <- c("BP", "ID", "lower_age", "upper_age", "ancestral/derived", "upstream", "downstream")
    }else{
      mut <- mut[,c("pos_of_snp", "rs-id", "age_begin", "age_end", "ancestral_allele/alternative_allele")]
      colnames(mut) <- c("BP", "ID", "lower_age", "upper_age", "ancestral/derived")
    }

    sele <- sele[,c("pos", "when_mutation_has_freq2")]
    colnames(sele) <- c("BP", "pvalue")

    freq <- freq[,c("pos", "TreeFreq")]
    colnames(freq) <- c("BP", "DAF")

    allele_ages <- merge(mut, freq, by = c("BP"), all.x = T)
    allele_ages <- merge(allele_ages, sele, by = c("BP"), all.x = T)

  }

  allele_ages$pvalue[allele_ages$pvalue == 1.0] <- NA
  return(allele_ages)

}

#' Filter allele_ages by quality of tree
#'
#' This function filters out selection pvalues on bad trees using summary statistics on tree quality stored in qual.
#'
#' @param allele_ages data.table. Obtained from get.allele_ages()
#' @param qual data.frame. qual file
#' @return Returns a data table.
#' @examples
#' Filter(allele_ages, qual)
#'
Filter <- function(allele_ages, qual){

  if(all(colnames(qual) != "BP")){
    stop("qual has no column named BP.")
  }

  before <- sum(!is.na(allele_ages$pvalue))

  if(any(colnames(allele_ages) == "CHR") & any(colnames(qual) == "CHR")){
    allele_ages <- merge(allele_ages, qual[,c(colnames(qual) != "ID")], by = c("CHR", "BP"))
  }else if(any(colnames(qual) == "CHR")){
    stop("qual has CHR column but allele_ages doesn't.")
  }else if(any(colnames(allele_ages) == "CHR")){
    stop("allele_ages has CHR column but qual doesn't.")
  }else{
    warning("No CHR column, merging allele_ages and qual only using BP.")
    allele_ages <- merge(allele_ages, qual[,c(colnames(qual) != "ID")], by = c("BP"))
  }

  threshold   <- c(quantile(qual$frac_branches_with_snp, probs = c(0.05)), quantile(qual$num_snps_on_tree, probs = c(0.05)))
  allele_ages <- subset(allele_ages, frac_branches_with_snp >= threshold[1] & num_snps_on_tree >= threshold[2])
  allele_ages <- allele_ages[,which(!(colnames(allele_ages) %in% c("frac_branches_with_snp", "num_snps_on_tree", "fraction_snps_not_mapping"))), with = FALSE]

  after <- sum(!is.na(allele_ages$pvalue))
  warning(paste0("Removed selection p-values from ", round(1.0-after/before,digits = 4)*100, "% of SNPs."))

  return(allele_ages)

}

