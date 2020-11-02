
object=GRset
snps = c("CpG", "SBE")
maf = 0
snpAnno = NULL

dropLociWithSnps
function (object, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)
{
    minfi:::.isGenomicOrStop(object)
    maf_cols <- paste0(snps, "_maf")
    snpDF <- getSnpInfo(object, snpAnno = snpAnno)
    choices <- c("Probe_maf", "CpG_maf", "SBE_maf")
    if (!all(choices %in% colnames(snpDF))) {
        stop("The specificed 'snpAnno' is not supported by this function")
    }
    if (sum(!(maf_cols %in% choices)) > 0) {
        stop("snps vector argument must be a combination of  \"Probe\", ",
            "\"CpG\" and \"SBE\"")
    }
    if (!is.numeric(maf) || maf < 0 || maf > 1) {
        stop("maf argument must be a numeric value between 0 and 1")
    }
    wh <- Reduce(union, lapply(maf_cols, function(xx) {
        which(snpDF[, xx] >= maf)
    }))
    wh <- sort(wh)
    if (length(wh) == 0)
        return(object)
    object[-wh, ]
}

object=GRset
snpAnno = NULL

getSnpInfo
function (object, snpAnno = NULL)
{
    av <- minfi:::.availableAnnotation(object)
    if (is.null(snpAnno)) {
        snpAnno <- grep(pattern = "^SNPs\\.", x = getAnnotationObject(object)@defaults,
            value = TRUE)
    } else {
        snpAnno <- sub("^SNPs\\.", "", snpAnno)
        if (!snpAnno %in% av$annoClassesChoices) {
            stop(sprintf("snpAnno '%s' is not part of the annotation",
                snpAnno))
        } else {
            snpAnno <- sprintf("SNPs.%s", snpAnno)
        }
    }
    snps <- getAnnotation(object, what = snpAnno)
    snps
}


