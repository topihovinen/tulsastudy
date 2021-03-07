library(DescTools)
# -------------------------------------------------------------------------------------------------------- #
# BPH Study (University of Turku & University of Helsinki) Statistical Analysis                            #
# Basic parameters of the data: 10 participants, 18 variables that were measured at 5 timepoints           #
# (Baseline, 3, 6, 9 and 12 months postop).                                                                #
#                                                                                                          #
# The "grand" null hypothesis is that the treatment has no effect to any parameter measured. Significant   #
# changes between baseline and 12 months postop are chosen to represent this change. To adjust for         #
# multiple statistical tests, we chose to use permutation-based adjustment method. We use all 256 possible #
# permutations (N=8 patients) between baseline and 12 month timepoint data to generate null distribution   #
# of p-values. Assuming no effect, using this method, significance level should be chosen so that there is #
# a 5% probability to observe at least one significant difference by chance in this dataset.               #
#                                                                                                          #
# Code originally written by: Topi Hovinen & Antti Viitala                                                 #
# -------------------------------------------------------------------------------------------------------- #

#create 256 permutations
pair1 <- DescTools::CombSet(c(1,5), 8, repl=TRUE, ord=TRUE) # 1 and 5 denoting timepoints 1 (baseline) and 5 (12 months postop) that are evaluated
pair2 <-matrix(6, nrow=256, ncol=8)
pair2 <- pair2 - pair1 #pair 1 and 2 are compliments of each other

# 0. Import data
{
  bphdata <- read.csv("bphdata.txt", header = T, dec = ',', sep = '\t', na.strings = "NA")
  bphdata_ncol <- length(bphdata[1,])
  bphdata[,3:bphdata_ncol]<-as.numeric(unlist(bphdata[,3:bphdata_ncol]))
}

# 1. Generate null distribution of significance levels
# We use wilcox.test with parameter paired = TRUE to evaluate statistical significance with Wilcoxon signed
# rank exact test. Rationale for choosing Wilcoxon signed rank exact test is discussed in more detail in
# the publication.
# We will use permutation-based method to correct for multiple testing to avoid "false discoveries" that
# occur at probability alpha = 5% with every independent statistical test.
# The logic and advantages of permutation tests are discussed in the Appendix of the publication

{
  perm_n <- 256 # Number of permutation rounds
  ids <- c(1,2,3,5,6,7,9,10) #Patients to consider. No data for patients 4 and 8 as their follow-up was discontinued prior to 12 months
  id_n <- length(ids) #Number of patients 
  
  # Initiate a variable for p values of each round for each variable. 18 different metrics were collected
  perm_pvalues <- data.frame(PSA = as.numeric(rep(NA,perm_n)),
                             ProstateVolume = as.numeric(rep(NA,perm_n)),
                             ResidualVolume = as.numeric(rep(NA,perm_n)),
                             AverageFlow = as.numeric(rep(NA,perm_n)),
                             MaxFlow = as.numeric(rep(NA,perm_n)),
                             VoidedVolume = as.numeric(rep(NA,perm_n)),
                             IPSSSympt = as.numeric(rep(NA,perm_n)),
                             IPSSQoL = as.numeric(rep(NA,perm_n)),
                             DANPSSV = as.numeric(rep(NA, perm_n)),
                             DANPSSS = as.numeric(rep(NA, perm_n)),
                             ICIQ = as.numeric(rep(NA, perm_n)),
                             IIEF = as.numeric(rep(NA,perm_n)),
                             IIEFQ2 = as.numeric(rep(NA,perm_n)),
                             EPIC26UI = as.numeric(rep(NA,perm_n)),
                             EPIC26UO = as.numeric(rep(NA,perm_n)),
                             EPIC26B = as.numeric(rep(NA,perm_n)),
                             EPIC26S = as.numeric(rep(NA,perm_n)),
                             EPIC26H = as.numeric(rep(NA,perm_n))
  )
  
  rpairs <- data.frame(ID = rep(bphdata$ID[1:8],2), #for each metric, pair of data at two timepoints
                       TimePoint = c(rep(1,8),rep(2,8)),
                       PSA = as.numeric(rep(NA, 16)),
                       ProstateVolume = as.numeric(rep(NA, 16)),
                       ResidualVolume = as.numeric(rep(NA, 16)),
                       AverageFlow = as.numeric(rep(NA, 16)),
                       MaxFlow = as.numeric(rep(NA, 16)),
                       VoidedVolume = as.numeric(rep(NA, 16)),
                       IPSSSympt = as.numeric(rep(NA, 16)),
                       IPSSQoL = as.numeric(rep(NA, 16)),
                       DANPSSV= as.numeric(rep(NA, 16)),
                       DANPSSS= as.numeric(rep(NA, 16)),
                       ICIQ = as.numeric(rep(NA, 16)),
                       IIEF = as.numeric(rep(NA, 16)),
                       IIEFQ2 = as.numeric(rep(NA, 16)),
                       EPIC26UI = as.numeric(rep(NA, 16)),
                       EPIC26UO = as.numeric(rep(NA, 16)),
                       EPIC26B= as.numeric(rep(NA, 16)),
                       EPIC26S = as.numeric(rep(NA, 16)),
                       EPIC26H = as.numeric(rep(NA, 16))
  )
  
  for(i in 1:perm_n){
    if(i == 1) {
      message("Starting permutations for generating null distribution of p-values in BPH Study\nTime is now: ", format(Sys.time(), "%H:%M.%S"))
    }

    for(j in 1:id_n){
      bphdata_j <- bphdata[bphdata$ID == ids[j],]
      rpairs[rpairs$ID == j,3:bphdata_ncol] <- rbind(bphdata_j[pair1[i,j],3:bphdata_ncol], bphdata_j[pair2[i,j],3:bphdata_ncol])
    }
    
    # For this permutation of data, calculate the p-values
    for(meas in 1:18) {
      col <- which(colnames(rpairs) == colnames(perm_pvalues[meas])) # Check the column number this way to avoid bugs
      if(meas %in% 1:6) { #Continuity correction off for continuous variables PSA, Prostate volume, Residual volume, Average flow, Max flow and Voided volume
        perm_pvalues[i,meas] = wilcox.test(rpairs[1:8,col], rpairs[9:16,col], paired = TRUE, exact = TRUE, correct = FALSE)$p.value
      } 
      else if(meas %in% 7:18) { # Continuity correction on for variables that stem from questionnaires (and are non-continuous by nature)
        perm_pvalues[i,meas] = wilcox.test(rpairs[1:8,col], rpairs[9:16,col], paired = TRUE, exact = TRUE, correct = TRUE)$p.value
      }
    }
  }
  perm_pvalues <- rbind(perm_pvalues)
  # Then, find a p-value limit alpha at which approximately 5% of these randomized datasets had at least one significantly different measurement, such that p < alpha.
  {
    # This can be done by checking through the smallest p-value of each round of 13 p-values, and selecting the 5th percentile
    min_pvalues <- apply(perm_pvalues, 1, min, na.rm = TRUE) # In fact, this also corresponds to cumulative distribution function of "corrected" p-values
    alpha <- quantile(min_pvalues, 0.05)
  }
  message("Permutations completed. The corrected significant p-value level for multiple comparisons at this dataset is\nalpha = ", alpha)
}

# 2. Next, calculate real p-values based on the actual timepoint-labels
{
  bph_pvalues <- data.frame(measure = as.character(),
                            comparison = as.character(),
                            test.statistic = as.numeric(),
                            p.value = as.numeric())
  for(meas in 1:18) { #18 variables
      
    # This will compare baseline to 12 months (Probably the most descriptive on )
    if(meas %in% 1:6) {
      testresult <- wilcox.test(bphdata[bphdata$TimePoint == "Baseline",(meas+2)], bphdata[bphdata$TimePoint == "12 Months",(meas+2)], exact = TRUE, paired = TRUE, correct=FALSE)
    } else if(meas %in% 7:18) {
      testresult <- wilcox.test(bphdata[bphdata$TimePoint == "Baseline",(meas+2)], bphdata[bphdata$TimePoint == "12 Months",(meas+2)], exact = TRUE, paired = TRUE, correct=TRUE)
    }
    bph_pvalues <- rbind(bph_pvalues,data.frame(measure = colnames(bphdata)[meas+2],
                                  comparison = "0-12",
                                  test.statistic = testresult$statistic[[1]],
                                  p.value = testresult$p.value))
  }
}

# Finally, calculate adjusted p-values based on the cumulative distribution function made from permutation tests
{
  pvalue_ecdf <- ecdf(min_pvalues)  #This creates an empirical cumulative distribution function of p-values
  bph_pvalues$adj.p.value <- NA
  for(meas in 1:18) {
    bph_pvalues$adj.p.value[meas] = pvalue_ecdf(bph_pvalues$p.value[meas])
  }    
}


