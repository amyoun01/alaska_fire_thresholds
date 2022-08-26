# This function calculates the AUC value. It is based off and inspired by the 
# "evaluate" function in the "dismo" package (Robert J. Hijmans et al. 2016)
#
# Robert J. Hijmans, Steven Phillips, John Leathwick and Jane Elith (2016). dismo: 
# Species Distribution Modeling. R package version 1.0-15. 
# http://CRAN.R-project.org/package=dismo

aucroc <- function(pres,abs,n){
  
  # ALLOCATE SPACE TO STORE RESULTS
  TPR <- matrix(NA,n,1) 
  FPR <- matrix(NA,n,1) 
  TNR <- matrix(NA,n,1) 
  FNR <- matrix(NA,n,1) 
  
  thresh <- seq(1,0,-(1/(n-1))) # SET PROBABILITY THRESHOLDS FOR CALCULATION OF 
                                # ROC VALUES
  
  for (i in 1:length(thresh)){ # FOR EACH THRESHOLD VALUE...
    
    p_i <- pres # SET PROBABILITY VALUES FOR OBSERVED PRESENCES TO THE VARIALE 'PRES'
    p_i[p_i >= thresh[i]] <- 1 # FOR ALL PROBABILITY VALUES GREATER THAN CURRENT THRESHOLD SET EQUAL TO 1
    p_i[p_i < thresh[i]] <- 0 # FOR ALL PROBABILITY VALUES LESS THAN CURRENT THRESHOLD SET EQUALT TO 0
    
    a_i <- abs
    a_i[a_i<=thresh[i]] <- 1.0001 # FOR PROBABILITY VALUES OF OBSERVED 
                                  # ABSENCES LESS THAN CURRENT THRESHOLD
                                  # SET EQUAL TO 1.0001. 
    a_i[a_i > thresh[i] & a_i < 1.0001] <- 0
    a_i[a_i == 1.0001] <- 1
    
    TPR[i,1] <- sum(p_i)/length(pres)     # CALCULATE TRUE POSITIVE RATE (SPECIFICITY)
    FPR[i,1] <- 1-(sum(a_i)/length(abs))  # CALCULATE FALSE POSITIVE RATE (OMMISSION)
    FNR[i,1] <- 1-(sum(p_i)/length(pres)) # CALCULATE FALSE NEGATIVE RATE(COMMISSION)
    TNR[i,1] <- sum(a_i)/length(abs)      # CALCULATE TRUE NEGATIVE RATE (SENSITIVITY)
  }
  
  # Calculate Area Under the Curve using the trapezoid rule
  idx <- 2:n # START ON SECOND VALUE, CAN'T CALCULATE AREA UNDER THE CURVE
             # FOR 1 VALUE.
  AUC <- as.double((FPR[idx] - FPR[idx - 1]) %*% (TPR[idx] + TPR[idx - 1]))/2

  # RETURN VALUES
  return(list(TPR=TPR,FPR=FPR,FNR=FNR,TNR=TNR,AUC=AUC,thresh=thresh))
  
}