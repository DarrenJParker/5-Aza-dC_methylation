#### data -  

### using command line args
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### aborts if no command line args provided
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

input_file_name <- trimws(args[1])

dat1 = read.csv(input_file_name, check.names=FALSE, stringsAsFactors=FALSE) 
# head(dat1)
### response:
# - meth_c   (successfully methylated sites)
# - unmeth_c (unsuccessfully methylated sites)
# factors
# - treat - three levels: Sucrose_24h 5aza_dC_24h 5aza_dC_48h 
# - collection - 3 levels: 0h, 24h and 48h
# - tiss_type - 2 levels: Bodies, Heads 


### set as factors + add ref level
dat1$treat <- as.factor(dat1$treat )
dat1$collection <- as.factor(dat1$collection)
dat1$tiss_type <- as.factor(dat1$tiss_type)

dat1$treat  <- relevel(dat1$treat, ref = "Sucrose_24h")
dat1$collection <- relevel(dat1$collection, ref = "0h")
dat1$tiss_type <- relevel(dat1$tiss_type, ref = "Bodies")

#### fit full_model (binomial glm)

m4 <- glm(cbind(dat1$meth_c, dat1$unmeth_c) ~ dat1$treat *  dat1$collection * dat1$tiss_type, family = "binomial")
#summary(m4)

#### use a log-L chi-sq to test significance 

# drop1(m4,~.,test="Chisq")

## same results as using 
# library(car)
# Anova(m4, type = 3 ) 

##### export 

### coeffs

write.csv(summary(m4)$coefficients, file= paste(input_file_name, "_coeff.csv", sep = ""), row.names=TRUE)

## pvals (overall effs)

write.csv(drop1(m4,~.,test="Chisq"), file= paste(input_file_name, "_LLChisq.csv", sep = ""), row.names=TRUE)
