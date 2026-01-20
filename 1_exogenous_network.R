#' Identifying Peer Effects with Unobserved Effort and Isolated Students
#' Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos
#' 
#' This file replicates the peer effect model estimation assuming that the network is exogenous
rm(list = ls()) 
library(PartialNetwork)
library(AER)
library(dplyr)
library(fastDummies)
library(Rcpp)

PATH_DATA_OUT <- "~/Dropbox/Data/AHdata/PEEffort/" # Path to the folder where the prepared data will be saved. Note: the trailing "/" is required.
PATH_RESULTS  <- "~/Dropbox/Academy/1.Papers/EffortGPA/Code-EffortGPA/_output/" # Path to the output folder. Note: the trailing "/" is required.
PATH_CODE     <- "~/Dropbox/Academy/1.Papers/EffortGPA/Code-EffortGPA/codefiles/" # Path to the code folder. Note: the trailing "/" is required.

#########################################################################################################
#########################################################################################################
############################################ Main estimations ###########################################
#########################################################################################################
#########################################################################################################

# Load source functions
sourceCpp(paste0(PATH_CODE, "SourceCpp.cpp"))
source(paste0(PATH_CODE, "SourceR.R"))

# Load data
load(file = paste0(PATH_DATA_OUT, "Gpa.rda"))
rm("Xlogit")
gc()

va.exo        <- head(va.names, -1)

### Prepare data
G             <- norm.network(G)
XY            <- mydata[,va.names]
GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
hasfriends    <- unlist(lapply(G, rowSums))
Ghasfriends   <- peer.avg(G, hasfriends)

J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
J3            <- lapply(1:nsch, function(x) {
  onehat      <- rowSums(G[[x]])
  onecheck    <- 1 - onehat
  sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
  sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
  diag(sch.size[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
})
F1            <- fdataFs(J1)
F3            <- fdataFs(J3)

F1XY          <- peer.avg(F1, cbind(XY, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(XY), "hasfriends"))
F3XY          <- peer.avg(F3, XY); colnames(F3XY) <- paste0("F3_", colnames(XY))
F1GXY         <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
F3GXY         <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
F3GGX         <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
mydataFE3     <- data.frame(F3XY, F3GXY, F3GGX)

rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F3XY", "F1GXY", "F3GXY", "F1GGX", "F3GGX", "J1", "J3"))
gc()

## Formula
### With fixed effects and without dummy for isolated students (Formula for Model 1 in Table 3)
form1      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
### With fixed effects and with dummy for isolated students (Formula for Model 2 in Table 3)
form2      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr2     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
### With fixed effects and with dummy for isolated students per school (Formula for Model 3 in Table 3)
form3      <- as.formula(paste0("F3_gpa ~", paste(c(-1, "F3G_gpa", paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
instr3     <- as.formula(paste0("~", paste(c(-1, paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))

## Estimation
if (!dir.exists(paste0(PATH_RESULTS, "Table 3"))) {
  dir.create(paste0(PATH_RESULTS, "Table 3"), recursive = TRUE)
}

### With fixed effects and without dummy for isolated students (Model 1 in Table 3)
print("######")
print("With fixed effects and without dummy for isolated students (Model 1 in Table 3)")
iv.FE1       <- ivreg(formula = form1, instruments = instr1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
mean(siv.FE1$residuals^2)
saveRDS(siv.FE1, file = paste0(PATH_RESULTS, "Table 3/Model1.RDS"))
write.csv(siv.FE1$coefficients, file = paste0(PATH_RESULTS, "Table 3/Model1.csv"))

### With fixed effects and with dummy for isolated students (Model 2 in Table 3)
print("######")
print("With fixed effects and with dummy for isolated students (Model 2 in Table 3)")
iv.FE2       <- ivreg(formula = form2, instruments = instr2, data = mydataFE1)
(siv.FE2     <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients["F1G_gpa"],
                       G, fixed.effects = TRUE, F1, start = c(2.491167, 0.6053381))
va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2, "gpa")
va.iv.FE2$output
va.iv.FE2$sargan.stat
va.iv.FE2$sargan.pvalue
va.iv.FE2$sigma2eta
va.iv.FE2$sigma2epsilon
va.iv.FE2$rho
saveRDS(va.iv.FE2, file = paste0(PATH_RESULTS, "Table 3/Model2.RDS"))
va.iv.FE2    <- readRDS(file = paste0(PATH_RESULTS, "Table 3/Model2.RDS"))
write.csv(va.iv.FE2$output, file = paste0(PATH_RESULTS, "Table 3/Model2.csv"))

### With fixed effects and with dummy for isolated students per school (Model 3 in Table 3)
print("######")
print("With fixed effects and with dummy for isolated students per school (Model 3 in Table 3)")
iv.FE3       <- ivreg(formula = form3, instruments = instr3, data = mydataFE3)
(siv.FE3     <- summary(iv.FE3, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F3), function(x) rep(x, nrow(F3[[x]]))))))
ml.iv.FE3    <- foptim(iv.FE3$residuals, iv.FE3$coefficients["F3G_gpa"],
                       G, fixed.effects = TRUE, F3, start = c(2.491167, 0.6053381))
va.iv.FE3    <- fvariance.iv(iv.FE3, ml.iv.FE3, "gpa")
va.iv.FE3$output
va.iv.FE3$sargan.stat
va.iv.FE3$sargan.pvalue
va.iv.FE3$sigma2eta
va.iv.FE3$sigma2epsilon
va.iv.FE3$rho
saveRDS(va.iv.FE3, file = paste0(PATH_RESULTS, "Table 3/Model3.RDS"))
write.csv(va.iv.FE3$output, file = paste0(PATH_RESULTS, "Table 3/Model3.csv"))

## Bootstrap to compare the models
## Bootstrap to compare the models
fboot <- function(k) {
  if (k %in% seq(1, 1000, 100)) {
    cat("Boot: ", k, "\n", sep = "")
  }
  scID <- unique(mydata$sschlcde)
  ID   <- lapply(scID, \(x) which(mydata$sschlcde == x))
  sesc <- sample(scID, length(scID), replace = TRUE)
  seID <- unlist(ID[sapply(sesc, \(x) which(scID == x))])
  iv1  <- ivreg(formula = form1, instruments = instr1, data = mydataFE1[seID, ])
  iv2  <- ivreg(formula = form2, instruments = instr2, data = mydataFE1[seID, ])
  iv3  <- ivreg(formula = form3, instruments = instr3, data = mydataFE3[seID, ])
  names(iv2$coefficients) <- paste0("F2", substr(names(iv2$coefficients), 3, nchar(names(iv2$coefficients))))
  c(iv1$coefficients, iv2$coefficients, iv3$coefficients)
}

set.seed(15950)
boot   <- sapply(1:100, fboot)
saveRDS(boot, file = paste0(PATH_RESULTS, "Table 3/boot.RDS"))
(boot["F3G_gpa",] - boot["F2G_gpa",]) / sd(boot["F3G_gpa",] - boot["F2G_gpa",])

#########################################################################################################
#########################################################################################################
######################## Estimation when removing completely isolated students ##########################
#########################################################################################################
#########################################################################################################
rm(list = ls()[!(ls() %in% c("PATH_DATA_OUT", "PATH_RESULTS", "PATH_CODE"))])

# load source functions
sourceCpp(paste0(PATH_CODE, "SourceCpp.cpp"))
source(paste0(PATH_CODE, "SourceR.R"))

# Load data
load(file = paste0(PATH_DATA_OUT, "Gpa.rda"))
rm("Xlogit")
gc()

va.exo        <- head(va.names, -1)

### We remove completely isolated students
nhafr         <- lapply(G, rowSums) # number of friends
nisfr         <- lapply(G, colSums) # number of times the student is a friend
keep          <- lapply(1:nsch, function(x) (nhafr[[x]] > 0) | (nisfr[[x]] > 0)) #
mydata        <- mydata %>% filter(unlist(keep)) %>% group_by(sschlcde) %>% mutate(nstudent = n()) %>% ungroup() %>% filter(nstudent > 2)
G             <- norm.network(lapply(1:nsch, function(x) G[[x]][keep[[x]], keep[[x]]]))
Gnrow         <- sapply(G, nrow); G <- G[Gnrow > 2] # We remove school with completely isolated students
nsch          <- length(G)
sch.size      <- Gnrow[Gnrow > 2]

### Prepare data
XY            <- mydata[,va.names]
GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
hasfriends    <- unlist(lapply(G, rowSums))
Ghasfriends   <- peer.avg(G, hasfriends)

J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
J3            <- lapply(1:nsch, function(x) {
  onehat      <- rowSums(G[[x]])
  onecheck    <- 1 - onehat
  sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
  sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
  diag(sch.size[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
})
F1            <- fdataFs(J1)
F3            <- fdataFs(J3)

F1XY          <- peer.avg(F1, cbind(XY, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(XY), "hasfriends"))
F3XY          <- peer.avg(F3, XY); colnames(F3XY) <- paste0("F3_", colnames(XY))
F1GXY         <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
F3GXY         <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
F3GGX         <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
mydataFE3     <- data.frame(F3XY, F3GXY, F3GGX)

rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F3XY", "F1GXY", "F3GXY", "F1GGX", "F3GGX", "J1", "J3"))
gc()

## Formula
### With fixed effects and without dummy for isolated students (Formula for Model 1 in Table 4)
form1      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
### With fixed effects and with dummy for isolated students (Formula for Model 2 in Table 4)
form2      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr2     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
### With fixed effects and with dummy for isolated students per school (Formula for Model 3 in Table 4)
form3      <- as.formula(paste0("F3_gpa ~", paste(c(-1, "F3G_gpa", paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
instr3     <- as.formula(paste0("~", paste(c(-1, paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))

## Estimation
if (!dir.exists(paste0(PATH_RESULTS, "Table 4"))) {
  dir.create(paste0(PATH_RESULTS, "Table 4"), recursive = TRUE)
}

### With fixed effects and without dummy for isolated students (Model 1 in Table 4)
print("######")
print("With fixed effects and without dummy for isolated students (Model 1 in Table 4)")
iv.FE1       <- ivreg(formula = form1, instruments = instr1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
mean(siv.FE1$residuals^2)
saveRDS(siv.FE1, file = paste0(PATH_RESULTS, "Table 4/Model1.RDS"))
write.csv(siv.FE1$coefficients, file = paste0(PATH_RESULTS, "Table 4/Model1.csv"))

### With fixed effects and with dummy for isolated students (Model 2 in Table 4)
print("######")
print("With fixed effects and with dummy for isolated students (Model 2 in Table 4)")
iv.FE2       <- ivreg(formula = form2, instruments = instr2, data = mydataFE1)
(siv.FE2     <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients["F1G_gpa"],
                       G, fixed.effects = TRUE, F1, start = c(1.73344, 0.1682905))
va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2, "gpa")
va.iv.FE2$output
va.iv.FE2$sargan.stat
va.iv.FE2$sargan.pvalue
va.iv.FE2$sigma2eta
va.iv.FE2$sigma2epsilon
va.iv.FE2$rho
saveRDS(va.iv.FE2, file = paste0(PATH_RESULTS, "Table 4/Model2.RDS"))
va.iv.FE2    <- readRDS(file = paste0(PATH_RESULTS, "Table 4/Model2.RDS"))
write.csv(va.iv.FE2$output, file = paste0(PATH_RESULTS, "Table 4/Model2.csv"))

### With fixed effects and with dummy for isolated students per school (Model 3 in Table 4)
print("######")
print("With fixed effects and with dummy for isolated students per school (Model 3 in Table 4)")
iv.FE3       <- ivreg(formula = form3, instruments = instr3, data = mydataFE3)
(siv.FE3     <- summary(iv.FE3, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F3), function(x) rep(x, nrow(F3[[x]]))))))
ml.iv.FE3    <- foptim(iv.FE3$residuals, iv.FE3$coefficients["F3G_gpa"],
                       G, fixed.effects = TRUE, F3, start = c(1.73344, 0.1682905))
va.iv.FE3    <- fvariance.iv(iv.FE3, ml.iv.FE3, "gpa")
va.iv.FE3$output
va.iv.FE3$sargan.stat
va.iv.FE3$sargan.pvalue
va.iv.FE3$sigma2eta
va.iv.FE3$sigma2epsilon
va.iv.FE3$rho
saveRDS(va.iv.FE3, file = paste0(PATH_RESULTS, "Table 4/Model3.RDS"))
va.iv.FE3    <- readRDS(file = paste0(PATH_RESULTS, "Table 4/Model3.RDS"))
write.csv(va.iv.FE3$output, file = paste0(PATH_RESULTS, "Table 4/Model3.csv"))


#########################################################################################################
#########################################################################################################
############################# Estimation when removing any isolated students ############################
#########################################################################################################
#########################################################################################################
rm(list = ls()[!(ls() %in% c("PATH_DATA_OUT", "PATH_RESULTS", "PATH_CODE"))])

# load source functions
sourceCpp(paste0(PATH_CODE, "SourceCpp.cpp"))
source(paste0(PATH_CODE, "SourceR.R"))

# Load data
load(file = paste0(PATH_DATA_OUT, "Gpa.rda"))
rm("Xlogit")
gc()

va.exo        <- head(va.names, -1)

### We remove any isolated students
nhafr         <- lapply(G, rowSums) # number of friends
keep          <- lapply(1:nsch, function(x) nhafr[[x]] > 0)
mydata        <- mydata %>% filter(unlist(keep)) %>% group_by(sschlcde) %>% mutate(nstudent = n()) %>% ungroup() %>% filter(nstudent > 2)
G             <- lapply(1:nsch, function(x) G[[x]][keep[[x]], keep[[x]], drop = FALSE])
G             <- norm.network(G[sapply(G, is.matrix)])
Gnrow         <- sapply(G, nrow); G <- G[Gnrow > 2] # We remove school with completely isolated students
nsch          <- length(G)
sch.size      <- sapply(G, nrow)

### Prepare data
XY            <- mydata[,va.names]
GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))

J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
F1            <- fdataFs(J1)

F1XY          <- peer.avg(F1, XY); colnames(F1XY) <- paste0("F1_", c(colnames(XY)))
F1GXY         <- peer.avg(F1, GXY); colnames(F1GXY) <- paste0("F1", c(colnames(GXY)))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)

rm(list = c("GXY", "GGX", "F1XY", "F1GXY", "F1GGX", "J1"))
gc()

## Formula
### With fixed effects and without dummy for isolated students (Formula for Model 1' in Table 4)
form1      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))

## Estimation
### With fixed effects and without dummy for isolated students (Model 1' in Table 4)
print("######")
print("With fixed effects and without dummy for isolated students (Model 1' in Table 4)")
iv.FE1       <- ivreg(formula = form1, instruments = instr1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
mean(siv.FE1$residuals^2)
saveRDS(siv.FE1, file = paste0(PATH_RESULTS, "Table 4/Model1'.RDS"))
write.csv(siv.FE1$coefficients, file = paste0(PATH_RESULTS, "Table 4/Model1'.csv"))


#########################################################################################################
#########################################################################################################
################### Estimation when removing links from students with 1 or 2 friends ####################
#########################################################################################################
#########################################################################################################
if (!dir.exists(paste0(PATH_RESULTS, "Table 5"))) {
  dir.create(paste0(PATH_RESULTS, "Table 5"), recursive = TRUE)
}

lprem         <- c(0.5, 1) # share of individual with one or two friends whose links should be deleted
ltdeg         <- list(1, c(1, 2)) # Degree to remove. For a share `lprem` of students with a degree in ltdeg, links will be removed

print("######")
print("Table 5")

for (prem in lprem) {
  for (tdeg in ltdeg) {
    rm(list = ls()[!(ls() %in% c("PATH_DATA_OUT", "PATH_RESULTS", "PATH_CODE", "lprem", "prem", "ltdeg", "tdeg"))])
    gc()
    # load source functions
    sourceCpp(paste0(PATH_CODE, "SourceCpp.cpp"))
    source(paste0(PATH_CODE, "SourceR.R"))
    
    # Load data
    load(file = paste0(PATH_DATA_OUT, "Gpa.rda"))
    rm("Xlogit")
    gc()
    
    va.exo        <- head(va.names, -1)
    
    ### Prepare data
    ### Remove links
    set.seed(15950)
    degree        <- lapply(G, rowSums)
    cat("\nBefore removing links\n")
    print(proportions(table(unlist(degree))))
    rem           <- lapply(degree, function(s) {
      id1         <- which(s %in% tdeg)
      id1[runif(length(id1)) <= prem]
    })
    G             <- lapply(1:length(G), function(s) {
      out         <- G[[s]]
      out[rem[[s]],]<- 0
      out
    })
    cat("\nAfter removing links\n")
    print(proportions(table(unlist(lapply(G, rowSums)))))
    
    G             <- norm.network(G)
    XY            <- mydata[,va.names]
    GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
    GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
    hasfriends    <- unlist(lapply(G, rowSums))
    Ghasfriends   <- peer.avg(G, hasfriends)
    
    J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
    J3            <- lapply(1:nsch, function(x) {
      onehat      <- rowSums(G[[x]])
      onecheck    <- 1 - onehat
      sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
      sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
      diag(sch.size[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
    })
    F1            <- fdataFs(J1)
    F3            <- fdataFs(J3)
    
    F1XY          <- peer.avg(F1, cbind(XY, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(XY), "hasfriends"))
    F3XY          <- peer.avg(F3, XY); colnames(F3XY) <- paste0("F3_", colnames(XY))
    F1GXY         <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
    F3GXY         <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
    F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
    F3GGX         <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
    mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
    mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
    mydataFE3     <- data.frame(F3XY, F3GXY, F3GGX)
    
    rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F3XY", "F1GXY", "F3GXY", "F1GGX", "F3GGX", "J1", "J3"))
    gc()
    
    ## Formula
    ### With fixed effects and without dummy for isolated students (Formula for Model 1 in Table 3)
    form1      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
    instr1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
    ### With fixed effects and with dummy for isolated students (Formula for Model 2 in Table 3)
    form2      <- as.formula(paste0("F1_gpa ~", paste(c(-1, "F1G_gpa", "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
    instr2     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
    ### With fixed effects and with dummy for isolated students per school (Formula for Model 3 in Table 3)
    form3      <- as.formula(paste0("F3_gpa ~", paste(c(-1, "F3G_gpa", paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
    instr3     <- as.formula(paste0("~", paste(c(-1, paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))
    
    ## Estimation
    ### With fixed effects and without dummy for isolated students (Model 1 in Table 5)
    iv.FE1       <- ivreg(formula = form1, instruments = instr1, data = mydataFE1)
    siv.FE1      <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]])))))
    print(siv.FE1)
    print(mean(siv.FE1$residuals^2))
    saveRDS(siv.FE1, file = paste0(PATH_RESULTS, "Table 5/Model1.p=", prem * 100, "%.deg=",
                                   paste0(tdeg, collapse = "."), ".RDS"))
    write.csv(siv.FE1$coefficients, file = paste0(PATH_RESULTS, "Table 5/Model1.p=", prem * 100, "%.deg=",
                                                  paste0(tdeg, collapse = "."), ".csv"))
    
    ### With fixed effects and with dummy for isolated students (Model 2 in Table 5)
    iv.FE2       <- ivreg(formula = form2, instruments = instr2, data = mydataFE1)
    siv.FE2      <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]])))))
    print(siv.FE2)
    ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients["F1G_gpa"],
                           G, fixed.effects = TRUE, F1, start = c(2.491167, 0.6053381))
    va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2, "gpa")
    print(va.iv.FE2$output)
    print(va.iv.FE2$sargan.stat)
    print(va.iv.FE2$sargan.pvalue)
    print(va.iv.FE2$sigma2eta)
    print(va.iv.FE2$sigma2epsilon)
    print(va.iv.FE2$rho)
    saveRDS(va.iv.FE2, file = paste0(PATH_RESULTS, "Table 5/Model2.p=", prem * 100, "%.deg=",
                                     paste0(tdeg, collapse = "."), ".RDS"))
    va.iv.FE2    <- readRDS(file = paste0(PATH_RESULTS, "Table 5/Model2.p=", prem * 100, "%.deg=",
                                          paste0(tdeg, collapse = "."), ".RDS"))
    write.csv(va.iv.FE2$output, file = paste0(PATH_RESULTS, "Table 5/Model2.p=", prem * 100, "%.deg=",
                                              paste0(tdeg, collapse = "."), ".csv"))
    
    ### With fixed effects and with dummy for isolated students per school (Model 3 in Table 5)
    iv.FE3       <- ivreg(formula = form3, instruments = instr3, data = mydataFE3)
    siv.FE3     <- summary(iv.FE3, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F3), function(x) rep(x, nrow(F3[[x]])))))
    print(siv.FE3)
    ml.iv.FE3    <- foptim(iv.FE3$residuals, iv.FE3$coefficients["F3G_gpa"],
                           G, fixed.effects = TRUE, F3, start = c(2.491167, 0))
    va.iv.FE3    <- fvariance.iv(iv.FE3, ml.iv.FE3, "gpa")
    print(va.iv.FE3$output)
    print(va.iv.FE3$sargan.stat)
    print(va.iv.FE3$sargan.pvalue)
    print(va.iv.FE3$sigma2eta)
    print(va.iv.FE3$sigma2epsilon)
    print(va.iv.FE3$rho)
    saveRDS(va.iv.FE3, file = paste0(PATH_RESULTS, "Table 5/Model3.p=", prem * 100, "%.deg=",
                                     paste0(tdeg, collapse = "."), ".RDS"))
    va.iv.FE3    <- readRDS(file = paste0(PATH_RESULTS, "Table 5/Model3.p=", prem * 100, "%.deg=",
                                          paste0(tdeg, collapse = "."), ".RDS"))
    write.csv(va.iv.FE3$output, file = paste0(PATH_RESULTS, "Table 5/Model3.p=", prem * 100, "%.deg=",
                                              paste0(tdeg, collapse = "."), ".csv"))
  }
}

#########################################################################################################
#########################################################################################################
############################################## Tobit model ##############################################
#########################################################################################################
#########################################################################################################
rm(list = ls()[!(ls() %in% c("PATH_DATA_OUT", "PATH_RESULTS", "PATH_CODE"))])

# load source functions
sourceCpp(paste0(PATH_CODE, "SourceCpp.cpp"))
source(paste0(PATH_CODE, "SourceR.R"))

# Load data
load(file = paste0(PATH_DATA_OUT, "Gpa.rda"))
rm("Xlogit")
gc()

### Prepare data
G             <- norm.network(G)
hasfriends    <- unlist(lapply(G, rowSums))
FE            <- dummy_cols(data.frame(scid = mydata$sschlcde))[,-1]
fe            <- colnames(FE)
IsFE          <- FE * (1 - hasfriends)
nIsFE         <- FE * hasfriends
colnames(IsFE)  <- paste0("Is", fe)
colnames(nIsFE) <- paste0("nIs", fe)
feis          <- paste0("Is", fe)
fenis         <- paste0("nIs", fe)
va.exo        <- head(va.names, -1)
X             <- mydata[,va.exo]
GX            <- peer.avg(G, X); colnames(GX) <- paste0("G_", colnames(X))

mydata        <- cbind(mydata, GX, "hasfriends" = hasfriends, FE, nIsFE, IsFE)

## X variables
### With fixed effects and with dummy for isolated students (Model 1 in Table B3)
Xvar1 <- c(va.exo, paste0("G_", va.exo), fe)
### With fixed effects and with dummy for isolated students (Model 2 in Table B3)
Xvar2 <- c(va.exo, paste0("G_", va.exo), "hasfriends", fe)
### With fixed effects and with dummy for isolated students per school (Model 3 in Table B3)
Xvar3 <- c(va.exo, paste0("G_", va.exo), feis, fenis)

## Estimation
if (!dir.exists(paste0(PATH_RESULTS, "Table B3"))) {
  dir.create(paste0(PATH_RESULTS, "Table B3"), recursive = TRUE)
}

### With fixed effects and without dummy for isolated students (Model 1 in Table B3)
parm0 <- c(0.386022, rep(0, length(Xvar1)), 0.968185, 1.12965, -0.738497)
post1 <- fTobit(mydata$gpa, V = as.matrix(mydata[,Xvar1]), G = G, sim = 5e4,
                nthreads = parallel::detectCores() - 1, lby = 1, uby = 4,
                target = 0.4, jumpmin = 0.001, jumpmax = 1, parm0 = parm0,
                seed = 15950)
colnames(post1$parms) <- c("lambda", Xvar1, "seta", "sepsilon", "rho")
write.csv(data.frame(variables = colnames(post1$parms), t(apply(post1$parms, 2, \(x){
  c(coef = mean(tail(x, 25e3)), sd = sd(tail(x, 25e3)))
}))), file = paste0(PATH_RESULTS, "Table B3/Model1.csv"))
saveRDS(post1, file = paste0(PATH_RESULTS, "Table B3/Model1.RDS"))

### With fixed effects and with dummy for isolated students (Model 2 in Table B3)
parm0 <- c(0.386022, rep(0, length(Xvar2)), 0.968185, 1.12965, -0.738497)
post2 <- fTobit(mydata$gpa, V = as.matrix(mydata[,Xvar2]), G = G, sim = 5e4,
                nthreads = parallel::detectCores() - 1, lby = 1, uby = 4,
                target = 0.4, jumpmin = 0.001, jumpmax = 1, parm0 = parm0,
                seed = 15950)
colnames(post2$parms) <- c("lambda", Xvar2, "seta", "sepsilon", "rho")
write.csv(data.frame(variables = colnames(post2$parms), t(apply(post2$parms, 2, \(x){
  c(coef = mean(tail(x, 25e3)), sd = sd(tail(x, 25e3)))
}))), file = paste0(PATH_RESULTS, "Table B3/Model2.csv"))
saveRDS(post2, file = paste0(PATH_RESULTS, "Table B3/Model2.RDS"))

### With fixed effects and with dummy for isolated students per school (Model 3 in Table B3)
parm0 <- c(0.489059, rep(0, length(Xvar3)), 0.767328, 0.823438, -0.534227)
post3 <- fTobit(mydata$gpa, V = as.matrix(mydata[,Xvar3]), G = G, sim = 5e4,
                nthreads = parallel::detectCores() - 1, lby = 1, uby = 4,
                target = 0.4, jumpmin = 0.001, jumpmax = 1, parm0 = parm0,
                seed = 15950)
colnames(post3$parms) <- c("lambda", Xvar3, "seta", "sepsilon", "rho")
write.csv(data.frame(variables = colnames(post3$parms), t(apply(post3$parms, 2, \(x){
  c(coef = mean(tail(x, 25e3)), sd = sd(tail(x, 25e3)))
}))), file = paste0(PATH_RESULTS, "Table B3/Model3.csv"))
saveRDS(post3, file = paste0(PATH_RESULTS, "Table B3/Model3.RDS"))
