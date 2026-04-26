############# 1.Load libraries and custom functions ###########

#install.packages("readxl")
#install.packages("dplyr")
#install.packages("Metrics")
#install.packages("data.table")
#install.packages("mice")
#install.packages("glmnet")
#install.packages("lubridate")
#install.packages("writexl")
#install.packages("FactoMineR")
#install.packages("forecast")
#install.packages("VIM")

library("readxl")
library("dplyr")
library("Metrics")
library("data.table")
library("mice")
library("glmnet")
library("lubridate")
library("writexl")
library("FactoMineR")
library("forecast")
library("VIM")

source("R/prep_spatial_ts.R") #Used in step 2
source("R/compute_spatial_window.R") #Used in step 3
source("R/vectors_s(q).R") #Used in step 4
source("R/vectors_p(r).R") #Used in step 5
source("R/tsfolds.R") #Used in step 5
source("R/lasso_ts.R") #One of two methods that can be used in step 5
source("R/enet_ts_strict.R") #One of two methods that can be used in step 5
source("R/pooledMASE.R") #Used in step 6
source("R/longest_na_run.R") #Used in step 7
source("R/repair.R") #Used in step 8
source("R/rv_similarity.R") #Used in step 9

############# 2.Parameter specification and data loading ##############

file_path <- "data/Additional file 2.xlsx" #File path of original input data
g <- 1 #Sheet number of working variable
d <- 17 #Column number of series displaying irregular pattern (in original data frame)
k_start <- 295   # ← row where irregular pattern STARTS
k_end   <- 588   # ← row where irregular pattern ENDS
f<- 4 #Column number of chosen reference series for comparing data similarity

#Adjust the following default model parameters as necessary

theta <- 24 # Periodicity of the data. For hourly data use theta = 24
k_t <- 5*theta # Desired degree of shift forwards and backwards. Algorithm will reduce k_t accordingly if there are less than k_t time periods on either side of the irregular pattern.
nfolds <- 10 #Chosen number of time folds for contiguous block time-series cross-validation
p <- 0 #Lagged/leading features included in imputation model from lag -p:p
m <- 5 #Number of imputations
step <- 1 #Change to 24 if you wish to use a seasonal step for hourly data
method="enet_ts_strict" #Either "enet_ts_strict" or "lasso_ts" 

#Load univariate spatial data for the working variable

my_spatial_ts_orig <- read_excel(file_path, sheet = g ) 

#Prepare working data frame 

my_spatial_ts=my_spatial_ts_orig %>%  relocate(d) 

#Extract time covariates from date variable and ensure they are placed in last columns as factor variables
#For example...
my_spatial_ts=my_spatial_ts %>% relocate(Date, .after = last_col()) 
my_spatial_ts=my_spatial_ts %>% relocate(Hour, .after = last_col()) 
my_spatial_ts$Hour=as.factor(my_spatial_ts$Hour)

#If irregular period spans more than 1 month, then include month covariate too
my_spatial_ts$Month=as.factor(month(my_spatial_ts$Date))

my_spatial_ts <- prep_spatial_ts(my_spatial_ts)

########## 3.Determine parameters for spatial window #######

window_params <- compute_spatial_window(
  my_spatial_ts = my_spatial_ts,
  k_t           = k_t,
  k_start       = k_start,
  k_end         = k_end
)

l       <- window_params$l
l_start <- window_params$l_start
l_end   <- window_params$l_end
k       <- window_params$k
k_t     <- window_params$k_t   # possibly reduced

################ 4.Generate interpolated s(q) vectors ##############

s0_shifts <- generate_s0_shifts(
  my_spatial_ts = my_spatial_ts,
  l_start       = l_start,
  l_end         = l_end,
  l             = l,
  k             = k,
  k_t           = k_t,
  theta         = theta
)

Series_irregular <- s0_shifts$Series_irregular
shift            <- s0_shifts$shift


############## 5.Generate p(r) vectors for r=1:m ################ 

#Predictor construction and multiple imputation

result <- generate_pr_imputed(
  my_spatial_ts = my_spatial_ts,   
  l_start       = l_start,
  l_end         = l_end,
  l             = l,
  k             = k,
  p             = p,
  m             = m,
  method        = method          
)

list2env(result, envir = .GlobalEnv)

############## 6.Generate pooled MASE for each value of q ############

mase_results <- compute_pooled_mase(
  Pr_imputed = Pr_imputed,
  shift      = shift,
  step       = step,
  m          = m,      
  k_t        = k_t     
)

MASE       <- mase_results$MASE
PooledMASE <- mase_results$PooledMASE


########### 7.Identify the optimal shift #################

#Extract, save and display results

shifts <- -k_t:k_t
results <- data.frame(Shift = shifts, PooledMASE = PooledMASE)
write.csv(results, "output/results.csv",row.names = FALSE)

opt_shift_idx <- which.min(PooledMASE)

sink("output/repair_summary.txt")
cat("=== REPAIR SUMMARY ===\n\n")
cat("Optimal shift:                               ", shifts[opt_shift_idx], "\n")
cat("Min pooled MASE:                             ", round(min(PooledMASE), 4), "\n")
cat("Complete missing periods (count):            ", total_missing_periods, "\n")
cat("Complete missing periods (%):                ", pct_missing_periods, "\n")
cat("Longest NA gap in irregular series:          ", longest_na_run(my_spatial_ts_orig[[d]][k_start:k_end]), "\n")
cat("Missing % in irregular series:               ", round(sum(is.na(my_spatial_ts_orig[(k_start:k_end), d])) / (k_end - k_start + 1) * 100, 2), "%\n")
cat("Missing % in window:                         ", round(percent_missing, 2), "%\n")
cat("Avg between-series correlation in window:    ", round(avg_cor, 4), "\n")
sink()

#Generate plot of pooled MASE for different degrees of shift 

ylim=c(0,9) #Adjust as necessary
x=c(-k_t:k_t)

png(filename = "output/pooledMASE_plot.png", width = 800, height = 600)
plot(x, PooledMASE, type = "l", xlab = "", ylab="", main="", xaxt="n", yaxs="i", ylim=ylim,cex.axis=1.5)
mytitle=expression("Evaluation of pooled MASE by degree of shift in s"^(q))
mtext(side=3, line=1.5, cex=1.5, mytitle)
title(ylab="Pooled MASE", line=2.5,cex.lab =1.5 )
title(xlab="Periods shifted (q)", line=3,cex.lab=1.5)
axis(1, at=c(-k_t,shifts[opt_shift_idx],k_t),cex.axis=1.5)
abline(h = 1, lty = 3)
abline(v=0, lty=3)
dev.off()

#Note: If you are unhappy with the results (e.g. if optimal shift is associated with a pooled MASE significantly greater than 1), then re-run the algorithm using a higher value of k_t.  Otherwise proceed with repair of full multivariate set.

##### 8.Repair irregular pattern in multivariate data ########

sheetnames <- readxl::excel_sheets(file_path)
j=length(sheetnames)

for (i in 1:j) {
  
  my_spatial_ts_orig <- read_excel(file_path, 
    sheet = sheetnames[i]          
  )
  
  repaired <- repair_spatial_ts(
    my_spatial_ts_orig = my_spatial_ts_orig,
    k_start = k_start,
    k_end = k_end,
    d = d,
    shifts = shifts,
    opt_shift_idx = opt_shift_idx
  )
  
  assign(sheetnames[i], repaired)
}

repaired_list <- mget(sheetnames, envir = .GlobalEnv)

writexl::write_xlsx(
  x = repaired_list,
  path = "output/my_spatial_ts_repaired.xlsx"
)

## 9.Assess improvement in data similarity between spatial series ####

# Empty list to store columns
col_list <- list()

#Extract original multivariate data for affected spatial location
for (i in 1:j) {
  df <- read_excel(file_path, sheet = i)
  col_list[[i]] <- df[[d]]
}

target_orig <- as.matrix(as.data.frame(col_list))

#Extract original multivariate data for chosen reference location
for (i in 1:j) {
  df <- read_excel(file_path, sheet = i)
  col_list[[i]] <- df[[f]]
}

ref_orig <- as.matrix(as.data.frame((col_list)))

#Extract repaired multivariate data for affected spatial location
target_repaired <- do.call(cbind, lapply(repaired_list, `[[`, d))

#Compute RV similarity versus reference location before and after repair
similarity=compute_rv_similarity(target_orig, target_repaired, ref_orig, j, k_start, k_end) 

write.csv(similarity, "output/similarity.csv", row.names = TRUE)


