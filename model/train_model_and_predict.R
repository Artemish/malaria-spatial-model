library(dplyr)
library(reshape2)
library(ggplot2)
library(MASS)
library(tidyr)
library(CARBayesST)
library(splines)

load("../data/full_m_district_PRCP.rda")
load("../data/malaria_related.rda")
load("../data/map_mat.rda")
load("../data/storm_district_final.rda")
load("../data/id.rda")
load("../data/temperature_district_final.rda")

# > load('../data/id.rda')
# > length(South_id) / 60
# [1] 30
# > length(North_id) / 60
# [1] 65
# > length(Central_id) / 60
# [1] 33

STUDY_MONTHS = 60
month_cols <- c(1,26:26+STUDY_MONTHS-1)

# PRCP
full_m_district_PRCP = full_m_district_PRCP[,1:4]
full_m_district_PRCP_wide_mean = dcast(full_m_district_PRCP, area~YEARMONTH, value.var = 'PRCP_timetotal_spatialmean_district')
full_m_district_PRCP_wide_total = dcast(full_m_district_PRCP, area~YEARMONTH, value.var = 'PRCP_timemean_spatialtotal_district')

full_m_district_PRCP_wide_mean_new = full_m_district_PRCP_wide_mean[match( marlaria_dst_case_final_new$district, full_m_district_PRCP_wide_mean$area),]
full_m_district_PRCP_wide_total_new = full_m_district_PRCP_wide_total[match( marlaria_dst_case_final_new$district, full_m_district_PRCP_wide_total$area),]

full_m_district_PRCP_wide_mean_new = full_m_district_PRCP_wide_mean_new[,month_cols]
full_m_district_PRCP_wide_total_new = full_m_district_PRCP_wide_total_new[,month_cols]

# TEMP
temperature_district_final_wide_max = dcast(temperature_district_final, area~YEARMONTH, value.var = 'lst_day_c_max')
temperature_district_final_wide_min = dcast(temperature_district_final, area~YEARMONTH, value.var = 'lst_day_c_min')
temperature_district_final_wide_mean = dcast(temperature_district_final, area~YEARMONTH, value.var = 'lst_day_c_mean')

temperature_district_final_wide_max_new = temperature_district_final_wide_max[match( marlaria_dst_case_final_new$district, temperature_district_final_wide_max$area),]
temperature_district_final_wide_min_new = temperature_district_final_wide_min[match( marlaria_dst_case_final_new$district, temperature_district_final_wide_min$area),]
temperature_district_final_wide_mean_new = temperature_district_final_wide_mean[match( marlaria_dst_case_final_new$district, temperature_district_final_wide_min$area),]

temperature_district_final_wide_max_new = temperature_district_final_wide_max_new[,month_cols]
temperature_district_final_wide_min_new = temperature_district_final_wide_min_new[,month_cols]
temperature_district_final_wide_mean_new = temperature_district_final_wide_mean_new[,month_cols]

marlaria_dst_case_final_new_tmp = cbind(marlaria_dst_case_final_new,
                                        spaceid=1:130)

tmp = cbind(marlaria_dst_population_final_new[,1],
            data.frame(rep(marlaria_dst_population_final_new[,-1], each=12))[,1:STUDY_MONTHS])

colnames(tmp) = c('district',1:STUDY_MONTHS)

# tom = melt(marlaria_dst_case_final_new_tmp,id.vars = 'spaceid')
# 
# tom = reshape(marlaria_dst_case_final_new_tmp, 
#               idvar = "spaceid",  direction = "wide")

marlaria_case_long = melt(marlaria_dst_case_final_new, id.vars = 'district', variable.name = "timeid")
PRCP_mean_long = melt(full_m_district_PRCP_wide_total_new, id.vars = 'area', variable.name = "timeid")
population_long = melt(tmp, id.vars = 'district', variable.name = "timeid")
TEMP_max_long = melt(temperature_district_final_wide_max_new, id.vars = 'area', variable.name = "timeid")
TEMP_min_long = melt(temperature_district_final_wide_min_new, id.vars = 'area', variable.name = "timeid")
TEMP_mean_long = melt(temperature_district_final_wide_mean_new, id.vars = 'area', variable.name = "timeid")

colnames(marlaria_case_long) = c('distrcit', 'timeid', 'marlaria')
colnames(PRCP_mean_long) = c('distrcit', 'timeid', 'PRCP')
colnames(population_long) = c('distrcit', 'timeid', 'population')
colnames(TEMP_max_long) = c('distrcit', 'timeid', 'TEMPmax')
colnames(TEMP_min_long) = c('distrcit', 'timeid', 'TEMPmin')
colnames(TEMP_mean_long) = c('distrcit', 'timeid', 'TEMPmean')

storm_full = PRCP_mean_long[,2:1]
storm_full$exists <- as.numeric(do.call(paste0, storm_full) %in% do.call(paste0, storm_district_final))

# final_data = cbind(spaceid = rep(1:130, STUDY_MONTHS),
#                    timeid = marlaria_case_long$timeid,
#                    marlaria = marlaria_case_long$marlaria,
#                    PRCP = PRCP_mean_long$PRCP,
#                    population = round(population_long$population))

South_district = c('Boane', 'Magude', 'Manhica', 'Marracuene', 'Matutuine', 'Moamba', 'Namaacha', 'MaputoCity', 'Funhalouro', 'Guvuro', 'Homoine', 'Inharime', 'Inhassoro', 'Jangamo', 'Mabote', 'Massinga', 'Morrumbene', 'Panda', 'Vilankulo', 'Zavala', 'Bilene', 'Chokwe', 'Chibuto', 'Chicualacuala', 'Chigubo', 'Guija', 'Mabalane', 'Mandlacaze', 'Massangena', 'Massingir', 'CIDADE DE XAI-XAI', 'Maputo')
Central_district = c('Barue', 'Gondola', 'Guro', 'Machaze', 'Macossa', 'Manica', 'Mossurize', 'Sussundenga', 'Tambara', "Buzi", "Caia", "Chemba", "Cheringoma", "Chibabava", "Dondo", "Gorongoza", "Gorongosa", "Machanga", "Maringue", "Marromeu", "Muanza", "Nhamatanda", "Angonia", "Cahora Bassa", "Changara", "Chifunde", "Chiuta", "Macanga", "Magoe", "Maravia", "Moatize", "Mutarara", "Tsangano", "Zumbo")
North_prov = c('Zambezia', 'Nampula', 'Nassa', 'Cabo Delgado')

# 
# cccc = map %>%
#   filter(NAME_1 %in% North_prov)

# North_district = cccc$NAME_2
# 
# PRCP_mean_long_new = PRCP_mean_long %>%
#   mutate(pos_district = case_when(distrcit %in% South_district ~ 'South',
#                                   distrcit %in% Central_district ~ 'Central',
#                                   distrcit %in% North_district ~ 'North'))
# 
# South_id = which(PRCP_mean_long_new$pos_district == 'South')
# Central_id = which(PRCP_mean_long_new$pos_district == 'Central')
# North_id = which(PRCP_mean_long_new$pos_district == 'North')
# 
# save(South_id, Central_id, North_id, file='id.RData')

N_SOUTH = length(South_id) / STUDY_MONTHS
N_CENTRAL = length(Central_id) / STUDY_MONTHS
N_NORTH = length(North_id) / STUDY_MONTHS


# center PRCP version
# final_data = cbind(spaceid = rep(1:130, 60),
#                    timeid = marlaria_case_long$timeid,
#                    marlaria = marlaria_case_long$marlaria,
#                    PRCP = c(scale(PRCP_mean_long$PRCP, scale=FALSE)),
#                    population = round(population_long$population),
#                    storm = storm_full$exists)
# final_data = cbind(spaceid = rep(1:130, 60),
#                   timeid = marlaria_case_long$timeid,
#                   marlaria = marlaria_case_long$marlaria,
#                   PRCP = c(scale(PRCP_mean_long$PRCP, scale=FALSE)),
#                   population = round(population_long$population),
#                   TEMPmax = c(scale(TEMP_max_long$TEMPmax, scale=FALSE)),
#                   TEMPmin = c(scale(TEMP_min_long$TEMPmin, scale=FALSE)),
#                   TEMPmean = c(scale(TEMP_mean_long$TEMPmean, scale=FALSE)),
#                   storm = storm_full$exists)

# South
# final_data = cbind(spaceid = rep(1:N_SOUTH, 60),
#                     timeid = marlaria_case_long$timeid[South_id],
#                     marlaria = marlaria_case_long$marlaria[South_id],
#                     PRCP = c(scale(PRCP_mean_long$PRCP[South_id], scale=FALSE)),
#                     population = round(population_long$population[South_id]),
#                     TEMPmax = c(scale(TEMP_max_long$TEMPmax[South_id], scale=FALSE)),
#                     storm = storm_full$exists[South_id])
# map_mat_south = map_mat[South_id[1:N_SOUTH],South_id[1:N_SOUTH]]
 
# # Central
# final_data = cbind(spaceid = rep(1:N_CENTRAL, 60),
#                     timeid = marlaria_case_long$timeid[Central_id],
#                     marlaria = marlaria_case_long$marlaria[Central_id],
#                     PRCP = c(scale(PRCP_mean_long$PRCP[Central_id], scale=FALSE)),
#                     population = round(population_long$population[Central_id]),
#                     TEMPmax = c(scale(TEMP_max_long$TEMPmax[Central_id], scale=FALSE)),
#                     storm = storm_full$exists[Central_id])
# map_mat_central = map_mat[Central_id[1:N_CENTRAL], Central_id[1:N_CENTRAL]]


# North
final_data = cbind(spaceid = rep(1:N_NORTH, 60),
                  timeid = marlaria_case_long$timeid[North_id],
                  marlaria = marlaria_case_long$marlaria[North_id],
                  PRCP = c(scale(PRCP_mean_long$PRCP[North_id], scale=FALSE)),
                  population = round(population_long$population[North_id]),
                  TEMPmax = c(scale(TEMP_max_long$TEMPmax[North_id], scale=FALSE)),
                  storm = storm_full$exists[North_id])

map_mat_north = map_mat[North_id[1:N_NORTH],North_id[1:N_NORTH]]


final_data = data.frame(final_data)
#final_data$add_PRCP = as.integer(final_data$PRCP>100)


# lagged PRCP
lagged1_PRCP = final_data[final_data$timeid<STUDY_MONTHS,]$PRCP
lagged2_PRCP = final_data[final_data$timeid<STUDY_MONTHS-1,]$PRCP
lagged3_PRCP = final_data[final_data$timeid<STUDY_MONTHS-2,]$PRCP
lagged4_PRCP = final_data[final_data$timeid<STUDY_MONTHS-3,]$PRCP
curr_col = ncol(final_data)

final_data$lag1_PRCP = final_data$PRCP
final_data$lag2_PRCP = final_data$PRCP
final_data$lag3_PRCP = final_data$PRCP
final_data$lag4_PRCP = final_data$PRCP
final_data[final_data$timeid==1,curr_col+1] = NA
final_data[final_data$timeid<=2,curr_col+2] = NA
final_data[final_data$timeid<=3,curr_col+3] = NA
final_data[final_data$timeid<=4,curr_col+4] = NA
final_data[final_data$timeid>1,curr_col+1] = lagged1_PRCP
final_data[final_data$timeid>2,curr_col+2] = lagged2_PRCP
final_data[final_data$timeid>3,curr_col+3] = lagged3_PRCP
final_data[final_data$timeid>4,curr_col+4] = lagged4_PRCP
#final_data_lagged4$storm = as.integer(final_data_lagged4$storm)

# lagged TEMP_k_max
lagged1_TEMPmax = final_data[final_data$timeid<STUDY_MONTHS,]$TEMPmax
lagged2_TEMPmax = final_data[final_data$timeid<STUDY_MONTHS-1,]$TEMPmax
lagged3_TEMPmax = final_data[final_data$timeid<STUDY_MONTHS-2,]$TEMPmax
lagged4_TEMPmax = final_data[final_data$timeid<STUDY_MONTHS-3,]$TEMPmax
curr_col = ncol(final_data)

final_data$lag1_TEMPmax = final_data$TEMPmax
final_data$lag2_TEMPmax = final_data$TEMPmax
final_data$lag3_TEMPmax = final_data$TEMPmax
final_data$lag4_TEMPmax = final_data$TEMPmax
final_data[final_data$timeid==1,curr_col+1] = NA
final_data[final_data$timeid<=2,curr_col+2] = NA
final_data[final_data$timeid<=3,curr_col+3] = NA
final_data[final_data$timeid<=4,curr_col+4] = NA
final_data[final_data$timeid>1,curr_col+1] = lagged1_TEMPmax
final_data[final_data$timeid>2,curr_col+2] = lagged2_TEMPmax
final_data[final_data$timeid>3,curr_col+3] = lagged3_TEMPmax
final_data[final_data$timeid>4,curr_col+4] = lagged4_TEMPmax


# # lagged TEMP_k_min
# lagged1_TEMPmin = final_data[final_data$timeid<60,]$TEMPmin
# lagged2_TEMPmin = final_data[final_data$timeid<59,]$TEMPmin
# lagged3_TEMPmin = final_data[final_data$timeid<58,]$TEMPmin
# lagged4_TEMPmin = final_data[final_data$timeid<57,]$TEMPmin
# curr_col = ncol(final_data)

# final_data$lag1_TEMPmin = final_data$TEMPmin
# final_data$lag2_TEMPmin = final_data$TEMPmin
# final_data$lag3_TEMPmin = final_data$TEMPmin
# final_data$lag4_TEMPmin = final_data$TEMPmin
# final_data[final_data$timeid==1,curr_col+1] = NA
# final_data[final_data$timeid<=2,curr_col+2] = NA
# final_data[final_data$timeid<=3,curr_col+3] = NA
# final_data[final_data$timeid<=4,curr_col+4] = NA
# final_data[final_data$timeid>1,curr_col+1] = lagged1_TEMPmin
# final_data[final_data$timeid>2,curr_col+2] = lagged2_TEMPmin
# final_data[final_data$timeid>3,curr_col+3] = lagged3_TEMPmin
# final_data[final_data$timeid>4,curr_col+4] = lagged4_TEMPmin

# # lagged TEMP_k_mean
# lagged1_TEMPmean = final_data[final_data$timeid<60,]$TEMPmean
# lagged2_TEMPmean = final_data[final_data$timeid<59,]$TEMPmean
# lagged3_TEMPmean = final_data[final_data$timeid<58,]$TEMPmean
# lagged4_TEMPmean = final_data[final_data$timeid<57,]$TEMPmean
# curr_col = ncol(final_data)

# final_data$lag1_TEMPmean = final_data$TEMPmean
# final_data$lag2_TEMPmean = final_data$TEMPmean
# final_data$lag3_TEMPmean = final_data$TEMPmean
# final_data$lag4_TEMPmean = final_data$TEMPmean
# final_data[final_data$timeid==1,curr_col+1] = NA
# final_data[final_data$timeid<=2,curr_col+2] = NA
# final_data[final_data$timeid<=3,curr_col+3] = NA
# final_data[final_data$timeid<=4,curr_col+4] = NA
# final_data[final_data$timeid>1,curr_col+1] = lagged1_TEMPmean
# final_data[final_data$timeid>2,curr_col+2] = lagged2_TEMPmean
# final_data[final_data$timeid>3,curr_col+3] = lagged3_TEMPmean
# final_data[final_data$timeid>4,curr_col+4] = lagged4_TEMPmean



final_data_lagged4 = final_data[final_data$timeid>4,]
final_data_lagged4$log_malaria = log(final_data_lagged4$marlaria+0.001)
final_data_lagged4$log_rate = log(final_data_lagged4$marlaria/final_data_lagged4$population+0.001)

#Ncar <- 500000
Ncar <- 20000
#burn.in.car <- 100000
burn.in.car <- 4000
thinning <- 40

lag1_PRCP_q1 = summary(final_data_lagged4$lag1_PRCP)[2]
lag1_PRCP_q2 = summary(final_data_lagged4$lag1_PRCP)[3]
lag1_PRCP_q3 = summary(final_data_lagged4$lag1_PRCP)[4]
lag2_PRCP_q1 = summary(final_data_lagged4$lag2_PRCP)[2]
lag2_PRCP_q2 = summary(final_data_lagged4$lag2_PRCP)[3]
lag2_PRCP_q3 = summary(final_data_lagged4$lag2_PRCP)[4]

TEMPmax_q1 = summary(final_data_lagged4$TEMPmax)[2]
TEMPmax_q2 = summary(final_data_lagged4$TEMPmax)[3]
TEMPmax_q3 = summary(final_data_lagged4$TEMPmax)[4]
lag1_TEMPmax_q1 = summary(final_data_lagged4$lag1_TEMPmax)[2]
lag1_TEMPmax_q2 = summary(final_data_lagged4$lag1_TEMPmax)[3]
lag1_TEMPmax_q3 = summary(final_data_lagged4$lag1_TEMPmax)[4]
lag2_TEMPmax_q1 = summary(final_data_lagged4$lag2_TEMPmax)[2]
lag2_TEMPmax_q2 = summary(final_data_lagged4$lag2_TEMPmax)[3]
lag2_TEMPmax_q3 = summary(final_data_lagged4$lag2_TEMPmax)[4]
lag3_TEMPmax_q1 = summary(final_data_lagged4$lag3_TEMPmax)[2]
lag3_TEMPmax_q2 = summary(final_data_lagged4$lag3_TEMPmax)[3]
lag3_TEMPmax_q3 = summary(final_data_lagged4$lag3_TEMPmax)[4]
# f = marlaria ~ offset(log(population)) + lag1_PRCP + lag2_PRCP + 
#   storm:lag1_PRCP + storm:lag2_PRCP +
#   TEMPmax + lag1_TEMPmax + lag2_TEMPmax + lag3_TEMPmax 


f = marlaria ~ offset(log(population)) + ns(lag1_PRCP, knots = c(lag1_PRCP_q1, lag1_PRCP_q2, lag1_PRCP_q3)) + ns(lag2_PRCP, knots = c(lag2_PRCP_q1, lag2_PRCP_q2, lag2_PRCP_q3)) + 
  ns(TEMPmax, knots = c(TEMPmax_q1, TEMPmax_q2, TEMPmax_q3)) +
  ns(lag1_TEMPmax, knots = c(lag1_TEMPmax_q1, lag1_TEMPmax_q2, lag1_TEMPmax_q3)) + ns(lag2_TEMPmax, knots = c(lag2_TEMPmax_q1, lag2_TEMPmax_q2, lag2_TEMPmax_q3)) + ns(lag3_TEMPmax, knots = c(lag3_TEMPmax_q1, lag3_TEMPmax_q2, lag3_TEMPmax_q3))
  
#f = marlaria ~ offset(log(population)) + lag1_PRCP + lag2_PRCP 
# f = marlaria ~ offset(log(population)) + TEMPmax + lag1_TEMPmax + lag2_TEMPmax + lag3_TEMPmax + lag4_TEMPmax +
#   TEMPmin + lag1_TEMPmin + lag2_TEMPmin + lag3_TEMPmin + lag4_TEMPmin +
#   TEMPmean + lag1_TEMPmean + lag2_TEMPmean + lag3_TEMPmean + lag4_TEMPmean 

# f = marlaria ~ offset(log(population)) + TEMPmax + lag1_TEMPmax + lag2_TEMPmax + lag3_TEMPmax + lag4_TEMPmax 
  
#f = marlaria ~ lag2_PRCP
#f = log_malaria ~ lag2_PRCP
# f = log_rate ~ lag2_PRCP
#f = marlaria ~ offset(log(population)) + storm
# M20_lagged2 <- Bcartime(formula=f, data=final_data_lagged2, 
#                         scol="spaceid", tcol= "timeid",  
#                W=map_mat, model="ar", AR=2, family="poisson", package="CARBayesST",
#                N=Ncar, burn.in=burn.in.car)
# fitts2 = M20_lagged2$fitteds



# CarBayesSt_lagged2 = ST.CARar(formula=f, data=final_data_lagged2, 
#                               family='poisson', W=map_mat, burnin=burn.in.car,
#                               n.sample=Ncar, thin=10, n.chains=2, n.cores=2, AR=2)


# CarBayesSt_lagged3 = ST.CARar(formula=f, data=final_data_lagged4, 
#                               family='poisson', W=map_mat, burnin=burn.in.car,
#                               n.sample=Ncar, thin=thinning, n.chains=3, n.cores=3, AR=2,
#                               rho.S=0.52, rho.T=c(0.68, 0.26))

CarBayesSt_lagged4 = ST.CARar(formula=f, data=final_data_lagged4, 
                              family='poisson', W=map_mat_north, burnin=burn.in.car,
                              n.sample=Ncar, thin=thinning, n.chains=3, n.cores=3, AR=2)

fit = CarBayesSt_lagged4
print(fit)

# save(f, final_data_lagged4, file='prediction_data.Rdata')

summary=fit$summary.results
save(summary, file='fit_beta_lagged3_max_summary_south_spline_full.RData')

fit_beta = CarBayesSt_lagged4$samples$beta
save(fit_beta, file='fit_beta_lagged3_max_south_spline_full.RData')

modelfit = CarBayesSt_lagged4$modelfit
save(modelfit, file='modelfit_lagged3_max_south_spline_full.RData')

## Start prediction/inference

source('carbayes_helpers.R')

final_data_lagged4[final_data_lagged4$timeid==60,'marlaria']=NA
final_data_lagged4[(nrow(final_data_lagged4)-1):nrow(final_data_lagged4),'marlaria']=NA
final_data_lagged4[nrow(final_data_lagged4),'marlaria']=NA

final_data_lagged4_aug = final_data_lagged4

frame.results <- common.frame(f, final_data_lagged4_aug, "poisson")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
Y.DA <- Y     




#### Check on the rho arguments
rho <- summary['rho.S','Mean']
fix.rho.S <- FALSE   

alpha <- c(summary['rho1.T','Mean'], summary['rho2.T','Mean'])
fix.rho.T <- FALSE   

#### CAR quantities
W.quants <- common.Wcheckformat.leroux(map_mat_north)
K <- W.quants$n
N <- N.all / K
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin



## Compute the blocking structure for beta     
block.temp <- common.betablock(p)
beta.beg  <- block.temp[[1]]
beta.fin <- block.temp[[2]]
n.beta.block <- block.temp[[3]]
list.block <- as.list(rep(NA, n.beta.block*2))
for(r in 1:n.beta.block)
{
  list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
  list.block[[r+n.beta.block]] <- length(list.block[[r]])
}




#############################
#### Initial parameter values
#############################
BETA_N_NORTH <- 27
BETA_N_SOUTH <- 25

beta <- summary[1:BETA_N_SOUTH,'Mean']

print('Beta: ')
print(dim(beta))
print('X.standardised: ')
print(dim(X.standardised))

log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd = res.sd)
tau2 <- summary['tau2','Mean']


#### Specify matrix quantities
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   
fitted <- exp(as.numeric(offset.mat + regression.mat + phi.mat))


###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
burnin=burn.in.car
n.sample=Ncar
thin=thinning
n.keep <- floor((n.sample - burnin)/thin)
samples.fitted <- array(NA, c(n.keep, N.all))
if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))


#### Specify the Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.05
proposal.sd.beta <- 0.01



#############################
#### Specify spatial elements
#############################
#### Spatial determinant
if(!fix.rho.S) 
{
  Wstar <- diag(apply(W,1,sum)) - W
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
}else
{}



accept.all <- rep(0,6)
accept <- accept.all

#use last month for Y
if(n.miss>0)
{
  Y.DA[which.miss==0] <- tail(Y.DA[which.miss!=0],65)   
}else
{}
Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)

#### Create the MCMC samples
for(j in 1:n.sample) {
  
  ####################
  ## Sample from phi
  ####################
  phi.offset <- offset.mat + regression.mat
  den.offset <- rho * W.triplet.sum + 1 - rho
  temp1 <- poissonar2carupdateRW(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, alpha[1], alpha[2], rho, Y.DA.mat, proposal.sd.phi, phi.offset, den.offset)      
  phi.temp <- temp1[[1]]
  phi <- as.numeric(phi.temp)  - mean(as.numeric(phi.temp))
  phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
  
  
  lp <- as.numeric(offset.mat + regression.mat + phi.mat)
  fitted <- exp(lp)
  
  ###################
  ## Save the results
  ###################
  if(j > burnin & (j-burnin)%%thin==0)
  {
    ele <- (j - burnin) / thin
    samples.fitted[ele, ] <- fitted
  }else
  {}
  
  print(j)
}

fitted.values <- apply(samples.fitted, 2, mean)
fitted.values_mat = matrix(fitted.values, nrow=K, ncol=N, byrow=FALSE)
fitted_pred = fitted.values_mat[,56]

save(fitted_pred, file = "fitted_pred_south")
