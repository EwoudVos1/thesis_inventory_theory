#Cite R
citation()
#Cite RStudio
RStudio.Version()



#Run the functions
Lambda <- function(t){
  if (t<=0){
    b*(t-tmin)
  } else if (t>0) {
    b*(t-tmin)+.5*a*t^2
  }
}

lambda1 <- function(t){
  if (t>=0){ 
    b+a*t #linear, so mean can be obtained halfway in period
  } else if (t<0){
    b #mean lambda over period
  } 
}


#aggregate demand distribution
aggr_dem <- function(x, arrivalparameter){
  N <- c()
  D <- c()
  for (k in 0:20) {
    N <- c(N, dpois(k, arrivalparameter))
    D <- c(D, dpois(x, k*mu))
  }
  return(N %*% D)
}



#calc predictive holding, backorder, and writeoff costs
exp_cost <- function(t, s, S) {
  if (s + t > tmax) {
    s_adj <- tmax - t
  } else {
    s_adj <- s
  }
  #calculate demand probabilities for zero up to maximum demand of interest S-1
  arrivalparam <- integrate(Vectorize(lambda1), t, t + s_adj)$value
  
  f_d <- seq(0, max(0, S - 1))
  f_d <- mapply(function(i, Lambda) aggr_dem(i, Lambda), i=f_d, Lambda=arrivalparam)
  
  ILplus <- 0
  for (i in seq(0, max(0, S - 1))) {
    ILplus <- ILplus + (S - i) * f_d[(i + 1)]
  }
  return((H + B) * ILplus - B * (S - arrivalparam * mu) + K * S * (s + t > tmax))
}


#calc predictive holding, backorder, and writeoff costs
exp_cost1 <- function(t, s, S) {
  if (s + t > tmax) {
    s_adj <- tmax - t
  } else {
    s_adj <- s
  }
  #calculate demand probabilities for zero up to maximum demand of interest S-1
  arrivalparam <- integrate(Vectorize(lambda1), t, t + s_adj)$value
  
  f_d <- seq(0, max(0, S - 1))
  f_d <- mapply(function(i, Lambda) aggr_dem(i, Lambda), i=f_d, Lambda=arrivalparam)
  
  ILplus <- 0
  for (i in seq(0, max(0, S - 1))) {
    ILplus <- ILplus + (S - i) * f_d[(i + 1)]
  }
  return((H + B) * ILplus - B * (S - arrivalparam * mu))
}


#simulation arrival using thinning
get_nhpp_realization <- function(mu){
  t <- tmin
  lambda_star <- max(sapply(seq(tmin, tmax, length.out=1000), lambda1))
  
  X <- numeric(0)
  while (t <= tmax) {
    e <- rexp(1, lambda_star)
    u <- runif(1)
    t <- t + e
    if (t > tmax) {
      break
    } else if (u < lambda1(t)/lambda_star){
      X <- c(X, t)
    }
  }
  X <- round(X, 2)
  
  Xdem <- rpois(length(X), mu)
  
  if (length(X) > 0) {
    return(data.frame("t_arrival"=X, "demand"=Xdem))
  } else {
    return(0)
  }
}



#### ordering fixed framework/S = 0 if D_t(t_max - t) \leq S ####
orderpolicy_fixed <- function(get_nhpp){
  #setup IL satisfies lead time demand in stationary D(tmin,tmin+L) = lambda(tmin, L)*L*mu
  #setup initial IL and IP
  IL <- integrate(Vectorize(lambda1), tmin, tmin+L)$value*mu 
  IP <- IL
  # outstanding_orders <- 0
  
  arrival_data <- rbind(c(tmin,0), get_nhpp, c(tmax,0))
  orders <- data.frame()
  row_indicator <- 1
  
  #iterate for optimal S
  j=0
  S_opt <- exp_cost(tmin, L, 0)
  S_next <- exp_cost(tmin, L, 5)
  while (S_opt > S_next) { #iteration with 5 steps
    j = j + 5
    S_opt <- exp_cost(tmin, L, j)
    S_next <- exp_cost(tmin, L, j+5) 
  }
  
  if (j==0){
    S_opt <- exp_cost(tmin, L, j)
    S_next <- exp_cost(tmin, L, j+1)
    while (S_opt > S_next) { #iterate upwards with 1 step
      j = j + 1
      S_opt <- exp_cost(tmin, L, j)
      S_next <- exp_cost(tmin, L, j+1) 
    }
  } else {
    index <- which.min(c(exp_cost(tmin, L, j-1), S_opt, exp_cost(tmin, L, j+1)))
    if (index==1){
      S_opt <- exp_cost(tmin, L, j)
      S_next <- exp_cost(tmin, L, j-1)
      while (S_opt > S_next) { #iterate to below with 1 step
        j = j - 1
        S_opt <- exp_cost(tmin, L, j)
        S_next <- exp_cost(tmin, L, j-1) 
      }
    } else if (index==3){
      S_opt <- exp_cost(tmin, L, j)
      S_next <- exp_cost(tmin, L, j+1)
      while (S_opt > S_next) { #iterate upwards with 1 step
        j = j + 1
        S_opt <- exp_cost(tmin, L, j)
        S_next <- exp_cost(tmin, L, j+1) 
      }
    }
  }
  
  for (i in arrival_data[,1]) {
    #update IL/IP/outstanding orders if order arrives
    arr <- orders$arrives_before %in% i
    if (any(arr==T)) {
      IL <- IL + sum(orders$ordersize[arr])
      # outstanding_orders <- outstanding_orders - sum(orders$ordersize[arr])
    }
    
    #subtract demand from IP and IL
    IL <- IL - arrival_data$demand[row_indicator]
    IP <- IP - arrival_data$demand[row_indicator]
    
    if (integrate(Vectorize(lambda1), i, tmax)$value*mu <= j){
      j <- 0
    }
    
    #update IP and store ordersize and ordertime
    if (IP < j) {
      orders <- rbind(orders, c(j-IP, 
                                i, 
                                min(Find(function(x) x>=(i+L), arrival_data$t_arrival), tmax),
                                j,
                                S_opt,
                                IL,
                                IP))
      # outstanding_orders <- outstanding_orders + j - IP
      IP <- j
    } else {
      orders <- rbind(orders, c(0, 
                                i, 
                                min(Find(function(x) x>=(i+L), arrival_data$t_arrival), tmax),
                                j,
                                S_opt,
                                IL,
                                IP))
    }
    
    if (i==tmin) {
      colnames(orders) <- c("ordersize", "ordertime", "arrives_before", 
                            "S_opt(t+L)", "E[c(S, t+L)]", "IL", "IP")
    }
    row_indicator <- row_indicator + 1
  }
  #orders$demand <- arrival_data$demand
  
  #cost evaluation part
  #create sorted arrivals of demand and orders
  time_changes_merged <- c(arrival_data$t_arrival, orders$ordertime + L)
  time_changes <- sort(time_changes_merged) #at this times
  
  inv_changes <- c(-arrival_data$demand, orders$ordersize)
  inv_changes <- inv_changes[order(time_changes_merged)] #this increases will be in IL
  
  #recalculate initial demand
  IL <- integrate(Vectorize(lambda1), tmin, tmin+L)$value*mu 
  
  IL_t <- c()
  for (i in 1:length(inv_changes)) {
    IL <- IL + inv_changes[i]
    IL_t <- c(IL_t, IL)
  }
  
  dt <- c(diff(time_changes),0)
  
  cost_realization <- data.frame(time_changes, inv_changes, IL_t, "backorder_amount"=0, dt, 
                                 "holding_cost"=0, "backorder_cost"=0, "write_off_cost"=0)
  
  for (i in 1:nrow(cost_realization)) {
    if (cost_realization$IL_t[[i]] > 0) {
      cost_realization$holding_cost[[i]] <- H*cost_realization$IL_t[[i]]*cost_realization$dt[[i]]
    } else if (cost_realization$IL_t[[i]] < 0) {
      cost_realization$backorder_cost[[i]] <- -B*cost_realization$IL_t[[i]]*cost_realization$dt[[i]]
      if (cost_realization$inv_changes[[i]]< 0) {
        cost_realization$backorder_amount[[i]] <- min(-cost_realization$inv_changes[[i]],-cost_realization$IL_t[[i]])
      } 
    }
    
    if (i==nrow(cost_realization)){
      end_inv <- cost_realization$IL_t[[i]]
      cost_realization$write_off_cost[[i]] <- K*end_inv
    }
  }
  cost_realization$total_cost <- cost_realization$holding_cost + 
    cost_realization$backorder_cost + 
    cost_realization$write_off_cost
  
  #achieved fill rate calculation
  #fill_rate <- (sum(arrival_data$demand) - sum(cost_realization$backorder_amount))/sum(arrival_data$demand)
  
  return(c(sum(arrival_data$demand), colSums(cost_realization[6:9]), end_inv, sum(cost_realization$backorder_amount))) 
}

#### ordering dynamic framework ####

#usable for dynamic system with write off forecast 
orderpolicy_dynamic1 <- function(get_nhpp){
  #inventory management part
  #setup IL satisfies lead time demand in stationary D(tmin,tmin+L) = lambda(tmin, L)*L*mu
  #setup initial IL and IP
  IL <- integrate(Vectorize(lambda1), tmin, tmin+L)$value*mu 
  IP <- IL
  # outstanding_orders <- 0
  
  arrival_data <- rbind(c(tmin,0), get_nhpp, c(tmax,0))
  orders <- data.frame()
  row_indicator <- 1
  
  for (i in arrival_data[,1]) {
    #update IL/IP/outstanding orders if order arrives
    arr <- orders$arrives_before %in% i
    if (any(arr==T)) {
      IL <- IL + sum(orders$ordersize[arr])
      # outstanding_orders <- outstanding_orders - sum(orders$ordersize[arr])
    }
    
    #subtract demand from IP and IL
    IL <- IL - arrival_data$demand[row_indicator]
    IP <- IP - arrival_data$demand[row_indicator]
    
    #iterate for optimal S
    j<-0
    S_opt <- exp_cost(i, L, 0)
    S_next <- exp_cost(i, L, 5)
    while (S_opt > S_next) { #iteration with 5 steps
      j = j + 5
      S_opt <- exp_cost(i, L, j)
      S_next <- exp_cost(i, L, j+5)
    }
    
    if (j==0){
      S_opt <- exp_cost(i, L, j)
      S_next <- exp_cost(i, L, j+1)
      while (S_opt > S_next) { #iterate upwards with 1 step
        j <- j + 1
        S_opt <- exp_cost(i, L, j)
        S_next <- exp_cost(i, L, j+1)
      }
    } else {
      index <- which.min(c(exp_cost(i, L, j-1), S_opt, exp_cost(i, L, j+1)))
      if (index==1){
        S_opt <- exp_cost(i, L, j)
        S_next <- exp_cost(i, L, j-1)
        while (S_opt > S_next) { #iterate to below with 1 step
          j <- j - 1
          S_opt <- exp_cost(i, L, j)
          S_next <- exp_cost(i, L, j-1)
        }
      } else if (index==3){
        S_opt <- exp_cost(i, L, j)
        S_next <- exp_cost(i, L, j+1)
        while (S_opt > S_next) { #iterate upwards with 1 step
          j = j + 1
          S_opt <- exp_cost(i, L, j)
          S_next <- exp_cost(i, L, j+1)
        }
      }
    }
    
    
    #update IP and store ordersize and ordertime
    if (IP < j) {
      orders <- rbind(orders, c(j-IP, 
                                i, 
                                min(Find(function(x) x>=(i+L), arrival_data$t_arrival), tmax),
                                j,
                                S_opt,
                                IL,
                                IP))
      # outstanding_orders <- outstanding_orders + j - IP
      IP <- j
    } else {
      orders <- rbind(orders, c(0, 
                                i, 
                                min(Find(function(x) x>=(i+L), arrival_data$t_arrival), tmax),
                                j,
                                S_opt,
                                IL,
                                IP))
    }
    
    if (i==tmin) {
      colnames(orders) <- c("ordersize", "ordertime", "arrives_before", 
                            "S_opt(t+L)", "E[c(S, t+L)]", "IL", "IP")
    }
    row_indicator <- row_indicator + 1
  }
  
  #cost evaluation part
  #create sorted arrivals of demand and orders
  time_changes_merged <- c(arrival_data$t_arrival, orders$ordertime + L)
  time_changes <- sort(time_changes_merged) #at this times
  
  inv_changes <- c(-arrival_data$demand, orders$ordersize)
  inv_changes <- inv_changes[order(time_changes_merged)] #this increases will be in IL
  
  #recalculate initial demand
  IL <- integrate(Vectorize(lambda1), tmin, tmin+L)$value*mu 
  
  IL_t <- c()
  for (i in 1:length(inv_changes)) {
    IL <- IL + inv_changes[i]
    IL_t <- c(IL_t, IL)
  }
  
  dt <- c(diff(time_changes),0)
  
  cost_realization <- data.frame(time_changes, inv_changes, IL_t, "backorder_amount"=0, dt, 
                                 "holding_cost"=0, "backorder_cost"=0, "write_off_cost"=0)
  
  for (i in 1:nrow(cost_realization)) {
    if (cost_realization$IL_t[[i]] > 0) {
      cost_realization$holding_cost[[i]] <- H*cost_realization$IL_t[[i]]*cost_realization$dt[[i]]
    } else if (cost_realization$IL_t[[i]] < 0) {
      cost_realization$backorder_cost[[i]] <- -B*cost_realization$IL_t[[i]]*cost_realization$dt[[i]]
      if (cost_realization$inv_changes[[i]]< 0) {
        cost_realization$backorder_amount[[i]] <- min(-cost_realization$inv_changes[[i]],-cost_realization$IL_t[[i]])
      } 
    }
    
    if (i==nrow(cost_realization)){
      end_inv <- cost_realization$IL_t[[i]]
      cost_realization$write_off_cost[[i]] <- K*end_inv
    }
  }
  cost_realization$total_cost <- cost_realization$holding_cost + 
    cost_realization$backorder_cost + 
    cost_realization$write_off_cost
  
  #achieved fill rate calculation
  #fill_rate <- (sum(arrival_data$demand) - sum(cost_realization$backorder_amount))/sum(arrival_data$demand)
  
  return(c(sum(arrival_data$demand), colSums(cost_realization[6:9]), end_inv, sum(cost_realization$backorder_amount)))
}


#dynamic system cutoff S = 0 if D_t(t_max - t) \leq S(t)
orderpolicy_dynamic4 <- function(get_nhpp){
  #inventory management part
  #setup IL satisfies lead time demand in stationary D(tmin,tmin+L) = lambda(tmin, L)*L*mu
  #setup initial IL and IP
  IL <- integrate(Vectorize(lambda1), tmin, tmin+L)$value*mu 
  IP <- IL
  # outstanding_orders <- 0
  
  arrival_data <- rbind(c(tmin,0), get_nhpp, c(tmax,0))
  orders <- data.frame()
  row_indicator <- 1
  stop_iterate <- F
  
  for (i in arrival_data[,1]) {
    #update IL/IP/outstanding orders if order arrives
    arr <- orders$arrives_before %in% i
    if (any(arr==T)) {
      IL <- IL + sum(orders$ordersize[arr])
      # outstanding_orders <- outstanding_orders - sum(orders$ordersize[arr])
    }
    
    #subtract demand from IP and IL
    IL <- IL - arrival_data$demand[row_indicator]
    IP <- IP - arrival_data$demand[row_indicator]
    
    if (stop_iterate==F) {
      #iterate for optimal S
      j <- 0
      S_opt <- exp_cost(i, L, 0)
      S_next <- exp_cost(i, L, 5)
      while (S_opt > S_next) {
        #iteration with 5 steps
        j = j + 5
        S_opt <- exp_cost(i, L, j)
        S_next <- exp_cost(i, L, j + 5)
      }
      
      if (j == 0) {
        S_opt <- exp_cost(i, L, j)
        S_next <- exp_cost(i, L, j + 1)
        while (S_opt > S_next) {
          #iterate upwards with 1 step
          j <- j + 1
          S_opt <- exp_cost(i, L, j)
          S_next <- exp_cost(i, L, j + 1)
        }
      } else {
        index <-
          which.min(c(exp_cost(i, L, j - 1), S_opt, exp_cost(i, L, j + 1)))
        if (index == 1) {
          S_opt <- exp_cost(i, L, j)
          S_next <- exp_cost(i, L, j - 1)
          while (S_opt > S_next) {
            #iterate to below with 1 step
            j <- j - 1
            S_opt <- exp_cost(i, L, j)
            S_next <- exp_cost(i, L, j - 1)
          }
        } else if (index == 3) {
          S_opt <- exp_cost(i, L, j)
          S_next <- exp_cost(i, L, j + 1)
          while (S_opt > S_next) {
            #iterate upwards with 1 step
            j = j + 1
            S_opt <- exp_cost(i, L, j)
            S_next <- exp_cost(i, L, j + 1)
          }
        }
      }
      
      #check if level of S should be forced to zero
      if (integrate(Vectorize(lambda1), i, tmax)$value*mu <= j){
        j <- 0
        stop_iterate <- T
      }
    }
    
    #update IP and store ordersize and ordertime
    if (IP < j) {
      orders <- rbind(orders, c(j-IP, 
                                i, 
                                min(Find(function(x) x>=(i+L), arrival_data$t_arrival), tmax),
                                j,
                                S_opt,
                                IL,
                                IP))
      # outstanding_orders <- outstanding_orders + j - IP
      IP <- j
    } else {
      orders <- rbind(orders, c(0, 
                                i, 
                                min(Find(function(x) x>=(i+L), arrival_data$t_arrival), tmax),
                                j,
                                S_opt,
                                IL,
                                IP))
    }
    
    if (i==tmin) {
      colnames(orders) <- c("ordersize", "ordertime", "arrives_before", 
                            "S_opt(t+L)", "E[c(S, t+L)]", "IL", "IP")
    }
    row_indicator <- row_indicator + 1
  }
  
  #cost evaluation part
  #create sorted arrivals of demand and orders
  time_changes_merged <- c(arrival_data$t_arrival, orders$ordertime + L)
  time_changes <- sort(time_changes_merged) #at this times
  
  inv_changes <- c(-arrival_data$demand, orders$ordersize)
  inv_changes <- inv_changes[order(time_changes_merged)] #this increases will be in IL
  
  #recalculate initial demand
  IL <- integrate(Vectorize(lambda1), tmin, tmin+L)$value*mu 
  
  IL_t <- c()
  for (i in 1:length(inv_changes)) {
    IL <- IL + inv_changes[i]
    IL_t <- c(IL_t, IL)
  }
  
  dt <- c(diff(time_changes),0)
  
  cost_realization <- data.frame(time_changes, inv_changes, IL_t, "backorder_amount"=0, dt, 
                                 "holding_cost"=0, "backorder_cost"=0, "write_off_cost"=0)
  
  for (i in 1:nrow(cost_realization)) {
    if (cost_realization$IL_t[[i]] > 0) {
      cost_realization$holding_cost[[i]] <- H*cost_realization$IL_t[[i]]*cost_realization$dt[[i]]
    } else if (cost_realization$IL_t[[i]] < 0) {
      cost_realization$backorder_cost[[i]] <- -B*cost_realization$IL_t[[i]]*cost_realization$dt[[i]]
      if (cost_realization$inv_changes[[i]]< 0) {
        cost_realization$backorder_amount[[i]] <- min(-cost_realization$inv_changes[[i]],-cost_realization$IL_t[[i]])
      } 
    }
    
    if (i==nrow(cost_realization)){
      end_inv <- cost_realization$IL_t[[i]]
      cost_realization$write_off_cost[[i]] <- K*end_inv
    }
  }
  cost_realization$total_cost <- cost_realization$holding_cost + 
    cost_realization$backorder_cost + 
    cost_realization$write_off_cost
  
  #achieved fill rate calculation
  #fill_rate <- (sum(arrival_data$demand) - sum(cost_realization$backorder_amount))/sum(arrival_data$demand)
  
  return(c(sum(arrival_data$demand), colSums(cost_realization[6:9]), end_inv, sum(cost_realization$backorder_amount)))
}

#set different scenarios lead times L
L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30




##### generate arrivaldata in list ######
mu <- 5;a<- -.03;b<- 1
nsim <- 300
summary_costs <- data.frame()
set.seed(123)
data_m5_ll <- list()
i <- 1
while (i<(nsim+1)) {
  data_m5_ll[[i]] <- get_nhpp_realization(mu)
  i <- i+ 1
}

mu <- 2;a<- -.03;b<- 1
nsim <- 300
summary_costs <- data.frame()
set.seed(123)
data_m2_ll <- list()
i <- 1
while (i<(nsim+1)) {
  data_m2_ll[[i]] <- get_nhpp_realization(mu)
  i <- i+ 1
}

# mu <- 5;a<- -.06;b<- 2
# nsim <- 300
# summary_costs <- data.frame()
# set.seed(123)
# data_m5_lh <- list()
# i <- 1
# while (i<(nsim+1)) {
#   data_m5_lh[[i]] <- get_nhpp_realization(mu)
#   i <- i+ 1
# }
# 
# mu <- 2;a<- -.06;b<- 2
# nsim <- 300
# summary_costs <- data.frame()
# set.seed(123)
# data_m2_lh <- list()
# i <- 1
# while (i<(nsim+1)) {
#   data_m2_lh[[i]] <- get_nhpp_realization(mu)
#   i <- i+ 1
# }
# 


####check convergence in costs#####
nsim <- 1000
sim_demand <- data.frame()
summary_demand <- data.frame()

set.seed(123)

for (j in 1:nsim) {
  sim_demand <- rbind(sim_demand, orderpolicy_fixed(data_m5_ll[[j]])[2:5])
  summary_demand <- rbind(summary_demand, c(sapply(sim_demand, mean), sapply(sim_demand, var)))
  print(j)
}


colnames(sim_demand) <- c("holding_cost", "backorder_cost", "write_off_cost",  "total_cost")
colnames(summary_demand) <- 
  c(paste0("m_", c("holding_cost", "backorder_cost", "write_off_cost",  "total_cost"), sep=""),
    paste0("v_", c("holding_cost", "backorder_cost", "write_off_cost",  "total_cost"), sep=""))

summary_demand$nsim <- 1:nsim

library(ggplot2)
library(gridExtra)
meanplot <- ggplot(summary_demand) + 
  geom_line(aes(x=nsim, y=m_total_cost)) + 
  theme_bw() +
  labs(x="Number of simulations", y="Estimated mean of total costs") +
  ylim(c(1500,2100))

varplot <- ggplot(summary_demand) + 
  geom_line(aes(x=nsim, y=v_total_cost)) +
  theme_bw() +
  labs(x="Number of simulations", y="Estimated variance of total costs") +
  ylim(c(600000,1200000))

grid.arrange(meanplot, varplot, nrow=1)





###### sim mu=5, lambda=low ####

L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30


library(writexl)
mu <- 5;a<- -.03;b<- 1
B<-2
nsim <- 300
summary_costs <- data.frame()
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  for (i in 1:nsim) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_fixed(data_m5_ll[[i]]), orderpolicy_dynamic1(data_m5_ll[[i]])[2:7],
                                      orderpolicy_dynamic4(data_m5_ll[[i]])[2:7]))
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu5_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

names <- c("holding_cost", "backorder_cost", "write_off_cost",
           "total_cost", "end_season_inv", "backorder")
colnames(summary_costs) <- c("L", "b", "mu", "tot_demand", 
                             paste("m", rep(1:3, times=1, each=6), rep(names,3), sep=""),
                             paste("v", rep(1:3, times=1, each=6), rep(names,3), sep=""))

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu5_lambdalow1.xlsx")



mu <- 5;a<- -.03;b<- 1
B<-5
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m5_ll[[i]]), orderpolicy_fixed(data_m5_ll[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m5_ll[[i]])[2:7], orderpolicy_dynamic2(data_m5_ll[[i]])[2:7]))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu5_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/mu5_lambdalow2.xlsx")

mu <- 5;a<- -.03;b<- 1
B<-10
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m5_ll[[i]]), orderpolicy_fixed(data_m5_ll[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m5_ll[[i]])[2:7], orderpolicy_dynamic2(data_m5_ll[[i]])[2:7]))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu5_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/mu5_lambdalow3.xlsx")





###### sim mu=2, lambda=low ####

L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30


mu <- 2;a<- -.03;b<- 1
B<-2
nsim <- 300
summary_costs <- data.frame()
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_fixed(data_m2_ll[[i]]), orderpolicy_dynamic1(data_m2_ll[[i]])[2:7],
                                      orderpolicy_dynamic4(data_m2_ll[[i]])[2:7]))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], sd)))
  print(c("L is ", j))
}

names <- c("holding_cost", "backorder_cost", "write_off_cost",
           "total_cost", "end_season_inv", "backorder")
colnames(summary_costs) <- c("L", "b", "mu", "tot_demand", 
                             paste("m", rep(1:3, times=1, each=6), rep(names,3), sep=""),
                             paste("v", rep(1:3, times=1, each=6), rep(names,3), sep=""))

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow1.xlsx")

mu <- 2;a<- -.03;b<- 1
B<-5
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m2_ll[[i]]), orderpolicy_fixed(data_m2_ll[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m2_ll[[i]])[2:7], orderpolicy_dynamic2(data_m2_ll[[i]])[2:7]))
    i <- i+1
    print(i)
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/mu2_lambdalow2.xlsx")

mu <- 2;a<- -.03;b<- 1
B<-10
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m2_ll[[i]]), orderpolicy_fixed(data_m2_ll[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m2_ll[[i]])[2:7], orderpolicy_dynamic2(data_m2_ll[[i]])[2:7]))
    i <- i+1
    print(i)
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}


write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/mu2_lambdalow3.xlsx")




###### sim mu=5, lambda=high ####

L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30


library(writexl)
mu <- 5;a<- -.06;b<- 2
B<-1
nsim <- 300
summary_costs <- data.frame()
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  for (i in 1:nsim) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m5_lh[[i]]), orderpolicy_fixed(data_m5_lh[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m5_lh[[i]])[2:7], orderpolicy_dynamic2(data_m5_lh[[i]])[2:7]))
  }
  path <- paste("sim_results/mu5_lambdahigh_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

names <- c("holding_cost", "backorder_cost", "write_off_cost",
           "total_cost", "end_season_inv", "backorder")
colnames(summary_costs) <- c("L", "b", "mu", "tot_demand", 
                             paste("m", rep(1:4, times=1, each=6), rep(names,4), sep=""),
                             paste("v", rep(1:4, times=1, each=6), rep(names,4), sep=""))

write_xlsx(summary_costs, "mu5_lambdahigh1.xlsx")



mu <- 5;a<- -.06;b<- 2
B<-5
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m5_lh[[i]]), orderpolicy_fixed(data_m5_lh[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m5_lh[[i]])[2:7], orderpolicy_dynamic2(data_m5_lh[[i]])[2:7]))
    i <- i+1
  }
  path <- paste("sim_results/mu5_lambdahigh_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "mu5_lambdahigh2.xlsx")

mu <- 5;a<- -.06;b<- 2
B<-10
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m5_lh[[i]]), orderpolicy_fixed(data_m5_lh[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m5_lh[[i]])[2:7], orderpolicy_dynamic2(data_m5_lh[[i]])[2:7]))
    i <- i+1
  }
  path <- paste("sim_results/mu5_lambdahigh_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "mu5_lambdahigh3.xlsx")





###### sim mu=2, lambda=high ####

L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30


mu <- 2;a<- -.06;b<- 2
B<-1
nsim <- 300
summary_costs <- data.frame()
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m2_lh[[i]]), orderpolicy_fixed(data_m2_lh[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m2_lh[[i]])[2:7], orderpolicy_dynamic2(data_m2_lh[[i]])[2:7]))
    i <- i+1
  }
  path <- paste("sim_results/mu2_lambdahigh_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

names <- c("holding_cost", "backorder_cost", "write_off_cost",
           "total_cost", "end_season_inv", "backorder")
colnames(summary_costs) <- c("L", "b", "mu", "tot_demand", 
                             paste("m", rep(1:4, times=1, each=6), rep(names,4), sep=""),
                             paste("v", rep(1:4, times=1, each=6), rep(names,4), sep=""))

write_xlsx(summary_costs, "mu2_lambdahigh1.xlsx")

mu <- 2;a<- -.06;b<- 2
B<-5
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m2_lh[[i]]), orderpolicy_fixed(data_m2_lh[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m2_lh[[i]])[2:7], orderpolicy_dynamic2(data_m2_lh[[i]])[2:7]))
    i <- i+1
    print(i)
  }
  path <- paste("sim_results/mu2_lambdahigh_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "mu2_lambdahigh2.xlsx")

mu <- 2;a<- -.06;b<- 2
B<-10
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_base(data_m2_lh[[i]]), orderpolicy_fixed(data_m2_lh[[i]])[2:7],
                                      orderpolicy_dynamic1(data_m2_lh[[i]])[2:7], orderpolicy_dynamic2(data_m2_lh[[i]])[2:7]))
    i <- i+1
    print(i)
  }
  path <- paste("sim_results/mu2_lambdahigh_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}


write_xlsx(summary_costs, "mu2_lambdahigh3.xlsx")




##### sim third method separately #####

#mu=5, lambda=low
L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30


library(writexl)
mu <- 5;a<- -.03;b<- 1
B<-1
nsim <- 300
summary_costs <- data.frame()
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_dynamic4(data_m5_ll[[i]])))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu5_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

names <- c("holding_cost", "backorder_cost", "write_off_cost",
           "total_cost", "end_season_inv", "backorder")
colnames(summary_costs) <- c("L", "b", "mu", "tot_demand", 
                             paste("m", names, sep=""),
                             paste("v", names, sep=""))

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu5_low1.xlsx")



mu <- 5;a<- -.03;b<- 1
B<-5
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_dynamic4(data_m5_ll[[i]])))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu5_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu5_low2.xlsx")

mu <- 5;a<- -.03;b<- 1
B<-10
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_dynamic4(data_m5_ll[[i]])))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu5_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu5_low3.xlsx")





###### sim mu=2, lambda=low ####

L0 <- c(1,5, 10, 15)
#fix holding
H <- 1; 
#cost of write off per unit
K=50
#end of period
tmin = -10
tmax = 30


mu <- 2;a<- -.03;b<- 1
B<-1
nsim <- 300
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_dynamic4(data_m2_ll[[i]])))
    i <- i+1
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

# names <- c("holding_cost", "backorder_cost", "write_off_cost",
#            "total_cost", "end_season_inv", "backorder")
# colnames(summary_costs) <- c("L", "b", "mu", "tot_demand", 
#                              paste("m", rep(1:4, times=1, each=6), rep(names,4), sep=""),
#                              paste("v", rep(1:4, times=1, each=6), rep(names,4), sep=""))

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu2_low1.xlsx")

mu <- 2;a<- -.03;b<- 1
B<-5
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_dynamic4(data_m2_ll[[i]])))
    i <- i+1
    print(i)
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}

write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu2_low2.xlsx")

mu <- 2;a<- -.03;b<- 1
B<-10
for (j in 1:4) {#change scenarios for lead time
  L <<- L0[j]
  sim_demand <- data.frame()
  i <- 1
  while (i<(nsim+1)) {
    sim_demand <- rbind(sim_demand, c(orderpolicy_dynamic4(data_m2_ll[[i]])))
    i <- i+1
    print(i)
  }
  path <- paste("/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/mu2_lambdalow_",B, L,".xlsx",sep = "")
  write_xlsx(sim_demand, path)
  summary_costs <- rbind(summary_costs, c(L,B,mu, sapply(sim_demand, mean), sapply(sim_demand[-1], var)))
  print(c("L is ", j))
}


write_xlsx(summary_costs, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/method3/mu2_low3.xlsx")
















###### Optimal S over time plots #####

#S opt max for L=15, B=10 and mu=5 is 106, for mu=2 is 45

#for mu=5
minS <- function(T) {
  which.min(sapply(0:110, exp_cost, t=T, s=L))
}

minS1 <- function(T) {
  which.min(sapply(0:50, exp_cost, t=T, s=L))
}
Sopt5 <- data.frame(time=seq(tmin, tmax, by=0.1))

mu=5
L=15;B=10
Sopt5$L15_B10 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=10;B=10
Sopt5$L10_B10 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=5;B=10
Sopt5$L5_B10 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=1;B=10
Sopt5$L1_B10 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=15;B=5
Sopt5$L15_B5 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=10;B=5
Sopt5$L10_B5 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=5;B=5
Sopt5$L5_B5 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=1;B=5
Sopt5$L1_B5 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=15;B=2
Sopt5$L15_B2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=10;B=2
Sopt5$L10_B2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=5;B=2
Sopt5$L5_B2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=1;B=2
Sopt5$L1_B2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

write_xlsx(Sopt5, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/Sopt_mu5.xlsx")
library(writexl)



mu=2
L=15;B=10
Sopt5$L15_B10_2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=10;B=10
Sopt5$L10_B10_2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=5;B=10
Sopt5$L5_B10_2<- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=1;B=10
Sopt5$L1_B10_2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=15;B=5
Sopt5$L15_B5_2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=10;B=5
Sopt5$L10_B5_2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=5;B=5
Sopt5$L5_B5_2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=1;B=5
Sopt5$L1_B5_2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=15;B=2
Sopt5$L15_B2_2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=10;B=2
Sopt5$L10_B2_2 <- sapply(seq(tmin, tmax, by=0.1), minS)-1

L=5;B=2
Sopt5$L5_B2_2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

L=1;B=2
Sopt5$L1_B2_2 <- sapply(seq(tmin, tmax, by=0.1), minS1)-1

library(writexl)
write_xlsx(Sopt5, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/Sopt_mu5.xlsx")




##### Determine S for base strategy code #####

cutoff <- function(s, S){
  mu*integrate(Vectorize(lambda1), s, tmax)$value - S
}

#if function cutoff is \leq zero, then expected future demand is less then S

Sinput <- Sopt5$L5_B5[1]
timing <- seq(tmin, tmax, by=0.1)
ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)


###### Determine S for strategy 3 code #####
Sinput <- Sopt5$L5_B5
timing <- seq(tmin, tmax, by=0.1)
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
S3_opt
Sinput


##### Graphical representation S strategy 2 #####
#In previous chapter dataset called Sopt5 is created with optimal S over time for scenarios where mu=5
library(ggplot2)
Sopt5$time <- seq(0,40, by=0.1)

#plot for L=5, mu=5 strategy 2
ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=L5_B10, colour="gr1")) +
  geom_step(aes(y=L5_B5, colour="gr2")) +
  geom_step(aes(y=L5_B2, colour ="gr3")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue"),
    labels =  c("b = 10", "b = 5", "b = 2")
  ) 



#plot for L=5, b=5 strategy 2
ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=L5_B5, colour="gr1")) +
  geom_step(aes(y=L5_B5_2, colour="gr2")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, b = 5",
    values = c("gr1"="navy", "gr2"="magenta"),
    labels = c("μ = 5", "μ = 2")
  )


#plot for b=5, mu=5 strategy 2
ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=L15_B5, colour="gr1")) +
  geom_step(aes(y=L10_B5, colour="gr2")) +
  geom_step(aes(y=L5_B5, colour ="gr3")) +
  geom_step(aes(y=L1_B5, colour="gr4")) +
  theme_bw() +
  scale_colour_manual(
    name = "b=5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue", "gr4"="green"),
    labels = c("L = 15", "L = 10", "L = 5", "L = 1")
  )


#### Plot strategy 1, 2, 3 combined #####

#create data mu=5, L=5, b=5
Sinput <- Sopt5$L5_B5[1]
timing <- seq(tmin, tmax, by=0.1)
#strategy 1
strategy1 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#strategy 3
Sinput <- Sopt5$L5_B5
strategy3 <- c()
for (i in 1:length(Sinput)) {
  strategy3 <- c(strategy3, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}

Sopt5$str1_L5_B5 <- strategy1
Sopt5$str3_L5_B5 <- strategy3



#plot for L=5, mu=5 strategy 1,2,3
ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str1_L5_B5, colour="gr1"), size=0.5) +
  geom_step(aes(y=L5_B5, colour="gr2"), size=0.9) +
  geom_step(aes(y=str3_L5_B5, colour ="gr3"), size=0.5) +
  theme_bw() +
  scale_colour_manual(
    name = "",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue"),
    labels = c("Strategy 1", "Strategy 2", "Strategy 3")
  ) + labs(y="order-up-to level S")



#### Create S opt for strategy 1 #####
#create data for b=5, L=5,mu=2
mu=2
Sinput <- Sopt5$L5_B5_2[1]
timing <- seq(tmin, tmax, by=0.1)
Sopt5$str1_L5_B5_2 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#L=10
Sinput <- Sopt5$L10_B5_2[1]
Sopt5$str1_L10_B5_2 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)


#create data for L=5,mu=5 different values b
mu=5
#b=5
Sinput <- Sopt5$L5_B5[1]
Sopt5$str1_L5_B5 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#b=2
Sinput <- Sopt5$L5_B2[1]
Sopt5$str1_L5_B2 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#b=10
Sinput <- Sopt5$L5_B10[1]
Sopt5$str1_L5_B10 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#create data for b=5 , mu=5 different values L
#L=5 already exists
#L=1
Sinput <- Sopt5$L1_B5[1]
Sopt5$str1_L1_B5 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#L=10
Sinput <- Sopt5$L10_B5[1]
Sopt5$str1_L10_B5 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)

#L=15
Sinput <- Sopt5$L15_B5[1]
Sopt5$str1_L15_B5 <- ifelse(sapply(timing, cutoff, S=Sinput)>0, Sinput, 0)




###### Create S opt for strategy 3  #####
#L=5, b=5, mu change
mu=5
Sinput <- Sopt5$L5_B5
timing <- seq(tmin, tmax, by=0.1)
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}

Sopt5$str3_L5_B5 <- S3_opt

mu=2
Sinput <- Sopt5$L5_B5_2
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}
Sopt5$str3_L5_B5_2 <- S3_opt

#L=5, mu=5, change b, b=5 already exists
mu=5 #reset to mu=5 again

#b=2
Sinput <- Sopt5$L5_B2
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}

Sopt5$str3_L5_B2 <- S3_opt

#b=10
Sinput <- Sopt5$L5_B10
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}

Sopt5$str3_L5_B10 <- S3_opt


#mu=5, b=5, change L L=5 exists
#L=15
Sinput <- Sopt5$L15_B5
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}
Sopt5$str3_L15_B5 <- S3_opt

#L=10
Sinput <- Sopt5$L10_B5
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}

Sopt5$str3_L10_B5 <- S3_opt

#L=1
Sinput <- Sopt5$L1_B5
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}



Sopt5$str3_L1_B5 <- S3_opt


#mu =2 L=10
mu=2
Sinput <- Sopt5$L10_B5_2
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}

Sopt5$str3_L10_B5_2 <- S3_opt

#mu =2 L=10
mu=2
Sinput <- Sopt5$L15_B5_2
S3_opt <- c()
for (i in 1:length(Sinput)) {
  S3_opt <- c(S3_opt, ifelse(cutoff(timing[i],Sinput[i])>0, Sinput[i], 0))
}
for (i in 2:length(Sinput)) { #bruteforce S to zero when S has been made 0 before
  if (S3_opt[i-1]==0){
    S3_opt[i] <- 0
  }
}

Sopt5$str3_L15_B5_2 <- S3_opt


write_xlsx(Sopt5, "/Users/ewoudvos/Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/Sopt5_total.xlsx")





##### Creating referencing plots ####
library(ggplot2)
library(ggpubr)
#plot for L=5, b=5 
#strategy 1
p1a <- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str1_L5_B5, colour="gr1")) +
  geom_step(aes(y=str1_L5_B5_2, colour="gr2")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, b = 5",
    values = c("gr1"="navy", "gr2"="magenta"),
    labels = c("μ = 5", "μ = 2")
  ) +
  labs(y="order-up-to level S", title = "Strategy 1")
p1a


#strategy 2
p1b <- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=L5_B5, colour="gr1")) +
  geom_step(aes(y=L5_B5_2, colour="gr2")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, b = 5",
    values = c("gr1"="navy", "gr2"="magenta"),
    labels = c("μ = 5", "μ = 2")
  )+
  labs(y="order-up-to level S", title = "Strategy 2")
p1b

#strategy 3
p1c <- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str3_L5_B5, colour="gr1")) +
  geom_step(aes(y=str3_L5_B5_2, colour="gr2")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, b = 5",
    values = c("gr1"="navy", "gr2"="magenta"),
    labels = c("μ = 5", "μ = 2")
  )+
  labs(y="order-up-to level S", title = "Strategy 3")+ylim(c(0,85))
p1c


plot_mu <- ggarrange(p1a, p1b, p1c, ncol=3, common.legend=T, legend="right")



#plot for mu=5, L=5
#strategy 1
p2a<- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str1_L5_B10, colour="gr1")) +
  geom_step(aes(y=str1_L5_B5, colour="gr2")) +
  geom_step(aes(y=str1_L5_B2, colour ="gr3")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue"),
    labels =  c("b = 10", "b = 5", "b = 2")
  ) +
  labs(y="order-up-to level S", title = "Strategy 1")




#strategy 2
p2b <- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=L5_B10, colour="gr1")) +
  geom_step(aes(y=L5_B5, colour="gr2")) +
  geom_step(aes(y=L5_B2, colour ="gr3")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue"),
    labels =  c("b = 10", "b = 5", "b = 2")
  ) +
  labs(y="order-up-to level S", title = "Strategy 2")
p2b

#strategy 3
p2c<- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str3_L5_B10, colour="gr1")) +
  geom_step(aes(y=str3_L5_B5, colour="gr2")) +
  geom_step(aes(y=str3_L5_B2, colour ="gr3")) +
  theme_bw() +
  scale_colour_manual(
    name = "L = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue"),
    labels =  c("b = 10", "b = 5", "b = 2")
  ) +
  labs(y="order-up-to level S", title = "Strategy 3")

plot_b <- ggarrange(p2a, p2b, p2c, ncol=3, common.legend=T, legend="right")



#plot for b=5, mu=5 
#strategy 1
p3a <- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str1_L15_B5, colour="gr1")) +
  geom_step(aes(y=str1_L10_B5, colour="gr2")) +
  geom_step(aes(y=str1_L5_B5, colour ="gr3")) +
  geom_step(aes(y=str1_L1_B5, colour="gr4")) +
  theme_bw() +
  scale_colour_manual(
    name = "b = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue", "gr4"="red"),
    labels = c("L = 15", "L = 10", "L = 5", "L = 1")
  )+
  labs(y="order-up-to level S", title = "Strategy 1")



#strategy 2
p3b<- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=L15_B5, colour="gr1")) +
  geom_step(aes(y=L10_B5, colour="gr2")) +
  geom_step(aes(y=L5_B5, colour ="gr3")) +
  geom_step(aes(y=L1_B5, colour="gr4")) +
  theme_bw() +
  scale_colour_manual(
    name = "b = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue", "gr4"="red"),
    labels = c("L = 15", "L = 10", "L = 5", "L = 1")
  )+
  labs(y="order-up-to level S", title = "Strategy 2")

#strategy 3
p3c <- ggplot(Sopt5, aes(x=time)) + geom_step(aes(y=str3_L15_B5, colour="gr1")) +
  geom_step(aes(y=str3_L10_B5, colour="gr2")) +
  geom_step(aes(y=str3_L5_B5, colour ="gr3")) +
  geom_step(aes(y=str3_L1_B5, colour="gr4")) +
  theme_bw() +
  scale_colour_manual(
    name = "b = 5, μ = 5",
    values = c("gr1"="navy", "gr2"="magenta", "gr3"="blue", "gr4"="red"),
    labels = c("L = 15", "L = 10", "L = 5", "L = 1")
  )+
  labs(y="order-up-to level S", title = "Strategy 3")


plot_L <- ggarrange(p3a, p3b, p3c, ncol=3, common.legend=T, legend="right")

ggarrange(plot_mu, plot_b, plot_L, ncol=1)
plot_mu
plot_b
plot_L




##### Visualize end of season inventory ####
end_season <- read_excel("Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/end_season_inv.xlsx", 
                             sheet = "L15")
end_season <- cbind(end_season, read_excel("Dropbox/EOR 3/Thesis/R/simulation_data/sim_results/end_season_inv.xlsx", 
                             sheet = "L10"))
#etc

table(end_season$mu5_b2_L15_1)[1]



##### Table cost performance to latex #####
table_tekst <- read_excel("Dropbox/EOR 3/Thesis/R/simulation_data/filtered output.xlsx", 
                         sheet = "table_tekst")
table_tekst[4:9] <- round(table_tekst[4:9])
table_tekst[10:14] <- round(table_tekst[10:14], 1)
table1 <- as.data.frame(table_tekst[1:14])

library(knitr)
library(kableExtra)
kable(table1, "latex", booktabs=T) %>%
  add_header_above(c(" "=3, "Total cost" = 6, "Cost decrease in %" = 2, "Achieved fill rate in %"))



table_appendix<- read_excel("Dropbox/EOR 3/Thesis/R/simulation_data/filtered output.xlsx", 
                          sheet = "table_appendix")
table_appendix[4:9] <- round(table_appendix[4:9])
table_appendix[12:17] <- round(table_appendix[12:17])
table_appendix[10:11] <- round(table_appendix[10:11],1)
table_appendix[18:19] <- round(table_appendix[18:19],1)


kable(table_appendix, "latex", booktabs=T) %>%
  add_header_above(c(" "=3, "Backorder cost" = 8, "Write off cost" = 8))


