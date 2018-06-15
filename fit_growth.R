
fit.growth <-  function(data, fit.all = FALSE, time_var ="time_sowing", only.MSE  = FALSE, delta_A_max = 0.2, min_points = 3){
  
  #data: data, with values concerning just one variable (ex: height_crop)
  #fit.all: should the function fit all points of the curve or just the points before the max is reached ?
  #time_var: "time_sowing" or "time_thermal"
  #only.MSE: should the function just return the MSE or also the parameters of the richards curve
  
  #delta_A_max: After having fitted the Richards curve, we get a value for A which is supposed to 
  #be the asymptotical value of the phenomenon. If the difference between A and the max of the values given is higher than
  #delta_A_max, we consider that the  nls algorithm has not been able to distinguish the asymptot
  #of the phenomenon correctly, so we perform a new fitting by adding artificial points to the originial data
  # --> see the code below for application
  
  #min_points: minimum number of points in order to try nls
  
  #order the data by time_var
  data <- data.frame(data)
  data <- data[order(data[,time_var]),]
  
  
  #Don't try the regression if not enough points
  if(nrow(data[1:which.max(data[,"value"]),]) < min_points)
  {return(list( asymp = NA, mu = NA, lambda = NA, 
                nu = NA, x_infl =  NA, y_infl = NA, MSE = NA, retry.nls = NA))}
  
  
  
  x_min <- min(data[,time_var])
  x_max <- max(data[,time_var])
  x <- seq(x_min,x_max, length=1000)
  loc_reg <- loess(value  ~ get(time_var), data = data)
  
  #Get the approximated values  
  pred <- predict(loc_reg,x)
  df_loess <- data.frame(x = x, y = pred)
  
  #We just keep the first part of the curve (before the lodging of the pea)
  if(!fit.all) {df_loess <- df_loess[1:which.max(df_loess[,2]),]}
  
  fit.control <- grofit.control(model.type = "richards")
  model.growth <- gcFitModel(time = df_loess[,1],data = df_loess[,2], control = fit.control)
  
  #Boolean indicating if we should retry the fitting of the model by adding new points
  retry.nls <- FALSE
  if(!is.atomic(model.growth) & !is.atomic(model.growth$nls)){
    
    coeff <- model.growth$nls %>% coefficients()
    A <- coeff["A"]
    mu <- coeff["mu"]
    lambda <- coeff["lambda"]
    nu <- coeff["addpar"]
    
    #Test coherence between A and the maximum of the smoothed curve
    if(abs(A - max(df_loess[,2])) > delta_A_max){
      retry.nls <- TRUE
    }
  }
  else {
    retry.nls <- TRUE
  }
  
  if(!retry.nls){
    
    #compute coordinates of the inflection point
    x_infl <- (1/mu)*(1+nu)^(-1/nu)*(A+mu*lambda*(1+nu)^(1/nu))
    y_infl <- richards(x_infl, A = A, mu = mu, lambda = lambda, addpar = nu)
    x_infl <- x_infl %>% as.double() ; y_infl <- y_infl %>% as.double()
    
    #compute MSE
    x_fit <- data[1:which.max(data[,"value"]),time_var]
    fit <- richards(x_fit, A = A, mu = mu, lambda = lambda, addpar = nu)
    ref <- data[1:which.max(data[,"value"]),"value"]
    MSE <- sum((ref-fit)^2)/length(fit)
    
    return(list(asymp = as.double(A), mu = as.double(mu), 
                lambda = as.double(lambda), nu = as.double(nu),
                x_infl =  x_infl, y_infl = y_infl, MSE = MSE, retry.nls = retry.nls))
    
  } else {
    
    print("Retrying nls by adding new values")
    
    df_value <- data[1:which.max(data[,"value"]),]
    df_value <- cbind(df_value[,time_var], df_value[,"value"])
    colnames(df_value) <- c(time_var, "value")
    
    
    time_max <- max(df_value[,time_var])
    val_max <- max(df_value[,"value"])
    
    #new values will be added at time_max + lag_time
    lag_time <- ifelse(time_var == "time_sowing", 15, 200)
    
    
    if(data$variable[1]== "height_crop"){
      lag_val <- 0.025
    }else if(data$variable[1]== "biomass_shoot") {
      lag_val <- 0.2
    }
    
    #Add new values
    df_value <-rbind(df_value, c(time_max + lag_time,val_max + lag_val))
    df_value <-rbind(df_value, c(time_max + 2*lag_time,val_max + lag_val + 0.001))
    df_value <- df_value %>% as.data.frame()
    
    
    #Redo the former procedure
    loc_reg <- loess(value  ~ get(time_var), data = df_value)
    
    
    x_min <- min(df_value[,time_var])
    x_max <- max(df_value[,time_var])
    x <- seq(x_min,x_max, length=1000)
    #Get the approximated values  
    pred <- predict(loc_reg,x)
    df_loess <- data.frame(x = x, y = pred)
    
    
    #We just keep the first part of the curve (before the lodging of the pea)
    if(!fit.all) {df_loess <- df_loess[1:which.max(df_loess[,2]),]}
    
    model.growth <- gcFitModel(time = df_loess[,1],data = df_loess[,2],
                               control = fit.control)
    
    if(!is.atomic(model.growth) & !is.atomic(model.growth$nls)){
      
      coeff <- model.growth$nls %>% coefficients()
      A <- coeff["A"]
      mu <- coeff["mu"]
      lambda <- coeff["lambda"]
      nu <- coeff["addpar"]
      
      x_infl <- (1/mu)*(1+nu)^(-1/nu)*(A+mu*lambda*(1+nu)^(1/nu))
      y_infl <- richards(x_infl,  A = A, mu = mu, lambda = lambda, addpar = nu)
      x_infl <- x_infl %>% as.double() ; y_infl <- y_infl %>% as.double()
      
      x_fit <- df_value[1:(nrow(df_value)-2),1]
      fit <- richards(x_fit, A = A, mu = mu, lambda = lambda, addpar = nu)
      ref <- df_value[1:(nrow(df_value)-2),2]
      
      MSE <- sum((ref-fit)^2)/length(fit)
      
      return(list(asymp = as.double(A), mu = as.double(mu), lambda = as.double(lambda),
                  nu = as.double(nu), x_infl =  x_infl, y_infl = y_infl, MSE = MSE,
                  retry.nls = retry.nls))
      
    }
    else
    {
      print("nls still not working")
      return(list( asymp = NA, mu = NA, lambda = NA, 
                   nu = NA, x_infl =  NA, y_infl = NA, MSE = NA, retry.nls = NA))
      
    }
  }
}


# ##Tests##

# #import data
data_traits <- read_rds("../data/phenotype/intercrop_traits.rds")
data_traits <- data_traits$crop
data_height <- data_traits %>%
  filter(variable == "height_crop")

data_test <- data_height%>%
  filter(microplot_id == "2013_TO_pea_geronimo_5_27")

#fit
model <- fit.growth(data_test)

  plot <- data_test %>% 
    ggplot(aes(time_sowing, value)) +
    geom_point(alpha=0.5) + geom_line(stat="smooth", method="loess", alpha=0.5, se=FALSE) +
    labs(x="Days after sowing", y="Crop height (m)")










