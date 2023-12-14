# load required library
library(survival)
library(flexsurv)
library(survminer)
library(tidyverse)
library(gt)

# build AFT Lognormal model based on file input
sLno<- flexsurvreg(Surv(heart$time,heart$DEATH_EVENT) ~ creatinine_phosphokinase+age+serum_creatinine+serum_sodium,dist='lnorm',data=heart)

# build test_df function that create a dataframe object using the input value as the data source.
test_df <- function(input1, input2, input3, input4) {
creatinine_phosphokinase <- input1
 age <- input2
serum_creatinine <- input3
 serum_sodium <- input4
test <-data.frame(creatinine_phosphokinase, age, serum_creatinine, serum_sodium)
return(test)
}

# Function to predict survival probabilities for new data and compare to reference population
# Arguments:
# - test: data frame of predictor variables for which to make predictions
# - time: time point at which to make predictions
pred <- function(test, time) {
  # Use pre-trained survival model to predict survival probabilities for test data
  new <- predict(sLno, newdata = test, type = 'survival', time = time)
  pred_test <- new$.pred_survival
  
  # Use pre-trained survival model to predict survival probabilities for reference population
  pred_obs <- predict(sLno, newdata = heart, type = 'survival', time = time)
  
  # Calculate mean predicted survival probability for reference population
  mean <- mean(pred_obs$.pred_survival)
  
  # Assign label of 1 or 0 to each observation in reference population based on survival probability
  pred_obs$label <- ifelse(pred_obs$.pred_survival < pred_test, 1, 0)
  
  # Create column indicating whether observation received label of 1
  pred_obs$cond <- pred_obs$label == 1
  
  # Calculate proportion of reference population that received label of 1
  prop <- mean(pred_obs$label)
  
  # Return list of results, including predicted survival probabilities for test data, reference population with labels, 
  # mean predicted survival probability for reference population, and proportion of reference population that received label of 1
  result <- list(pred_test, pred_obs, mean, prop)
  return(result)
}



# Function to predict individual survival at 1 month
# Inputs:
#   - test: data frame containing the predictor variables for the patient to be predicted
# Outputs:
#   - pred_test: predicted probability of survival at 1 month for the patient
pred1m<- function(test, time){
  # Use the trained survival model 'sLno' to predict the survival probability of the patient
  # at the time of 1 month (time=30)
  new <- predict(sLno,newdata=test, type='survival',time=time)
  
  # Extract the predicted survival probability for the patient
  pred_test <- new$.pred_survival
  
  # Return the predicted survival probability for the patient
  return(pred_test)
}


# Function to create a histogram of survival probabilities for a specific time point
# Inputs:
#   - test: data frame containing the predictor variables for the patient of interest
#   - time: time point (in days) to calculate survival probabilities at
# Returns:
#   - plotly object of the histogram
MakeHist<-function(test,time){
  # Calculate observed survival probabilities for the population at the specified time point
  obs <- predict(sLno,newdata=heart, type='survival',time=time)
  # Calculate predicted survival probabilities for the patient of interest at the specified time point
  new <- predict(sLno,newdata=test, type='survival',time=time)
  # Create the histogram plot using ggplot2
  histplot<-ggplot(data=obs,aes(x=.pred_survival,fill=obs$cond)) +
    geom_histogram(aes(y = after_stat(count / sum(count))),binwidth =0.05)+
    scale_x_continuous(name = "Month Survivals") +
    scale_fill_manual(values = c("orange", "purple"))+
    # Add vertical lines indicating the predicted survival probability for the patient of interest and the population mean
    geom_vline(aes(xintercept = new$.pred_survival,hoverinfo = "x+text",text=paste("Patient",round(new$.pred_survival,2))), color = "red", linetype = "dashed")+ 
    geom_vline(aes(xintercept = mean(obs$.pred_survival),hoverinfo = "x+text",text=paste("Population Mean",round(mean(obs$.pred_survival),2))), color = "blue", linetype = "dashed") + 
    # Set plot theme and axis labels
    theme_bw() + theme(axis.line = element_line(size=1, colour = "black"),
                       panel.grid.major = element_line(colour = "#d3d3d3"),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(), panel.background = element_blank(),
                       plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
                       text=element_text(family="Tahoma"),
                       axis.text.x = element_text(colour = "black", size = 9),
                       axis.text.y = element_text(colour = "black", size = 9),
                       legend.text = element_text(size = 10),
                       legend.position = "bottom") +
    ggtitle("Percentile of Individual Survival in Selected Month") +
    ylab("Proportion of Survivals") +
    scale_fill_manual(values = c("orange", "purple"), labels = c("Proportion of survivals >= patient", "Proportion of survivals < patient")) +
    labs(fill = "Proportion")
  # Convert plot to plotly object and return it
  plotly_plot <- plotly::ggplotly(histplot, tooltip = c("x", "text"))
  return(plotly_plot)
}

# This function creates a survival plot with two survival curves: one for the overall population and one for a patient/test dataset.
# The input 'test' is a survival dataset for a patient/test.
# The function returns a survival plot.
Makesurvplot <- function(test){
  # Create a new plot
  plot.new()
  # Set the x-axis and y-axis limits
  plot.window(xlim = c(0,365), ylim = c(0,1))
  # Add a light gray background
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "#ebebeb")
  # Add white grid lines
  grid(nx = NULL, ny = NULL,
       col = "white", lwd = 2)
  # Add x-axis and y-axis labels, adjust their size
  axis(1, cex.axis = 1.5); axis(2, cex.axis = 1.5); box(); title(ylab="Survival", xlab="Time (days)", cex.lab = 1.5)
  # Add a blue survival curve for the overall population
  lines(sLno, type='survival', col="blue")
  # Add a red survival curve for the patient/test dataset
  lines(sLno,newdata=test, type='survival',col = "red")
  # Add a legend to distinguish the two curves
  legend("topright", legend = c("Patient","Overall"),
         lty = 1, col = c('red',"blue"))
}

# This function creates a table that summarizes the survival statistics of a patient/test dataset compared to the overall population.
# The input 'pred' is a list that contains the following items:
#   - pred[[1]]: the survival rate for the patient/test dataset
#   - pred[[3]]: the mean survival rate for the overall population
#   - pred[[4]]: the proportion of patients with a survival rate lower than the patient/test dataset
# The function returns a data frame that contains the statistics in a readable format.
Makesurvtable <- function(pred){
  # Extract the survival statistics from the input list
  meansurv <- pred[[3]]
  patsurv <- pred[[1]]
  prop <- pred[[4]]
  # Create a data frame that contains the survival statistics
  data <- data.frame(
    Statistic = c('Mean survival of population', 'Patient survival','Proprotion of survivals (< patient)'),
    Value = c(meansurv, patsurv, prop)
  )
  # Add formatting to the data frame using the 'gt' package (commented out)
  # data <-  gt(data) %>% 
  #   cols_label(
  #     Statistic = md("<b>Statistic</b>"),
  #     Value = md("<b>Value</b>"))
  # Return the data frame
  return(data)
}




