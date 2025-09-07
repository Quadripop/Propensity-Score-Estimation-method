load("~/Downloads/STAT 530/STAT 530/results/logitR1000 size500 parATE.Rdata")
load("~/Downloads/STAT 530/STAT 530/results/R1000 size500 parATE.Rdata")
load("~/Downloads/STAT 530/STAT 530/results/rfR1000 size500 parATE.Rdata")
load("~/Downloads/STAT 530/STAT 530/results/treeR1000 size500 parATE.Rdata")
load("~/Downloads/STAT 530/STAT 530/results/bagR1000 size500 parATE.Rdata")

colnames(lgresults)<-colnames(bagresults)<-colnames(rfresults)<-colnames(treeresults)<-colnames(twresults)<-scenarios

source("Project propensity.R")
# Define the scenarios and methods
scenarios <- c("A", "B", "C", "D", "E", "F", "G")
methods <- c( "logit", "tree")

# Store results
results <- data.frame()

# Simulate data for each scenario


  for (method in methods) {
    for(scenario in scenarios){
      for(j in 1:100){
        set.seed(j)
        covfun(200,scenario)
        res <- funsim(sim_data_df, psmethod = method, par = "ATE")
        temp <- data.frame(
          Scenario = scenario,
          Method = method,
          Mean_ASAM = res$masb_aw,
          Rel_Bias = res$absrbias,
          CI_Width = 2 * res$hatgsew # Approximate 95% CI width
        )
        results <- rbind(results, temp)
    }
   
     
    }
    
   
   
  }


# Convert Method and Scenario to factors for ordered plotting
results$Method <- factor(results$Method, levels = methods)
results$Scenario <- factor(results$Scenario, levels = scenarios)

# Define a common theme
custom_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

# Function to create and save separate plots for each scenario
plot_metric <- function(data, metric, xlab, title_prefix) {
  for (scenario in unique(data$Scenario)) {
    p <- ggplot(subset(data, Scenario == scenario), aes_string(x = metric, y = "Method", fill = "Method")) +
      geom_boxplot() +
      labs(title = paste(title_prefix, " - Scenario", scenario), x = xlab, y = "Method") +
      custom_theme
    print(p)  # Print each plot separately
  }
}

# Generate and print the plots
plot_metric(results, "Mean_ASAM", "Mean ASAM", "Mean ASAM across Methods")
plot_metric(results, "Rel_Bias", "Relative Bias (%)", "Relative Bias across Methods")
plot_metric(results, "CI_Width", "CI Width", "95% CI Width across Methods")