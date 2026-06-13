# rm(list=ls())
source("./funs.R")

##### Fix data summary #####
load("logitfix0424.RData")

summarylist = list()
betatrue=rep(0, 1000); betatrue[1:10]=0.5
MC = length(LogitFix); M = dim( LogitFix[[1]][[1]] )[2]

library(stringr)

for( i in 1:7 ){ # i-th method
  methodi = matrix(0, 5, M)
  row.names(methodi)= c( "l_2", "l_sigma", "l_inf", "mcc", "FDRFNR" )
  colnames(methodi) = str_c(1:M, "-batch")
  
  for(j in 1:MC){
    methodi = methodi + apply( LogitFix[[j]][[i]], 2, mysummary, real = betatrue )
    cat(j, "/", MC, "repetition, ", i, "/", 7, "method", "\r")
  }
  summarylist[[i]] = methodi/MC
}

FixLogit = summarylist[-4]


##### Fix data plot #####
library(ggplot2)
library(dplyr)
library(tidyr)

method_labels <- c("AD-IHT", "AD-Lasso", "Renew-Lasso",
                   "Renew-SIM", "RADAR-GLM", "Oracle")

six_colors <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#4D4D4D")
six_linetypes <- c("solid", "42", "82", "4111", "8121",  "11")


plot_df <- lapply(1:length(FixLogit), function(i) {
  mat <- FixLogit[[i]]
  df <- as.data.frame(mat)
  
  df$Metric <- rownames(mat)
  df$Method <- method_labels[i] 
  
  df %>% pivot_longer(cols = contains("batch"), 
                 names_to = "Batch", 
                 values_to = "Value") %>%
    mutate(Batch = as.numeric(gsub("-batch", "", Batch)))
}) %>% bind_rows()


target_metrics <- c("l_2", "l_sigma", "mcc", "FDRFNR")
plot_df_final <- plot_df %>% 
  filter(Metric %in% target_metrics) %>%
  mutate(Metric = factor(Metric, levels = target_metrics,
                         labels = c("L2 error", "L_Sigma error", "MCC", "FDP+FNP")),
         Method = factor(Method, levels = method_labels))

base_font_size <- 22


LogitFixPlot = ggplot(plot_df_final, aes(x = Batch, y = Value, color = Method, linetype = Method, shape = Method, group = Method)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Metric, scales = "free", nrow = 2, axes = "all") +
  
  scale_color_manual(values = six_colors) +
  scale_linetype_manual(values = six_linetypes) +
  scale_x_continuous(breaks = c(1, (1:6)*5 )) +
  
  theme_bw() +
  labs(x = "Batch number",
       y = NULL, 
       color = "Method",
       linetype = "Method",
       shape = "Method") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = base_font_size + 6),
    axis.title.x = element_text(size = base_font_size),
    axis.text.x = element_text(size = base_font_size - 4),
    axis.text.y = element_text(size = base_font_size - 4),
    strip.text = element_text(size = base_font_size - 4),
    legend.text = element_text(size = base_font_size - 6),
    legend.title = element_text(size = base_font_size - 4 ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(5, "line"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(2, "lines"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

ggsave("LogitFixPlot24.pdf", LogitFixPlot, width = 12, height = 8, dpi = 600)




##### Increasing data summary #####
load("logitincrease0418.RData")
summarylistInc = list()
betatrue = rep(0, 1000); betatrue[1:10] = 0.5
MC = length(LogitIncrease); M = dim( LogitIncrease[[1]][[1]] )[2]

library(stringr)

for( i in 1:7 ){ # i-th method
  methodi = matrix(0, 5, M)
  row.names(methodi)=  c( "l_2", "l_sigma", "l_inf", "mcc", "FDRFNR" )
  colnames(methodi) = str_c(1:M, "-batch")
  
  for(j in 1:MC){
    methodi = methodi + apply( LogitIncrease[[j]][[i]], 2, mysummary, real = betatrue )
    cat(j, "/", MC, "repetition, ", i, "/", 7, "method", "\r")
  }
  summarylistInc[[i]] = methodi/MC
}

IncLogit = summarylistInc[-4]


##### Increase data plot #####
library(ggplot2)
library(dplyr)
library(tidyr)

method_labels <- c("AD-IHT", "AD-Lasso", "Renew-Lasso",
                   "Renew-SIM", "RADAR-GLM", "Oracle")
six_colors <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#4D4D4D")
six_linetypes <- c("solid", "42", "82", "4111", "8121",  "11")


plot_df <- lapply(1:length(IncLogit), function(i) {
  mat <- IncLogit[[i]]
  df <- as.data.frame(mat)
  
  df$Metric <- rownames(mat)
  df$Method <- method_labels[i] 
  
  df %>% 
    pivot_longer(cols = contains("batch"), 
                 names_to = "Batch", 
                 values_to = "Value") %>%
    mutate(Batch = as.numeric(gsub("-batch", "", Batch)))
}) %>% bind_rows()



target_metrics <- c("l_2", "l_sigma", "mcc", "FDRFNR")
plot_df_final <- plot_df %>% 
  filter(Metric %in% target_metrics) %>%
  mutate(Metric = factor(Metric, levels = target_metrics,
                         labels = c("L2 error", "L_Sigma error", "MCC", "FDP+FNP")),
         Method = factor(Method, levels = method_labels))




LogitIncPlot = ggplot(plot_df_final, aes(x = Batch, y = Value, color = Method, linetype = Method, shape = Method, group = Method)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Metric, scales = "free", nrow = 2, axes = "all") +
  scale_color_manual(values = six_colors) +
  scale_linetype_manual(values = six_linetypes) +
  scale_x_continuous(breaks = 1:12 ) +
  theme_bw() +
  labs(x = "Batch number",
       y = NULL, 
       color = "Method",
       linetype = "Method",
       shape = "Method") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = base_font_size + 6),
    axis.title.x = element_text(size = base_font_size),
    axis.text.x = element_text(size = base_font_size - 4),
    axis.text.y = element_text(size = base_font_size - 4),
    strip.text = element_text(size = base_font_size - 4),
    legend.text = element_text(size = base_font_size - 6),
    legend.title = element_text(size = base_font_size - 4 ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(5, "line"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(2, "lines"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

ggsave("LogitIncPlot.pdf", LogitIncPlot, width = 12, height = 8, dpi = 600)


