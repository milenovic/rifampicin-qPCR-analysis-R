library(ggplot2)
library(cowplot)
library(tidyr)
library(plyr)
library(multcompView)
library(MASS)
library(readr)
library(naniar)
library(gridExtra)

## locate and read data
raw_data_sample <- read_csv("Data_samples.csv")
#View(raw_data_sample)
raw_data_curves <- read_csv("Data_curves.csv")
#View(raw_data_curves)

#Clean up the data
raw_data_sample$Group <- gsub('-','_', raw_data_sample$Group)
raw_data_sample <- subset(raw_data_sample, Group != "NC" & Group != "PC_D1" & Group != "PC_D2" & Group != "PC_D3")

raw_data_sample <- raw_data_sample %>% replace_with_na(replace = list(Actin=40.00000))


#Function for plotting standard curves
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    ggtitle(paste("Standard curve for", names(fit$model)[1]),
            subtitle = substitute(paste("Adj. R"^2," = ",r,
                                  ", Intercept = ",i,
                                  ", Slope = ",s,
                                  ", p = ",p),
                                  list(r= signif(summary(fit)$adj.r.squared, 5),i= signif(fit$coef[[1]],5 ),s=signif(fit$coef[[2]], 5),p=signif(summary(fit)$coef[2,4], 5))
                                 )
           ) + xlab(paste(names(fit$model)[2], "(log10)")) + ylab(paste(names(fit$model)[1], "(Ct)"))
}


#Make standard curves, save models to "curves" and parameters to "curves_param", and plot the curves
curves_param <- data.frame()
curves <- list()
raw_data_curves_trans <- raw_data_curves
raw_data_curves_trans[, "Quantity"] = log10(raw_data_curves_trans[, "Quantity"])
for(i in names(raw_data_curves[, -c(2, 3)])[-1]){
  curves[[i]] <- lm(paste(i, " ~ Quantity"), raw_data_curves_trans[, -c(2, 3)])
  tmp <- data.frame(lapply(curves[[i]][["coefficients"]][1:2], type.convert), stringsAsFactors=FALSE)
  tmp[, "rsquared"] = signif(summary(curves[[i]])$adj.r.squared, 5)
  tmp[, "p"] = signif(summary(curves[[i]])$coef[2,4], 5)
  names(tmp)[1] <- "intercept"
  names(tmp)[2] <- "slope"
  curves_param <-  rbind(curves_param, cbind(Gene = i, tmp))
  rm(tmp)
  print(ggplotRegression(curves[[i]]))
  }
rm(raw_data_curves_trans)

pcr_amount <- function(vec, a, b) {
  res <- 10 ^ ((vec - a)/b)
  return(res)
}

get_amounts <- function(df, intercept, slope) {
  
  amounts <- mapply(function(d, a, b) pcr_amount(d, a, b),
                    d = df, a = intercept, b = slope)
  amounts <- as.data.frame(amounts)
  
  return(amounts)
}


relative_quantities_perSample <- get_amounts(
  raw_data_sample[3:7],
  intercept = curves_param$intercept,
  slope = curves_param$slope
)

not_normalized_perSample <- relative_quantities_perSample
not_normalized_perSample <- cbind(not_normalized_perSample, raw_data_sample[c("ID", "Group")])

powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1
  )
}

# run the box-cox transformation to determine optimal lambda
bc_nn <- boxcox(not_normalized_perSample$Actin ~ not_normalized_perSample$Group)

#save and print lambda
(lambda_nn <- bc_nn$x[which.max(bc_nn$y)])

#BoxCox transform the data
not_normalized_perSample <- mutate(not_normalized_perSample, bcActin = powerTransform(not_normalized_perSample$Actin, lambda_nn))

#ANOVA model using NOT normalized, transformed data
ANOVA_model_nn <- aov(not_normalized_perSample$bcActin~Group, data=not_normalized_perSample)
par(mfrow = c(2, 2))
plot(ANOVA_model_nn)
summary(ANOVA_model_nn)
ANOVA_model_nn

#Post-hoc Tukey test
tHSD_nn <- TukeyHSD(ANOVA_model_nn, ordered = FALSE, conf.level = 0.95)
tHSD_nn

#Function to generate letter groupings for graphing
generate_label_df <- function(d, HSD, flev, colNo){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  names(d)[colNo] <- "transformed"
  boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$transformed)) + 0.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}

not_normalized_perSample <- mutate(not_normalized_perSample, OsActin = not_normalized_perSample$Actin)

library(tidyverse)
not_normalized_perSample <- not_normalized_perSample %>% 
  mutate(sex = ifelse(grepl("M$", Group), "Male", "Female"))

         
#Box plot NOT normalized data with Tukey letters, on transformed scale
p3<-ggplot(not_normalized_perSample, aes(x=Group, y=bcActin)) + 
  geom_boxplot(show.legend = FALSE , aes(fill=sex)) + scale_fill_manual(values=c("#bdbdbd", "#636363")) + 
  geom_text(data = generate_label_df(not_normalized_perSample, tHSD_nn, 'Group', 8), aes(x = plot.labels, y = V1, label = labels)) +
  ggtitle("Actin", subtitle = paste("Not normalized, transformed scale, Box Cox lambda = ", round(lambda_nn, digits=3))) +
  xlab("Treatments") + ylab("Relative quantity (transformed)") +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(face="bold")) +
  scale_x_discrete(labels = str_wrap(c("Antibiotic female", "Antibiotic male", "Control female", "Control male", "Colony female", "Colony male"), width = 10))
p3
ggsave(filename = "Actin-NN-Ts.tiff", path="Paper/Graphs", width = 185, height = 185, units="mm", device='tiff', dpi=300, compression = "lzw")



#Box plot NOT normalized data with Tukey letters, on original scale
p4<-ggplot(not_normalized_perSample, aes(x=Group, y=OsActin)) + 
  geom_boxplot(show.legend = FALSE , aes(fill=sex)) + scale_fill_manual(values=c("#bdbdbd", "#636363")) + 
  geom_text(data = generate_label_df(not_normalized_perSample, tHSD_nn, 'Group', 9), aes(x = plot.labels, y = V1, label = labels)) +
  ggtitle("Actin", subtitle = paste("Not normalized, original scale, Box Cox lambda = ", round(lambda_nn, digits=3))) +
  xlab("Treatments") + ylab("Relative quantity") +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(face="bold")) +
  scale_x_discrete(labels = str_wrap(c("Antibiotic female", "Antibiotic male", "Control female", "Control male", "Colony female", "Colony male"), width = 10))
p4
ggsave(filename = "Actin-NN-Os.tiff", path="Paper/Graphs", width = 185, height = 185, units="mm", device='tiff', dpi=300, compression = "lzw")


#write.csv2(normalized_perSample,'Normalized per sampleP.csv')