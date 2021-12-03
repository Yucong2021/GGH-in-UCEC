

library(tidyverse)
library(ggplot2)
library(survival)
library(DESeq2)
library(clusterProfiler)
library(flextable)

gene_id <- "ENSG00000107864"
gene_name <- "ENSG00000237489"
my_color <- c("#5abdc7","#f99821","#bcc74d","#f8e54f")

load("./data/expr.Rdata") 
load("./data/clinical.Rdata") 

###
rownames(dat_tpm) <- substr(rownames(dat_tpm), 1, 12)
dat_clincal <- dat_tpm %>%  
  dplyr::select(UQ(gene)) %>%
  rownames_to_column("Samples") %>% 
  left_join(dat_clincal, ., by = "Samples") %>% 
  dplyr::rename("Expression"=UQ(gene)) %>% 
  dplyr::filter(!is.na(Expression))   
dat_clincal$Group <- ifelse(dat_clincal$Expression >= median(dat_clincal$Expression), "High", "Low")
dat_clincal$Group <- factor(dat_clincal$Group, levels = c("Low","High"))

dat_clincal[,gene_name] <- ifelse(dat_clincal$Expression >= median(dat_clincal$Expression), "High", "Low")
dat_clincal[,gene_name] <- factor(dat_clincal[,gene_name], levels = c("Low","High"))
dat_clincal$OS.event <- factor(dat_clincal$OS.event, levels = c(0,1), labels = c("Alive","Dead"))


library(tableone)
tmp <- CreateTableOne(vars =  colnames(dat_clincal)[c(14:18,20:26)],
                      strata = gene_name,
                      data = dat_clincal)
output <- print(tmp,contDigits = 1, showAllLevels = TRUE) %>% as.data.frame()



### 
table.logi <- function(dat_colnames, data, group){
  var1 <- NULL
  
  for(i in dat_colnames) {
    num <- sum(!is.na(data[,i])) 
    fml <- as.formula(paste0(i," ~ ", group))
    mod1 <- glm(data = data, formula = fml, family = 'binomial')
    pvalue <- ifelse(summary(mod1)$coefficients[2,"Pr(>|z|)"] < 0.001, 
                     "<0.001", 
                     sprintf("%.3f",summary(mod1)$coefficients[2,"Pr(>|z|)"]))
    CI <- paste0(sprintf("%.2f",exp(confint(mod1))[2,]),collapse = '-' )
    OR <- paste0(sprintf("%.2f",exp(coef(mod1))[2]), '(', CI, ')')
    name <- gsub("."," ", i,fixed = T)
    char <- paste0(name,' (',paste0(levels(data[,i]), collapse = ' vs. '),')')
    var2 <- data.frame('Characteristics' = char, 'Total(N)' = num, 'OR' = OR, 'pvalue' = pvalue)
    var1 <- rbind(var1, var2)
  }
  var1$Characteristics <- gsub("_"," ",var1$Characteristics)
  var1$Characteristics <- Hmisc::capitalize(var1$Characteristics) ## 首字母大写
  tmp <- paste0("Odds ratio in ", gene_name, " expression")
  colnames(var1)[2:4] <- c("Total (N)", tmp, "P value") 
  
  return(var1)
}

output <- table.logi(data = dat_clincal, group = 'Group',
                     dat_colnames = colnames(dat_clincal)[c(14:18,20:26)])




Cox.subgroup <- function(colname, subgroup, time = "OS.time", event = "OS.event", data){  
  data <- data[!is.na(data[,time]),]
  data <- data[!is.na(data[,event]),]
  tmp_filter <- data[data[, colname] %in% subgroup,]
  mysur <- Surv(time = tmp_filter[,time]/30, event = tmp_filter[,event])
  fml <- as.formula(paste0("mysur ~", "group"))
  fit2 <- coxph(fml, data = tmp_filter)
  sum_fit2 <- summary(fit2)
  percent <- round(sum_fit2$n / sum(!is.na(data[, colname]))*100,0)
  number <- paste0(sum_fit2$n, " (", percent, ")")
  # HR <- round(sum_fit2$coefficients[,2],3)
  HR <- sprintf("%0.3f", sum_fit2$coefficients[,2])
  CI <- paste0(sprintf("%0.3f", sum_fit2$conf.int[,3:4]),
               collapse = '-')
  HR <- paste0(HR, " (", CI, ")")  
  pvalue <- ifelse(sum_fit2$coefficients[,5] < 0.001,  "<0.001",
                   sprintf("%0.3f",sum_fit2$coefficients[,5]))
  
  data.frame("Name" = paste0(colname, ".", subgroup, collapse = "+"),
             'Characteristics' = paste0(subgroup, collapse = "+"),
             "N (%)" = number,
             'Hazard ratio' = HR,
             'P value' = pvalue, check.names = F)
}



Cox.singlegroup <- function(colname, subgroup="all",time="OS.time", event="OS.event",data){
  
  if(subgroup[1] == "all"){
    subgroup = levels(data[,colname])
    tmp <- data.frame("Name" = colname,
                      'Characteristics' = gsub("."," ",colname,fixed = T),
                      "N (%)" = NA,
                      'Hazard ratio' = NA,
                      'P value' = NA, check.names = F)
    for(i in subgroup){
      tmp1 <- Cox.subgroup(colname=colname, subgroup = i, time = time,event = event, data=data)
      tmp <- rbind(tmp, tmp1)
    }
  } else {
    tmp <- Cox.subgroup(colname = colname, subgroup = subgroup, time = time, event = event, data=data)
  }
  return(tmp)
}



Cox.singlegroup(colname = colnames(dat_clincal)[14],subgroup = levels(dat_clincal[,13]),
                data = dat_clincal)


Cox.somegroup <- function(colVar, time="OS.time", event="OS.event", data){
  tmp <- list()
  for(i in 1:length(colVar)){
    tmp[[i]] <- Cox.singlegroup(colname = colVar[i], subgroup = "all", time = time, event = event, data = data)
  }
  return(tmp)
}

output <- Cox.somegroup(colVar = colnames(dat_clincal)[c(14:18,20:26)], 
                        time  = "OS.time",
                        event = "OS.event",
                        data  = dat_clincal)



Cox.Univariate <- function(colname ,time="OS.time", event="OS.event", data){
  dd <- NULL
  for(i in colname){
    mysur <- Surv(time = data[,time]/30, event = data[,event])
    fit2 <- coxph(mysur ~ data[,i])
    sum_fit2 <- summary(fit2)
    HR <- round(sum_fit2$coefficients[,"exp(coef)"],3)
    CI <- paste0(round(sum_fit2$conf.int[,3:4],3), collapse = '-')
    HR <- paste0(HR, "(", CI, ")")
    pvalue <- ifelse(sum_fit2$coefficients[,"Pr(>|z|)"] < 0.001, 
                     "<0.001",
                     round(sum_fit2$coefficients[,"Pr(>|z|)"],3))
    
    name <- gsub("."," ",i,fixed = T)
    char <- paste0(name,' (',paste0(levels(data[,i]), collapse = ' vs. '),')')
    
    d <- data.frame("Column"=i,
                    'Characteristics' = char,
                    'Univariate analysis.Hazard ratio(95% CI)' = HR,
                    'Univariate analysis.P value' = pvalue,
                    check.names = F)
    dd <- rbind(dd, d)
  }
  return(dd)
}

output <- Cox.Univariate(colname = colnames(dat_clincal)[c(14:18,20:26)], 
                         time="OS.time",
                         event="OS.event",
                         data = dat_clincal)

Cox.Multivariate <- function(uniResult, time="OS.time", event="OS.event", data, threshold=0.1, variable.show = F){
  
  variables <- uniResult$Column[uniResult$`Univariate analysis.P value` < 0.1]
  mysur <- Surv(time = data[,time]/30, event = data[,event])
  fml <- as.formula(paste0('mysur~',paste0(variables, collapse = '+')))
  multicox <- coxph(fml,data = data)
  multisum <- summary(multicox)
  cindex <- multisum$concordance
  c.ci <- sprintf("%.3f",c(cindex[1] - cindex[2], cindex[1] + cindex[2]))
  cindex <- sprintf("%.3f",cindex[1])
  cindex <- paste0("C-index:",cindex,"(",c.ci[1], "-", c.ci[2],")")
  print(cindex)
  
  MHR <- round(multisum$coefficients[,"exp(coef)"],3) 
  MCIL <- round(multisum$conf.int[,3],3)
  MCIR <- round(multisum$conf.int[,4],3)
  MCI <- paste0(MCIL, '-', MCIR)
  MHR <- paste0(MHR, "(", MCI ,")")
  Mpvalue <-ifelse(multisum$coefficients[,"Pr(>|z|)"] < 0.001, 
                   "<0.001",
                   round(multisum$coefficients[,"Pr(>|z|)"],3)) 
  dd <- data.frame('Column' = variables,
                   'Multivariate analysis.Hazard ratio(95% CI)' = MHR,
                   'Multivariate analysis.P value' = Mpvalue,
                   check.names = F)
  
  dd <- dplyr::left_join(uniResult, dd, by = "Column")
  if(!variable.show){dd <- dd[,-1]}
  
  return(dd)
}

output <- Cox.Multivariate(uniResult = output.coxU,
                           time="OS.time",
                           event="OS.event",
                           data = dat_clincal)




f1.1 <- function(data, colour){
  
  data$expr <- log2(data$expr + 1)
  data$type2 <- factor(data$type2)
  
  ## 手动分析差异
  for(i in 1:31){
    tmp <- data[data$tissue_type %in% cancerlist[c(2*i-1,2*i)],]
    data$pvalue[data$tissue_type %in% cancerlist[c(2*i-1,2*i)]] <- 
      wilcox.test(expr~tissue_type, data = tmp)$p.value
  }
  data$pvalue <- ifelse(data$pvalue > 0.05, "ns",
                        ifelse(data$pvalue > 0.01, "*",
                               ifelse(data$pvalue > 0.001, "**","***")))
  data$y <- max(data$expr)
  
  
  
  p1.1 <- ggplot(data, aes(x = tissue, y = expr)) +
    geom_point(mapping = aes(colour = type2),
               position = position_jitterdodge(jitter.width = .5),size = .6,alpha = .2) +
    geom_boxplot(aes(fill = type2),outlier.shape= NA, 
                 notch = F, alpha = .6, size = 0.35,width = .8) +
    
    labs(x = "",y = paste0("The expression of ",gene_name,"\nlog2(TPM + 1)")) +
    
    theme(axis.text.x = element_text(angle = 45,hjust=1)) +  ## X轴的text旋转45度
    theme(axis.title.y = element_text(size = 7),
          axis.text.y = element_text(size = 6)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7)) +
    theme(axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35)) +
    
    theme(panel.grid = element_blank(),
          panel.background = element_blank()) +
    
    theme(legend.key.size = unit(0.4,'cm'),
          legend.position = "top",
          legend.margin = unit(c(0,0,0,0), "cm"), ## 白底的边距
          legend.background = element_blank(), ## 去掉白底
          legend.box.margin = unit(c(-0.2,0,-0.3,0), "cm"),  ## 两边的空白
          legend.text = element_text(size = 7),
          legend.title = element_blank()) +
    
    scale_fill_manual(values = colour,labels=c("TCGA&GTEx_Normal","TCGA_Tumor")) +
    scale_color_manual(values = colour) +
    guides(colour = F) +
    geom_text(aes(label = pvalue,y = y),check_overlap = T,size = 2.3)
  
  
  return(p1.1)
}

p <- f1.1(dat_clincal, colour = my_color[1:2])


f1.4 <- function(data, colour="#f8766d"){
  library(pROC)
  
  roc1 <- roc(response = data$type2, predictor = data$expr)
  
  auc <- auc(roc1) %>% sprintf("%.3f",.)
  auc <- paste0("AUC: ",auc)
  ci <- ci.auc(roc1) %>% sprintf("%.3f",.)
  ci <- paste0("CI: ",min(ci),"-",max(ci))
  anno <- paste0(auc, "/n",ci)
  
  
  tmp <-  data.frame(x = 1-roc1$specificities,
                     y = roc1$sensitivities)
  tmp <- tmp[order(tmp$x, tmp$y),]
  
  library(plotROC)
  
  p1.4 <- ggplot() +
    # geom_roc(data = data, aes(d = group, m = expr), colour = colour,
    #          n.cuts = 0,show.legend = F, size = 0.465) +
    geom_area(data = tmp,aes(x = x, y = y),position = "identity",
              fill = colour,alpha = .2) +
    geom_line(data = tmp, aes(x = x, y = y), colour = colour,size =0.35) +
    geom_line(aes(x = c(0,1), y = c(0,1)), colour = "grey", size = 0.35)+
    annotate("text", x=0.55, y=0.19, size = 2.5, hjust = 0,
             label = gene_name) +
    annotate("text", x=0.55, y=0.13, size = 2.5, hjust = 0,
             label = auc) +
    annotate("text", x=0.55, y=0.07, size = 2.5, hjust = 0,
             label = ci) +
    scale_x_continuous(expand = c(0, 0),breaks = 0:5*0.2, limits = c(-0.01,1.01)) +
    scale_y_continuous(expand = c(0, 0),breaks = 0:5*0.2, limits = c(-0.01,1.01)) +
    labs(y = "Sensitivity (TPR)", x = "1-Specificity (FPR)") +
    theme_bw() +
    # theme(panel.grid = element_blank()) +
    theme(axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 7)) +
    theme(axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 7)) +
    theme(axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35))
  return(p1.4)
}

p <- f1.4(dat_clincal,colour = my_color[2])



boxplot.1Clinical <- function(colname, subgroup = "all", data, level = NULL, colour){
  
  
  data <- data[!is.na(data$Expression),]
  data <- data[!is.na(data[,colname]),]
  
  if(subgroup[1] == "all"){
    subgroup = levels(data[,colname])}
  data <- data[data[,colname] %in% subgroup, ] 
  data[,colname] <- as.character(data[,colname])
  data[,colname] <- factor(data[,colname], levels = subgroup)
  if(length(setdiff(level, levels(data[,colname]))) != 0){
    combineLevel <- setdiff(levels(data[,colname]),level)
    combineName <- setdiff(level, levels(data[,colname]))
    
    data[,colname] <- as.character(data[,colname])
    data[,colname][data[,colname] %in% combineLevel] <- combineName
    data[,colname] <- factor(data[,colname], levels = level)
  }
  colour = colour[1:length(levels(data[,colname]))]
  
  ## 加个名字
  num <- table(data[,colname]) %>% as.data.frame()
  x_name <- paste0(num$Var1, "\n(n=",num$Freq, ")")
  
  data$Expression <- log2(data$Expression + 1)
  
  p <- ggplot() +  
    geom_violin(aes_string(x = colname, y = "Expression", fill = colname),data,
                alpha = .5, size = 0.35) +
    geom_boxplot(aes_string(x = colname, y = "Expression", fill = colname),data,
                 outlier.shape= NA, notch = F, alpha = .8, size = 0.35) +
    # geom_point(aes_string(x = colname, y = "Expression", colour = colname),data,
    #            position=position_jitter(width=.2,height=0), size = 1.5, shape= 21, alpha = .8,fill = "white") +
    stat_summary(aes_string(x = colname, y = "Expression"),data,
                 fun.y=mean,geom="point",fill='white',shape=21,size=2) +
    scale_y_continuous(limits = c(0, max(data$Expression)+0.5),
                       breaks = 0:5
    ) +
    labs(x = gsub("."," ",colname,fixed = T),
         y = paste0("The expression of ",gene_name,"\n log2(TPM+1)")) +
    guides(fill = F, colour = F) +
    theme_classic() +
    theme(panel.background = element_blank()) +
    theme(axis.title.y = element_text(size = 7),
          axis.text.y = element_text(size = 7)) +
    theme(axis.text.x = element_text(size = 7),
          axis.title.x = element_text(size = 7)) +
    theme(axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35)) +
    scale_fill_manual(values = colour) +
    scale_color_manual(aesthetics = "colour",values = colour)
  
  
  if(length(levels(data[,colname])) > 2){ 
    pvalue <- kruskal.test(data$Expression~data[,colname])$p.value 
  } else {
    pvalue <- wilcox.test(data$Expression~data[,colname])$p.value 
  } 
  if(pvalue > 0.001){
    anno <- paste0("p = ",round(pvalue, 3))
  } else {
    anno <- "p < 0.001"
  }
  
  p + geom_line(aes(x = c(1,length(levels(data[,colname]))),
                    y = rep(max(data$Expression)+0.3, 2)),size = 0.35)  +
    geom_line(aes(x = c(1,1),
                  y = c(max(data$Expression)+0.2,max(data$Expression)+0.4)),size = 0.35)  +
    geom_line(aes(x = c(rep(length(levels(data[,colname])),2)),
                  y = c(max(data$Expression)+0.2,max(data$Expression)+0.4)),size = 0.35)  +
    geom_text(x = length(levels(data[,colname]))/2 + 0.5,
              y = max(data$Expression)+0.5,
              aes(label = anno), size = 2.3)
  
}

boxplot.1Clinical(colname = colnames(dat_clincal)[13], data = dat_clincal,colour = my_color)


plot.KM <- function(dat=dat_clincal, time="OS.time", event="OS.event", colour, 
                    best_seperation = F, cut=0.25){
  dat <- dat[!is.na(dat[,time]),]
  dat <- dat[!is.na(dat[,event]),]
  
  if(time %in% "OS.time"){type = "Overall Survival"}
  if(time %in% "PFI.time"){type = "Progression-free Interval"}
  if(time %in% "DSS.time"){type = "Disease-specific Survival "}
  
  
  if(best_seperation){
    res.cut <- survminer::surv_cutpoint(dat, 
                                        time = time,
                                        event = event,
                                        variables = "Expression",
                                        minprop = cut) 
    dat$Group <- ifelse(dat$Expression > res.cut$cutpoint[1,1], "High","Low")
    dat$Group <- factor(dat$Group, levels = c("Low","High"))}
  
  mysur <-  Surv(time = dat[,time]/30, event = dat[,event])
  fml <- as.formula(paste0('mysur~', 'Group')) 
  fit <- survfit(fml, data = dat)
  # return(fit)
  
  Low <- paste0("Low ",gene_name, " (", fit$n[1], ")")
  High <- paste0("High ",gene_name, " (", fit$n[2], ")")
  
  legend.lab = c(Low, High)
  
  p1 <- survminer::ggsurvplot(fit = fit, data = dat, 
                              size = 0.35, censor.size = 2,
                              # surv.scale = "percent",
                              # # conf.int = F, ## 置信区间
                              # risk.table = F,
                              # # ncensor.plot = F,
                              legend.title = gene_name,
                              legend.labs = legend.lab,
                              xlab = "Time (months)",
                              # surv.median.line = "hv",
                              # # tables.theme = theme_cleantable(),
                              palette  = colour ##这个调整颜色的参数需要注意
  )
  
  fit2 <- coxph(fml, data = dat)
  sum_fit2 <- summary(fit2)
  HR <- sprintf("%0.2f",sum_fit2$coefficients[,"exp(coef)"])
  CI <- paste0(sprintf("%0.2f",sum_fit2$conf.int[,3:4]), collapse = '-')
  HR <- paste0("HR=",HR, "(", CI, ")")
  pvalue <- ifelse(sum_fit2$coefficients[,"Pr(>|z|)"] < 0.001, 
                   "p<0.001",
                   paste0("p=",sprintf("%0.3f",sum_fit2$coefficients[,"Pr(>|z|)"])))
  
  p1$plot +
    # geom_text(aes(label = paste0(type,"\n",HR, "\n", pvalue),), 
    #           x = max(dat[,time],na.rm = T)/30-80, y = .85, check_overlap = T, size = 2.5) +
    theme(axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 7)) +
    theme(axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 7)) +
    theme(axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35)) +
    
    theme(legend.key.size = unit(0.4,'cm'),
          legend.key.height = unit(0.12,'cm'), 
          legend.position = c(0.95,0.85), 
          legend.justification = c(1,.5),  
          legend.margin = unit(c(0,0,0,0), "cm"),
          legend.background = element_blank(), 
          legend.box.margin = unit(c(-0.2,0,-0.3,0), "cm"), 
          legend.text = element_text(size = 7),
          legend.title = element_blank()) +
    geom_text(aes(label = type), hjust = 1,
              x = (max(dat[,time],na.rm = T)/30)*1.1, y = .95, check_overlap = T, size = 2.5) +
    geom_text(aes(label = HR), hjust = 1,
              x = (max(dat[,time],na.rm = T)/30)*1.1, y = .74, check_overlap = T, size = 2.5) +
    geom_text(aes(label = pvalue), hjust = 1,
              x = (max(dat[,time],na.rm = T)/30)*1.1, y = .67, check_overlap = T, size = 2.5)
}

plot.KM(dat = dat_clincal,time = "OS.time",event = "OS.event",colour = my_color[1:2])
