#' Summary of clinical performances
#'
#' Evaluate biomarker individual performances using best classification metrics
#' 
#' @param data provide a data.frame
#' @param group group by a categorical variable
#' @param direction Define the positive case
#' 
#' @import dplyr
#' @import pROC
#' @import tidyr
#' @importFrom reshape2 melt
#' @importFrom nnet multinom
#' @importFrom MLmetrics MultiLogLoss
#' 
#' @details 
#' `Logistic regression` was used for binary class problem. \cr
#' `Multinomial logistic regression` was used for multi class problem.
#' 
#' 
#' @export
geneseng_summary_class_metrics <- function(data, group, direction = c("auto", ">", "<")){
  
  direction <- match.arg(direction)
  
  if(length(unique(data[,group])) == 2){
    
    geneseng_summary_binary_class_metrics(data = data, group = group, direction =  direction)
    
  } else {
    
    geneseng_summary_multi_class_metrics(data = data, group = group, direction =  direction)
    
  }
  
}


#' @rdname geneseng_summary_class_metrics
#' @export
geneseng_summary_binary_class_metrics <- function(data, group, direction = c("auto", ">", "<")){
  
  directon <- match.arg(direction)
  
  # Compute logLoss
  logLoss <- function(pred, actual){
    -mean(actual * log(pred) + (1 - actual) * log(1 - pred))
  }
  
  df <- data %>%
    pivot_longer(cols = !group, names_to = "biomarker") %>%
    rename(group2 = names(.)[1]) %>%
    as.data.frame() %>%
    mutate_at(vars(group2), list(factor))
  
  
  results <- lapply(1:length(unique(df$biomarker)), function(x){
    
    tmp <- df[df$biomarker == unique(df$biomarker)[x], ]
    model <- glm(group2 ~ value, data = tmp, family = "binomial")
    rocs <- roc(predictor = tmp$value, response = tmp$group2)
    metrics <- coords(roc = rocs, x = "best", ret = "all", transpose = FALSE)
    
    cbind.data.frame(
      biomarker = unique(df$biomarker)[x],
      class1 = ifelse(rocs[["direction"]] == ">", yes = levels(df$group2)[1], no = levels(df$group2)[2]),
      class2 = ifelse(rocs[["direction"]] == ">", yes = levels(df$group2)[2], no = levels(df$group2)[1]),
      model = "logistic reg.",
      logLoss = logLoss(pred = model$fitted.values, actual = model$y),
      best.method = ifelse(metrics$youden > metrics$closest.topleft, yes = "youden", no = "closest.topleft"),
      threshold = metrics$threshold,
      auc = rocs$auc[1],
      sens = metrics$sensitivity,
      spe = metrics$specificity,
      PPV = metrics$ppv,
      NPV = metrics$npv,
      accuracy = metrics$accuracy,
      no.info.rate = max(
        sum(metrics[,"tp"], metrics[,"fn"]) / sum(metrics[,c("tp", "fp", "tn", "fn")]), 
        sum(metrics[,"fp"], metrics[,"fn"]) / sum(metrics[,c("tp", "fp", "tn", "fn")])
      ),
      balanced.accuracy = (metrics$sensitivity + metrics$specificity) / 2,
      precision = metrics$precision,
      f1 = 2 * ((metrics$recall * metrics$precision) / (metrics$recall + metrics$precision)),
      TP = metrics$tp,
      FN = metrics$fn,
      FP = metrics$fp,
      TN = metrics$tn
    )
    
  })
  
  newdata <- do.call(rbind, results) %>%
    mutate_at(7, list(~ round(x = ., digits = 2))) %>%
    mutate_at(c(8:17), list(~ round(x = ., digits = 3)))
  
  rownames(newdata) <- NULL
  
  return(newdata)
  
}

#' @rdname geneseng_summary_class_metrics
#' @export
geneseng_summary_multi_class_metrics <- function(data, group, direction = c("auto", ">", "<")){
  
  direction <- match.arg(direction)
  
  data <- data.frame(data)
  data[,group] <- as.factor(data[,group])
  newdf <- melt(data, id.vars = group)
  names(newdf)[1:2] <- c("group", "biomarker")
  
  # Multi classification - One vs All
  res01 <- lapply(1:length(levels(newdf$biomarker)), function(x){
    
    # MAKE SUMMARY METRICS FOR ALL BMKs
    dataset <- na.omit(newdf[newdf$biomarker == levels(newdf$biomarker)[x], ])
    
    # Define target class
    target <- as.character(unique(dataset$group))[x]
    
    # Confusion Matrix
    model <- multinom(group ~ value, data = dataset, trace = FALSE)
    pred <- levels(dataset$group)[apply(model$fitted.values, 1, function(x) which.max(x))]
    
    if(length(pred) != length(dataset$group)){
      idx <- which(is.na(dataset$value))
      tbl <- table(dataset$group[-c(idx)], factor(pred, levels = levels(dataset$group)))
    } else {
      tbl <- table(dataset$group, factor(pred, levels = levels(dataset$group)))
    }
    
    mtx <- caret::confusionMatrix(tbl)
    
    tbl01 <- lapply(1:length(unique(dataset$group)), function(x){
      data.frame(
        TP =  tbl[x,x],
        FN =  sum(tbl[x,-x]),
        FP =  sum(tbl[-x,x]),
        TN =  sum(tbl[-x,-x])
      )
    })
    
    tbl01 <- do.call(rbind, tbl01)
    
    # Compute Accuracy
    acc <- lapply(1:nrow(tbl01), function(x){
      (tbl01[x,1] + tbl01[x,4]) / sum(tbl01[x,])
    })
    
    no.info <- lapply(1:nrow(tbl01), function(x){
      max(
        (tbl01[x,1] + tbl01[x,2]) / sum(tbl01[x,]),
        (tbl01[x,3] + tbl01[x,4]) / sum(tbl01[x,])
      )
    })
    
    # Metrics
    if(length(pred) != length(dataset$group)){
      idx <- which(is.na(dataset$value))
      multi_logLoss <- MultiLogLoss(model$fitted.values, dataset$group[-c(idx)])
    } else {
      multi_logLoss <- MultiLogLoss(model$fitted.values, dataset$group)
    }
    

    
    metrics <- data.frame(
      
      # Biomarker & groups
      biomarker = rep(unique(dataset$biomarker), length(unique(dataset$group))),
      class1 = unique(dataset$group),
      class2 = rep("other", length(unique(dataset$group))),
      
      # Model
      model = "multinomial log. reg.",
      
      # Metrics
      multiLogLoss = multi_logLoss,
      sens = mtx$byClass[,1],
      spe = mtx$byClass[,2],
      PPV = mtx$byClass[,3],
      NPV = mtx$byClass[,4],
      accuracy = do.call(rbind, acc),
      no.info.rate = do.call(rbind, no.info),
      precision = mtx$byClass[,5],
      balanced.accuracy = mtx$byClass[,11],
      f1 = mtx$byClass[,7],
      
      # Table
      TP = tbl01[,1],
      FN = tbl01[,2],
      FP = tbl01[,3],
      TN = tbl01[,4]
    )
    
  })
  
  res01 <- do.call(rbind, res01)
  rownames(res01) <- NULL
  res01 <- res01 %>% mutate_at(c(6:14), list(~ round(x = ., digits = 3)))
  
  
  # Multi classification - One vs One
  comb <- combn(x = length(levels(data[,group])), m = 2)
  
  res02 <- lapply(1:ncol(comb), function(y){
    tmp <- data[data[,group] %in% levels(newdf$group)[comb[,y]], ]
    tmp[,group] <- factor(tmp[,group])
    geneseng_summary_binary_class_metrics(data = tmp, group = group, direction = direction)
    
  })
  
  res02 <- do.call(rbind, res02)
  res02 <- res02[order(res02$biomarker), ]
  rownames(res02) <- NULL
  
  list(
    `One vs All` = res01,
    `One vs One` = res02
  )
  
}
