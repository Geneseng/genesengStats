#' Identify best features using a logistic/linear reg. with stepwise algorithm
#' 
#' @param data provide a data.frame
#' @param group group by a categorical variable
#' @param direction direction of the stepwise algorithm
#' 
#' @import dplyr
#' @import pROC
#' @importFrom tidyr pivot_longer
#' @importFrom stats as.formula predict step
#' 
#' @export
geneseng_best_model <- function(data, group, direction =  "backward"){
  
  if(length(unique(data[,group])) == 2){
    
    geneseng_best_binary_log_model(data, group, direction =  "backward")
    
  }  else {
    
    geneseng_best_linear_reg_model(data, group, direction =  "backward")
  }
  
}


#' @rdname geneseng_best_model
#' @export
geneseng_best_binary_log_model <- function(data, group, direction =  "backward"){
  
  df <- data
  df[,group] <-  as.integer(factor(df[,group]))
  
##################### Best model performance ###################################
  
  best_model <- as.formula(paste(group, "~", ".")) %>%
    glm(data = na.omit(df)) %>%
    step(direction = direction, trace = 0)
  
  best_model_roc <- roc(
    response = best_model$y , 
    predictor = best_model$fitted.values, 
    ci = TRUE
  )
  
  model_metrics <- data.frame(
    biomarker = toString(names(best_model[["model"]])[-1]),
    class1 = ifelse(best_model_roc[["direction"]] == ">", yes = levels(factor(data[,group]))[1], no = levels(factor(data[,group]))[2]),
    class2 = ifelse(best_model_roc[["direction"]] == ">", yes = levels(factor(data[,group]))[2], no = levels(factor(data[,group]))[1]),
    AIC = round(best_model[["aic"]], digits = 2),
    best.method = "youden",
    threshold = as.numeric(
      coords(best_model_roc, "best", best.method = "youden", transpose = TRUE)[1]
    ),
    auc = best_model_roc$auc[1],
    auc_Lower =  best_model_roc$ci[1],
    auc_Upper =  best_model_roc$ci[3]
  )

  
##################### Singular BMKs performances ###############################
  
  singular_metrics <- df[,names(best_model[["model"]])] %>%
    pivot_longer(cols = !group, names_to = "biomarker") %>%
    rename(group = names(.)[1])  %>%
    group_by(biomarker) %>%
    do(
      glm_roc = roc(group ~ value, data = ., ci = TRUE, direction = best_model_roc[["direction"]])
    ) %>%
    summarise(
      biomarker = biomarker,
      class1 = model_metrics$class1[1],
      class2 = model_metrics$class2[1],
      auc = glm_roc$ci[2],
      `auc_Lower` = glm_roc$ci[1],
      `auc_Upper` = glm_roc$ci[3],
      `DeLong's test` = roc.test(best_model_roc, glm_roc, paired = FALSE, method = "delong")$p.value
    )

######################### Final performances ###################################
  
  combined_metrics <-  model_metrics %>%
    bind_rows(singular_metrics) %>%
    mutate_at(6, list(~ round(x = ., digits = 2))) %>%
    mutate_at(7:9, list(~ round(x = ., digits = 3))) %>%
    mutate_at(10, list(~ formatC(x = ., format = "e", digits = 2)))
  
  return(combined_metrics)
  
}

#' @rdname geneseng_best_model
#' @export
geneseng_best_linear_reg_model <- function(data, group, direction =  "backward"){
  
##################### Best model performance ###################################
  
  best_model <- as.formula(paste(group, "~", ".")) %>%
    lm(data = na.omit(data)) %>%
    step(direction = direction, trace = 0)

  best_aic  <- best_model[["anova"]][["AIC"]]
  
  model_metrics <- data.frame(
    biomarker = toString(names(best_model[["model"]])[-1]),
    target = group,
    model = "linear reg.",
    AIC = round(min(best_aic), digits = 2),
    rmse = sqrt(mean((data[,group] - best_model$fitted.values) ^ 2, na.rm = TRUE)),
    mse = mean((data[,group] - best_model$fitted.values) ^ 2, na.rm = TRUE),
    mae = mean(abs(data[,group] - best_model$fitted.values), na.rm = TRUE)
  )
  
  
##################### Singular BMKs performances ###############################
  
  singular_metrics <- data[,names(best_model[["model"]])] %>%
    geneseng_summary_reg_metrics(group = group)
  
######################### Final performances ###################################
  
  combined_metrics <- model_metrics %>%
    bind_rows(singular_metrics)
  
  return(combined_metrics)
  
}

