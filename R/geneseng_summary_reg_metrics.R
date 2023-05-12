#' Summary of regression metrics 
#'
#' Evaluate biomarker individual performances using best regression metrics
#' 
#' @param data provide a data.frame
#' @param group group by a continuous variable
#'
#' @examples
#' library(genesengStats)
#' geneseng_summary_reg_metrics(data = mtcars, group = "mpg")
#' geneseng_summary_reg_metrics(data = iris, group = "Sepal.Length")
#' 
#' @import dplyr
#' @importFrom stats lm
#' 
#' @details 
#' `Linear regression` was used for regression problem.
#' 
#' @export
geneseng_summary_reg_metrics <- function(data, group){
  
  nb <- data[,names(data)[!(names(data) %in% group)]]
  
  res <- lapply(1:ncol(nb), function(x){
    
    # Rearrange columns
    newdata <- cbind.data.frame(resp = data[,group], nb)
    newdf <- newdata[, c("resp", names(nb)[x])]
    
    # Create model
    model <- lm(resp ~ ., data = newdf)
    pred <- model$fitted.values
    
    # Metrics
    df <- data.frame(
      biomarker = names(nb)[x],
      target = group,
      model = "linear reg.",
      rmse = sqrt(mean((data[,group] - pred) ^ 2, na.rm = TRUE)),
      mse = mean((data[,group] - pred) ^ 2, na.rm = TRUE),
      mae = mean(abs(data[,group] - pred), na.rm = TRUE)
    )
    
    return(df)
    
  })
  
  do.call(rbind, res)

}


