#' Summary of quantiative or categorical variables
#'
#' Univariate or bivariate description of variables
#' 
#' @param data provide a data.frame
#' @param group group by one variable
#' 
#' 
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats sd
#' 
#' 
#' @export
geneseng_summary_stats <- function(data, group = NULL){
  
  options(dplyr.summarise.inform = FALSE)
  data <- data.frame(data)
  
 if(missing(group) | is.null(group)){
   
   make_summary(data = data)
   
 } else {
   
    make_summary_by_group(data = data, group = group)
   
  }
  
}


#' @rdname geneseng_summary_stats
#' @export
make_summary <- function(data){
  
  n <- apply(data, 2, function(x) length(unique(x)))
  
  if(sum(n > 7) == 0){
    
    res01 <- NULL
    
  } else {
    
    res01 <- data[, which(n > 7), drop = FALSE] %>%
      pivot_longer(cols = names(.), names_to = "biomarker") %>%
      group_by(biomarker) %>% 
      summarise(
        n = n(),
        `n distinct` = n_distinct(value, na.rm = TRUE),
        min = min(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        iqr = IQR(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        `NA's` = sum(is.na(value)),
        `Shapiro's test` = shapiro.test(value)$p.value,
        `normality` = if_else(`Shapiro's test` >= 0.5, "yes", "no")
      ) %>%
      mutate_at(4:9, list(~round(., digits = 2))) %>%
      mutate_at(11, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      as.data.frame()
    
  }
  
  if(sum(n <= 7) == 0){
    
    res02 <- NULL
    
  } else {
    
    res02 <- data[, which(n <= 7), drop = FALSE] %>%
      pivot_longer(cols = names(.), names_to = "biomarker") %>%
      group_by(biomarker, value) %>%
      summarise(
        n = n()
      ) %>%
      mutate(percent = n / sum(n)) %>%
      mutate_at(4, list(~round(., digits = 2))) %>%
      as.data.frame()
    
  }
  
  list(
    quantitative = res01,
    categorical = res02
  )
  
}

#' @rdname geneseng_summary_stats
#' @export
make_summary_by_group <- function(data, group){
  
  n <- apply(data, 2, function(x) length(unique(x)))
  
  if(sum(n > 7) == 0){
    
    res01 <- NULL
    
  } else {
    
    res01 <- data[, which(n > 7 | names(n) %in% group), drop = FALSE] %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group = names(.)[1]) %>%
      group_by(biomarker, group) %>% 
      summarise(
        n = n(),
        `n distinct` = n_distinct(value, na.rm = TRUE),
        min = min(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        iqr = IQR(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        `NA's` = sum(is.na(value)),
        `Shapiro's test` = shapiro.test(value)$p.value,
        `normality` = if_else(`Shapiro's test` >= 0.5, "yes", "no")
      ) %>%
      mutate_at(5:10, list(~round(., digits = 2))) %>%
      mutate_at(12, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      as.data.frame()
    
  }
  
  if(sum(n <= 7) == 1){
    
    res02 <- NULL
    
  } else {
    
    res02 <- data[, which(n <= 7), drop = FALSE] %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group = names(.)[1]) %>%
      group_by(biomarker, group, value) %>%
      summarise(
        n = n()
      ) %>%
      mutate(percent = n / sum(n)) %>%
      mutate_at(5, list(~round(., digits = 2))) %>%
      as.data.frame()
    
  }
  
  list(
    quantitative = res01,
    categorical = res02
  )
  
}
