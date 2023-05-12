#' Evaluate correlation between all kind of variables
#'
#' @param data a data.frame
#'
#' @import dplyr 
#' @import tidyr
#' @importFrom reshape2 dcast
#' @importFrom utils combn
#' @importFrom stats cor
#' @importFrom stats cor.test
#' @importFrom purrr map2
#' @importFrom broom tidy
#' 
#' @export
geneseng_summary_corr <- function(data){
  
  data <- data.frame(data)
  n <- apply(data, 2, function(x) length(unique(x)))
  
  res01 <- continuous_to_continuous(data)
  res02 <- continuous_to_categ(data)
  res03 <- categ_to_categ(data = data)
  
  list(
    `Continuous vs Continuous` = res01,
    `Continuous vs Categorical` = res02,
    `Categorical vs Categorical` = res03
  )
  
}

#' @rdname geneseng_summary_corr
#' @export
continuous_to_continuous <- function(data){
  
  n <- apply(data, 2, function(x) length(unique(x)))
  
  if(sum(n > 7) < 2){
    
    res <- NULL
    
  } else {
    
    suppressMessages(
      
      suppressWarnings(
        
        res <- data[, which(n > 7), drop = FALSE] %>%
          names() %>%
          combn(m = 2) %>%
          t() %>%
          data.frame(stringsAsFactors = FALSE) %>%
          mutate(
            map2(X1, X2, ~tidy(spearmanCI(data[,.x], data[,.y])))
          ) %>%
          unnest() %>%
          dcast(X1+X2 ~ names) %>% 
          rename("Var1" = "X1", "Var2" = "X2") %>%
          mutate(
            method = rep("Spearman", times = nrow(.))
          ) %>%
          select(1,2,6,3,4,5) %>%
          mutate(
            pvalue = unlist(map2(Var1, Var2, ~tidy(cor.test(data[,.x], data[,.y], method = "spearman")$p.value)))
          ) %>%
          mutate_at(7, list(~ formatC(x = ., format = "e", digits = 2)))
        
      )
      
    )
    
  }
  
  return(res)
  
}

#' @rdname geneseng_summary_corr
#' @export
continuous_to_categ <- function(data){
  
  n <- apply(data, 2, function(x) length(unique(x)))
  continuous <- data[, which(n > 7), drop = FALSE]
  categorical <- data[, which(n <= 7), drop = FALSE]
  
  if(ncol(continuous) == 0 | ncol(categorical) == 0){
    
    res <- NULL
    
  } else {
    
    lst <- lapply(1:ncol(categorical), function(x){
      
      res <- lapply(1:ncol(continuous), function(y){
        
        df <- data.frame(
          Var1 = names(categorical)[x],
          Var2 = names(continuous)[y],
          method = "PBC",
          t(
            spearmanCI(
              var1 = as.numeric(as.factor(categorical[, names(categorical)[x]])),
              var2 = continuous[, names(continuous)[y]]
            )
          ),
          pvalue = cor.test(
            as.numeric(as.factor(categorical[, names(categorical)[x]])),
            continuous[, names(continuous)[y]],
            method = "spearman"
          )$p.value
        )
        
        df$pvalue <- formatC(x = df$pvalue, format = "e", digits = 2)
        
        return(df)
        
      })
      
      do.call(rbind, res)
      
    })
    
    res <- do.call(rbind, lst)
    rownames(res) <- NULL
    
  }
  
  return(res)
  
}


#' @rdname geneseng_summary_corr
#' @export
categ_to_categ <- function(data){
  
  n <- apply(data, 2, function(x) length(unique(x)))
  
  if(sum(n <= 7) < 2){
    
    res <- NULL
    
  } else {
    
    res <- suppressWarnings({
      data[, which(n <= 7), drop = FALSE] %>%
        names() %>%
        combn(m = 2) %>%
        t() %>%
        data.frame(stringsAsFactors = FALSE) %>%
        mutate(
          map2(X1, X2, ~tidy(chisq.test(data[,.x], data[,.y])))
        ) %>%
        unnest() %>%
        mutate_at(3, list(~round(., digits = 2))) %>%
        mutate_at(4, list(~ formatC(x = ., format = "e", digits = 2))) %>%
        mutate_at(6, list(~ ifelse(
          test = . == "Pearson's Chi-squared test", 
          yes = "Chi-squared test",  
          no = "Chi-squared test with Yates's correction"))
        ) %>%
        select(1,2,6,5,3,4) %>%
        rename("Var1" = "X1", "Var2" = "X2", "pvalue" = "p.value") %>%
      as.data.frame()
    })
    
  }
  
  return(res)
  
}
