#' Summary of statistical tests 
#'
#' Statistical tests at 95% CI
#' 
#' @param data provide a data.frame
#' @param group group by a categorical variable
#'
#' 
#' @import dplyr
#' @import tidyr
#' @import stats  
#' 
#' @export
geneseng_summary_tests <- function(data, group){
  
  options(dplyr.summarise.inform = FALSE)
  data <- data.frame(data)
  
  suppressMessages(
    
    suppressWarnings(
      
      if(length(unique(data[,group])) == 2){
        
        make_two_samples_tests(data = data, group = group)
        
      } else {
        
        make_two_or_more_samples_tests(data = data, group = group)
        
      }
      
    )
    
  )
  
}

#' @rdname geneseng_summary_tests
#' @export
make_two_samples_tests <- function(data, group){
  
############################### Continuous #####################################
  
  n <- apply(data, 2, function(x) length(unique(x)))
  
  if(sum(n > 7) == 0){
    
    res01 <- NULL
    
  } else {
    
    suppressWarnings({
      res01 <- data[, which(n > 7 | names(n) %in% group), drop = FALSE] %>%
        pivot_longer(cols = !group, names_to = "biomarker") %>%
        rename(group = names(.)[1]) %>%
        group_by(biomarker) %>%
        do(
          mwtest = wilcox.test(value ~ group, data = ., paired = FALSE, correct = FALSE),
          ftest = var.test(value ~ group, data = .),
          ttest = t.test(value ~ group, data = ., paired = FALSE, var.equal = TRUE),
          welchtest = t.test(value ~ group, data = ., paired = FALSE, var.equal = FALSE)
        ) %>% 
        summarise(
          biomarker = biomarker,
          group = group,
          `Mann-Whitney's test` = mwtest$p.value,
          `F-test` = ftest$p.value,
          `T-test` = ttest$p.value,
          `Welch's t-test` = welchtest$p.value
        ) %>%
        mutate_at(3:6, list(~ formatC(x = ., format = "e", digits = 2))) %>%
        as.data.frame()
    })
    

    tbl <- table(data[,group])
    
    if(all(sapply(tbl, function(x) x == tbl))) {
      
      suppressWarnings({
      res01 <- data[, which(n > 7 | names(n) %in% group), drop = FALSE] %>%
        pivot_wider(names_from = group, values_from = value) %>%
        unnest() %>%
        na.omit() %>%
        pivot_longer(cols = !biomarker, names_to = "group") %>%
        group_by(biomarker) %>%
        do(
          pairedttest = t.test(value ~ group, data = ., paired = TRUE, var.equal = TRUE),
          wilcoxtest = wilcox.test(value ~ group, data = ., paired = TRUE)
        ) %>%
        summarise(
          `Paired t-test` = pairedttest$p.value,
          `Wilcox's test` = wilcoxtest$p.value
        ) %>%
        mutate_at(1:2, list(~ formatC(x = ., format = "e", digits = 2))) %>%
        cbind(res01, .) %>%
        as.data.frame()
      })
      
    }
    
  }
  
############################### Categorical ####################################
  
  m <- data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
    apply(2, function(x) length(unique(x)))
  
  if(sum(m == 2) == 1 & sum(m > 2) == 0){
    
    res02 <- NULL
    
  } else if(sum(m == 2) > 1 & sum(m > 2) == 0){
    
    res02 <- data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
      select(which(m == 2)) %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group2 = names(.)[1]) %>%
      group_by(biomarker) %>%
      summarise(
        biomarker = unique(biomarker),
        group = group,
        `Chi-squared test` = chisq.test(group2, value)$p.value,
        `Fisher's exact test` = fisher.test(group2, value)$p.value,
        `McNemar's test` = mcnemar.test(group2, value)$p.value
      ) %>%
      mutate_at(3:5, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      as.data.frame()
    
  } else if(sum(m == 2) == 1 & sum(m > 2) >= 1){
    
    res02 <- data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
      select(which(m > 2 | names(m) %in% group)) %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group2 = names(.)[1]) %>%
      group_by(biomarker) %>%
      summarise(
        biomarker = unique(biomarker),
        group = group,
        `Chi-squared test` = chisq.test(group2, value)$p.value,
      ) %>%
      mutate_at(3, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      as.data.frame()
    
  } else if(sum(m == 2) >= 1 & sum(m > 2) >= 1){
    
    res02 <- data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
      select(which(m == 2)) %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group2 = names(.)[1]) %>%
      group_by(biomarker) %>%
      summarise(
        biomarker = unique(biomarker),
        group = group,
        `Chi-squared test` = chisq.test(group2, value)$p.value,
        `Fisher's exact test` = fisher.test(group2, value)$p.value,
        `McNemar's test` = mcnemar.test(group2, value)$p.value
      ) %>%
      mutate_at(3:5, list(~ formatC(x = ., format = "e", digits = 2)))

    
    res03 <- data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
      select(which(m > 2 | names(m) %in% group)) %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group2 = names(.)[1]) %>%
      group_by(biomarker) %>%
      summarise(
        biomarker = unique(biomarker),
        group = group,
        `Chi-squared test` = chisq.test(group2, value)$p.value
      ) %>%
      mutate_at(3, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      full_join(res02, .) %>%
      as.data.frame()
    
  }
  
  list(
    quantitative = res01,
    categorical = res02
  )
  
}

#' @rdname geneseng_summary_tests
#' @export
make_two_or_more_samples_tests <- function(data, group){
  
  n <- apply(data, 2, function(x) length(unique(x)))
  
############################## Continuous ######################################
  
  if(sum(n > 7) == 1){
    
    res01 <- NULL
    
  } else {
    
    res01 <- data[, which(n > 7 | names(n) %in% group), drop = FALSE] %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group = names(.)[1]) %>%
      group_by(biomarker) %>%
      do(
        kw = kruskal.test(value ~ group, data = .),
        bart = bartlett.test(value ~ group, data = .),
        anova = summary(aov(value ~ group, data = .))[[1]][["Pr(>F)"]]
      ) %>% 
      summarise(
        biomarker = biomarker,
        group = group,
        `Kruskal-Wallis's test` = kw$p.value,
        `Bartlett's test` = bart$p.value,
        `One-way Anova` = anova[1]
      ) %>%
      mutate_at(3:5, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      as.data.frame()
    
    
    idx <- table(data[,group])
    
    if(all(sapply(idx, function(x) x == idx))){
      
      res01 <- data[, which(n > 7 | names(n) %in% group), drop = FALSE] %>%
        mutate(subject = rep(x = paste0("S_", 1:idx[1]), times = length(idx))) %>%
        pivot_longer(cols = !c(group, subject), names_to = "biomarker") %>%
        pivot_wider(names_from = group, values_from = value) %>%
        unnest() %>%
        na.omit() %>%
        pivot_longer(cols = !c(subject, biomarker), names_to = "group") %>%
        group_by(biomarker) %>%
        do(
          friedman = friedman.test(value ~ group | subject, data = .),
          repeated_anova = summary(aov(formula = value ~ group + Error(subject/group), data = .))
        ) %>%
        summarise(
          `Friedman's test` = friedman$p.value,
          `One-way Anova with RM` = repeated_anova[[2]][[1]][["Pr(>F)"]][1]
        ) %>%
        mutate_at(1:2, list(~ formatC(x = ., format = "e", digits = 2))) %>%
        cbind(res01, .)
      
    }
    
  }
  
############################### Categorical ####################################

  if(sum(n <= 7) == 1){
    
    res02 <- NULL
    
  } else {
    
    res02 <- data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
      pivot_longer(cols = !group, names_to = "biomarker") %>%
      rename(group2 = names(.)[1]) %>%
      group_by(biomarker)  %>% 
      summarise(
        biomarker = unique(biomarker),
        group = group,
        `Chi-squared test` = chisq.test(group2, value)$p.value
      ) %>%
      mutate_at(3, list(~ formatC(x = ., format = "e", digits = 2))) %>%
      as.data.frame()
    
    
    idx <- table(data[,group])
    
    
    if(all(sapply(idx, function(x) x == idx))){
      
      data[, which(n <= 7 | names(n) %in% group), drop = FALSE] %>%
        mutate(subject = rep(x = paste0("S_", 1:idx[1]), times = length(idx))) %>%
        pivot_longer(cols = !c(group, subject), names_to = "biomarker") %>%
        group_by(biomarker) %>%
        do(
          cochran = friedman.test(value ~ group | subject, data = .)
        ) %>% 
        summarise(
          `Cochran's Q test` = cochran$p.value,
        ) %>%
        mutate_at(1, list(~ formatC(x = ., format = "e", digits = 2))) %>%
        cbind(res02, .) %>%
        as.data.frame()
      
    }
    
  }
  
  list(
    quantitative = res01,
    categorical = res02
  )
}