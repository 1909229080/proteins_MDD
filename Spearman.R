#Packages
library(ggplot2)    
library(dplyr)    
library(purrr)      
library(gridExtra)  
library(grid)       
library(tidyr)      
library(openxlsx)   
library(boot)       
library(stringr)  
# Expand UniProt IDs separated by semicolons into multiple lines
spread_uniprot <- function(dat,
                           col_uniprot = "UniProt",
                           sep = ";") {
  dat_output <- tibble()
  for (i in 1:nrow(dat)) {
    uniprot_i <- dat[i,][[col_uniprot]] %>% 
      str_split(sep) %>% 
      unlist()
    for (each_uniprot in uniprot_i) {
      dat_output_i <- dat[i,]
      dat_output_i[[col_uniprot]] <- each_uniprot
      dat_output <- dat_output %>% 
        bind_rows(dat_output_i)
    }
  }
  dat_output
}
# Format P-value display
format_threshold <- function(p) {
  case_when(
    p < 0.001 ~ format(p, scientific = TRUE, digits = 1), 
    p == 1 ~ "1.00",                                      
    TRUE ~ as.character(round(p, 2))                       
  )
}
# Correlation analysis + Bootstrap confidence interval
local_cor_test <- function(data, 
                           ..., 
                           vars = NULL, 
                           vars2 = NULL,
                           method = "spearman",
                           conf.level = 0.95,
                           n_rep = 1000) {
  # Extract variables
  data$var1 <- data[[as.character(substitute(vars))]]
  data$var2 <- data[[as.character(substitute(vars2))]]
  # Traditional correlation test
  cor_test <- cor.test(data$var1, data$var2, 
                       method = method, 
                       exact = FALSE)
  # Bootstrap confidence interval calculation function
  boot_ci <- function(x, y, method, nrep, conf.level) {
    cor_fun <- function(d, i) cor(d$x[i], d$y[i], method = method)
    boot_res <- boot::boot(data = data.frame(x = x, y = y), 
                           statistic = cor_fun, 
                           R = nrep)
    boot_ci <- boot::boot.ci(boot_res, 
                             conf = conf.level, 
                             type = "bca")  # Use the BCA method
    c(boot_ci$bca[4], boot_ci$bca[5])
  }
  # Calculate the confidence interval
  ci <- boot_ci(data$var1, data$var2, 
                method = method, 
                nrep = n_rep, 
                conf.level = conf.level)
  # Return results
  tibble(
    var1 = "b.csf", 
    var2 = "b.plasma",
    cor = cor_test$estimate,    # correlation coefficient
    conf.low = ci[1],           # lower limit of the confidence interval
    conf.high = ci[2],          # upper confidence limit
    p.value = cor_test$p.value, 
    method = method,            
    n = sum(complete.cases(data$var1, data$var2))  
  )
}
# Generate Chart A (correlation at specified threshold)
correlation_bbb <- function(dat1 = cor_plasma, 
                            dat2 = cor_csf,
                            p_threshold = 1, 
                            title_text = "(A)",
                            method = "spearman",
                            n_rep = 1000) {
  # Step 1: Data Filtering ----
  cor_plasma_filtered <- dat1 %>% 
    filter(pval < p_threshold)  
  cor_csf_filtered <- dat2 %>% 
    filter(pval < p_threshold)  
  # Step 2: Obtain overlapping proteins ----
  intersect_uniprot <- intersect(cor_plasma_filtered$UniProt,
                                 cor_csf_filtered$UniProt)
  # Step 3: Merge data ----
  cor_both_protein <- bind_rows(
    cor_plasma_filtered %>% 
      filter(UniProt %in% intersect_uniprot),  
    cor_csf_filtered %>% 
      filter(UniProt %in% intersect_uniprot)   
  ) %>%
    mutate(b_se = paste0(b, "_", se)) %>%      
    select(UniProt, Tissue, b_se) %>%
    distinct() %>%
    pivot_wider(
      names_from = "Tissue", 
      values_from = "b_se"
    ) %>%                                      
    mutate(
      b.csf = map_dbl(str_split(CSF, "_"), ~ as.numeric(.x[1])),    
      se.csf = map_dbl(str_split(CSF, "_"), ~ as.numeric(.x[2])),    
      b.plasma = map_dbl(str_split(Plasma, "_"), ~ as.numeric(.x[1])),
      se.plasma = map_dbl(str_split(Plasma, "_"), ~ as.numeric(.x[2])) 
    )
  # Step 4: Conduct correlation analysis ----
  if (nrow(cor_both_protein) >= 2) {
    cor_test <- local_cor_test(
      cor_both_protein, 
      vars = b.csf, 
      vars2 = b.plasma,
      method = method,
      n_rep = n_rep
    )
    # Generate annotated text
    cor_label <- sprintf(
      "Spearman's ρ: %.2f\n95%% CI: (%.2f, %.2f)\nP = %s\nn = %d",
      cor_test$cor,
      cor_test$conf.low,
      cor_test$conf.high,
      format_threshold(cor_test$p.value),
      cor_test$n
    )
  } else {
    cor_label <- "Insufficient data\n(n < 2)"
  }
  # Step 5: Visualization ----
  ggplot(cor_both_protein, aes(x = b.csf, y = b.plasma)) +
    # error bars
    geom_errorbar(
      aes(
        ymin = b.plasma - 1.96 * se.plasma,
        ymax = b.plasma + 1.96 * se.plasma
      ),
      color = "grey70",
      width = 0
    ) +
    geom_errorbarh(
      aes(
        xmin = b.csf - 1.96 * se.csf,
        xmax = b.csf + 1.96 * se.csf
      ),
      color = "grey70",
      height = 0
    ) +
    # Trend line
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      color = "lightgray",
      alpha = 0.5,
      se = FALSE,
      linewidth = 0.8
    ) +
    # Scatter plot
    geom_point(
      color = "#2E86C1",  # 蓝色系
      size = 2.5
    ) +
    # reference line
    geom_hline(
      yintercept = 0,
      linetype = 2,
      color = "#E74C3C"   # 红色系
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      color = "#E74C3C"
    ) +
    labs(
      title = paste("(A) Correlation at P-value <", format_threshold(p_threshold)),
      x = "CSF effect size (β)",
      y = "Plasma effect size (β)"
    ) +
    annotate(
      "text",
      x = -Inf,          
      y = Inf,          
      label = cor_label,
      hjust = -0.1,      
      vjust = 1.2,       
      size = 5           
    ) +
    theme_bw(base_size = 18) +  
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16,
                                margin = margin(b = 5)  
                                 ),  
      axis.title = element_text(size = 14),         
      axis.text = element_text(size = 14),   
      axis.title.x = element_text(
        margin = margin(t = 6),  
        vjust = -0.2              
      ),
      axis.title.y = element_text(
        margin = margin(r = 6),  
        vjust = 1.0               
      ),
      axis.text.x = element_text(
        vjust = 0.5               
      ),
      axis.text.y = element_text(
        hjust = 0.5               
      )
    )
}
# Generate Chart B (Multi-threshold Dynamic Analysis) ----
correlation_bbb_step <- function(dat1 = cor_plasma,
                                 dat2 = cor_csf,
                                 thresholds = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
                                                0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 
                                                0.9, 0.95, 1),
                                 title_text = "(B) Correlation across different P-value thresholds",
                                 method = "spearman",
                                 n_rep = 1000,
                                 precomputed_threshold = NULL,
                                 precomputed_result = NULL) {
  
  # Step 1: Calculate correlations for each threshold iteratively ----
 results <- map_df(thresholds, ~ {
    if (!is.null(precomputed_threshold) && .x == precomputed_threshold) {
      return(tibble(
        threshold = .x,
        rho = precomputed_result$cor,
        lower = precomputed_result$conf.low,
        upper = precomputed_result$conf.high,
        p.value = precomputed_result$p.value,
        n = precomputed_result$n
      ))
    }
    cor_plasma_filtered <- dat1 %>% filter(pval < .x)
    cor_csf_filtered <- dat2 %>% filter(pval < .x)
    intersect_uniprot <- intersect(cor_plasma_filtered$UniProt, cor_csf_filtered$UniProt)
    n_intersect <- length(intersect_uniprot)
    if (n_intersect >= 2) {
      cor_both_protein <- bind_rows(
        cor_plasma_filtered %>% filter(UniProt %in% intersect_uniprot),
        cor_csf_filtered %>% filter(UniProt %in% intersect_uniprot)
      ) %>%
        mutate(b_se = paste0(b, "_", se)) %>%
        select(UniProt, Tissue, b_se) %>%
        distinct() %>%
        pivot_wider(names_from = "Tissue", values_from = "b_se") %>%
        mutate(
          b.csf = map_dbl(str_split(CSF, "_"), ~ as.numeric(.x[1])),
          b.plasma = map_dbl(str_split(Plasma, "_"), ~ as.numeric(.x[1]))
        )
      cor_test <- local_cor_test(cor_both_protein, vars = b.csf, vars2 = b.plasma,
                                 method = method, n_rep = n_rep)
      tibble(
        threshold = .x,
        rho = cor_test$cor,
        lower = cor_test$conf.low,
        upper = cor_test$conf.high,
        p.value = cor_test$p.value,
        n = cor_test$n
      )
    } else {
      tibble(
        threshold = .x,
        rho = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        p.value = NA_real_,
        n = n_intersect
      )
    }
  }) %>% 
    mutate(
      row_num = row_number(),
      label_y = ifelse(row_number() %% 2 == 1, 
                       upper + 0.15, 
                       lower - 0.15)
    )
  
  # Step 2: Visualization ----
  ggplot(results, aes(x = factor(format_threshold(threshold),
                                 levels = format_threshold(sort(thresholds, decreasing = TRUE))),
                      y = rho)) +
    geom_hline(yintercept = 0, color = "#E74C3C", linetype = 2) +
    geom_pointrange(aes(ymin = lower, ymax = upper), 
                    color = "#3498DB", 
                    size = 0.7,
                    position = position_dodge(width = 0.5)) +
    geom_text(aes(label = ifelse(is.na(rho), "NA",
                                 sprintf("ρ=%.2f\n(%.2f, %.2f)", 
                                         rho, lower, upper)),
                  y = label_y),
              size = 4.5,                                 
              position = position_nudge(x = 0.15)) +     
    geom_text(aes(label = ifelse(is.na(rho), "",
                                 sprintf("P=%s\nn=%d", 
                                         format_threshold(p.value), n)),
                  y = label_y),
              size = 4.5,                                 
              position = position_nudge(x = 0.15, y = -0.09)) +  
    scale_x_discrete(breaks = format_threshold(c(0.01, 0.05, 
                                                 0.25, 0.5, 0.75, 1))) +
    scale_y_continuous(
      limits = c(-0.4, 1),
      breaks = seq(-0.4, 1, 0.2)
    ) +
    labs(title = title_text, 
         x = "P-value threshold", 
         y = expression(paste("Spearman's ", rho))) +
    theme_minimal(base_size = 18) +  
    theme(
      plot.title = element_text(
        hjust = 0.5, 
        size = 16,
        margin = margin(b = 8)  
      ),
      axis.title.x = element_text(
        margin = margin(t = 1.5), 
        vjust = 4             
      ),
      axis.title.y = element_text(
        margin = margin(r = 4),  
        vjust = -0.2             
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 2,            
        size = 14
      ),
      axis.text.y = element_text(size = 14), 
      axis.title = element_text(size = 16),  
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(r = 20, l = 10)
    )
}
# Data retrieval ----
#
mr_csf <- read.xlsx("mr_res_csf_MDD.xlsx")  
mr_csf <- mr_csf[,-12]
csf <- read.xlsx("CSF_protein_cis.xlsx")  
cor_csf <- mr_csf %>% 
  filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>% 
  left_join(csf, by = "exposure") %>%
  spread_uniprot() %>%  
  distinct(UniProt, .keep_all = TRUE)  
#
mr_plasma <- read.xlsx("mr_res_plasma_MDD.xlsx") 
mr_plasma <- mr_plasma[,-12]
plasma <- read.xlsx("Plasma_protein_cis.xlsx")  
cor_plasma <- mr_plasma %>% 
  filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>% 
  left_join(plasma, by = "exposure") %>% 
  spread_uniprot() %>%
  distinct(UniProt, .keep_all = TRUE)
# Primary analysis process ----
# Generate Image A and obtain the results
p_threshold_a <- 0.05
p1 <- correlation_bbb(p_threshold = p_threshold_a, title_text = "A", n_rep = 1000)
# Extract the calculation results from Diagram A
cor_plasma_filtered_a <- cor_plasma %>% filter(pval < p_threshold_a)
cor_csf_filtered_a <- cor_csf %>% filter(pval < p_threshold_a)
intersect_uniprot_a <- intersect(cor_plasma_filtered_a$UniProt, cor_csf_filtered_a$UniProt)
cor_both_protein_a <- bind_rows(
  cor_plasma_filtered_a %>% filter(UniProt %in% intersect_uniprot_a),
  cor_csf_filtered_a %>% filter(UniProt %in% intersect_uniprot_a)
) %>%
  mutate(b_se = paste0(b, "_", se)) %>%
  select(UniProt, Tissue,b_se) %>%
  distinct() %>%
  pivot_wider(names_from = "Tissue", values_from = "b_se") %>%
  mutate(
    b.csf = map_dbl(str_split(CSF, "_"), ~ as.numeric(.x[1])),
    b.plasma = map_dbl(str_split(Plasma, "_"), ~ as.numeric(.x[1]))
  )
cor_test_a <- local_cor_test(cor_both_protein_a, vars = b.csf, vars2 = b.plasma,
                             method = "spearman", n_rep = 1000)
# Pass precomputed results when generating Chart B
p2 <- correlation_bbb_step(precomputed_threshold = p_threshold_a, 
                           precomputed_result = cor_test_a,
                           n_rep = 1000)
# Composite graphic output
final_plot <- grid.arrange(
  p1, p2,
  ncol = 2,
  widths = c(1, 1.5),
  top = textGrob("Cross-tissue Protein Effect Size Correlation", 
                 gp = gpar(fontsize = 20, fontface = "bold",
                           lineheight = 1.2,
                           margin = margin(b = 10)  
                 ),
                 vjust = 0.5)              
)