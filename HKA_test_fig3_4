library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(purrr)

df <- read.csv("supp_3_030326.csv")

add_ratio_stems <- function(tbl) {
  tbl %>%
    mutate(
      stem_sin_div= stem_sin / stem_div,
      stem_pol_div=stem_pol/stem_div,
      stem_pd=stem_pol/stem_div,
      loop_pd=loop_pol/loop_div,
      stem_sin_pd=stem_sin/stem_div,
      loop_sin_pd=loop_sin/loop_div
    )
}
df <- add_ratio_stems(df)

hominini_df       <- filter(df, G4.group == "hominini")
ape_df            <- filter(df, G4.group == "ape")
human_spec_df     <- filter(df, G4.group == "human-specific")


add_ratio <- function(tbl) {
  tbl %>%
    mutate(
      stem_loop_pd=stem_pd/loop_pd,
      
      stem_loop_sin_pd=stem_sin_pd/loop_sin_pd,
      ngnr_total_sin = total_sin[annotation == "NFNR"][1],
      ngnr_total_div = total_div[annotation == "NFNR"][1],
      ngnr_total_pol = total_pol[annotation == "NFNR"][1],
      ngnr_stem_sin = stem_sin[annotation == "NFNR"][1],
      ngnr_stem_pol = stem_pol[annotation == "NFNR"][1],
      ngnr_stem_div = stem_div[annotation == "NFNR"][1],
      ngnr_loop_sin = loop_sin[annotation == "NFNR"][1],
      ngnr_loop_pol = loop_pol[annotation == "NFNR"][1],
      ngnr_loop_div = loop_div[annotation == "NFNR"][1],      
      ngnr_total_sin_div = total_sin_div[annotation == "NFNR"][1],   
      total_sin_div_ratio = total_sin_div / ngnr_total_sin_div,            
      ngnr_total_pol_div = total_pol_div[annotation == "NFNR"][1],   
      total_pol_div_ratio = total_pol_div / ngnr_total_pol_div,            
      ngnr_stem_sin_div = stem_sin_div[annotation == "NFNR"][1],   
      stem_sin_div_ratio = stem_sin_div / ngnr_stem_sin_div,            
      ngnr_stem_pol_div = stem_pol_div[annotation == "NFNR"][1],  
      stem_pol_div_ratio = stem_pol_div / ngnr_stem_pol_div
    )
}


hominini_df   <- add_ratio(hominini_df)
ape_df        <- add_ratio(ape_df)
human_spec_df <- add_ratio(human_spec_df)


result_df <- bind_rows(hominini_df, ape_df)


fish_df <- result_df %>%
  filter(annotation!="genome wide") %>%
  select(G4.group, annotation, stem_sin, stem_div, loop_sin, loop_div) %>%
  group_by(G4.group, annotation) %>%
  nest() %>%  # <-- create 'data' column
  mutate(
    fish_result = map(data, function(tbl) {
      # Construct the 2x2 contingency table
      mat <- matrix(c(
        tbl$stem_sin, tbl$stem_div,
        tbl$loop_sin, tbl$loop_div
      ), nrow = 2, byrow = TRUE)
      
      rownames(mat) <- c("stem", "loop")
      colnames(mat) <- c("polymorphic", "divergence")
      
      # Only test if no structural issues
      if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
        fisher.test(mat)
      } else {
        NA
      }
    }),
    p_value = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$p.value, NA)),
    odds_ratio = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$estimate, NA)),
    ci_lower = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$conf.int[1], NA)),
    ci_upper = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$conf.int[2], NA))
  ) %>%
  select(G4.group, annotation, p_value, odds_ratio, ci_lower, ci_upper)

run_fisher_test <- function(df, poly_1, div_1, poly_2, div_2, group_name) {
  df %>%
    select(G4.group, annotation,
           poly_1 = {{poly_1}},
           div_1  = {{div_1}},
           poly_2 = {{poly_2}},
           div_2  = {{div_2}}) %>%
    
    group_by(G4.group, annotation) %>%
    nest() %>%
    mutate(
      fish_result = map(data, function(tbl) {
        mat <- matrix(c(
          tbl$poly_1, tbl$div_1,
          tbl$poly_2, tbl$div_2
        ), nrow = 2, byrow = TRUE)
        
        rownames(mat) <- c("regions1", "region2")
        colnames(mat) <- c("polymorphic", "non_polymorphic")
        
        if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
          fisher.test(mat)
        } else {
          NA
        }
      }),
      p_value = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$p.value, NA)),
      odds_ratio = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$estimate, NA)),
      ci_lower = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$conf.int[1], NA)),
      ci_upper = map_dbl(fish_result, ~ ifelse(is.list(.x), .x$conf.int[2], NA)),
      test_label = group_name
    ) %>%
    ungroup() %>%
    select(G4.group, annotation, test_label, p_value, odds_ratio, ci_lower, ci_upper)
}


test_groups2 <-list(list("stem_pol", "loop_pol", "stem_length","loop_length","Polymorphism"),
                    list("stem_div", "loop_div", "stem_length","loop_length","Divergence"),
                    list("stem_sin", "loop_sin", "stem_length","loop_length","Singletons")
)


result_df <- bind_rows(hominini_df, ape_df, human_spec_df)

result_df[, c("stem_length", "loop_length")] <- 
  lapply(result_df[, c("stem_length", "loop_length")], as.numeric)
fish_df_fig3 <- map_dfr(test_groups2, ~ run_fisher_test(result_df, .x[[1]], .x[[2]], .x[[3]],.x[[4]],.x[[5]]))


annotation_order <- c("Genome wide","5’UTRs", "3’UTRs","CDS", "Introns", "Promoters", "Enhancers", "CpG islands", "Replication origins","Non-protein-coding genes","NFNR")




library(ggplot2)


fish_df_fig3 <- fish_df_fig3 %>% filter(test_label %in% c("Divergence","Polymorphism"))
#fish_df_fig3$Variant_Type <- fish_df_fig3$test_label
annotation_order <- c("Genome wide","5’UTRs", "3’UTRs","CDS", "Introns", "Promoters", "Enhancers", "CpG islands", "Replication origins","Non-protein-coding genes","NFNR")
fish_df_fig3$annotation <- factor(fish_df_fig3$annotation, levels = annotation_order)




library(ggplot2)
library(patchwork)
library(dplyr)


custom_colors <- c("#66C2A5", "#FFD92F")

pd <- position_dodge(width = 0.5)


plot_panel_custom <- function(data, group_name, y_label) {
  pd <- position_dodge(width = 0.5)
  custom_colors <- c("#66C2A5", "#FFD92F") 
  
  p <- data %>% 
    filter(G4.group == group_name) %>% 
    ggplot(aes(x = annotation, y = odds_ratio, color = test_label)) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper), 
      width = 0.2, 
      size = 0.8, 
      position = pd
    ) +
    geom_point(size = 3, position = pd) +
    #geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = custom_colors) +
    coord_cartesian(ylim = c(0, 2.5)) + 
    labs(
      title = group_name,
      y = y_label, 
      x = "Annotation", 
      color = "Variant Type"
    ) +
    theme_bw() +
    theme(

      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),

      panel.grid.major.y = element_line(color = "grey92"),
      panel.grid.minor.y = element_line(color = "grey95"),

      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.background = element_rect(fill = "white")
    )
  return(p)
}

p1_with_legend <- plot_panel_custom(fish_df_fig3, "ape", "Ratio of variant frequency between\nstems and loops at ape G4s") + theme(legend.position = "top")
p2_no_legend <- plot_panel_custom(fish_df_fig3, "hominini", "Ratio of variant frequency between\nstems and loops at Hominini G4s") + theme(legend.position = "none")
p3_no_legend <- plot_panel_custom(fish_df_fig3, "human-specific", "Ratio of variant frequency between\nstems and loops at human−specific G4s") + theme(legend.position = "none")

figure3_final <- plot_grid(
  p1_with_legend, 
  p2_no_legend, 
  p3_no_legend,
  labels = c('A', 'B', 'C'),
  ncol = 1,
  align = 'v',
  rel_heights = c(1, 1, 1)
)


figure3_final

ggsave("/Users/xzhang/Desktop/Figure3_cowplot2.pdf", figure3_final, width = 8, height = 12)


fish_df_all <- map_dfr(test_groups1, ~ run_fisher_test(result_df, .x[[1]], .x[[2]], .x[[3]],.x[[4]],.x[[5]]))



result_df <- bind_rows(hominini_df, ape_df, human_spec_df)
#result_df[, c("stem_length", "loop_length")] <- gsub(",","",result_df[, c("stem_length", "loop_length")]) 
result_df[, c("stem_length", "loop_length")] <- 
  lapply(result_df[, c("stem_length", "loop_length")], as.numeric)

panel_order <- c("Singletons to divergence between stems and loops (S/D)", "Polymorphism to divergence between stems and loops (P/D)", "S/D ratios for whole G4s in each functional region vs. for whole G4s in NFNR", "P/D raios for whole G4s in each functional region vs. for whole G4 in NFNR", "S/D ratios for G4 stems in each functional region vs. in NFNR", "P/D ratios for G4 stems in each functional region vs. in NFNR")
fish_df_all$test_label <- factor(fish_df_all$test_label, levels = panel_order)

fish_df_all <- fish_df_all %>%
  distinct(test_label) %>%
  arrange(test_label) %>%
  mutate(tagged_label = paste0(LETTERS[1:n()], ". ", test_label)) %>%
  right_join(fish_df_all, by = "test_label")
write.csv(fish_df_all,"fish_df_fig4.csv")


fish_df_all<- read.csv("fish_df_fig4.csv")
fig4_data <- fish_df_all %>% 
  filter(G4.group!="human-specific" & annotation!="genome wide" & !test_label %in% c("S/D ratios for whole G4s in each functional region vs. for whole G4s in NFNR", "P/D raios for whole G4s in each functional region vs. for whole G4 in NFNR"))

annotation_order <- c("5'UTR", "3'UTR","CDS", "Intron", "Promoter", "Enhancer", "CpG islands", "Replication origins","Non-protein-coding genes","NFNR")
fig4_data$annotation <- factor(fig4_data$annotation, levels = annotation_order)
p4_new <- ggplot(fig4_data, aes(x = annotation, y = odds_ratio, color = G4.group)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.6),
                width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~ tagged_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Odds Ratios of Polymorphism by Region and G4 Group",
    x = "Functional Annotation",
    y = "Odds Ratio (95% CI)",
    color = "G4 Group"
  ) +

  scale_color_manual(
    values = c(
      "ape" = "#E6A8B6",  # soft pink
      "hominini" = "#7A4F2C"   # deep brown
    )
  ) +
  theme(
    strip.text = element_text(face = "bold", hjust = 0),  # top-left corner
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major.y = element_line(color = "gray85"),
    panel.spacing = unit(1.2, "lines")
  )
p4_new
ggsave("p4.pdf",p4_new)
