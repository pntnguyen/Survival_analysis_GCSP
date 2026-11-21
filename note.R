library(tidyverse)
library(magrittr)

df %<>%
  select(1:31)

## Cancer cells are positive or negative for estrogen receptors

df_analysis <- df |>
  select(
    patient_id,
    overall_survival_months,
    overall_survival,
    death_from_cancer,
    chemotherapy,
    er_status,
    her2_status,
    pr_status,
    hormone_therapy,
    cohort
  ) |>
  mutate(event = if_else(overall_survival == 1, 0, 1))


surv_obj <- Surv(time = df_analysis$overall_survival_months, event = df_analysis$event)

library(survival)
library(ggsurvfit)

survfit2(Surv(overall_survival_months, event) ~ pr_status + her2_status, data = df_analysis) |>
  ggsurvfit() +
  add_confidence_interval() +
  facet_wrap(~her2_status) + # Facet by the generated strata
  theme_minimal()+
  theme(legend.position = "bottom")

library(gtsummary)

tbl_summary(
  df,
  by = pr,
  include = c(
    age_at_diagnosis,
    er_status,
    her2_status,
    tumor_size,
    tumor_stage,
    neoplasm_histologic_grade,
    lymph_nodes_examined_positive,
    chemotherapy,
    hormone_therapy,
    radio_therapy
  )
)

df |>
  tbl_summary()

df |>
  select(`pam50_+_claudin-low_subtype`) |>
  unique()



df |>
  count(pr, npi_cat) |>
  group_by(pr) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = npi_cat, y = pct)) +
  geom_col() +
  geom_text(aes(label = scales::percent(pct, accuracy = 1)), vjust = -0.3, size = 3) +
  facet_wrap(~pr) +
  labs(x = "NPI Category", y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df |>
  count(pr, `pam50_+_claudin-low_subtype`) |>
  group_by(pr) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = `pam50_+_claudin-low_subtype`, y = pct)) +
  geom_col() +
  geom_text(aes(label = scales::percent(pct, accuracy = 1)), vjust = -0.3, size = 3) +
  facet_wrap(~pr) +
  labs(x = "NPI Category", y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



df %<>%
  select(1:31) %>%
  mutate(
    pr = factor(pr_status, levels = c("Positive", "Negative")),
    npi_cat = cut(
      nottingham_prognostic_index,
      breaks = c(-Inf, 2.4, 3.4, 4.4, 5.4, 6.4, Inf),
      labels = c(
        "Excellent",
        "Good",
        "Moderate 1",
        "Moderate 2",
        "Poor",
        "Very poor"
      ),
      right = TRUE
    )
  )

install.packages("ggdist")
library(ggdist)
library(ggplot2)

df |>
  select(pr,
         age_at_diagnosis,
         lymph_nodes_examined_positive) |>
  pivot_longer(
    cols = c(
      age_at_diagnosis,
      lymph_nodes_examined_positive
    ),
    names_to = "variable",
    values_to = "value"
  ) |>
  ggplot(aes(x = pr, y = value, fill = pr)) +
  stat_halfeye(adjust = 0.5, width = 0.6, justification = -0.2, .width = 0, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  # stat_dots(side = "left", justification = 1.1, dotsize = 0.3) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "PR status", y = "Value") +
  theme_bw() +
  theme(legend.position = "none")


df$pr_status <- levels(c("Positive","Negative"))

## EDA for 5 cohorts

# Fit the survival object
# survival_months = time variable, overall_survival = event indicator (1 = event/death, 0 = censored)
surv_obj <- Surv(time = df_analysis$overall_survival_months, event = df_analysis$event)

# Fit KM curves stratified by er_status
fit <- survfit(surv_obj ~ cohort, data = df_analysis)

# Plot
ggsurvplot(
  fit,
  data = df_analysis,
  # pval = TRUE,
  # conf.int = TRUE,
  # risk.table = TRUE,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  legend.title = "ER Status",
  legend.labs = levels(df_analysis$cohort),
  palette = "jco",
  ggtheme = theme_minimal(),
  xlim = c(0,200)
)

## eda for each cohort


df_analysis |>
  group_by(cohort) |>
  count(er_status) |>
  ggplot(aes(x = er_status, y = n)) +
  geom_col() +
  facet_wrap(~cohort)

df_analysis |>
  group_by(cohort) |>
  count(hormone_therapy) |>
  ggplot(aes(x = hormone_therapy, y = n)) +
  geom_col() +
  facet_wrap(~cohort)

df |>
  ggplot(aes(x = tumor_size)) +
  geom_histogram(bins = 30) +
  facet_wrap(~cohort) +
  labs(x = "Tumor Size", y = "Count", title = "Tumor Size Distribution by Cohort")


df |>
  select(er_status,hormone_therapy,cohort) |>
  group_by(cohort) |>
  count()

## is hormone therapy associate with survival probability

fit2 <- survfit(surv_obj ~ hormone_therapy, data = df_analysis)

ggsurvplot(
  fit2,
  data = df_analysis,
  # pval = TRUE,
  # conf.int = TRUE,
  # risk.table = TRUE,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  legend.title = "ER Status",
  legend.labs = levels(df_analysis$hormone_therapy),
  palette = "jco",
  ggtheme = theme_minimal(),
  xlim = c(0,200)
)

df |>
  group_by(pr_status) |>
  count(hormone_therapy)

## explore by pr_status

fit3 <- survfit(surv_obj ~ pr_status, data = df)

ggsurvplot(
  fit3,
  data = df,
  # pval = TRUE,
  # conf.int = TRUE,
  # risk.table = TRUE,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  surv.median.line = "hv",
  # legend.title = "ER Status",
  legend.labs = levels(df$pr_status),
  palette = "jco",
  ggtheme = theme_minimal(),
  # xlim = c(0,200)
  {}
)


df |>
  group_by(pr_status) |>
  count(hormone_therapy) |>
  mutate(per = n/sum(n)) |>
  ggplot(aes(x = hormone_therapy, y = per))+
  geom_col()+
  scale_y_continuous(limits = c(0,1))+
  facet_wrap(~pr_status)+
  theme_bw()

lung








### week 5 exercise
library(haven)

## vets data
vets <- read_dta("data/vets.dta")

colnames(vets)

vets |>
  select(survt,ct1,status)


addicts <- read.table(
  "./data/addicts.dat",
  skip = 19,
  # skip comment/metadata lines above the data
  header = FALSE,
  col.names = c(
    "id",
    "clinic",
    "status",
    "survival_time",
    "prison",
    "methadone_dose"
  )
)

library(survival)
library(ggsurvfit)


df_km <- df |>
  select(overall_survival_months, overall_survival, death_from_cancer, brca1) |>
  mutate(
    surv = case_when(death_from_cancer == "Living" ~ 0, .default = 1),
    brca1_group = cut(
      brca1,
      breaks = quantile(brca1, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("Low", "Medium", "High"),
      include.lowest = TRUE
    )
  )

survfit2(Surv(overall_survival_months, surv) ~ brca1_group, data = df_km) |>
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  add_pvalue(caption = "Log-rank {p.value}") +
  labs(x = "Time (months)", y = "Survival probability", color = "BRCA1 expression", fill = "BRCA1 expression") +
  scale_ggsurvfit() +
  theme_minimal() +
  theme(legend.position = "bottom")

library(msigdbr)

library(dplyr)

# Pull Hallmark gene sets (older msigdbr versions use category = "H" instead of collection = "H")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

pathway_map <- tribble(
  ~gs_name,                           ~pathway,
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "Estrogen / hormone signaling",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",  "Estrogen / hormone signaling",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "PI3K-AKT-mTOR signaling",
  "HALLMARK_E2F_TARGETS",             "Cell cycle / proliferation",
  "HALLMARK_G2M_CHECKPOINT",          "Cell cycle / proliferation",
  "HALLMARK_DNA_REPAIR",              "DNA repair"
)

pathway_genes <- hallmark |>
  inner_join(pathway_map, by = "gs_name") |>
  distinct(pathway, gene_symbol) |>
  mutate(gene_lower = tolower(gene_symbol))   # your columns are lowercase (e.g. brca1)


# Use the FULL data frame, not the one trimmed with select(1:31)
data_genes <- tolower(names(df))

pathway_genes <- pathway_genes |>
  mutate(in_data = gene_lower %in% data_genes)

# Summary: how many genes per pathway are available
pathway_genes |>
  group_by(pathway) |>
  summarise(
    total_genes = n(),
    available   = sum(in_data),
    pct         = round(100 * available / total_genes, 1)
  )

# Which genes are available, per pathway
pathway_genes |>
  filter(in_data) |>
  group_by(pathway) |>
  summarise(genes = paste(gene_symbol, collapse = ", "))


library(dplyr)
library(purrr)
library(tibble)

# 1. Map each available gene to its actual column name in the data
col_lookup <- tibble(col = names(df), gene_lower = tolower(names(df)))

pathway_cols <- pathway_genes |>
  filter(in_data) |>
  inner_join(col_lookup, by = "gene_lower")   # adds the `col` column

# 2. Z-score each gene column (mean 0, SD 1 across patients)
gene_cols <- unique(pathway_cols$col)

z_mat <- df |>
  select(all_of(gene_cols)) |>
  mutate(across(everything(), ~ as.numeric(scale(.x))))

# 3. Pathway score = mean z-score of the genes in G_p, per patient
pathway_scores <- pathway_cols |>
  split(~pathway) |>
  map(\(x) rowMeans(z_mat[, x$col, drop = FALSE], na.rm = TRUE)) |>
  as_tibble() |>
  rename(
    score_hormone    = `Estrogen / hormone signaling`,
    score_pi3k_akt   = `PI3K-AKT-mTOR signaling`,
    score_cell_cycle = `Cell cycle / proliferation`,
    score_dna_repair = `DNA repair`
  )

# 4. Attach to your data
df_full <- bind_cols(df, pathway_scores)

summary(pathway_scores)

### clustering algorithm

library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(survival)
library(ggsurvfit)

# 1. Prepare the score matrix (complete cases only, rescaled so all 4 scores have equal weight)
score_df <- df_full |>
  mutate(patient_row = row_number()) |>
  select(patient_row, starts_with("score_")) |>
  drop_na()

score_mat <- score_df |>
  select(starts_with("score_")) |>
  scale()

# 2. Choose k BEFORE looking at survival (elbow + silhouette)
fviz_nbclust(score_mat, kmeans, method = "wss")
fviz_nbclust(score_mat, kmeans, method = "silhouette")

# 3. Run k-means for k = 2 and k = 3
set.seed(123)
km2 <- kmeans(score_mat, centers = 2, nstart = 50)
km3 <- kmeans(score_mat, centers = 3, nstart = 50)

# 4. Attach cluster labels back to the data
score_df <- score_df |>
  mutate(cluster_k2 = factor(km2$cluster),
         cluster_k3 = factor(km3$cluster))

df_clust <- df_full |>
  mutate(patient_row = row_number()) |>
  left_join(score_df |> select(patient_row, cluster_k2), by = "patient_row")

# 5. Profile each cluster (cluster numbers are arbitrary, so look at the mean scores)
df_clust |>
  filter(!is.na(cluster_k3)) |>
  group_by(cluster_k3) |>
  summarise(n = n(), across(starts_with("score_"), mean))

# Visual profile
df_clust |>
  filter(!is.na(cluster_k3)) |>
  select(cluster_k3, starts_with("score_")) |>
  pivot_longer(-cluster_k3, names_to = "pathway", values_to = "score") |>
  ggplot(aes(x = cluster_k3, y = score, fill = cluster_k3)) +
  geom_boxplot() +
  facet_wrap(~pathway) +
  theme_bw() +
  theme(legend.position = "none")

# 6. KM curves (use the event definition you decided on; overall survival shown here)
df_clust <- df_clust |>
  mutate(event = if_else(overall_survival == 0, 1, 0))

survfit2(Surv(overall_survival_months, event) ~ cluster_k2, data = df_clust) |>
  ggsurvfit() +
  add_confidence_interval() +
  # add_risktable() +
  add_pvalue(caption = "Log-rank {p.value}") +
  labs(x = "Time (months)", y = "Overall survival probability",
       color = "Cluster", fill = "Cluster") +
  scale_ggsurvfit() +
  theme_minimal() +
  theme(legend.position = "bottom")

df |>
  select(all_of(expression_cols)) |>
  select(1:10) |>
  skim()

## exploratory analysis

library(gtsummary)

df_clust |>
  mutate(npi_cat = cut(
    nottingham_prognostic_index,
    breaks = c(-Inf, 2.4, 3.4, 4.4, 5.4, 6.4, Inf),
    labels = c("Excellent", "Good", "Moderate 1", "Moderate 2", "Poor", "Very poor"),
    right = TRUE
  )) |>
  tbl_summary(
    by = cluster_k2,
    include = c(
      age_at_diagnosis,
      er_status,
      her2_status,
      pr_status,
      chemotherapy,
      hormone_therapy,
      radio_therapy,
      npi_cat
    ),
    label = list(
      age_at_diagnosis ~ "Age at diagnosis (years)",
      er_status        ~ "ER status",
      her2_status      ~ "HER2 status",
      pr_status        ~ "PR status",
      chemotherapy     ~ "Chemotherapy",
      hormone_therapy  ~ "Hormone therapy",
      radio_therapy    ~ "Radiotherapy",
      npi_cat          ~ "Nottingham Prognostic Index group"
    ),
    statistic = list(
      age_at_diagnosis ~ "{mean} ({sd})"
    ),
    digits = list(
      age_at_diagnosis ~ 1
    )
  )









