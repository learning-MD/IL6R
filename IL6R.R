library(mice)
library(epitools)
library(survival)
library(tidyverse)
library(gtsummary)
library(sjmisc)
library(survival)
library(survminer)
library(metafor)
library(meta)
library(dmetar)

# Loading data from VUMC for general analyses
data <- read_csv("BioVU_SD.csv")
data <- data %>% distinct(grid, .keep_all = TRUE)

# Adding Yes/No column for homozygous IL6R protection
data <- data %>%
  mutate(
    il6r = case_when(
      Count_IL6R_D358A == 2 ~ 1,
      TRUE ~ 0
    )
  )

# Multiple imputation
set.seed(1)
list <- data %>% dplyr::select(-c(cav_mention, last_fup, date_of_death, dna_collection,
                                  grid, htx_date, dob, age_at_death, age_at_fup,
                                  Count_IL6R_D358A))
imp2 <- mice(list, m = 16)

# Loading CUIMC data
cumc <- readxl::read_xlsx("CUMC_data.xlsx")
cumc$gender <- ifelse(cumc$gender=="M",1,0)

cumc <- cumc %>% mutate(
  il6r = case_when(
    Count_IL6R_D358A == 2 ~ 1,
    TRUE ~ 0
  )
)

cumc.table_1 <- cumc %>% select(il6r, age_at_htx, gender, ischemic, htn, diabetes, statin, TwoRplus, amr, dsa, cav, deceased)
cumc.table_1 <- cumc.table_1 %>% rename(acr = TwoRplus)

# Table 1 demographics
mat <- complete(imp2)
table_2 <- mat %>% select(il6r, age_at_htx, gender, ischemic, htn, diabetes, statin, acr, amr, dsa, cav, deceased)

vumc.imp.il6r <- tbl_summary(table_2,
                             by = il6r,
                             missing = "no",
                             label = list(age_at_htx ~ "Age at transplant (years)",
                                          gender ~ "Gender (male)",
                                          ischemic ~ "Ischemic etiology",
                                          htn ~ "Hypertension",
                                          diabetes ~ "Diabetes mellitus",
                                          statin ~ "Statin use",
                                          acr ~ "Acute cellular rejection",
                                          amr ~ "Antibody-mediated rejection",
                                          dsa ~ "Donor specific antibodies",
                                          cav ~ "Cardiac allograft vasculopathy",
                                          deceased ~ "Deceased")) %>%
  add_p() %>%
  modify_header(label = "**Demographics by IL6R signaling status**",
                stat_1 = "**Preserved IL6R Signaling (N = 535)**",
                stat_2 = "**Impaired IL6R Signaling (N = 71)**") %>% bold_labels()

vumc.imp.il6r

cuimc.il6r <- tbl_summary(cumc.table_1,
                          by = il6r,
                          missing = "no",
                          label = list(age_at_htx ~ "Age at transplant (years)",
                                       gender ~ "Gender (male)",
                                       ischemic ~ "Ischemic etiology",
                                       htn ~ "Hypertension",
                                       diabetes ~ "Diabetes mellitus",
                                       statin ~ "Statin use",
                                       acr ~ "Acute cellular rejection",
                                       amr ~ "Antibody-mediated rejection",
                                       dsa ~ "Donor specific antibodies",
                                       cav ~ "Cardiac allograft vasculopathy",
                                       deceased ~ "Deceased")) %>%
  add_p() %>%
  modify_header(label = "**Demographics by IL6R signaling status**",
                stat_1 = "**Preserved IL6R Signaling (N = 152)**",
                stat_2 = "**Impaired IL6R Signaling (N = 26)**") %>% bold_labels()

# Merge Table 1
tbl_merge(
  tbls = list(vumc.imp.il6r, cuimc.il6r),
  tab_spanner = c("**Vanderbilt University**", "**Columbia University**")
) %>%
  as_gt()

# CAV and IL6R
cav.vumc.logit <- imp2 %>%
  with(glm(cav ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa, family = "binomial"))

cav.vumc <- tbl_regression(cav.vumc.logit,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "Diabetes mellitus",
                                        acr ~ "History of ACR",
                                        amr ~ "History of AMR",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**")

cumc.cav.logit <- glm(cav ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                      data = cumc.table_1, family = "binomial")

cumc.cav <- tbl_regression(cumc.cav.logit,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "Diabetes mellitus",
                                        acr ~ "History of ACR",
                                        amr ~ "History of AMR",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**")

tbl_merge(
  tbls = list(cav.vumc, cumc.cav),
  tab_spanner = c("**Vanderbilt University**", "**Columbia University**")
) %>%
  as_gt()

# ACR
acr.vumc.logit <- imp2 %>%
  with(glm(acr ~ il6r + age_at_htx + ischemic + htn + diabetes + amr + dsa, family = "binomial"))

acr.vumc <- tbl_regression(acr.vumc.logit,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "Diabetes mellitus",
                                        amr ~ "History of AMR",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Acute Cellular Rejection**")

cumc.acr.logit <- glm(acr ~ il6r + age_at_htx + ischemic + htn + diabetes + amr + dsa,
                      data = cumc.table_1, family = "binomial")

acr.cumc <- tbl_regression(cumc.acr.logit,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "Diabetes mellitus",
                                        amr ~ "History of AMR",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Acute Cellular Rejection**")

tbl_merge(
  tbls = list(acr.vumc, acr.cumc),
  tab_spanner = c("**Vanderbilt University**", "**Columbia University**")
) %>%
  as_gt()

# AMR
amr.vumc.logit <- imp2 %>%
  with(glm(amr ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + dsa, family = "binomial"))

amr.vumc <- tbl_regression(amr.vumc.logit,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "Diabetes mellitus",
                                        acr ~ "History of ACR",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Antibody-Mediated Rejection**")

amr.cumc.logit <- glm(amr ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + dsa,
                      data = cumc.table_1, family = "binomial")

amr.cumc <- tbl_regression(amr.cumc.logit,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "Diabetes mellitus",
                                        acr ~ "History of ACR",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Antibody-Mediated Rejection**")

tbl_merge(
  tbls = list(amr.vumc, amr.cumc),
  tab_spanner = c("**Vanderbilt University**", "**Columbia University**")
) %>%
  as_gt()

# Mortality
mort.vumc.logit <- imp2 %>%
  with(glm(deceased ~ il6r + age_at_htx + ischemic + cav + acr + amr + dsa, family = "binomial"))

mort.vumc <- tbl_regression(mort.vumc.logit,
                            exponentiate = TRUE,
                            pvalue_fun = ~style_pvalue(.x, digits = 2),
                            label = list(il6r ~ "IL6R p.Asp358Ala",
                                         age_at_htx ~ "Age at transplant",
                                         ischemic ~ "Ischemic etiology",
                                         cav ~ "Cardiac allograft vasculopathy",
                                         acr ~ "History of ACR",
                                         amr ~ "History of AMR",
                                         dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Mortality**")

mort.cumc.logit <- glm(deceased ~ il6r + age_at_htx + ischemic + cav + acr + amr + dsa,
                       data = cumc.table_1, family = "binomial")

mort.cumc <- tbl_regression(mort.cumc.logit,
                            exponentiate = TRUE,
                            pvalue_fun = ~style_pvalue(.x, digits = 2),
                            label = list(il6r ~ "IL6R p.Asp358Ala",
                                         age_at_htx ~ "Age at transplant",
                                         ischemic ~ "Ischemic etiology",
                                         cav ~ "Cardiac allograft vasculopathy",
                                         acr ~ "History of ACR",
                                         amr ~ "History of AMR",
                                         dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Mortality**")

tbl_merge(
  tbls = list(mort.vumc, mort.cumc),
  tab_spanner = c("**Vanderbilt University**", "**Columbia University**")
) %>%
  as_gt()

### Fixed-effect meta-analyses
power.analysis(OR = 1.36, k = 2, n1 = 606, n2 = 178, p = 0.05, heterogeneity = "fixed")

meta_cav <- read_csv("meta_il6r_cav.csv")
meta_mortality <- read_csv("meta_il6r_mortality.csv")
meta_acr <- read_csv("meta_il6r_acr.csv")
meta_amr <- read_csv("meta_il6r_amr.csv")

# CAV
m.bin.cav <- metabin(Ee, Ne, Ec, Nc,
                     data = meta_cav,
                     studlab = paste(Group),
                     comb.fixed = T, comb.random = F,
                     method = "MH", sm = "OR")
m.bin.cav
forest(m.bin.cav, leftcols = c('studlab'))

# Mortality
m.bin.mortality <- metabin(Ee, Ne, Ec, Nc,
                           data = meta_mortality,
                           studlab = paste(Group),
                           comb.fixed = T, comb.random = F,
                           method = "MH", sm = "OR")
m.bin.mortality
forest(m.bin.mortality, leftcols = c('studlab'))

# ACR
m.bin.acr <- metabin(Ee, Ne, Ec, Nc,
                     data = meta_acr,
                     studlab = paste(Group),
                     comb.fixed = T, comb.random = F,
                     method = "MH", sm = "OR")
m.bin.acr
forest(m.bin.acr, leftcols = c('studlab'))

# AMR
m.bin.amr <- metabin(Ee, Ne, Ec, Nc,
                     data = meta_amr,
                     studlab = paste(Group),
                     comb.fixed = T, comb.random = F,
                     method = "MH", sm = "OR")
m.bin.amr
forest(m.bin.amr, leftcols = c('studlab'))

mace_data <- mat %>% mutate(
  mace = case_when(
    deceased == 1 | acr == 1 | amr == 1 | cav == 1 ~ 1,
    TRUE ~ 0
  )
)

mace.vumc.logit <- glm(mace ~ il6r + age_at_htx + ischemic + htn + diabetes + dsa,
                       data = mace_data, family = "binomial")

mace.vumc.regression <- tbl_regression(mace.vumc.logit,
                                       exponentiate = TRUE,
                                       pvalue_fun = ~style_pvalue(.x, digits = 2),
                                       label = list(il6r ~ "IL6R p.Asp358Ala",
                                                    age_at_htx ~ "Age at transplant",
                                                    ischemic ~ "Ischemic etiology",
                                                    htn ~ "Hypertension",
                                                    diabetes ~ "Diabetes mellitus",
                                                    dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**MACE**")

mace <- imp2 %>%
  with(glm(deceased | cav | acr | amr ~ il6r + age_at_htx + ischemic + htn + diabetes + dsa, family = "binomial"))

mace.reg <- tbl_regression(test,
                           exponentiate = TRUE,
                           pvalue_fun = ~style_pvalue(.x, digits = 2),
                           label = list(il6r ~ "IL6R p.Asp358Ala",
                                        age_at_htx ~ "Age at transplant",
                                        ischemic ~ "Ischemic etiology",
                                        htn ~ "Hypertension",
                                        diabetes ~ "History of diabetes mellitus",
                                        dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**MACE**")

# CUIMC MACE data
cumc.mace_data <- cumc.table_1 %>% mutate(
  mace = case_when(
    deceased == 1 | acr == 1 | amr == 1 | cav == 1 ~ 1,
    TRUE ~ 0
  )
)

cumc.reg.mace <- glm(mace ~ il6r + age_at_htx + ischemic + htn + diabetes + dsa,
                     data = cumc.mace_data, family = "binomial")

cuimc.mace.reg <- tbl_regression(cumc.reg.mace,
                                 exponentiate = TRUE,
                                 pvalue_fun = ~style_pvalue(.x, digits = 2),
                                 label = list(il6r ~ "IL6R p.Asp358Ala",
                                              age_at_htx ~ "Age at transplant",
                                              ischemic ~ "Ischemic etiology",
                                              htn ~ "Hypertension",
                                              diabetes ~ "History of diabetes mellitus",
                                              dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**MACE**")

tbl_merge(
  tbls = list(test.reg, cuimc.mace.reg),
  tab_spanner = c("**Vanderbilt University**", "**Columbia University**")
) %>%
  as_gt()

# MACE meta
meta_mace <- read_csv("meta_il6r_mace2.csv")

# MACE
m.bin.mace <- metabin(Ee, Ne, Ec, Nc,
                      data = meta_mace,
                      studlab = paste(Group),
                      comb.fixed = T, comb.random = F,
                      method = "MH", sm = "OR")
m.bin.mace
forest(m.bin.mace, leftcols = c('studlab'))

# Odds ratios for individual variables

program <- c("Variant", "Normal")
outcome <- c("Event", "No Event")

vumc.diabetes <- matrix(c(44, 27, 270, 268), nrow = 2, ncol = 2,
                        byrow = TRUE)
dimnames(vumc.diabetes) <- list('Program'=program, 'Outcome'=outcome)
vumc.diabetes
oddsratio(vumc.diabetes)

cumc.diabetes <- matrix(c(9, 17, 53, 99), nrow = 2, ncol = 2,
                        byrow = TRUE)
dimnames(cumc.diabetes) <- list('Program'=program, 'Outcome'=outcome)
cumc.diabetes
oddsratio(cumc.diabetes)

vumc.htn <- matrix(c(62, 9, 486, 52), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(vumc.htn) <- list('Program'=program, 'Outcome'=outcome)
vumc.htn
oddsratio(vumc.htn)

cumc.htn <- matrix(c(13, 13, 91, 61), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(cumc.htn) <- list('Program'=program, 'Outcome'=outcome)
cumc.htn
oddsratio(cumc.htn)

vumc.dsa <- matrix(c(10, 61, 122, 416), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(vumc.dsa) <- list('Program'=program, 'Outcome'=outcome)
vumc.dsa
oddsratio(vumc.dsa)

cumc.dsa <- matrix(c(4, 22, 47, 105), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(cumc.dsa) <- list('Program'=program, 'Outcome'=outcome)
cumc.dsa
oddsratio(cumc.dsa)

vumc.acr <- matrix(c(33, 38, 210, 328), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(vumc.acr) <- list('Program'=program, 'Outcome'=outcome)
vumc.acr
oddsratio(vumc.acr)

cumc.acr <- matrix(c(3, 23, 19, 133), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(cumc.acr) <- list('Program'=program, 'Outcome'=outcome)
cumc.acr
oddsratio(cumc.acr)

vumc.amr <- matrix(c(6, 65, 73, 465), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(vumc.amr) <- list('Program'=program, 'Outcome'=outcome)
vumc.amr
oddsratio(vumc.amr)

cumc.amr <- matrix(c(2, 24, 17, 135), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(cumc.amr) <- list('Program'=program, 'Outcome'=outcome)
cumc.amr
oddsratio(cumc.amr)

vumc.cav <- matrix(c(35, 36, 260, 278), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(vumc.cav) <- list('Program'=program, 'Outcome'=outcome)
vumc.cav
oddsratio(vumc.cav)

cumc.cav <- matrix(c(11, 15, 69, 83), nrow = 2, ncol = 2,
                   byrow = TRUE)
dimnames(cumc.cav) <- list('Program'=program, 'Outcome'=outcome)
cumc.cav
oddsratio(cumc.cav)

