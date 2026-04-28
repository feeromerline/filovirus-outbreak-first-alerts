library(dplyr)
library(ggplot2)
library(readr)

trim  <- function(x) gsub("^\\s+|\\s+$", "", x)
lower <- function(x) tolower(trim(x))

normalize_country <- function(x) {
  dplyr::case_when(
    x %in% c("DRC") ~ "Democratic Republic of the Congo",
    x %in% c("RoC", "ROC") ~ "Republic of the Congo",
    x %in% c("Equitorial Guinea", "Equitorial") ~ "Equatorial Guinea",
    TRUE ~ x
  )
}

# Load the authoritative outbreak-level dataset.
# Extractions_R.csv has one row per filovirus outbreak with hand-coded
# analytical variables (Single_or_Cluster, Healthcare_Worker_Infections,
# Days_to_Declaration, Sex, With_Haemorrhage, Exposure_Known).
# The per-source-document public CSV lacks these explicit codings and
# would require imprecise text mining to approximate them.
csv_path <- file.path(path.expand("~"), "Downloads", "Extractions_R.csv")
if (!file.exists(csv_path)) {
  stop("Missing analysis dataset: ", csv_path)
}

df <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(
    Start_Year    = suppressWarnings(as.integer(Start_Year)),
    Virus         = trim(Virus),
    Country       = normalize_country(trim(Country)),
    Report_Source = trim(Report_Source)
  )

####################################
###Reporting Source Analysis###
####################################

library(dplyr)
library(ggplot2)

#Decades
dec_vals   <- sort(unique(floor(as.integer(df$Start_Year) / 10) * 10))
dec_levels <- paste0(dec_vals, "s")

df <- df %>%
  mutate(
    Decade = ifelse(is.na(Start_Year), NA_character_,
                    paste0(floor(as.integer(Start_Year) / 10) * 10, "s")),
    Decade = factor(Decade, levels = dec_levels, ordered = TRUE),
    
    # First alert source categories
    First_Alert_Source = case_when(
      grepl(" and ", Report_Source, ignore.case = TRUE) ~ "Mixed",
      grepl("National Epidemic Alert|community|alert", Report_Source, ignore.case = TRUE) ~ "Community/CBS",
      grepl("physician|doctor|nurse", Report_Source, ignore.case = TRUE) ~ "Physician/Hospital Staff",
      grepl("laboratory|\\bMoH\\b", Report_Source, ignore.case = TRUE) ~ "Formal Indicator",
      TRUE ~ "Not Reported"
    )
  )

head(df)

#N per decade
q_df <- df %>%
  filter(!is.na(Decade)) %>%
  count(Decade) %>%
  rename(n_total = n)

q_df

#Counts
tab <- df %>%
  filter(!is.na(Decade)) %>%
  count(Decade, First_Alert_Source) %>%
  group_by(Decade) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

#Barchart
ggplot(tab, aes(x = Decade, y = pct, fill = First_Alert_Source)) +
  geom_col(color = "white", width = 0.7) +
  geom_text(
    aes(label = paste0(round(pct, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  # add n
  geom_text(
    data = q_df,
    aes(x = Decade, y = 105, label = paste0("n=", n_total)),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  coord_cartesian(ylim = c(0, 106)) +
  labs(x = "Decade", y = "Percent of Outbreaks", fill = "First alert source") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

#Table
tab_src_dec <- with(df, table(Decade, First_Alert_Source, useNA = "no"))
tab_src_dec
round(100 * prop.table(tab_src_dec, 1), 1)  # row-wise %

#Fisher's exact
if (all(dim(tab_src_dec) > 1)) {
  f_src_dec <- fisher.test(tab_src_dec)
  cat(sprintf("\nFisher's exact p = %.4f\n", f_src_dec$p.value))
}


####################################
###Clusters vs Single Index Cases###
####################################

df$Virus <- trim(df$Virus)
df$Single_or_Cluster <- trim(df$Single_or_Cluster)
df$Days_to_declaration <- trim(df$Days_to_Declaration)
df$Start_Year <- as.integer(df$Start_Year)

#Early Signals
signal <- rep(NA_character_, nrow(df))
signal[grepl("single", df$Single_or_Cluster, ignore.case = TRUE)] <- "Single"
signal[grepl("cluster", df$Single_or_Cluster, ignore.case = TRUE)] <- "Cluster"
df$SignalType <- factor(signal, levels = c("Single", "Cluster"))

# Viruses
virus_group <- rep(NA_character_, nrow(df))
virus_group[grepl("bundib", df$Virus, ignore.case = TRUE)] <- "BDBV"
virus_group[grepl("sudan", df$Virus, ignore.case = TRUE)] <- "SUDV"
virus_group[grepl("zaire", df$Virus, ignore.case = TRUE)] <- "EBOV"
virus_group[grepl("marburg", df$Virus, ignore.case = TRUE)] <- "MARV"
df$VirusGroup <- factor(virus_group,
                        levels = c("EBOV", "SUDV", "BDBV", "MARV")
)

#Counts for table
tab_signal <- table(df$SignalType, useNA = "no")
tab_signal

#Proportion calc
prop_signal <- prop.table(tab_signal)
round(100 * prop_signal, 1)

#Virus calc
tab_signal_virus <- with(df, table(VirusGroup, SignalType, useNA = "no"))
tab_signal_virus
round(100 * prop.table(tab_signal_virus, 1), 1)

#Summary
props <- df |>
  dplyr::filter(!is.na(VirusGroup), !is.na(SignalType)) |>
  dplyr::count(VirusGroup, SignalType, name = "n") |>
  dplyr::group_by(VirusGroup) |>
  dplyr::mutate(pct = n / sum(n),
                lbl = paste0(round(pct * 100), "%"))

#N per virus group
n_df <- props |>
  dplyr::summarise(n = sum(n), .groups = "drop")

#Plotting
ggplot(props, aes(x = VirusGroup, y = pct, fill = SignalType)) +
  geom_col(position = "fill", color = "white", width = 0.7) +
  # % labels inside stacks
  geom_text(aes(label = lbl),
            position = position_fill(vjust = 0.5),
            color = "white", size = 3.2) +
  # n labels above bars
  geom_text(data = n_df,
            aes(x = VirusGroup, y = 1.02, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3.2) +
  coord_cartesian(ylim = c(0, 1.06)) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Virus",
       y = "Percent of Index Cases",
    fill = "First signal") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 25, hjust = 1))

#Counts
tab <- with(df, table(VirusGroup, SignalType, useNA = "no"))
tab
round(100 * prop.table(tab, 1), 1)

#Association test
fisher.test(tab)

################################
###Time to Detection Analyses###
################################

#Establish iqr
iqr_bounds <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  list(p25 = unname(q[1]), med = unname(q[2]), p75 = unname(q[3]))
}

#Re-parse for clarity
parse_days <- function(x) {
  x <- lower(x)
  if (x %in% c("", "nr", "not declared", "na", "n/a")) return(c(NA_real_, FALSE))
  if (grepl("^>\\s*\\d+", x)) {
    val <- as.numeric(gsub("[^0-9]", "", x))
    return(c(val, TRUE))  # lower bound
  }
  if (grepl("^[0-9]+$", x)) return(c(as.numeric(x), FALSE))
  return(c(NA_real_, FALSE))
}

#Time to declaration
if (!("TimeDecl_days" %in% names(df))) {
  parsed <- t(vapply(df$Days_to_Declaration, parse_days, c(NA_real_, FALSE)))
  df$TimeDecl_days     <- as.numeric(parsed[, 1])
  df$TimeDecl_censored <- as.logical(parsed[, 2])
}

#Signal
if (!("SignalType" %in% names(df))) {
  sig_raw <- lower(trim(df$Single_or_Cluster))
  df$SignalType <- NA_character_
  df$SignalType[grepl("single", sig_raw)]  <- "Single"
  df$SignalType[grepl("cluster", sig_raw)] <- "Cluster"
  df$SignalType <- factor(df$SignalType, levels = c("Single","Cluster"))
}

# Virus
if (!("VirusGroup" %in% names(df))) {
  v <- lower(trim(df$Virus))
  df$VirusGroup <- dplyr::case_when(
    grepl("bundib|bdbv", v)                ~ "BDBV",
    grepl("sudan|sudv",  v)                ~ "SUDV",
    grepl("zaire",       v)                ~ "EBOV",
    grepl("\\bebov\\b|ebola( virus disease)?$", v) ~ "EBOV",
    grepl("marburg|marv",v)                ~ "MARV",
    TRUE ~ NA_character_
  )
  df$VirusGroup <- factor(df$VirusGroup, levels = c("EBOV","SUDV","BDBV","MARV"))
}

#HCW infections
col_hcw <- names(df)[grepl("health.*worker.*infect", tolower(names(df)))]
if (length(col_hcw) != 1) {
  stop("Couldn't uniquely identify the HCW infections column. Found: ",
       paste(col_hcw, collapse = ", "))
}

hcw_raw <- lower(trim(df[[col_hcw]]))
df$HCW_YN <- dplyr::case_when(
  hcw_raw %in% c("y","yes","1") ~ "Yes",
  hcw_raw %in% c("n","no","0")  ~ "No",
  hcw_raw %in% c("", "nr", "not reported", "na", "n/a") ~ NA_character_,
  TRUE ~ NA_character_
)
df$HCW_YN <- factor(df$HCW_YN, levels = c("No","Yes"))

#Time to declaration
sub_hcw <- df %>%
  filter(!is.na(TimeDecl_days), !is.na(HCW_YN))

n_hcw <- sub_hcw %>% count(HCW_YN, name = "n")
stats_hcw <- sub_hcw %>%
  group_by(HCW_YN) %>%
  summarise(
    n   = n(),
    p25 = quantile(TimeDecl_days, 0.25),
    med = median(TimeDecl_days),
    p75 = quantile(TimeDecl_days, 0.75),
    .groups = "drop"
  )

test_hcw <- wilcox.test(TimeDecl_days ~ HCW_YN, data = sub_hcw, exact = FALSE)

cat("\n--- Time to declaration by early HCW infections ---\n")
print(stats_hcw)
cat(sprintf("Mann–Whitney (Wilcoxon rank-sum) p = %.4f\n", test_hcw$p.value))


#Figure 5: Boxplot HCW infections vs Time to Declaration
x_labs_hcw <- sub_hcw %>% count(HCW_YN) %>%
  mutate(lab = paste0(HCW_YN, " (n=", n, ")")) %>%
  as.data.frame()

p_hcw <- ggplot(sub_hcw, aes(HCW_YN, TimeDecl_days, fill = HCW_YN)) +
  geom_boxplot(width = 0.65, alpha = 0.75, outlier.alpha = 0.3, color = "black") +
  geom_jitter(aes(shape = VirusGroup), width = 0.15, height = 0, size = 2.2,
              alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("No"="#afd1e7","Yes"="#f6f364"), guide = "none") +
  scale_shape_manual(values = c(21,24,22,23), name = "Virus") +
  scale_x_discrete(labels = setNames(x_labs_hcw$lab, x_labs_hcw$HCW_YN)) +
  labs(x = "Early HCW infections",
       y = "Days from earliest report to official declaration") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "bottom")
p_hcw
ggsave("fig6_time_by_HCW.png", p_hcw, width = 6, height = 4, dpi = 300)

#Figure 6: Cluster or Single vs Time to Declaration
sub_sig <- df %>%
  filter(!is.na(TimeDecl_days), !is.na(SignalType))

stats_sig <- sub_sig %>%
  group_by(SignalType) %>%
  summarise(
    n   = n(),
    p25 = quantile(TimeDecl_days, 0.25),
    med = median(TimeDecl_days),
    p75 = quantile(TimeDecl_days, 0.75),
    .groups = "drop"
  )

test_sig <- wilcox.test(TimeDecl_days ~ SignalType, data = sub_sig, exact = FALSE)

cat("\n--- Time to declaration: Single vs Cluster ---\n")
print(stats_sig)
cat(sprintf("Mann–Whitney (Wilcoxon rank-sum) p = %.4f\n", test_sig$p.value))

#Boxplot 
x_labs_sig <- sub_sig %>% count(SignalType) %>%
  mutate(lab = paste0(SignalType, " (n=", n, ")")) %>%
  as.data.frame()

p_sig <- ggplot(sub_sig, aes(SignalType, TimeDecl_days, fill = SignalType)) +
  geom_boxplot(width = 0.65, alpha = 0.75, outlier.alpha = 0.3, color = "black") +
  geom_jitter(aes(shape = VirusGroup), width = 0.15, height = 0, size = 2.2,
              alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("Single"="#b2df8a","Cluster"="#08306b"), guide = "none") +
  scale_shape_manual(values = c(21,24,22,23), name = "Virus") +
  scale_x_discrete(labels = setNames(x_labs_sig$lab, x_labs_sig$SignalType)) +
  labs(x = "First signal type",
       y = "Days from earliest report to official declaration") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "bottom")
p_sig
ggsave("fig7_time_by_signal.png", p_sig, width = 6, height = 4, dpi = 300)

fmt <- function(med, p25, p75) sprintf("%g (IQR %g–%g)", med, p25, p75)



######################################
###Time to Detection Across Decades###
######################################

#Subset
sub_dec <- df %>%
  filter(!is.na(TimeDecl_days), !is.na(Decade))

#N per decade
x_labs_dec <- sub_dec %>%
  count(Decade) %>%
  mutate(lab = paste0(as.character(Decade), " (n=", n, ")")) %>%
  as.data.frame()

#Style
shape_map <- c("EBOV" = 21, "SUDV" = 24, "BDBV" = 22, "MARV" = 23)

#Colors
dec_levels <- levels(sub_dec$Decade)
color_map  <- setNames(
  colorRampPalette(c("#C1E1C1", "#b2df8a", "#6B8E23", "#228B22"))(length(dec_levels)),
  dec_levels
)

p_dec <- ggplot(sub_dec, aes(Decade, TimeDecl_days, fill = Decade)) +
  geom_boxplot(width = 0.65, alpha = 0.85, outlier.alpha = 0.3, color = "black") +
  #shapes by virus, black outline, no color fill
  geom_jitter(aes(shape = VirusGroup),
              width = 0.12, height = 0, size = 2.2,
              alpha = 0.7, color = "black", show.legend = TRUE) +
  scale_fill_manual(values = color_map, guide = "none") +
  scale_shape_manual(values = shape_map, name = "Virus") +
  scale_x_discrete(labels = setNames(x_labs_dec$lab, x_labs_dec$Decade)) +
  labs(x = "Decade",
       y = "Days from earliest report to official declaration") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "bottom")

p_dec

#Medians & IQR
stats_dec <- sub_dec %>%
  group_by(Decade) %>%
  summarise(
    n   = n(),
    p25 = quantile(TimeDecl_days, 0.25),
    med = median(TimeDecl_days),
    p75 = quantile(TimeDecl_days, 0.75),
    .groups = "drop"
  )
cat("\n--- Time to declaration by decade ---\n")
print(stats_dec)

#Kruskal–Wallis
kw_dec <- kruskal.test(TimeDecl_days ~ Decade, data = sub_dec)
cat(sprintf("Kruskal–Wallis: chi^2 = %.3f, df = %d, p = %.4f\n",
            kw_dec$statistic, kw_dec$parameter, kw_dec$p.value))

####Overall trend graphs#####

#Cleanup
df$Virus <- trim(df$Virus)

df$VirusGroup <- dplyr::case_when(
  grepl("bundib|bdbv", tolower(df$Virus)) ~ "BDBV",
  grepl("sudan|sudv",  tolower(df$Virus)) ~ "SUDV",
  grepl("zaire",       tolower(df$Virus)) ~ "EBOV",         
  grepl("\\bebov\\b|ebola( virus disease)?$", tolower(df$Virus)) ~ "EBOV",
  grepl("marburg|marv",tolower(df$Virus)) ~ "MARV",
  TRUE ~ NA_character_
)

df$VirusGroup <- factor(df$VirusGroup, levels = c("EBOV","SUDV","BDBV","MARV"))


#Parse days
parse_days <- function(x) {
  x <- tolower(trim(x))
  if (x %in% c("", "nr", "not declared", "na", "n/a")) return(c(NA_real_, FALSE))
  if (grepl("^>\\s*\\d+", x)) {
    val <- as.numeric(gsub("[^0-9]", "", x))
    return(c(val, TRUE))
  }
  if (grepl("^[0-9]+$", x)) return(c(as.numeric(x), FALSE))
  return(c(NA_real_, FALSE))
}

parsed <- t(vapply(df$Days_to_Declaration, parse_days, c(NA_real_, FALSE)))
df$TimeDecl_days     <- as.numeric(parsed[, 1])
df$TimeDecl_censored <- as.logical(parsed[, 2])

#Cleanup
df$Country    <- trim(df$Country)
df$Virus      <- trim(df$Virus)
df$VirusGroup <- factor(df$VirusGroup,
                        levels = c("EBOV", "SUDV", "BDBV", "MARV")
)

#Plotting
plotdf <- df %>%
  dplyr::filter(!is.na(Start_Year),
                !is.na(TimeDecl_days),
                !is.na(VirusGroup),
                !is.na(Country))

plotdf$Country <- dplyr::recode(plotdf$Country,
                                "The Republic of the Congo" = "Republic of the Congo",
                                "Guinea/West Africa"        = "Guinea"
)

shape_map <- c("EBOV"=16, "SUDV"=17, "BDBV"=15, "MARV"=18)
color_map <- c("EBOV"="#afd1e7", "SUDV"="#08306b",
               "BDBV"="#f7fbff", "MARV"="#f6f364")

#Country limiter
keepers <- plotdf %>% dplyr::count(Country) %>% dplyr::filter(n >= 2) %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::pull(Country)
facetdf <- plotdf %>% dplyr::filter(Country %in% keepers)
facetdf$Country <- factor(facetdf$Country, levels = keepers)

topN <- 2
top_countries <- plotdf %>% count(Country, sort = TRUE) %>% slice_head(n = topN) %>% pull(Country)
facetdf2 <- plotdf %>% filter(Country %in% top_countries)

p_country <- ggplot(facetdf2, aes(Start_Year, TimeDecl_days)) +
  # overall (all viruses) per-country trend
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "grey55", linetype = "dashed", linewidth = 0.7) + 
  # black outline for points (slightly larger, underneath)
  geom_point(
    aes(shape = VirusGroup),
    color = "black",
    size = 3.5,
    alpha = 0.9
  ) +
  # virus-specific points and trends
  geom_point(aes(shape = VirusGroup, color = VirusGroup),
             size = 3, alpha = 0.9) +
  scale_shape_manual(values = shape_map, name = "Virus") +
  scale_color_manual(values = color_map, name = "Virus") +
  labs(x = "Outbreak start year",
       y = "Days from earliest report to official declaration") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ Country, ncol = 4) 
p_country


#Shapes
shape_map <- c("EBOV" = 21,   # circle
               "SUDV" = 24,   # triangle
               "BDBV" = 22,   # square
               "MARV" = 23)   # diamond

color_map <- c("EBOV"= "#afd1e7", "SUDV"="#08306b",
               "BDBV"="#f7fbff", "MARV"="#f6f364")

p_overall <- ggplot(plotdf, aes(Start_Year, TimeDecl_days)) +
  # black outline for points
  geom_point(
    aes(shape = VirusGroup, fill = VirusGroup),
    color = "black",   # outline
    size = 3, alpha = 0.9) +
  # global solid trend (all outbreaks)
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  # colored dashed lines by virus
  geom_smooth(
    aes(color = VirusGroup, group = VirusGroup),
    method = "lm", se = FALSE,
    linetype = "dashed", linewidth = 0.9) +
  # manual scales (only once each!)
  scale_shape_manual(values = shape_map, name = "Virus") +
  scale_fill_manual(values  = color_map, name = "Virus") +
  scale_color_manual(values = color_map, name = "Virus") +
  labs(
    x = "Outbreak start year",
    y = "Days from earliest report to official declaration") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom")

p_overall

#########################
###Additional Analyses###
#########################
trim <- function(x) gsub("^\\s+|\\s+$", "", x)
lower <- function(x) tolower(trim(x))

col_signal <- if ("Single_or_Cluster" %in% names(df)) "Single_or_Cluster" else names(df)[grepl("single|cluster", tolower(names(df)))]
col_sex    <- if ("Sex" %in% names(df)) "Sex" else names(df)[grepl("^sex$", tolower(names(df)))]
col_hemo   <- if ("With_haemorrhage" %in% names(df)) "With_haemorrhage" else names(df)[grepl("haem|hemorr", tolower(names(df)))]
col_expo   <- if ("Exposure_known" %in% names(df)) "Exposure_known" else names(df)[grepl("exposure", tolower(names(df)))]

if (length(col_signal) != 1) stop("Couldn't uniquely identify 'Single or Cluster' column")
if (length(col_sex)    != 1) stop("Couldn't uniquely identify 'Sex' column")
if (length(col_hemo)   != 1) stop("Couldn't uniquely identify 'With haemorrhage' column")
if (length(col_expo)   != 1) stop("Couldn't uniquely identify 'Exposure known' column")

#Cleanup
df$Signal_raw <- lower(df[[col_signal]])
df$Sex_raw    <- lower(df[[col_sex]])
df$Hemo_raw   <- lower(df[[col_hemo]])
df$Expo_raw   <- lower(df[[col_expo]])

#Standardize
df$SignalType <- dplyr::case_when(
  grepl("single",  df$Signal_raw)  ~ "Single",
  grepl("cluster", df$Signal_raw)  ~ "Cluster",
  df$Signal_raw %in% c("", "NR", "Not Reported", "NA", "N/A") ~ "NR",
  TRUE ~ "NR"
)

# 1) Index case sex 
df$Sex_cat <- dplyr::case_when(
  df$SignalType != "Single"                              ~ "Not applicable",
  df$Sex_raw %in% c("", "NR", "not reported")            ~ "Not reported",
  grepl("^m", df$Sex_raw)                                ~ "Male",
  grepl("^f", df$Sex_raw)                                ~ "Female",
  df$Sex_raw %in% c("na", "N/A")                         ~ "Not applicable",
  TRUE                                                   ~ "Other/unspecified"
)

sex_counts <- df %>%
  count(Sex_cat, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(match(Sex_cat,
                c("Male", "Female", "Other/unspecified",
                  "Not reported", "Not applicable")))
print(sex_counts)

#By Virus
 df$VirusGroup <- dplyr::case_when(
   grepl("bundib", tolower(df$Virus)) ~ "EBOV-BDBV",
   grepl("sudan",  tolower(df$Virus)) ~ "EBOV-SUDV",
   grepl("zaire",  tolower(df$Virus)) ~ "EBOV-Zaire",
   grepl("marburg",tolower(df$Virus)) ~ "MARV",
   TRUE ~ "Other/NR"
 )
 
 df %>% count(VirusGroup, Sex_cat) %>%
   group_by(VirusGroup) %>%
   mutate(pct = round(100 * n / sum(n), 1)) %>%
   arrange(VirusGroup, desc(n)) %>% print(n = 100)

#Exposure known
df$Expo_cat <- dplyr::case_when(
  df$Expo_raw %in% c("", "nr", "not reported", "na", "n/a")            ~ "Not reported or Cluster with unknown exposure",
  grepl("^y", df$Expo_raw)                                ~ "Yes",
  grepl("^n", df$Expo_raw)                                ~ "No",
  TRUE                                                    ~ "Not reported or Cluster with unknown exposure"
)

expo_counts <- df %>%
  count(Expo_cat, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(match(Expo_cat,
                c("Yes", "No", "Not reported or Cluster with unknown exposure")))

print(expo_counts)


# 3) With haemorrhage?
df$Hemo_cat <- dplyr::case_when(
  df$Hemo_raw %in% c("", "nr", "not reported")             ~ "Not reported",
  grepl("^y|^yes", df$Hemo_raw)                            ~ "Yes",
  grepl("^n|^no",  df$Hemo_raw)                            ~ "No",
  df$Hemo_raw %in% c("na", "n/a")                          ~ "Not applicable",
  TRUE                                                     ~ "Not reported"
)

hemo_counts <- df %>%
  count(Hemo_cat, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(match(Hemo_cat,
                c("Yes", "No", "Not reported", "Not applicable")))

print(hemo_counts)