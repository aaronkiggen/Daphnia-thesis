library(qqman)
library(dplyr)
library(ggplot2)

##### MARIA
#w15000_s2500

df1 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_1.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df2 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_2.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df3 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_3.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df4 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_4.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df5 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_5.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df6 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_6.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df7 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_7.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df8 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_8.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df9 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_9.stats_B_M_w15000_s2500_eggStats.stats",sep = ",", header = TRUE)
df10 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_10.stats_B_M_w15000_s2500_eggStats.stats", sep = ",",header = TRUE)
df11 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_11.stats_B_M_w15000_s2500_eggStats.stats", sep = ",",header = TRUE)
df12 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_12.stats_B_M_w15000_s2500_eggStats.stats", sep = ",",header = TRUE)
df13 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_13.stats_B_M_w15000_s2500_eggStats.stats", sep = ",",header = TRUE)


# Convert all columns to numeric before combining
convert_numeric <- function(df) {
  df %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
}

df1  <- convert_numeric(df1)
summary(df1) 
sapply(df1, function(x) sum(is.na(x)))  # Count NAs per column

df2  <- convert_numeric(df2)
df3  <- convert_numeric(df3)
df4  <- convert_numeric(df4)
df5  <- convert_numeric(df5)
df6  <- convert_numeric(df6)
df7  <- convert_numeric(df7)
df8  <- convert_numeric(df8)
df9  <- convert_numeric(df9)
df10 <- convert_numeric(df10)

# Combine all into one dataframe
all_stats_w15000_s2500 <- bind_rows(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)


fst_avg_per_chr <- all_stats_w15000_s2500 %>%
  group_by(scaffold) %>%
  summarise(mean_fst = mean(Fst_B_M, na.rm = TRUE)) %>%
  arrange(scaffold)

print(fst_avg_per_chr)
# ── (1) Compute chromosome/scaffold lengths and cumulative starts ─────────────
chr_info <- all_stats_w15000_s2500 %>%
  group_by(scaffold) %>%
  summarise(chr_len = max(end)) %>%
  # put them in the correct order (1,2,…,9,10,scaffold_11,12,13)
  arrange(as.integer(gsub("chr_scaffold_([0-9]+)|chr_([0-9]+)", "\\1\\2", scaffold))) %>%
  mutate(
    cum_start = lag(cumsum(chr_len), default = 0),
    center    = cum_start + chr_len/2
  )

# ── (2) Add cumulative midpoint to your main data ──────────────────────────────
all_stats_w15000_s2500 <- all_stats_w15000_s2500 %>%
  left_join(
    chr_info %>% select(scaffold, cum_start),
    by = "scaffold"
  ) %>%
  mutate(
    cum_mid = mid + cum_start
  )

# ── (3) (Optional) drop negative FST and flag outliers ≥0.15 ────────────────────
# 1) Zero‐floor negatives

str(all_stats_w15000_s2500)

all_stats_w15000_s2500 <- all_stats_w15000_s2500 %>%
  mutate(FstWC_B_M_egg = pmax(Fst_B_M, 0))

# 2) Compute 95th percentile cutoff
fst_thr <- quantile(all_stats_w15000_s2500$FstWC_B_M_egg, 0.99, na.rm = TRUE)

# 3) Flag outliers
all_stats_w15000_s2500 <- all_stats_w15000_s2500 %>%
  mutate(is_outlier = FstWC_B_M_egg >= fst_thr)

# Check the threshold
print(fst_thr)

# ── (4) Plot using cum_mid, with vertical separators & labels ──────────────────
ggplot(all_stats_w15000_s2500, aes(x = cum_mid, y = FstWC_B_M_egg)) +
  # vertical lines at the start of each chr
  geom_vline(
    data = chr_info,
    aes(xintercept = cum_start),
    color = "grey80", size = 0.3
  ) +
  # the points
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  # x axis: ticks & labels at the center of each chr block
  scale_x_continuous(
    breaks = chr_info$center,
    labels = chr_info$scaffold,
    expand = c(0, 0)
  ) +
  labs(
    x        = "Chromosome / Scaffold",
    y        = expression(F[ST]),
    title    = "Genome-wide FST Manhattan Plot w15000_s2500",
    subtitle = "Windows with top 5% FST highlighted in red"
  ) +
  theme_minimal() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )












# Subset to scaffold 1
scaf1_stats <- all_stats_w15000_s2500 %>%
  filter(scaffold == "1")    # keep only rows where chr equals "chr_1" :contentReference[oaicite:0]{index=0}

# Then plot using scaf1_stats
ggplot(scaf1_stats, aes(x = mid, y = FstWC_B_M_egg)) +
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  geom_hline(yintercept = fst_thr, linetype = "dashed", color = "red") +
  labs(
    x = "Position on Scaffold 1 (bp)",
    y = expression(F[ST]),
    title = "FST Along Scaffold 1 w15000_s2500",
    subtitle = sprintf("Outliers ≥ 95ᵗʰ percentile (≥ %.3f)", fst_thr)
  ) +
  theme_minimal() +
  theme(legend.position = "none")








##### MARIA
#w20000_s5000

df1 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_1.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df2 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_2.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df3 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_3.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df4 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_4.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df5 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_5.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df6 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_6.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df7 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_7.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df8 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_8.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df9 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_9.stats_B_M_w20000_s5000_eggStats.stats",sep = ",", header = TRUE)
df10 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_10.stats_B_M_w20000_s5000_eggStats.stats", sep = ",",header = TRUE)
df11 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_11.stats_B_M_w20000_s5000_eggStats.stats", sep = ",",header = TRUE)
df12 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_12.stats_B_M_w20000_s5000_eggStats.stats", sep = ",",header = TRUE)
df13 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_13.stats_B_M_w20000_s5000_eggStats.stats", sep = ",",header = TRUE)


# Convert all columns to numeric before combining
convert_numeric <- function(df) {
  df %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
}

df1  <- convert_numeric(df1)
summary(df1) 
sapply(df1, function(x) sum(is.na(x)))  # Count NAs per column

df2  <- convert_numeric(df2)
df3  <- convert_numeric(df3)
df4  <- convert_numeric(df4)
df5  <- convert_numeric(df5)
df6  <- convert_numeric(df6)
df7  <- convert_numeric(df7)
df8  <- convert_numeric(df8)
df9  <- convert_numeric(df9)
df10 <- convert_numeric(df10)

# Combine all into one dataframe
all_stats_w20000_s5000 <- bind_rows(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)


fst_avg_per_chr <- all_stats_w20000_s5000 %>%
  group_by(scaffold) %>%
  summarise(mean_fst = mean(Fst_B_M, na.rm = TRUE)) %>%
  arrange(scaffold)

print(fst_avg_per_chr)
# ── (1) Compute chromosome/scaffold lengths and cumulative starts ─────────────
chr_info <- all_stats_w20000_s5000 %>%
  group_by(scaffold) %>%
  summarise(chr_len = max(end)) %>%
  # put them in the correct order (1,2,…,9,10,scaffold_11,12,13)
  arrange(as.integer(gsub("chr_scaffold_([0-9]+)|chr_([0-9]+)", "\\1\\2", scaffold))) %>%
  mutate(
    cum_start = lag(cumsum(chr_len), default = 0),
    center    = cum_start + chr_len/2
  )

# ── (2) Add cumulative midpoint to your main data ──────────────────────────────
all_stats_w20000_s5000 <- all_stats_w20000_s5000 %>%
  left_join(
    chr_info %>% select(scaffold, cum_start),
    by = "scaffold"
  ) %>%
  mutate(
    cum_mid = mid + cum_start
  )

# ── (3) (Optional) drop negative FST and flag outliers ≥0.15 ────────────────────
# 1) Zero‐floor negatives

str(all_stats_w20000_s5000)

all_stats_w20000_s5000 <- all_stats_w20000_s5000 %>%
  mutate(FstWC_B_M_egg = pmax(Fst_B_M, 0))

# 2) Compute 95th percentile cutoff
fst_thr <- quantile(all_stats_w20000_s5000$FstWC_B_M_egg, 0.99, na.rm = TRUE)

# 3) Flag outliers
all_stats_w20000_s5000 <- all_stats_w20000_s5000 %>%
  mutate(is_outlier = FstWC_B_M_egg >= fst_thr)

# Check the threshold
print(fst_thr)

# ── (4) Plot using cum_mid, with vertical separators & labels ──────────────────
ggplot(all_stats_w20000_s5000, aes(x = cum_mid, y = FstWC_B_M_egg)) +
  # vertical lines at the start of each chr
  geom_vline(
    data = chr_info,
    aes(xintercept = cum_start),
    color = "grey80", size = 0.3
  ) +
  # the points
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  # threshold line
  # x axis: ticks & labels at the center of each chr block
  scale_x_continuous(
    breaks = chr_info$center,
    labels = chr_info$scaffold,
    expand = c(0, 0)
  ) +
  labs(
    x        = "Chromosome / Scaffold",
    y        = expression(F[ST]),
    title    = "Genome-wide FST Manhattan Plot w20000_s5000",
    subtitle = "Windows with top 5% FST  highlighted in red"
  ) +
  theme_minimal() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )












# Subset to scaffold 1
scaf1_stats <- all_stats %>%
  filter(scaffold == "1")    # keep only rows where chr equals "chr_1" :contentReference[oaicite:0]{index=0}

# Then plot using scaf1_stats
ggplot(scaf1_stats, aes(x = mid, y = FstWC_B_M_egg)) +
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  geom_hline(yintercept = fst_thr, linetype = "dashed", color = "red") +
  labs(
    x = "Position on Scaffold 1 (p) ",
    y = expression(F[ST]),
    title = "FST Along Scaffold 1 w20000_s5000",
    subtitle = sprintf("Outliers ≥ 95ᵗʰ percentile (≥ %.3f)", fst_thr)
  ) +
  theme_minimal() +
  theme(legend.position = "none")



##### MARIA
#w50000_s25000

df1 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_1.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df2 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_2.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df3 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_3.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df4 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_4.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df5 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_5.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df6 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_6.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df7 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_7.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df8 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_8.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df9 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_9.stats_B_M_w50000_s25000_eggStats.stats",sep = ",", header = TRUE)
df10 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_10.stats_B_M_w50000_s25000_eggStats.stats", sep = ",",header = TRUE)
df11 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_11.stats_B_M_w50000_s25000_eggStats.stats", sep = ",",header = TRUE)
df12 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_12.stats_B_M_w50000_s25000_eggStats.stats", sep = ",",header = TRUE)
df13 <- read.table("C:/Users/u0173289/Desktop/FST/final/chaturvedi_chr_scaffold_13.stats_B_M_w50000_s25000_eggStats.stats", sep = ",",header = TRUE)


# Convert all columns to numeric before combining
convert_numeric <- function(df) {
  df %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
}

df1  <- convert_numeric(df1)
summary(df1) 
sapply(df1, function(x) sum(is.na(x)))  # Count NAs per column

df2  <- convert_numeric(df2)
df3  <- convert_numeric(df3)
df4  <- convert_numeric(df4)
df5  <- convert_numeric(df5)
df6  <- convert_numeric(df6)
df7  <- convert_numeric(df7)
df8  <- convert_numeric(df8)
df9  <- convert_numeric(df9)
df10 <- convert_numeric(df10)

# Combine all into one dataframe
all_stats_w50000_s25000 <- bind_rows(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)


fst_avg_per_chr <- all_stats_w50000_s25000 %>%
  group_by(scaffold) %>%
  summarise(mean_fst = mean(Fst_B_M, na.rm = TRUE)) %>%
  arrange(scaffold)

print(fst_avg_per_chr)
# ── (1) Compute chromosome/scaffold lengths and cumulative starts ─────────────
chr_info <- all_stats_w50000_s25000 %>%
  group_by(scaffold) %>%
  summarise(chr_len = max(end)) %>%
  # put them in the correct order (1,2,…,9,10,scaffold_11,12,13)
  arrange(as.integer(gsub("chr_scaffold_([0-9]+)|chr_([0-9]+)", "\\1\\2", scaffold))) %>%
  mutate(
    cum_start = lag(cumsum(chr_len), default = 0),
    center    = cum_start + chr_len/2
  )

# ── (2) Add cumulative midpoint to your main data ──────────────────────────────
all_stats_w50000_s25000 <- all_stats_w50000_s25000 %>%
  left_join(
    chr_info %>% select(scaffold, cum_start),
    by = "scaffold"
  ) %>%
  mutate(
    cum_mid = mid + cum_start
  )

# ── (3) (Optional) drop negative FST and flag outliers ≥0.15 ────────────────────
# 1) Zero‐floor negatives

str(all_stats_w50000_s25000)

all_stats_w50000_s25000 <- all_stats_w50000_s25000 %>%
  mutate(FstWC_B_M_egg = pmax(Fst_B_M, 0))

# 2) Compute 95th percentile cutoff
fst_thr <- quantile(all_stats_w50000_s25000$FstWC_B_M_egg, 0.95, na.rm = TRUE)

# 3) Flag outliers
all_stats_w50000_s25000 <- all_stats_w50000_s25000 %>%
  mutate(is_outlier = FstWC_B_M_egg >= fst_thr)

# Check the threshold
print(fst_thr)

# ── (4) Plot using cum_mid, with vertical separators & labels ──────────────────
ggplot(all_stats_w50000_s25000, aes(x = cum_mid, y = FstWC_B_M_egg)) +
  # vertical lines at the start of each chr
  geom_vline(
    data = chr_info,
    aes(xintercept = cum_start),
    color = "grey80", size = 0.3
  ) +
  # the points
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  # threshold line
  geom_hline(yintercept = 0.15, linetype = "dashed", color = "red") +
  # x axis: ticks & labels at the center of each chr block
  scale_x_continuous(
    breaks = chr_info$center,
    labels = chr_info$scaffold,
    expand = c(0, 0)
  ) +
  labs(
    x        = "Chromosome / Scaffold",
    y        = expression(F[ST]),
    title    = "Genome-wide FST Manhattan Plot  w50000_s25000",
    subtitle = "Windows with top 5% FST  highlighted in red"
  ) +
  theme_minimal() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )












# Subset to scaffold 1
scaf1_stats <- all_stats %>%
  filter(scaffold == "1")    # keep only rows where chr equals "chr_1" :contentReference[oaicite:0]{index=0}

# Then plot using scaf1_stats
ggplot(scaf1_stats, aes(x = mid, y = FstWC_B_M_egg)) +
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  geom_hline(yintercept = fst_thr, linetype = "dashed", color = "red") +
  labs(
    x = "Position on Scaffold 1 (p) ",
    y = expression(F[ST]),
    title = "FST Along Scaffold 1 w50000_s25000",
    subtitle = sprintf("Outliers ≥ 95ᵗʰ percentile (≥ %.3f)", fst_thr)
  ) +
  theme_minimal() +
  theme(legend.position = "none")




str(all_stats_w15000_s2500)








##################### FST LOCATIONS
out15 <- all_stats_w15000_s2500 %>% filter(is_outlier) %>% select(scaffold, start, end,is_outlier,FstWC_B_M_egg)
out20 <- all_stats_w20000_s5000 %>% filter(is_outlier) %>% select(scaffold, start, end,is_outlier,FstWC_B_M_egg)
out50 <- all_stats_w50000_s25000 %>% filter(is_outlier) %>% select(scaffold, start, end,is_outlier,FstWC_B_M_egg)

head(out15)
head(out20)
head(out50)

library(dplyr)

# assume out15, out20, out50 all have columns: scaffold, start, end

# helper that checks overlap for one reference set
flag_overlap <- function(query, subject) {
  # for each row i in query, test if any row in subject overlaps fully
  mapply(function(chr, st, en) {
    any(
      subject$scaffold == chr &
        subject$start    <= st  &
        subject$end      >= en
    )
  },
  query$scaffold, query$start, query$end)
}

out15 <- out15 %>%
  mutate(
    overlaps20 = flag_overlap(out15, out20),
    overlaps50 = flag_overlap(out15, out50)
  )

# quick contingency
table(out15$overlaps20, out15$overlaps50)


library(dplyr)

# Helper: tests whether each query window overlaps any subject window
# by at least 50% of its own length, on the same scaffold
flag_overlap_50pct <- function(query, subject, pct = 0.5) {
  mapply(function(chr, st, en) {
    # filter to same scaffold once
    subj <- subject[subject$scaffold == chr, ]
    if (nrow(subj) == 0) return(FALSE)
    # compute intersection lengths
    inter <- pmin(subj$end,   en) - pmax(subj$start, st) + 1
    inter[inter < 0] <- 0       # negative → no overlap
    # proportion of query covered
    prop <- inter / (en - st + 1)
    any(prop >= pct)
  },
  query$scaffold, query$start, query$end)
}

# Apply to your 15kb outliers:
out15 <- out15 %>%
  mutate(
    overlaps20_50pct = flag_overlap_50pct(out15, out20, pct = 0.5),
    overlaps50_50pct = flag_overlap_50pct(out15, out50, pct = 0.5)
  )

# Check the new contingency
table(out15$overlaps20_50pct, out15$overlaps50_50pct)

out15
head(out15)
write_xlsx(out50, path = "out50_windows.xlsx")






library(ggplot2)
library(patchwork)

# 1) Define each plot (replace these with your actual ggplot calls)
p15 <- ggplot(all_stats_w15000_s2500, aes(x = cum_mid, y = FstWC_B_M_egg)) +
  geom_vline(data = chr_info, aes(xintercept = cum_start), color = "grey80", size = 0.3) +
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  scale_x_continuous(breaks = chr_info$center, labels = chr_info$scaffold, expand = c(0,0)) +
  labs(title = "15 kb/2.5 kb windows", y = expression(F[ST]), x = NULL) +
  theme_minimal() +   theme(
    legend.position      = "none",
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    # remove all vertical grid lines
    panel.grid.major.x   = element_blank(),
    panel.grid.minor.x   = element_blank(),
  )

p20 <- ggplot(all_stats_w20000_s5000, aes(x = cum_mid, y = FstWC_B_M_egg)) +
  geom_vline(data = chr_info, aes(xintercept = cum_start), color = "grey80", size = 0.3) +
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  scale_x_continuous(breaks = chr_info$center, labels = chr_info$scaffold, expand = c(0,0)) +
  labs(title = "20 kb/5 kb windows", y = expression(F[ST]), x = NULL) +
  theme_minimal() +   theme(
    legend.position      = "none",
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    # remove all vertical grid lines
    panel.grid.major.x   = element_blank(),
    panel.grid.minor.x   = element_blank())

p50 <- ggplot(all_stats_w50000_s25000, aes(x = cum_mid, y = FstWC_B_M_egg)) +
  geom_vline(data = chr_info, aes(xintercept = cum_start), color = "grey80", size = 0.3) +
  geom_point(aes(color = is_outlier), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  scale_x_continuous(breaks = chr_info$center, labels = chr_info$scaffold, expand = c(0,0)) +
  labs(title = "50 kb/25 kb windows", y = expression(F[ST]), x = "Chromosome / Scaffold") +
  theme_minimal() +   theme(
    legend.position      = "none",
    # remove all vertical grid lines
    panel.grid.major.x   = element_blank(),
    panel.grid.minor.x   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x       = element_text(size = 16))

# 2) Stack them
combined_plot <- p15 / p20 / p50 + plot_layout(ncol = 1, heights = c(1,1,1))

# 3) Display (and optionally save)
print(combined_plot)
ggsave("fst_manhattan_stacked.png", combined_plot, width = 10, height = 12, dpi = 300)








