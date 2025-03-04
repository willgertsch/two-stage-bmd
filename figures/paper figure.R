# figures for manuscript
library(ggplot2)
library(dplyr)
library(tidyr)

`%nin%` <- function(x, y) !(x %in% y)

# Naive designs ################################################################
eval_results_naive <- readRDS("~/Rprojects/two-stage-simulations/results/naive/Kociba et al. 1978-JBRC-Deguelin-Kociba log-logistic-Single stage uniform design-Single stage log-uniform design/Varying N-num_doses/eval_results.rds")
naive_results = eval_results_naive$`all eval`

# kociba
# N
naive_results %>%
  filter(.dgp_name == 'Kociba et al. 1978', num_doses == 4) %>%
  group_by(.method_name, N) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=N, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  labs(
    x = 'Total sample size',
    y = 'Median sim. CI length',
    color = ''
    ) +
  theme_bw() 

# number of doses
naive_results %>%
  filter(.dgp_name == 'Kociba et al. 1978', N == 117) %>%
  group_by(.method_name, num_doses) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=num_doses, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  labs(
    x = 'Number of doses',
    y = 'Median sim. CI length',
    color = ''
  ) +
  theme_bw()

# Degeulin
naive_results %>%
  filter(.dgp_name == 'Deguelin', num_doses == 4) %>%
  group_by(.method_name, N) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=N, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  labs(
    x = 'Total sample size',
    y = 'Median sim. CI length',
    color = ''
  ) +
  theme_bw()

naive_results %>%
  filter(.dgp_name == 'Deguelin', N == 146) %>%
  group_by(.method_name, num_doses) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=num_doses, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  labs(
    x = 'Number of doses',
    y = 'Median sim. CI length',
    color = ''
  ) +
  theme_bw()

# JBRC data
naive_results %>%
  filter(.dgp_name == 'JBRC', num_doses == 5) %>%
  group_by(.method_name, N) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=N, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  labs(
    x = 'Total sample size',
    y = 'Median sim. CI length',
    color = ''
  ) +
  theme_bw()

naive_results %>%
  filter(.dgp_name == 'JBRC', N == 99) %>%
  group_by(.method_name, num_doses) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=num_doses, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  labs(
    x = 'Number of doses',
    y = 'Median sim. CI length',
    color = ''
  ) +
  theme_bw()


# Kociba data scenario #########################################################
eval_results_kociba <- readRDS("~/Rprojects/two-stage-simulations/results/kociba_qlinear/Kociba et al. 1978/Varying N1/eval_results.rds")
all_eval_kociba = eval_results_kociba$`all eval`


# boxplots at equal allocation
all_eval_kociba %>%
  filter(N1 == 117) %>%
  mutate(method_name = reorder(as.factor(.method_name), CI_length, FUN = median)) %>%
  ggplot(aes(x = method_name, y = CI_length)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(
    y = 'Simulated CI length',
    x = ''
  )

# length by N1
all_eval_kociba %>% 
  group_by(.method_name, N1) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=N1, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  theme_bw() +
  labs(
    x = 'Stage 1 sample size',
    y = 'Median sim. CI length',
    color = ''
  )

# bias by N1
all_eval_kociba %>%
  mutate(bmd_diff = as.numeric(bmd_diff)) %>%
  group_by(.method_name, N1) %>%
  summarise(bias = mean(bmd_diff, na.rm = T)) %>%
  ggplot(aes(x=N1, y=bias, color = .method_name)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    y = 'Bias(BMD)',
    x = 'Stage 1 sample size',
    color = ''
  )

# coverage by N1
all_eval_kociba %>%
  mutate(coverage = as.numeric(coverage)) %>%
  group_by(.method_name, N1) %>%
  summarise(coverage_rate = mean(coverage, na.rm = T)) %>%
  ggplot(aes(x = N1, y = coverage_rate, color = .method_name)) +
  geom_point(size = 2) + 
  geom_line(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = 'dashed') +
  labs(
    x = 'Stage 1 sample size',
    y = 'BMD CI coverage rate',
    color = ''
  )

# dose placements
Nsim = 3000
Nsamples = 100
all_eval_kociba %>%
  filter(N1 == 117) %>%
  filter(.method_name %nin% c('Single stage uniform design', 'Single stage log-uniform design')) %>%
  #slice_sample(n = Nsamples) %>%
  select(.method_name, N1, bmd, stage2_x, stage2_n) %>%
  unnest(cols = c(stage2_x, stage2_n)) %>%
  ggplot(aes(x = stage2_x, y = stage2_n/(N1*Nsim))) +
  geom_col(fill = 'blue', color = 'black') +
  facet_wrap(~.method_name, scales = 'free_y') +
  theme_bw() +
  labs(
    y = 'Prop. stage 2 dose placements',
    x = 'Dose'
  )


# Deguelin data scenario #######################################################
eval_results_deguelin <- readRDS("~/Rprojects/two-stage-simulations/results/deguelin_loglogistic/Deguelin/Varying N1/eval_results.rds")
all_eval_deguelin = eval_results_deguelin$`all eval`

# boxplots at equal allocation
all_eval_deguelin %>%
  filter(N1 == 146) %>%
  filter(.method_name != 'Single stage log-uniform design') %>%
  mutate(method_name = reorder(as.factor(.method_name), CI_length, FUN = median)) %>%
  ggplot(aes(x = method_name, y = CI_length)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(
    y = 'Simulated CI length',
    x = ''
  )

# length by N1
all_eval_deguelin %>% 
  filter(.method_name != 'Single stage log-uniform design') %>%
  group_by(.method_name, N1) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=N1, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  theme_bw() +
  labs(
    x = 'Stage 1 sample size',
    y = 'Median sim. CI length',
    color = ''
  )

# bias by N1
all_eval_deguelin %>%
  filter(.method_name != 'Single stage log-uniform design') %>%
  mutate(bmd_diff = as.numeric(bmd_diff)) %>%
  group_by(.method_name, N1) %>%
  summarise(bias = mean(bmd_diff, na.rm = T)) %>%
  ggplot(aes(x=N1, y=bias, color = .method_name)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    y = 'Bias(BMD)',
    x = 'Stage 1 sample size',
    color = ''
  )

# coverage by N1
all_eval_deguelin %>%
  filter(.method_name != 'Single stage log-uniform design') %>%
  mutate(coverage = as.numeric(coverage)) %>%
  group_by(.method_name, N1) %>%
  summarise(coverage_rate = mean(coverage, na.rm = T)) %>%
  ggplot(aes(x = N1, y = coverage_rate, color = .method_name)) +
  geom_point(size = 2) + 
  geom_line(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = 'dashed') +
  labs(
    x = 'Stage 1 sample size',
    y = 'BMD CI coverage rate',
    color = ''
  )

# dose placements at equal allocation
Nsim = 3000
Nsamples = 200
all_eval_deguelin %>%
  filter(N1 == 146) %>%
  filter(.method_name %nin% c('Single stage uniform design', 'Single stage log-uniform design')) %>%
  #slice_sample(n = Nsamples) %>%
  select(.method_name, N1, bmd, stage2_x, stage2_n) %>%
  unnest(cols = c(stage2_x, stage2_n)) %>%
  ggplot(aes(x = stage2_x, y = stage2_n/(N1*Nsim))) +
  geom_col(fill = 'blue', color = 'black') +
  facet_wrap(~.method_name, scales = 'free_y') +
  theme_bw() +
  labs(
    y = 'Prop. stage 2 dose placements',
    x = 'Dose'
  )

# JBRC data scenario ###########################################################
eval_results_jbrc <- readRDS("~/Rprojects/two-stage-simulations/results/jbrc_hill_logunif/JBRC-Locally D-optimal design-Augmented D-optimal design-Locally c-optimal design-Augmented c-optimal design-Single stage uniform design-Single stage log-uniform design/Varying N1-first_stage-num_doses/eval_results.rds")
all_eval_jbrc = eval_results_jbrc$`all eval`

# boxplots at equal allocation
all_eval_jbrc %>%
  filter(N1 == 99) %>%
  filter(.method_name != 'Single stage uniform design') %>%
  mutate(method_name = reorder(as.factor(.method_name), CI_length, FUN = median)) %>%
  ggplot(aes(x = method_name, y = CI_length)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(
    y = 'Simulated CI length',
    x = ''
  )

# length by N1
all_eval_jbrc %>% 
  filter(.method_name != 'Single stage uniform design') %>%
  group_by(.method_name, N1) %>%
  summarise(med_ci = median(CI_length, na.rm = T)) %>%
  ggplot(aes(x=N1, y=med_ci, color = .method_name)) +
  geom_point() + geom_line() +
  theme_bw() +
  labs(
    x = 'Stage 1 sample size',
    y = 'Median sim. CI length',
    color = ''
  )

# bias by N1
all_eval_jbrc %>%
  filter(.method_name != 'Single stage uniform design') %>%
  mutate(bmd_diff = as.numeric(bmd_diff)) %>%
  group_by(.method_name, N1) %>%
  summarise(bias = mean(bmd_diff, na.rm = T)) %>%
  ggplot(aes(x=N1, y=bias, color = .method_name)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    y = 'Bias(BMD)',
    x = 'Stage 1 sample size',
    color = ''
  )

# coverage by N1
all_eval_jbrc %>%
  filter(.method_name != 'Single stage uniform design') %>%
  mutate(coverage = as.numeric(coverage)) %>%
  group_by(.method_name, N1) %>%
  summarise(coverage_rate = mean(coverage, na.rm = T)) %>%
  ggplot(aes(x = N1, y = coverage_rate, color = .method_name)) +
  geom_point(size = 2) + 
  geom_line(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = 'dashed') +
  labs(
    x = 'Stage 1 sample size',
    y = 'BMD CI coverage rate',
    color = ''
  )

# dose placements at equal allocation
Nsim = 3000
all_eval_jbrc %>%
  filter(N1 == 99) %>%
  filter(.method_name %nin% c('Single stage uniform design', 'Single stage log-uniform design')) %>%
  select(.method_name, N1, bmd, stage2_x, stage2_n) %>%
  unnest(cols = c(stage2_x, stage2_n)) %>%
  ggplot(aes(x = stage2_x, y = stage2_n/(N1*Nsim))) +
  geom_col(fill = 'blue', color = 'black') +
  facet_wrap(~.method_name, scales = 'free_y') +
  theme_bw() +
  labs(
    y = 'Prop. stage 2 dose placements',
    x = 'Dose'
  )
