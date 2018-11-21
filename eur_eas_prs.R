# Set up data/directories -------------------------------------------------

lapply(c('plyr', 'dplyr', 'ggplot2', 'tidyr', 'rcompanion', 'forcats', 'plotrix', 'cowplot'), require, character.only = TRUE)

setwd('/Users/alicia/daly_lab/pgc_scz/eur_eas')
options(stringsAsFactors = FALSE)

prefix <- c('mix_cno1_asn_ml-qc.hg19.ch.fl.bg', 'ms.scz_tai1_asn_hh-qc.hg19.ch.fl.bg', 'ms.scz_tai2_asn_hh-qc.hg19.ch.fl.bg', 'scz_bix1_asn_hh-qc.hg19.ch.fl.bg', 'scz_bix2_asn_hh-qc.hg19.ch.fl.bg', 'scz_bix3_asn_hh-qc.hg19.ch.fl.bg', 'scz_hnk1_asn_ml-qc.hg19.ch.fl.bg', 'scz_imh1_asn_ml-qc.hg19.ch.fl.bg', 'scz_imh2_asn_ml-qc.hg19.ch.fl.bg', 'scz_jpn1_asn_ml-qc.hg19.ch.fl.bg', 'scz_umc1_asn_hh-qc.hg19.ch.fl.bg', 'scz_uwa1_asn_hh-qc.hg19.ch.fl.bg', 'scz_xju1_asn_ml-qc.hg19.ch.fl.bg')
cohorts = separate(data.frame(name=prefix), name, into=c('first', 'cohort'), sep='_')$cohort
p <- rev(c('all', 'p5', 'p2', 'p1', 'p05', 'pe2', 'pe3', 'pe4', 'pe6', 'p5e8'))
n_snps <- data.frame(code=paste('n', 1:10, sep=''),
                     num=c(100, 1500, 3000, 5000, 10000, 15000, 25000, 35000, 50000, 65000))
ps <- data.frame(code=paste('s', 1:10, sep=''),
                 num=c(5e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.5, 1))
covar <- read.table('../merge14_pgc_scz_asia.menv.mds_cov', header=T)
covar <- covar %>%
  select(-c(st1:st13)) %>%
  separate(FID, c('cc', 'pheno', 'name'), sep='_', remove=F)

# Read/process PRS --------------------------------------------------------

read_files <- function(pop_cutoff, cohort_name) {
  my_file <- read.table(paste0(cohort_name, '.', pop_cutoff, '.profile'), header=T)
  my_file$PHENO <- ifelse(grepl('cas', my_file$FID), 2, my_file$PHENO) # add back in phenotypes when they're lost (e.g. missing sexes)
  my_file$PHENO <- ifelse(grepl('con', my_file$FID), 1, my_file$PHENO)
  my_file$pop_cutoff <- pop_cutoff
  my_file$cohort_name <- cohort_name
  return(my_file)
}

# all_scores <- ldply(p, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })
# all_scores <- separate(all_scores, cohort_name, into=c('first', 'cohort'), sep='_')
# all_scores$PHENO <- all_scores$PHENO - 1
# all_scores$PHENO <- factor(all_scores$PHENO)
# all_scores <- subset(all_scores, PHENO!='-10')
# all_scores$p_cutoff <- factor(all_scores$p_cutoff, levels=p)
# 

# all_score_covar <- merge(all_scores, covar, by=c('FID', 'IID'))

# Same num SNPs -----------------------------------------------------------

setwd('/Users/alicia/daly_lab/pgc_scz/eur_eas/same_snp_num5')

n_snps_eas <- paste('EAS.n', 1:10, sep='')
n_snps_eas2 <- paste('looeas.n', 1:10, sep='')
n_snps_eur <- paste('eur.n', 1:10, sep='')
all_scores_eas <- ldply(n_snps_eas, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })
all_scores_eas2 <- ldply(n_snps_eas2, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })
all_scores_eur <- ldply(n_snps_eur, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })

bind_eas_eur <- function(eur, eas) {
  all_scores_bound <- bind_rows(eur, eas)
  all_scores_bound <- separate(all_scores_bound, cohort_name, into=c('first', 'cohort'), sep='_')
  all_scores_bound$PHENO <- all_scores_bound$PHENO - 1
  all_scores_bound$PHENO <- factor(all_scores_bound$PHENO)
  all_scores_bound <- subset(all_scores_bound, PHENO!='-10')
  all_scores_bound_covar <- all_scores_bound %>%
    inner_join(covar, by=c('FID'='FID', 'IID'='IID', 'cohort'='name'))
  return(all_scores_bound_covar)
}

bind_n_snps <- function(eur, eas) {
  all_scores_bound_covar <- bind_eas_eur(eur, eas)
  all_scores_bound <- separate(all_scores_bound_covar, pop_cutoff, into=c('study_pop', 'num_snps'))
  return(all_scores_bound)
}

all_scores_eureas_covar <- bind_n_snps(all_scores_eur, all_scores_eas)
all_scores_eureas2_covar <- bind_n_snps(all_scores_eur, all_scores_eas2)

h2l_R2N <- function(k, r2n, p) {
  # k = population prevalence
  # r2n = Nagelkerke's attributable to genomic profile risk score
  # p = proportion of cases in the case-control sample
  # calculates proportion of variance explained on the liability scale
  # from ABC at http://www.complextraitgenomics.com/software/
  # Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  h2l_R2N <- cc * e * r2n / (1 + cc * e * theta * r2n) 
}

compute_r2 <- function(my_subset, coh, study) {
  #print(paste(coh, study, n_snps, sep=', '))
  model1 <- glm(PHENO~SCORE+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10, data=my_subset, family='binomial')
  model0 <- glm(PHENO~C1+C2+C3+C4+C5+C6+C7+C8+C9+C10, data=my_subset, family='binomial')
  my_r2 <- nagelkerke(model1, null=model0)
  r2n <- my_r2$Pseudo.R.squared.for.model.vs.null[3]
  r2l <- h2l_R2N(k=0.01, r2n=r2n, p=sum(my_subset$PHENO==1)/sum(my_subset$PHENO %in% c(0,1)))
  df <- data.frame(cohort=coh,
                   study_pop=study,
                   r2_nagelkerke=r2n,
                   r2_liability=r2l,
                   p=my_r2$Likelihood.ratio.test[4])
  return(df)
}

make_subset_n_snps <- function(dataset, coh, n_snps, study) {
  my_subset <- subset(dataset, cohort==coh&num_snps==n_snps&study_pop==study)
  my_df <- compute_r2(my_subset, coh, study)
  my_df$num_snps=n_snps
  return(my_df)
}

study_pops <- c('EAS', 'eur')
r2_cohorts_eureas <- ldply(cohorts, function(x) { ldply(n_snps$code, function(y) { ldply(study_pops, function(z) { make_subset_n_snps(all_scores_eureas_covar, x, y, z) } ) } ) } )
study_pops <- c('looeas', 'eur')
r2_cohorts_eureas2 <- ldply(cohorts, function(x) { ldply(n_snps$code, function(y) { ldply(study_pops, function(z) { make_subset_n_snps(all_scores_eureas2_covar, x, y, z) } ) } ) } )

r2_summarize <- function(dataset) {
  r2_summary <- dataset %>%
    group_by(study_pop, num_snps) %>%
    dplyr::summarize(r2n_mean=mean(r2_nagelkerke), semn = std.error(r2_nagelkerke),
                     r2l_mean=mean(r2_liability), seml = std.error(r2_liability))
  
  r2_summary <- merge(r2_summary, n_snps, by.x='num_snps', by.y='code')
  #r2_summary$num <- factor(r2_summary$num)
  r2_summary <- r2_summary %>% 
    mutate(study_pop = recode(study_pop, eur = "European", looeas = "East Asian", EAS = 'East Asian'))
  return(r2_summary)
}

#left_join(mean_snps, by=c('cohort', 'p_cutoff')) %>%

r2_summarize_tgp <- r2_summarize(r2_cohorts_eureas)
r2_summarize_bg <- r2_summarize(r2_cohorts_eureas2)


plot_prs <- function(dataset, r2_mean, sem, r2_type, adj_x=F, xlab='P-value threshold') {
  dataset$ymin <- dataset[,r2_mean] - dataset[,sem]
  dataset$ymax <- dataset[,r2_mean] + dataset[,sem]
  p_eas_eur <- ggplot(dataset, aes_string(x='num', y=r2_mean, color='study_pop', ymin='ymin', ymax='ymin')) +
    geom_point(position=pd) +
    geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd)  +
    labs(y=bquote(East~Asian~R[.(r2_type)]^2), x=xlab) +
    scale_color_brewer(palette='Set1', name='Training data') +
    theme_classic() +
    theme(legend.justification = c(1, 0), 
          legend.position = c(1, 0),
          text = element_text(size=18),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill="transparent"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          axis.text.x = element_text(angle = 45, hjust = 1))
  if(adj_x) {
    p_eas_eur <- p_eas_eur +
      scale_x_sqrt(breaks = c(100, 1500, 5000, 15000, 25000, 35000, 50000, 65000),
                   labels = c(100, 1500, 5000, 15000, 25000, 35000, 50000, 'all'))
  }
  return(p_eas_eur)
}

pd <- position_dodge(10) # move them .05 to the left and right
p_eas_eur_tgp_snp_n <- plot_prs(r2_summarize_tgp, 'r2n_mean', 'semn', 'Nagelkerke', adj_x=T, xlab='Number of SNPs in PRS')
p_eas_eur_tgp_snp_l <- plot_prs(r2_summarize_tgp, 'r2l_mean', 'seml', 'Liability', adj_x=T, xlab='Number of SNPs in PRS')
p_eas_eur_bg_snp_n <- plot_prs(r2_summarize_bg, 'r2n_mean', 'semn', 'Nagelkerke', adj_x=T, xlab='Number of SNPs in PRS')
p_eas_eur_bg_snp_l <- plot_prs(r2_summarize_bg, 'r2l_mean', 'seml', 'Liability', adj_x=T, xlab='Number of SNPs in PRS')
ggsave('p_eas_eur_pval_n.tgp.pdf', p_eas_eur_tgp_snp_n, height=6, width=5)
ggsave('p_eas_eur_pval_l.tgp.pdf', p_eas_eur_tgp_snp_l, height=6, width=5)
ggsave('p_eas_eur_pval_n.bg.pdf', p_eas_eur_bg_snp_n, height=6, width=5)
ggsave('p_eas_eur_pval_l.bg.pdf', p_eas_eur_bg_snp_l, height=6, width=5)

# to fix: legend names, make size correspond to max(-log10(p))

# Same p-val threshold (updated) ------------------------------------------

setwd('/Users/alicia/daly_lab/pgc_scz/eur_eas/same_p_threshold5')

ps_eas <- paste('EAS.s', 1:10, sep='')
ps_eas2 <- paste('looeas.s', 1:10, sep='')
ps_eur <- paste('eur.s', 1:10, sep='')
all_scores_eas <- ldply(ps_eas, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })
all_scores_eas2 <- ldply(ps_eas2, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })
all_scores_eur <- ldply(ps_eur, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })

bind_p_val <- function(eur, eas) {
  all_scores_bound_covar <- bind_eas_eur(eur, eas)
  all_scores_bound <- separate(all_scores_bound_covar, pop_cutoff, into=c('study_pop', 'ps'))
  return(all_scores_bound)
}

all_scores_eureas_covar <- bind_p_val(all_scores_eur, all_scores_eas)
all_scores_eureas2_covar <- bind_p_val(all_scores_eur, all_scores_eas2)

make_subset_p_val <- function(dataset, coh, p_threshold, study) {
  my_subset <- subset(dataset, cohort==coh&ps==p_threshold&study_pop==study)
  my_df <- compute_r2(my_subset, coh, study) #something isn't getting passed through here
  my_df$p_val=p_threshold
  return(my_df)
}


# my_r2 <- nagelkerke(model1, null=model0)
# r2_observed <- my_r2$Pseudo.R.squared.for.model.vs.null[2]
# r2_nagelkerke <- my_r2$Pseudo.R.squared.for.model.vs.null[3]
# # K = population prevalence
# # P = cohort prevalence
# P = sum(my_subset$PHENO==1) / ( sum(my_subset$PHENO==1) + sum(my_subset$PHENO==0) )
# z = dnorm(qnorm(K))
# m = z/K
# t_threshold = qnorm(K, lower=F)
# theta = m*(P-K)/(1-K)*(m*(P-K)/(1-K)-t_threshold)
# C = (K*(1-K)/z^2)*(K*(1-K)/(P*(1-P)))
# 
# r2_liability = 
#   df <- data.frame(p_val=p_threshold,
#                    cohort=coh,
#                    study_pop=study,
#                    r2_observed=r2_observed,
#                    r2_nagelkerke=r2_nagelkerke,
#                    #r2_liability=r2_liability,
#                    p=my_r2$Likelihood.ratio.test[4])

study_pops <- c('EAS', 'eur')
r2_cohorts_eureas <- ldply(cohorts, function(x) { ldply(ps$code, function(y) { ldply(study_pops, function(z) { make_subset_p_val(all_scores_eureas_covar, x, y, z) } ) } ) } )
study_pops <- c('looeas', 'eur')
r2_cohorts_eureas2 <- ldply(cohorts, function(x) { ldply(ps$code, function(y) { ldply(study_pops, function(z) { make_subset_p_val(all_scores_eureas2_covar, x, y, z) } ) } ) } )

r2_summarize <- function(dataset) {
  r2_summary <- dataset %>%
    group_by(study_pop, p_val) %>%
    dplyr::summarize(r2n_mean=mean(r2_nagelkerke), semn = std.error(r2_nagelkerke),
                     r2l_mean=mean(r2_liability), seml = std.error(r2_liability))
  
  r2_summary <- merge(r2_summary, ps, by.x='p_val', by.y='code')
  r2_summary$num <- factor(r2_summary$num)
  r2_summary <- r2_summary %>% 
    mutate(study_pop = recode(study_pop, eur = "European", looeas = "East Asian", EAS = 'East Asian'))
  return(r2_summary)
}

#left_join(mean_snps, by=c('cohort', 'p_cutoff')) %>%

r2_summarize_tgp <- r2_summarize(r2_cohorts_eureas)
r2_summarize_bg <- r2_summarize(r2_cohorts_eureas2)

pd <- position_dodge(0.25) # move them .05 to the left and right
p_eas_eur_tgp_p_n <- plot_prs(r2_summarize_tgp, 'r2n_mean', 'semn', 'Nagelkerke')
p_eas_eur_tgp_p_l <- plot_prs(r2_summarize_tgp, 'r2l_mean', 'seml', 'Liability')
p_eas_eur_bg_p_n <- plot_prs(r2_summarize_bg, 'r2n_mean', 'semn', 'Nagelkerke')
p_eas_eur_bg_p_l <- plot_prs(r2_summarize_bg, 'r2l_mean', 'seml', 'Liability')
ggsave('p_eas_eur_pval_n.tgp.pdf', p_eas_eur_tgp_p_n, height=6, width=5)
ggsave('p_eas_eur_pval_l.tgp.pdf', p_eas_eur_tgp_p_l, height=6, width=5)
ggsave('p_eas_eur_pval_n.bg.pdf', p_eas_eur_bg_p_n, height=6, width=5)
ggsave('p_eas_eur_pval_l.bg.pdf', p_eas_eur_bg_p_l, height=6, width=5)

main_prs <- plot_grid(p_eas_eur_tgp_p_l, p_eas_eur_tgp_snp_l, labels=c('A', 'B'), ncol=2)
supp_prs <- plot_grid(p_eas_eur_bg_p_l, p_eas_eur_bg_snp_l,
                      p_eas_eur_tgp_p_n, p_eas_eur_tgp_snp_n,
                      p_eas_eur_bg_p_n, p_eas_eur_bg_snp_n,
                      labels=c('A','B','C','D','E','F'),
                      ncol=2)
save_plot('../p_eas_eur_l.tgp.pdf', main_prs, base_width=8, base_height=4)
save_plot('../p_eas_eur.supp.pdf', supp_prs, base_width=8, base_height=12)

# Compute R^2 --------------------------------------------------------------

compute_r2 <- function(dataset, coh, pval) {
  #print(paste(coh, pval, sep=', '))
  my_subset <- subset(dataset, cohort==coh&pop_cutoff==pval)
  model1 <- glm(PHENO~SCORE+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10, data=my_subset, family='binomial')
  model0 <- glm(PHENO~C1+C2+C3+C4+C5+C6+C7+C8+C9+C10, data=my_subset, family='binomial')
  my_r2 <- nagelkerke(model1, null=model0)
  df <- data.frame(p_cutoff=pval,
                   cohort=coh,
                   r2=my_r2$Pseudo.R.squared.for.model.vs.null[3],
                   p=my_r2$Likelihood.ratio.test[4])
  return(df)
}

r2_cohorts <- ldply(cohorts, function(x) { ldply(p, function(y) { compute_r2(all_score_covar, x, y) } ) })
#current_r2 <- (1-exp((logLik(model0)-logLik(model1))[1])^(2/n))/(1-(exp(logLik(model0)[1])^(2/n))) #Cox and Snell
#current_r2 <- (1-exp(2*(logLik(model0)-logLik(model1))[1]/n))/(1-(exp(2*logLik(model0)[1]/n))) #Nagelkerke
r2_cohorts$p_cutoff <- factor(r2_cohorts$p_cutoff, levels=p)

# Aggregate data summaries ------------------------------------------------

mean_snps <- all_scores %>%
  group_by(cohort, p_cutoff) %>%
  dplyr::summarize(CNT=mean(CNT)) %>%
  ungroup()

mean_cc_score <- all_scores %>%
  group_by(PHENO, cohort, p_cutoff) %>%
  dplyr::summarize(prs=mean(SCORE)) %>%
  mutate(phenotype = fct_recode(PHENO, case = '1', control = '0')) %>%
  ungroup() %>%
  select(-PHENO) %>%
  spread(phenotype, prs) %>%
  mutate(case_greater=case>control)

spread(PHENO, value=SCORE)

r2_cohorts <- r2_cohorts %>%
  left_join(mean_snps, by=c('cohort', 'p_cutoff')) %>%
  left_join(mean_cc_score, by=c('cohort', 'p_cutoff')) %>%
  mutate(p_string = fct_recode(p_cutoff, '5e-8' = 'p5e8', '1e-6' = 'pe6', '1e-4' = 'pe4', '1e-3' = 'pe3', '1e-2' = 'pe2',
                               '0.05' = 'p05', '0.1' = 'p1', '0.2' = 'p2', '0.5' = 'p5', '1' = 'all'))

p_snp_r2 <- ggplot(r2_cohorts, aes(x=CNT, y=r2)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~p_cutoff, scales='free', ncol=5) +
  theme_bw() +
  xlab('Number of alleles') +
  ylab(expression("Nagelkerke's"~R^2))
ggsave('snp_r2.pdf', p_snp_r2, width=12, height=10)

write.table(t(spread(r2_cohorts[,1:3], p_cutoff, r2)), 'eur_eas_r2.txt', col.names=F, quote=F)
write.table(t(spread(r2_cohorts[,c(1,2,4)], p_cutoff, p)), 'eur_eas_p.txt', col.names=F, quote=F)
write.table(t(spread(r2_cohorts[,c(1,2,5)], p_cutoff, CNT)), 'eur_eas_mean_snps.txt', col.names=F, quote=F)

eas_eas <- read.table('/Users/alicia/daly_lab/pgc_scz/eas_eas/eas_eas_r2.txt', header=T, sep='\t')
eas_eas_cnt <- read.table('/Users/alicia/daly_lab/pgc_scz/eas_eas/eas_eas_cnt.txt', header=T, sep='\t')
r2_cohorts$P.value.threshold <- as.numeric(as.character(r2_cohorts$p_string))
r2_cohorts$pop='European'
r2_eas_eas <- eas_eas %>%
  gather(cohort, r2, -P.value.threshold) %>%
  mutate(pop='East Asian')
r2_eas_eas_cnt <- eas_eas_cnt %>%
  gather(cohort, CNT, -P.value.threshold) %>%
  mutate(pop='East Asian')
r2_eas_eas$CNT <- r2_eas_eas_cnt$CNT


r2_all <- bind_rows(r2_cohorts, r2_eas_eas)
r2_all$P.value.threshold <- factor(r2_all$P.value.threshold)
r2_summary <- subset(r2_all, !P.value.threshold %in% c('1e-06', '5e-08')) %>%
  group_by(pop, P.value.threshold) %>%
  dplyr::summarize(r2_mean=mean(r2), sem = std.error(r2), snps=mean(CNT)/2)

pd <- position_dodge(0.1) # move them .05 to the left and right
#num_snps <- read.table('num_snps_p_val.txt', header=T)
p_eas_eur <- ggplot(r2_summary, aes(x=P.value.threshold, y=r2_mean, color=pop)) +
  geom_point(position=pd, aes(size=snps)) +
  geom_errorbar(aes(ymin=r2_mean-sem, ymax=r2_mean+sem), width=.1, position=pd)  +
  labs(y=expression(East~Asian~Prediction~R^2), x='p-value threshold') +
  scale_color_brewer(palette='Set1', name='Training data') +
  scale_size_area(trans='log10', breaks=c(100,1000,10000), name='# SNPs') +
  theme_classic() +
  theme(legend.position = c(0.2, 0.8),
        text = element_text(size=18),
        panel.background = element_rect(fill = "transparent",colour = NA),
        #panel.grid.minor = element_blank(), 
        #panel.grid.major = element_blank(),
        legend.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('p_eas_eur.pdf', p_eas_eur, height=6, width=5)

r2_01 <- subset(r2_all, P.value.threshold==1)
p_r2 <- ggplot(r2_01, aes(x=pop, y=r2, fill=pop)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill='white') +
  scale_fill_brewer(palette='Set1') +
  ylim(c(0,max(r2_01$r2))) +
  theme_classic() +
  labs(x='Training population', y=expression(East~Asian~Prediction~R^2)) +
  theme(legend.position = "none",
        text = element_text(size=18))
ggsave('predict_r2_eur_eas.pdf', p_r2, height=5,width=5)

p_r2_cohort <- ggplot(r2_cohorts) +
  geom_bar(aes(x=cohort, y=r2, fill=p_string), stat='identity', position='dodge') +
  scale_fill_brewer(palette='Spectral', name='P value') +
  ylab(expression("Nagelkerke's"~R^2)) +
  xlab('Cohort') +
  ylim(0,0.15) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=18))

ggsave('all_prs_cohort_r2.pdf', p_r2_cohort, width=12, height=10)

compute_r2_combined <- function(pval) {
  #print(paste(coh, pval, sep=', '))
  my_subset <- subset(all_score_covar, p_cutoff==pval)
  model1 <- glm(PHENO~SCORE+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+cohort, data=my_subset, family='binomial')
  model0 <- glm(PHENO~C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+cohort, data=my_subset, family='binomial')
  my_r2 <- nagelkerke(model1, null=model0)
  df <- data.frame(p_cutoff=pval,
                   r2=my_r2$Pseudo.R.squared.for.model.vs.null[3],
                   p=my_r2$Likelihood.ratio.test[4])
  return(df)
}

mean_snps <- all_scores %>%
  group_by(p_cutoff) %>%
  summarize(CNT=mean(CNT)) %>%
  ungroup()
#spread(cohort, value=CNT)

r2_comb <- ldply(p, function(x) { compute_r2_combined(x) })
r2_comb <- r2_comb %>%
  left_join(mean_snps, by=c('p_cutoff'))


# Case/control scores -----------------------------------------------------

cc_scores <- all_scores %>%
  group_by(cohort, p_cutoff, PHENO) %>%
  summarize(mean=mean(SCORE)) %>%
  spread(PHENO, value=mean)
write.table(cc_scores, 'cc_scores.txt', row.names=F, quote=F)

cc_scores <- all_scores %>%
  group_by(cohort, p_cutoff, PHENO) %>%
  summarize(mean=mean(SCORE))
p_cc <- ggplot(cc_scores, aes(x=PHENO, y=mean)) +
  facet_wrap(~p_cutoff, ncol=5, scales='free') +
  geom_violin() +
  theme_bw()
ggsave('cc_p.pdf', p_cc, width=12, height=10)

p1 <- ggplot(all_scores, aes(x=SCORE, color=PHENO)) +
  facet_grid(p_cutoff~cohort, scales='free') +
  #facet_wrap(~cohort) +
  geom_density(alpha=0.5) +
  theme_bw()

p2 <- ggplot(all_scores, aes(x=SCORE, color=PHENO)) +
  facet_wrap(~p_cutoff) +
  geom_density(alpha=0.75) +
  theme_bw()

ggsave('all_prs_cohort.pdf', p1, width=12, height=10)
ggsave('all_prs_comb.pdf', p2)

permute_r2 <- function(p, cohort) {
  my_subset <- subset(all_scores, p_cutoff==p & cohort_name==cohort)
  my_subset$PHENO <- sample(my_subset$PHENO)
  
  prs <- merge(my_subset, covar, by=c('FID', 'IID')) %>%
    subset(PHENO %in% c(1, 2))
  #prs$study <- factor(prs$study)
  prs$PHENO <- prs$PHENO - 1
  
  model0 <- glm(PHENO~C1+C2+C3+C4+C5+C6+C7+C8+C9+C10, data=prs, family='binomial')
  model1 <- glm(PHENO~SCORE+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10, data=prs, family='binomial')
  
  #my_anova <- anova(model0, model1)
  n <- nrow(prs)
  current_r2 <- (1-exp((logLik(model0)-logLik(model1))[1])^(2/n))/(1-(exp(logLik(model0)[1])^(2/n)))
  my_data <- data.frame(p_cutoff=p, cohort_name=cohort, r2=current_r2)
  return(my_data)
}

test <- permute_r2('p1', 'mix_cno1_asn_ml-qc')
test <- permute_r2('p2', 'mix_cno1_asn_ml-qc')


all_permuted <- function() {
  print(date())
  return(ldply(p, function(x) { ldply(prefix, function(y) { permute_r2(x, y) } ) }))
}

my_permutations <- replicate(100, all_permuted(), simplify = FALSE) %>% bind_rows()

summary_permute <- my_permutations %>%
  group_by(p_cutoff, cohort_name) %>%
  dplyr::summarize(mean=mean(r2), lower_ci=quantile(r2, probs=c(0.025))[[1]], upper_ci=quantile(r2, probs=c(0.975))[[1]])

my_table <- summary_permute %>%
  select(p_cutoff, cohort_name, mean) %>%
  dcast(p_cutoff~cohort_name)

lower_ci <- summary_permute %>%
  select(p_cutoff, cohort_name, lower_ci) %>%
  dcast(p_cutoff~cohort_name)

upper_ci <- summary_permute %>%
  select(p_cutoff, cohort_name, upper_ci) %>%
  dcast(p_cutoff~cohort_name)

# 10 p-value cutoffs
# 12 prefixes
# 100 permutations each

## questions for hailiang:
# convert con_ IDs with PHENO=-9 to controls? 


# European SNPs, East Asian effect sizes ----------------------------------

setwd('/Users/alicia/daly_lab/pgc_scz/eur_snps_eas_effects/')

all_scores_eurasn <- ldply(p, function(x) { ldply(prefix, function(y) { read_files(x, y) } ) })
all_scores_eurasn <- separate(all_scores_eurasn, cohort_name, into=c('first', 'cohort'), sep='_')
all_scores_eurasn$PHENO <- all_scores_eurasn$PHENO - 1
all_scores_eurasn$PHENO <- factor(all_scores_eurasn$PHENO)
all_scores_eurasn <- subset(all_scores_eurasn, PHENO!='-10')
all_scores_eurasn$p_cutoff <- factor(all_scores_eurasn$p_cutoff, levels=p)
all_scores_eurasn_covar <- merge(all_scores_eurasn, covar, by=c('FID', 'IID'))

r2_cohorts_eurasn <- ldply(cohorts, function(x) { ldply(p, function(y) { compute_r2(all_scores_eurasn_covar, x, y) } ) })

mean_snps <- all_scores_eurasn %>%
  group_by(cohort, p_cutoff) %>%
  dplyr::summarize(CNT=mean(CNT)) %>%
  ungroup()

mean_cc_score <- all_scores_eurasn %>%
  group_by(PHENO, cohort, p_cutoff) %>%
  dplyr::summarize(prs=mean(SCORE)) %>%
  mutate(phenotype = fct_recode(PHENO, case = '1', control = '0')) %>%
  ungroup() %>%
  select(-PHENO) %>%
  spread(phenotype, prs) %>%
  mutate(case_greater=case>control)

r2_cohorts_eurasn <- r2_cohorts_eurasn %>%
  left_join(mean_snps, by=c('cohort', 'p_cutoff')) %>%
  left_join(mean_cc_score, by=c('cohort', 'p_cutoff')) %>%
  mutate(p_string = fct_recode(p_cutoff, '5e-8' = 'p5e8', '1e-6' = 'pe6', '1e-4' = 'pe4', '1e-3' = 'pe3', '1e-2' = 'pe2',
                               '0.05' = 'p05', '0.1' = 'p1', '0.2' = 'p2', '0.5' = 'p5', '1' = 'all'))

# write.table(t(spread(r2_cohorts_eurasn[,1:3], p_cutoff, r2)), 'eur_eas_r2.txt', col.names=F, quote=F)
# write.table(t(spread(r2_cohorts[,c(1,2,4)], p_cutoff, p)), 'eur_eas_p.txt', col.names=F, quote=F)
# write.table(t(spread(r2_cohorts[,c(1,2,5)], p_cutoff, CNT)), 'eur_eas_mean_snps.txt', col.names=F, quote=F)

eas_eas <- read.table('/Users/alicia/daly_lab/pgc_scz/eas_eas/eas_eas_r2.txt', header=T, sep='\t')
r2_cohorts_eurasn$P.value.threshold <- as.numeric(as.character(r2_cohorts_eurasn$p_string))
r2_cohorts_eurasn$pop='European SNPs, East Asian effects'
r2_eas_eas <- eas_eas %>%
  gather(cohort, r2, -P.value.threshold) %>%
  mutate(pop='East Asian')

r2_all <- bind_rows(r2_cohorts, r2_eas_eas, r2_cohorts_eurasn)
r2_all$P.value.threshold <- factor(r2_all$P.value.threshold)

library(plotrix)
r2_summary <- subset(r2_all, !P.value.threshold %in% c('1e-06', '5e-08')) %>%
  group_by(pop, P.value.threshold) %>%
  dplyr::summarize(r2_mean=mean(r2), sem = std.error(r2))

pd <- position_dodge(0.1) # move them .05 to the left and right
p_eas_eur <- ggplot(r2_summary, aes(x=P.value.threshold, y=r2_mean, color=pop)) +
  geom_point(position=pd, size=3) +
  geom_errorbar(aes(ymin=r2_mean-sem, ymax=r2_mean+sem), width=.1, position=pd)  +
  labs(y=expression(East~Asian~Prediction~R^2), x='p-value threshold') +
  scale_color_brewer(palette='Set1', name='Training data') +
  theme_classic() +
  theme(legend.position = c(0.2, 0.9),
        text = element_text(size=18),
        panel.background = element_rect(fill = "transparent",colour = NA),
        #panel.grid.minor = element_blank(), 
        #panel.grid.major = element_blank(),
        legend.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA))
ggsave('p_eas_eur.pdf', p_eas_eur, height=6, width=5)
