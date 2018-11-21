library(snp.plotter)

setwd('/Users/alicia/daly_lab/Travel/23andMe')
setwd('/Users/alicia/daly_lab/pgc_scz/ld_fig')

date(); snp.plotter(config.file = 'EUR.config'); date()
date(); snp.plotter(config.file = 'EAS.config'); date()

date(); snp.plotter(config.file = 'EUR_10kb.config'); date()
date(); snp.plotter(config.file = 'EAS_10kb.config'); date()
