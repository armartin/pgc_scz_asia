"""
take combinations of cohorts and write metal scripts to perform inverse variance-weighted meta-analysis
"""

import argparse
import gzip

def write_cohort_info(out, path):
    out.write('MARKER\tSNP\n')
    #out.write('WEIGHT\tN\n')
    out.write('ALLELE\tA1 A2\n')
    #out.write('FREQ\tFRQ\n')
    out.write('EFFECT\tlog(OR)\n')
    out.write('STDERR\tSE\n')
    out.write('PVAL\tP\n')
    out.write('PROCESS ' + path + '\n\n')
        
        
def main(args):
    paths = open(args.paths)
    file_shorthand = []
    file_paths = []
    for line in paths:
        file_shorthand.append(line.strip().split('/')[-1].split()[0].split('scz_')[1].split('_eur')[0])
        file_paths.append(line.strip())
    print(file_shorthand)
    print(file_paths)
    for i in range(len(file_paths)):
        print(file_paths[i])
        out = open('/Users/alicia/daly_lab/pgc_scz/cross_pop/loo_' + file_shorthand[i], 'w')
        out.write('SCHEME STDERR\n\n')
        for j in range(len(file_paths)):
            if i != j:
                write_cohort_info(out, file_paths[j])

        out.write('OUTFILE /home/armartin/pgc_scz_asia/imp_collect/prs_score/armartin/cross_pop/eur_loometa/loo_' + file_shorthand[i] + '.meta .tbl\n')
        out.write('ANALYZE HETEROGENEITY\n')
        out.close()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paths')
    parser.add_argument('--dirname')
    args = parser.parse_args()
    main(args)
