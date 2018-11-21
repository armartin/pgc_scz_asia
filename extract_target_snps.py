import argparse
import gzip
import numpy as np
from datetime import datetime


def open_file(filename):
    if filename.endswith('.gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return(my_file)


def current_time():
    return(' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']')


def main(args):
    sumstats = open_file(args.sumstats)
    header = sumstats.readline().strip().split()
    freq_col = header.index(args.freq)
    info_col = header.index('INFO')
    a1_col = header.index('A1')
    a2_col = header.index('A2')
    snp_col = header.index('SNP')

    # make list of snps to keep with maf > 0.05, info > 0.9, no indels in original summary stats files
    print "QC original sumstats" + current_time()
    keep_snps = set()
    for line in sumstats:
        line = line.strip().split()
        try:
            if 0.05 < float(line[freq_col]) < 0.95 and float(line[info_col]) > 0.9 and len(line[a1_col]) == 1 and len(line[a2_col]) == 1:
                keep_snps.add(line[snp_col])
        except IndexError:
            print [freq_col, info_col, a1_col, a2_col, line]
    print "Num SNPs: " + str(len(keep_snps))

    # filter to snps present in all bim files
    print "QC all bim files" + current_time()
    all_bim_snps = set()
    bims = open(args.bim_list)
    for f in bims:
        current_set = set()
        my_file = open(f.strip())
        for line in my_file:
            line = line.strip().split()
            current_set.add(line[1])
        if len(all_bim_snps) > 0:
            all_bim_snps = all_bim_snps.intersection(current_set)
        else:
            all_bim_snps = current_set
        print [f.strip(), len(all_bim_snps), len(current_set)]
    print "Num SNPs: " + str(len(all_bim_snps))

    # get all snps in mtag file
    print "QC mtag sumstats" + current_time()
    mtag = open_file(args.mtag)
    header = mtag.readline().strip().split()
    mtag_snp_col = header.index('SNP')
    mtag_snps = set()
    for line in mtag:
        line = line.strip().split()
        mtag_snps.add(line[mtag_snp_col])
    print "Num SNPs: " + str(len(mtag_snps))

    print [len(keep_snps), len(all_bim_snps), len(mtag_snps)]
    final_snps = keep_snps.intersection(all_bim_snps)
    final_snps = final_snps.intersection(mtag_snps)
    print "Num SNPs: " + str(len(final_snps))

    print "Writing output"
    out = open(args.out, 'w')
    for snp in final_snps:
        out.write(snp + '\n')
    out.close()
    print 'Done!' + current_time()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mtag')
    parser.add_argument('--sumstats')
    parser.add_argument('--freq')
    parser.add_argument('--bim_list')
    parser.add_argument('--out')

    args = parser.parse_args()
    main(args)
