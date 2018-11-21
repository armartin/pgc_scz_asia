zcat /home/gwas/pgc_scz_asia/meta_analysis/daner_pgc_scz_asia_meta_run27_replic_run3.meta.gz | \
awk '{ if ( (length($4)==1 && length($5)==1 && $6>0.05 && $7>0.05 && $6<0.95 && $7<0.95 && $9>0.9) )
	{ print $2 }
}' > daner_pgc_scz_asia_meta_run27_replic_run3.maf05.info90.snps2

#	($1==6 && $3 < 25000000 || $1==6 && $3 > 35000000 || $1!=6) )

