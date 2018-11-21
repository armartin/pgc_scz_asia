zcat daner_PGC_SCZ52.p3.0215a.euro | \
awk '{ if ( (length($4)==1 && length($5)==1 && $6>0.05 && $7>0.05 && $6<0.95 && $7<0.95 && $9>0.9) )
	{ print $2 }
}' > daner_PGC_SCZ52.p3.0215a.euro.maf05.info90.snps2

#	($1==6 && $3 < 25000000 || $1==6 && $3 > 35000000 || $1!=6) )

