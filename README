FastQTL
---------

Modified version of FastQTL (Ongen, H. et al., 2016) that handles structural variants

[Ongen, H., Buil, A., Brown, A. A., Dermitzakis, E. T. & Delaneau, O. Fast and efficient QTL mapper for thousands of molecular phenotypes. Bioinformatics 32, 1479–1485 (2016).](http://bioinformatics.oxfordjournals.org/content/32/10/1479)

Documentation for FastQTL can be found here: [http://fastqtl.sourceforge.net](http://fastqtl.sourceforge.net)

#### Description

A structural variant (SV) lies within the cis window if:
* any part of the spanned region is within the cis window *except* inversions, for which one (or both) of the breakpoints must fall are within the cis window

The spanned region for deletions, duplications, and inversions is inferred from the END field in the 8th column of the VCF. If the END field is absent for a variant, it assumes that the END is the same as the POS field (identical to standard FastQTL behavior).

The variant type is inferred from the SVTYPE field of the 8th column of the VCF. Inversions must be designated as SVTYPE=INV for proper behavior
