# VCF file transformation

## Overview

Simple transformation of VCF file to BED-like format with the addutional `patient` field

Data reshaping using `pandas.DataFrame`

## Similar scripts

* [`pyvcf`](https://github.com/jamescasbon/PyVCF) seems to have a [similar VCF melting script](https://github.com/jamescasbon/PyVCF/blob/master/scripts/vcf_melt). However, it does not exclude samples without calls.
* BEDOPS's [`vcf2bed`](http://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/vcf2bed.html) converts VCF to BED format
* [gtamazian/bioformats](https://github.com/gtamazian/bioformats) has the function [`variants.vcf2bed`](https://github.com/gtamazian/bioformats/blob/master/bioformats/variants.py#L156)
* [erscott/pandasVCF](https://github.com/erscott/pandasVCF) parses a VCF file into a `pandas.DataFrame` object

## Other links

* http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it
* https://github.com/jamescasbon/PyVCF/issues/89