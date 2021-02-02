from subprocess import call


def get_options():
    import argparse

    description = 'converts multi-allele to biallelic variants and annotate'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--vcf',
                        help='Complete path  to the vcf file ')
    parser.add_argument('--out',
                        help='Complete path to output (sims.vcf)')

    return parser.parse_args()

options = get_options()

                                
def vcfReformater(vcf,out): #bcftools annotate -x removes the column you do not need (refer to manual for full detail)
            CallString=' bcftools norm -Ou -m -any  %s | \
             bcftools annotate -x ID -I %s -Ov -o %s' % \
            (vcf,'+%CHROM:%POS:%REF:%ALT',out)
            call(CallString,shell=True)
            
vcfReformater(options.vcf,options.out)