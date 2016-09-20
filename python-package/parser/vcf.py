import numpy as np

import vcfnp

def vcf2snp(vcf_filepath, missing=3, cache=True):
    """
        Return a SNP matrix based on a VCF file.

        This take a VCF file and create a SNP matrix. It will keep only the SNP with 2
        variants. This function is based on the vcfnp package.

        :param vcf_filepath: The path of the VCF file
        :param missing: How to annotate missing data
        :param cache: Use cache
        :type vcf_filepath: string
        :type missing: np.int8
        :type cache: boolean
        :return: The genotype matrix containing SNP
        :rtype: np.array of np.int8

        :Example:

        >>> G = vcf2snp('file.vcf')

        .. warnings:: This function is not maintain.
    """

    c = vcfnp.calldata_2d(vcf_filepath, cache=cache).view(np.recarray)
    G = c.genotype

    ## create a mask to keep only 0/0, 1/0, 0/1, 1/1 and missing datas
    mask = np.logical_and.reduce(np.logical_and(G >= -1, G <= 1), axis = 2)
    mask = np.logical_and.reduce(mask, axis=1)

    G = G[mask, :]
    G = np.sum(G.T, axis=0, dtype=np.int8)

    G[G == -2] = missing

    return G