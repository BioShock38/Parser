import numpy as np

import vcfnp

def vcf2snp(filename, missing=3, cache=True):
    """
    Return a SNP matrix based on a VCF file.

    This take a VCF file and create a SNP matrix. It will keep only the SNP with 2
    variants. This function is based on the vcfnp package.

    :param filename: The path of the VCF file
    :param missing: How to annotate missing data
    :param cache: Use cache
    :type filename: string
    :type missing: np.int8
    :type cache: boolean
    :return: The genotype matrix containing SNP
    :rtype: np.array of np.int8

    :Example:

    >>> G = vcf2snp('file.vcf')

    ... warnings:: This function is not maintain.
    """

    c = vcfnp.calldata_2d(filename, cache=cache).view(np.recarray)
    G = c.genotype

    ## create a mask to keep only 0/0, 1/0, 0/1, 1/1 and missing datas
    mask = np.logical_and.reduce(np.logical_and(G >= -1, G <= 1), axis = 2)
    mask = np.logical_and.reduce(mask, axis=1)

    G = G[mask, :]
    mask_missing = np.logical_and.reduce(G == -1, axis=2)
    G = np.sum(G.T, axis=0, dtype=np.int8)

    G[mask_missing.T] = missing

    return G