/*
 *  AlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

/**
 * Provides DNA read depth data for a genotype table.
 *
 * @author Terry Casstevens
 */
public interface AlleleDepth {

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     * The array of depths is sized as determined by NUMBER_NUCLEOTIDE_ALLELES
     * (6), and it is ordered as in NUCLEOTIDE_ALLELES (A, C, G, T,
     * +(insertion), - (deletion))
     *
     * The depths from 0 to 127 are recorded exactly. Depths above 127 are
     * stored as approximate logs (equation), and they are represented by
     * negative numbers.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return array of counts
     */
    public int[] depthForAlleles(int taxon, int site);

    /**
     * Returns allele depth for given taxon, site, and allele.
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     *
     * @return allele depth
     */
    public int depthForAllele(int taxon, int site, int allele);

    /**
     * Returns allele depth for given taxon, site, and allele in byte format
     * (negative values need to be decoded for real value).
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     *
     * @return allele depth byte encoded.
     */
    public byte depthForAlleleByte(int taxon, int site, int allele);

    /**
     * Returns allele depths for all alleles and sites for a taxon in byte
     * format (negative values need to be decoded for real value).
     *
     * @param taxon taxon
     *
     * @return allele depths byte encoded
     */
    public byte[][] depthAllSitesByte(int taxon);

    /**
     * Returns total depth for given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return total depth
     */
    public int depth(int taxon, int site);

    /**
     * Returns total depth of all alleles and sites for given taxon.
     *
     * @param taxon taxon
     *
     * @return total depth
     */
    public int depthForTaxon(int taxon);

    /**
     * Returns total depth of all alleles and taxa for given site.
     *
     * @param site site
     *
     * @return depth
     */
    public int depthForSite(int site);

}
