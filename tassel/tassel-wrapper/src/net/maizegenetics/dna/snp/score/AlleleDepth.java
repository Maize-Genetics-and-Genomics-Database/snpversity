/*
 *  AlleleDepth
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleDepth extends SiteScore {

    public static final SiteScore.SITE_SCORE_TYPE[] ALLELE_DEPTH_TYPES = new SiteScore.SITE_SCORE_TYPE[]{
        SiteScore.SITE_SCORE_TYPE.DepthA, SiteScore.SITE_SCORE_TYPE.DepthC,
        SiteScore.SITE_SCORE_TYPE.DepthG, SiteScore.SITE_SCORE_TYPE.DepthT,
        SiteScore.SITE_SCORE_TYPE.DepthGap, SiteScore.SITE_SCORE_TYPE.DepthInsertion};

    AlleleDepth(Byte2D[] values) {
        super(values);
    }

    /**
     * Returns the depth of nucleotide (scoreType) at given taxon and site.
     * Depth values are stored in bytes and translated to integer using
     * AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     * @param scoreType nucleotide
     *
     * @return depth
     */
    public int value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        return AlleleDepthUtil.depthByteToInt(byteStorage(scoreType).valueForAllele(taxon, site));
    }

    /**
     * Returns the depth values of all nucleotides at given taxon and site.
     * Depth values are stored in bytes and translated to integer using
     * AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     *
     * @return depths
     */
    public int[] values(int taxon, int site) {
        int[] result = new int[ALLELE_DEPTH_TYPES.length];
        int count = 0;
        for (SITE_SCORE_TYPE current : ALLELE_DEPTH_TYPES) {
            result[count++] = value(taxon, site, current);
        }
        return result;
    }

    /**
     * Returns depth values of all nucleotides and sites for given taxon. The
     * first dimension of returned array is nucleotides (ALLELE_DEPTH_TYPES) and
     * second dimension is sites.
     *
     * @param taxon taxon
     *
     * @return depths
     */
    public int[][] values(int taxon) {
        int[][] result = new int[ALLELE_DEPTH_TYPES.length][numSites()];
        int count = 0;
        for (SITE_SCORE_TYPE current : ALLELE_DEPTH_TYPES) {
            for (int site = 0; site < numSites(); site++) {
                result[count][site] = value(taxon, site, current);
            }
            count++;
        }
        return result;
    }

    /**
     * Returns the depth (byte representation) of nucleotide (scoreType) at
     * given taxon and site. Depth values are stored in bytes and translated to
     * integer using AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     * @param scoreType nucleotide
     *
     * @return depth
     */
    public byte valueByte(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        return byteStorage(scoreType).valueForAllele(taxon, site);
    }

    /**
     * Returns the depth values (byte representation) of all nucleotides at
     * given taxon and site. Depth values are stored in bytes and translated to
     * integer using AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     *
     * @return depths
     */
    public byte[] valuesByte(int taxon, int site) {
        byte[] result = new byte[ALLELE_DEPTH_TYPES.length];
        int count = 0;
        for (SITE_SCORE_TYPE current : ALLELE_DEPTH_TYPES) {
            result[count++] = valueByte(taxon, site, current);
        }
        return result;
    }

    /**
     * Returns depth values (byte representation) of all nucleotides and sites
     * for given taxon. The first dimension of returned array is nucleotides
     * (ALLELE_DEPTH_TYPES) and second dimension is sites.
     *
     * @param taxon taxon
     *
     * @return depths
     */
    public byte[][] valuesForTaxonByte(int taxon) {
        byte[][] result = new byte[ALLELE_DEPTH_TYPES.length][numSites()];
        int count = 0;
        for (SITE_SCORE_TYPE current : ALLELE_DEPTH_TYPES) {
            for (int site = 0; site < numSites(); site++) {
                result[count][site] = valueByte(taxon, site, current);
            }
            count++;
        }
        return result;
    }

    /**
     * Returns depth values (byte representation) of all nucleotides and taxa
     * for given site. The first dimension of returned array is nucleotides
     * (ALLELE_DEPTH_TYPES) and second dimension is taxa.
     *
     * @param site site
     *
     * @return depths
     */
    public byte[][] valuesForSiteByte(int site) {
        byte[][] result = new byte[ALLELE_DEPTH_TYPES.length][numTaxa()];
        int count = 0;
        for (SITE_SCORE_TYPE current : ALLELE_DEPTH_TYPES) {
            for (int taxon = 0; taxon < numTaxa(); taxon++) {
                result[count][taxon] = valueByte(taxon, site, current);
            }
            count++;
        }
        return result;
    }

    /**
     * Returns sum of all nucleotide depths at given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return sum of depths
     */
    public int depth(int taxon, int site) {
        int result = 0;
        for (SITE_SCORE_TYPE current : ALLELE_DEPTH_TYPES) {
            result += value(taxon, site, current);
        }
        return result;
    }

    /**
     * Returns sum of all nucleotide depths across all sites at given taxon.
     *
     * @param taxon taxon
     *
     * @return sum of depths
     */
    public int depthForTaxon(int taxon) {
        int result = 0;
        for (int site = 0; site < numSites(); site++) {
            result += depth(taxon, site);
        }
        return result;
    }

    /**
     * Returns sum of all nucleotide depths across all taxa at given site.
     *
     * @param site site
     *
     * @return sum of depths
     */
    public int depthForSite(int site) {
        int result = 0;
        for (int taxon = 0; taxon < numTaxa(); taxon++) {
            result += depth(taxon, site);
        }
        return result;
    }

}
