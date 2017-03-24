/*
 *  AlleleProbability
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleProbability extends SiteScore {

    public static final SiteScore.SITE_SCORE_TYPE[] ALLELE_PROBABILITY_TYPES = new SiteScore.SITE_SCORE_TYPE[]{
        SiteScore.SITE_SCORE_TYPE.ProbA, SiteScore.SITE_SCORE_TYPE.ProbC,
        SiteScore.SITE_SCORE_TYPE.ProbG, SiteScore.SITE_SCORE_TYPE.ProbT,
        SiteScore.SITE_SCORE_TYPE.ProbGap, SiteScore.SITE_SCORE_TYPE.ProbInsertion};

    AlleleProbability(Byte2D[] values) {
        super(values);
    }

    public float value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        return SiteScoreUtil.byteToFloatPercentage(byteStorage(scoreType).valueForAllele(taxon, site));
    }

}
