/*
 *  SiteScore
 */
package net.maizegenetics.dna.snp.score;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public abstract class SiteScore {

    public static enum SITE_SCORE_TYPE {

        None, QualityScore, ReferenceProbablity, Dosage,
        DepthA, DepthC, DepthG, DepthT, DepthGap, DepthInsertion,
        ProbA, ProbC, ProbG, ProbT, ProbGap, ProbInsertion
    };

    private final Map<SITE_SCORE_TYPE, Byte2D> myValues = new HashMap<>();
    protected final SITE_SCORE_TYPE myOnlyScoreType;
    private final int myNumTaxa;
    private final int myNumSites;

    SiteScore(Byte2D[] values) {
        if (values.length == 0) {
            throw new IllegalArgumentException("SiteScore: init: no values provided.");
        }
        int numTaxa = values[0].numTaxa();
        int numSites = values[0].numSites();
        for (int i = 0; i < values.length; i++) {
            if ((numTaxa != values[i].numTaxa()) || (numSites != values[i].numSites())) {
                throw new IllegalArgumentException("SiteScore: init: number of taxa or sites don't match for all values.");
            }
            myValues.put(values[i].siteScoreType(), values[i]);
        }
        myNumTaxa = numTaxa;
        myNumSites = numSites;
        if (values.length == 1) {
            myOnlyScoreType = values[0].siteScoreType();
        } else {
            myOnlyScoreType = null;
        }
    }

    /**
     * Return the site scores types.
     *
     * @return site score types.
     */
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        return myValues.keySet();
    }

    public int numTaxa() {
        return myNumTaxa;
    }

    public int numSites() {
        return myNumSites;
    }

    public int numAlleles() {
        return myValues.size();
    }

    protected Byte2D byteStorage(SITE_SCORE_TYPE type) {
        return myValues.get(type);
    }

    Collection<Byte2D> byteStorage() {
        return myValues.values();
    }
}
