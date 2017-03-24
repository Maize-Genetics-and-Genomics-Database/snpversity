/*
 *  ReferenceProbability
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class ReferenceProbability extends SiteScore {

    private final Byte2D myStorage;

    ReferenceProbability(Byte2D value) {
        super(new Byte2D[]{value});
        myStorage = value;
    }

    public float value(int taxon, int site) {
        return SiteScoreUtil.byteToFloatPercentage(myStorage.valueForAllele(taxon, site));
    }

}
