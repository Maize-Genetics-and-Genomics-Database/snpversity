/*
 *  Dosage
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class Dosage extends SiteScore {

    private final Byte2D myStorage;

    Dosage(Byte2D value) {
        super(new Byte2D[]{value});
        myStorage = value;
    }

    public byte value(int taxon, int site) {
        return myStorage.valueForAllele(taxon, site);
    }

}
