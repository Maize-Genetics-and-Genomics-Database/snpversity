/*
 *  MaskGenotypeTableBuilder
 * 
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.genotypecall.MaskGenotypeCallTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskGenotypeTableBuilder {
    
    private final byte DEFAULT_MASK_VALUE = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;

    private final GenotypeTable myBase;
    private final BitSet[] myBitSets;
    private byte myMaskValue = DEFAULT_MASK_VALUE;

    public MaskGenotypeTableBuilder(GenotypeTable base) {
        myBase = base;
        int numSites = myBase.numberOfSites();
        int numTaxa = myBase.numberOfTaxa();
        myBitSets = new BitSet[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            myBitSets[t] = new OpenBitSet(numSites);
        }
    }
    
    public void setMaskValue (byte value) {
        myMaskValue = value;
    }
    
    public void set (int taxon, int site) {
        myBitSets[taxon].fastSet(site);
    }
    
    public GenotypeTable build() {
        return new CoreGenotypeTable(
            new MaskGenotypeCallTable(myBase.genotypeMatrix(), myMaskValue, myBitSets), 
            myBase.positions(),
            myBase.taxa(),
            myBase.depth(),
            myBase.alleleProbability(),
            myBase.referenceProbability(),
            myBase.dosage(),
            myBase.annotations()
        );
    }

}
