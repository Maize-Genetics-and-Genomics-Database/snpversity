/*
 *  MaskGenotypeCallTable
 * 
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.util.BitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable myBase;
    private final byte myMaskValue;
    private final BitSet[] myBitSets;

    public MaskGenotypeCallTable(GenotypeCallTable base, byte maskValue, BitSet[] bitSets) {
        super(base.numberOfTaxa(), base.numberOfSites(), base.isPhased(), base.alleleDefinitions());
        myBase = base;
        myMaskValue = maskValue;
        myBitSets = bitSets;
    }

    @Override
    public byte genotype(int taxon, int site) {
        if (myBitSets[taxon].fastGet(site)) {
            return myMaskValue;
        } else {
            return myBase.genotype(taxon, site);
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myBase.diploidAsString(site, genotype(taxon, site));
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        return myBase.genotypeAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myBase.diploidAsString(site, value);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBase.transposeData(siteInnerLoop);
    }

}
