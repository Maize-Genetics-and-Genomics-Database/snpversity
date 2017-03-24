/*
 *  FilterByte2D
 */
package net.maizegenetics.dna.snp.byte2d;

import net.maizegenetics.dna.snp.FilterGenotypeTable;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class FilterByte2D extends AbstractByte2D {

    private final FilterGenotypeTable myFilterGenotypeTable;
    private final Byte2D myBase;

    FilterByte2D(Byte2D base, FilterGenotypeTable filterGenotypeTable) {
        super(base.siteScoreType(), filterGenotypeTable.numberOfTaxa(), filterGenotypeTable.numberOfSites());
        myBase = base;
        myFilterGenotypeTable = filterGenotypeTable;
    }

    @Override
    public byte valueForAllele(int taxon, int site) {
        return myBase.valueForAllele(myFilterGenotypeTable.translateTaxon(taxon),
                myFilterGenotypeTable.translateSite(site));
    }
}
