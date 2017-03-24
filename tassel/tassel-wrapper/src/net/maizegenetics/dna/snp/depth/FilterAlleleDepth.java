/*
 *  MemoryAlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

import net.maizegenetics.dna.snp.FilterGenotypeTable;

/**
 * In memory allele depth class. In memory allele depth occupies 6-fold more
 * memory than the allele calls.
 *
 * @author Ed Buckler
 */
public class FilterAlleleDepth extends AbstractAlleleDepth {

    private final FilterGenotypeTable myFilterGenotypeTable;
    private final AlleleDepth myBaseAlleleDepth;

    public FilterAlleleDepth(AlleleDepth baseAlleleDepth, FilterGenotypeTable filterGenotypeTable) {
        super(filterGenotypeTable.maxNumAlleles(), filterGenotypeTable.numberOfTaxa(), filterGenotypeTable.numberOfSites());
        myBaseAlleleDepth = baseAlleleDepth;
        myFilterGenotypeTable = filterGenotypeTable;
    }

    @Override
    public byte depthForAlleleByte(int taxon, int site, int allele) {
        return myBaseAlleleDepth.depthForAlleleByte(myFilterGenotypeTable.translateTaxon(taxon),
                myFilterGenotypeTable.translateSite(site), allele);
    }
}
