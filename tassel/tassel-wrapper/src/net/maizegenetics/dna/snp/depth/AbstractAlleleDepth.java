/*
 *  AbstractAlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractAlleleDepth implements AlleleDepth {

    private final int myNumAlleles;
    private final int myNumTaxa;
    private final int myNumSites;

    AbstractAlleleDepth(int numAlleles, int numTaxa, int numSites) {
        myNumAlleles = numAlleles;
        myNumTaxa = numTaxa;
        myNumSites = numSites;
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        int[] result = new int[myNumAlleles];
        for (int a = 0; a < myNumAlleles; a++) {
            result[a] = depthForAllele(taxon, site, a);
        }
        return result;
    }

    @Override
    public int depthForAllele(int taxon, int site, int allele) {
        return AlleleDepthUtil.depthByteToInt(depthForAlleleByte(taxon, site, allele));
    }

    @Override
    public int depth(int taxon, int site) {
        int result = 0;
        for (int a = 0; a < myNumAlleles; a++) {
            result += depthForAllele(taxon, site, a);
        }
        return result;
    }

    @Override
    public int depthForTaxon(int taxon) {
        int result = 0;
        for (int site = 0; site < myNumSites; site++) {
            result += depth(taxon, site);
        }
        return result;
    }

    @Override
    public int depthForSite(int site) {
        int result = 0;
        for (int taxon = 0; taxon < myNumTaxa; taxon++) {
            result += depth(taxon, site);
        }
        return result;
    }

    @Override
    public byte[][] depthAllSitesByte(int taxon) {
        byte[][] result = new byte[myNumAlleles][myNumSites];
        for (int allele = 0; allele < myNumAlleles; allele++) {
            for (int site = 0; site < myNumSites; site++) {
                result[allele][site] += depthForAllele(taxon, site, allele);
            }
        }
        return result;
    }
}
