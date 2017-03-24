package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;
import java.util.List;

public class RefProbAdditiveSite extends AbstractAdditiveSite {

    private static final long serialVersionUID = 2040665024409852166L;
    private int ntaxa;
    private float[] cov;
    private int[] taxaIndex = null;

    public RefProbAdditiveSite(int site, String chr, int pos, String id,
            CRITERION selectionCriteria, float[] covariate) {
        super(site, chr, pos, id, selectionCriteria);
        cov = covariate;
        ntaxa = cov.length;
    }

    @Override
    public double[] getCovariate() {
        double[] dcov = new double[ntaxa];
        if (taxaIndex == null) {
            for (int i = 0; i < ntaxa; i++)
                dcov[i] = cov[i];
            return dcov;
        } else {
            for (int i = 0; i < ntaxa; i++)
                dcov[i] = cov[taxaIndex[i]];
            return dcov;
        }
    }

    @Override
    public double[] getCovariate(int[] subset) {
        if (taxaIndex == null)
            return Arrays.stream(subset).mapToDouble(i -> cov[i]).toArray();
        return Arrays.stream(subset).mapToDouble(i -> cov[taxaIndex[i]]).toArray();
    }

    @Override
    public void reindexTaxa(int[] taxaIndex, List<Integer> uniqueTaxa) {
        this.taxaIndex = taxaIndex;
    }

}
