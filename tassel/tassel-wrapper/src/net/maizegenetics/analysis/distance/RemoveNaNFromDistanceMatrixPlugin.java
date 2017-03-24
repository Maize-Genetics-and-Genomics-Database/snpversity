/*
 *  RemoveNaNFromDistanceMatrixPlugin
 * 
 *  Created on Oct 5, 2015
 */
package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class RemoveNaNFromDistanceMatrixPlugin extends AbstractPlugin {

    public RemoveNaNFromDistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> temp = input.getDataOfType(DistanceMatrix.class);
        if (temp.size() != 1) {
            throw new IllegalArgumentException("RemoveNaNFromDistanceMatrixPlugin: processData: Must select one Distance Matrix");
        }

        DistanceMatrix origMatrix = (DistanceMatrix) temp.get(0).getData();
        String origName = temp.get(0).getName();

        TaxaList taxa = origMatrix.getTaxaList();
        int numTaxa = taxa.numberOfTaxa();
        BitSet[] whichNaN = new BitSet[numTaxa];

        for (int t = 0; t < numTaxa; t++) {
            whichNaN[t] = new OpenBitSet(numTaxa);
        }

        for (int x = 0; x < numTaxa; x++) {
            for (int y = x; y < numTaxa; y++) {
                if (Double.isNaN(origMatrix.getDistance(x, y))) {
                    whichNaN[x].fastSet(y);
                    whichNaN[y].fastSet(x);
                }
            }
        }

        List<Integer> taxaToRemove = new ArrayList<>();
        long highestCount = -1;
        while (highestCount != 0) {
            highestCount = 0;
            int highestTaxon = -1;
            for (int t = 0; t < numTaxa; t++) {
                long currentCount = whichNaN[t].cardinality();
                if (currentCount > highestCount) {
                    highestCount = currentCount;
                    highestTaxon = t;
                }
            }
            if (highestCount != 0) {
                taxaToRemove.add(highestTaxon);
                for (int t = 0; t < numTaxa; t++) {
                    whichNaN[t].fastClear(highestTaxon);
                    whichNaN[highestTaxon].fastClear(t);
                }
            }
        }

        if (!taxaToRemove.isEmpty()) {

            int newNumTaxa = numTaxa - taxaToRemove.size();
            int[] origIndex = new int[newNumTaxa];
            TaxaListBuilder builder = new TaxaListBuilder();
            int count = 0;
            for (int t = 0; t < numTaxa; t++) {
                if (!taxaToRemove.contains(t)) {
                    builder.add(taxa.get(t));
                    origIndex[count++] = t;
                }
            }
            TaxaList newTaxaList = builder.build();

            double[][] result = new double[newNumTaxa][newNumTaxa];
            for (int x = 0; x < newNumTaxa; x++) {
                for (int y = x; y < newNumTaxa; y++) {
                    result[x][y] = result[y][x] = origMatrix.getDistance(origIndex[x], origIndex[y]);
                }
            }

            DistanceMatrix resultMatrix = new DistanceMatrix(result, newTaxaList);
            return new DataSet(new Datum(origName + " with no NaN", resultMatrix, null), this);

        } else {
            return input;
        }

    }

    public static DataSet runPlugin(DataSet input) {
        return new RemoveNaNFromDistanceMatrixPlugin(null, false).processData(input);
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Remove NaN From Distance Matrix";
    }

    @Override
    public String getToolTipText() {
        return "Remove NaN From Distance Matrix";
    }

}
