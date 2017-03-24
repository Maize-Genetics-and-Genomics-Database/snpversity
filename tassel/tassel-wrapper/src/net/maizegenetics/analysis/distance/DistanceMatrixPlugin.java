/*
 * DistanceMatrixPlugin
 */
package net.maizegenetics.analysis.distance;

import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import javax.swing.*;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.ProgressListener;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class DistanceMatrixPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(DistanceMatrixPlugin.class);

    public DistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
            if (alignInList.size() < 1) {
                String message = "Invalid selection.  Please select sequence alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error(message);
                }
                return null;
            }

            List<DataSet> result = new ArrayList<>();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                Datum current = itr.next();
                DataSet tds = processDatum(current);
                result.add(tds);
                if (tds != null) {
                    fireDataSetReturned(new PluginEvent(tds, DistanceMatrixPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public DataSet processDatum(Datum input) {
        GenotypeTable aa = (GenotypeTable) input.getData();
        DistanceMatrix adm = getDistanceMatrix(aa, this);
        return new DataSet(new Datum("Matrix:" + input.getName(), adm, "Distance Matrix"), this);
    }

    public static DistanceMatrix getDistanceMatrix(GenotypeTable alignment) {
        return IBSDistanceMatrix3Alleles.getInstance(alignment);
    }

    public static DistanceMatrix getDistanceMatrix(GenotypeTable alignment, ProgressListener listener) {
        return IBSDistanceMatrix3Alleles.getInstance(alignment, listener);
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = DistanceMatrixPlugin.class.getResource("/net/maizegenetics/analysis/images/DistanceMatrix.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Distance Matrix";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Create a distance matrix";
    }
}
