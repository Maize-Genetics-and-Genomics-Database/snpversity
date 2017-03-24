/*
 * FilterTaxaAlignmentPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.gui.SelectFromAvailableDialog;
import net.maizegenetics.gui.TaxaAvailableListModel;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.phenotype.FilterPhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class FilterTaxaAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterTaxaAlignmentPlugin.class);
    private TaxaList myIdsToKeep = null;
    private TaxaList myIdsToRemove = null;

    /**
     * Creates a new instance of FilterTaxaAlignmentPlugin
     */
    public FilterTaxaAlignmentPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List inputData = input.getDataSet();
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                } else {
                    myLogger.error("performFunction: Please input a single sequence or phenotype.");
                }
                return null;
            }

            Datum inputDatum = (Datum) inputData.get(0);

            if (!(inputDatum.getData() instanceof GenotypeTable) && !(inputDatum.getData() instanceof Phenotype)) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                } else {
                    myLogger.error("performFunction: Please input a single sequence or phenotype.");
                }
                return null;
            }

            Datum td = processDatum(inputDatum, isInteractive());
            if (td == null) {
                return null;
            }

            DataSet output = new DataSet(td, this);
            //I am setting the firing class as the metadata - so that the control panel know where the event is coming from
            fireDataSetReturned(new PluginEvent(output, FilterTaxaAlignmentPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }
    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {

        Object theData = inDatum.getData();

        if (isInteractive) {
            TaxaList taxaList = null;
            SelectFromAvailableDialog dialog = null;
            if (theData instanceof GenotypeTable) {
                final GenotypeTable alignment = (GenotypeTable) theData;
                taxaList = alignment.taxa();
                TaxaAvailableListModel listModel = new TaxaAvailableListModel(taxaList);
                dialog = new SelectFromAvailableDialog(getParentFrame(), "Taxa Filter", listModel);
            } else if (theData instanceof Phenotype) {
                final Phenotype phenotype = (Phenotype) theData;
                taxaList = phenotype.taxa();
                TaxaAvailableListModel listModel = new TaxaAvailableListModel(taxaList);
                dialog = new SelectFromAvailableDialog(getParentFrame(), "Taxa Filter", listModel);
            } else {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                return null;
            }
            dialog.setLocationRelativeTo(getParentFrame());
            dialog.setVisible(true);
            if (dialog.isCanceled()) {
                return null;
            }
            int[] indicesToKeep = dialog.getDesiredIndices();
            Taxon[] ids = new Taxon[indicesToKeep.length];
            for (int i = 0; i < indicesToKeep.length; i++) {
                ids[i] = taxaList.get(indicesToKeep[i]);
            }
            myIdsToKeep = new TaxaListBuilder().addAll(ids).build();
            dialog.dispose();
        }

        if (((myIdsToKeep == null) || (myIdsToKeep.numberOfTaxa() == 0))
                && ((myIdsToRemove == null) || (myIdsToRemove.numberOfTaxa() == 0))) {
            return null;
        }

        Object result = null;
        long count = 0;
        if (theData instanceof GenotypeTable) {
            if (myIdsToKeep != null) {
                result = FilterGenotypeTable.getInstance((GenotypeTable) theData, myIdsToKeep, false);
            } else if (myIdsToRemove != null) {
                result = FilterGenotypeTable.getInstanceRemoveIDs((GenotypeTable) theData, myIdsToRemove);
            }
            count = ((GenotypeTable) result).numberOfTaxa();
        } else if (theData instanceof Phenotype) {
            if (myIdsToKeep != null) {
                result = new PhenotypeBuilder().fromPhenotype((Phenotype) theData).keepTaxa(myIdsToKeep).build();
            } else if (myIdsToRemove != null) {
                result = new PhenotypeBuilder().fromPhenotype((Phenotype) theData).removeTaxa(myIdsToRemove).build();
            }
            count = ((FilterPhenotype) result).getRowCount();
        } else {
            myLogger.error("processDatum: Please input a single sequence or phenotype.  Unknown data type: " + theData.getClass().getName());
            return null;
        }

        String theName, theComment;
        theName = inDatum.getName() + "_" + count + " Rows";
        theComment = "Subset of " + count + " Taxa\n" + inDatum.getComment();
        return new Datum(theName, result, theComment);

    }

    public TaxaList getIdsToKeep() {
        return myIdsToKeep;
    }

    public void setIdsToKeep(TaxaList idsToKeep) {
        myIdsToKeep = idsToKeep;
    }

    public void setIdsToRemove(TaxaList idsToRemove) {
        myIdsToRemove = idsToRemove;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterTaxaAlignmentPlugin.class.getResource("/net/maizegenetics/analysis/images/Filter_horizontal.gif");
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
    public String getButtonName() {
        return "Taxa Names";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Select Taxa Within Dataset";
    }
}
