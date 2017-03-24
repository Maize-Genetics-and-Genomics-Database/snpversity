/*
 * GetTaxaListPlugin.java
 *
 * Created on June 16, 2014
 *
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;

import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class GetTaxaListPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(GetTaxaListPlugin.class);

    public GetTaxaListPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet processData(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

            if (alignInList.size() != 1) {
                String gpMessage = "Invalid selection.  Please select one genotype alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }

            Datum current = alignInList.get(0);
            GenotypeTable genotypeTable = (GenotypeTable) current.getData();
            String name = current.getName();

            Datum outputDatum = new Datum(name + "_TaxaList", genotypeTable.taxa(), "Taxa List from " + name);
            DataSet output = new DataSet(outputDatum, this);

            return output;

        } finally {
            fireProgress(100);
        }

    }

    public String getToolTipText() {
        return "Get Taxa List";
    }

    public ImageIcon getIcon() {
        URL imageURL = GetTaxaListPlugin.class.getResource("/net/maizegenetics/analysis/images/lists.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Get Taxa List";
    }

}
