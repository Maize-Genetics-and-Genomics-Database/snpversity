/*
 * XYMultipleYToolTipGenerator
 */
package net.maizegenetics.analysis.chart;

import java.text.DecimalFormat;
import org.jfree.chart.labels.AbstractXYItemLabelGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.data.xy.XYDataset;

/**
 *
 * @author yz79
 */
public class XYMultipleYToolTipGenerator extends AbstractXYItemLabelGenerator implements XYToolTipGenerator {

    @Override
    public String generateToolTip(XYDataset dataset, int series, int item) {
        TableReportManhattanDataset myDataset = (TableReportManhattanDataset) dataset;
        String[] markers = myDataset.getMarkers();
        String[] chroms = myDataset.getChroms();
        DecimalFormat df = new DecimalFormat("#0.000");
        StringBuilder sb = new StringBuilder("SNP ID: ");
        sb.append(markers[item]);
        sb.append(", Chromosome: ");
        sb.append(chroms[item]);
        sb.append(", -Log10(P-Value): ");
        sb.append(df.format(myDataset.getYValue(series, item)));
        sb.append(", Position: ");
        sb.append(myDataset.getXValue(series, item));
        return sb.toString();
    }

}
