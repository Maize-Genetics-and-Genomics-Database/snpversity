/*
 * XYScatterMultipleYPanel
 */
package net.maizegenetics.analysis.chart;

import java.awt.Component;
import java.awt.event.ActionEvent;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import net.maizegenetics.util.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;

/**
 *
 * @author yz79
 */
public class XYScatterMultipleYPanel extends BasicChartPanel {

    ManhattanDisplayPlugin myManhattanDisplayPlugin;
    ChartPanel myChartPanel;
    JButton saveButton = new JButton("save...");
    TableReportManhattanDataset dataset;
    TableReport myTableReport;

    public XYScatterMultipleYPanel(ManhattanDisplayPlugin plugin, TableReport theTable, int start, int end) {
        myManhattanDisplayPlugin = plugin;
        myTableReport = theTable;
        dataset = new TableReportManhattanDataset(theTable, start, end);
        chart = createChart(dataset);
        myChartPanel = new ChartPanel(chart);
        myChartPanel.setPreferredSize(new java.awt.Dimension(900, 500));
        myTableReport = theTable;
        jbInit();
    }

    private void jbInit() {
        this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
        this.add(myChartPanel);
        saveButton.setAlignmentX(Component.RIGHT_ALIGNMENT);
        saveButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveButton_actionPerformed(e);
            }
        });
        this.add(saveButton);
    }

    private void saveButton_actionPerformed(ActionEvent e) {
        myManhattanDisplayPlugin.saveDataToFile(myChartPanel);
    }

    public JFreeChart createChart(TableReportManhattanDataset dataset) {
        String name = "Please select numeric variables";
        String xName = "X";
        String y1Name = "Y";
        String y2Name = "Y2";
        if (dataset != null) {
            xName = dataset.getXName();
            y1Name = "-Log10(P-Value)";
            name = "P-Values by Chromosome for " + dataset.getTrait();
        }
        chart = ChartFactory.createScatterPlot(
                name,
                xName, y1Name,
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        chart.getXYPlot().setForegroundAlpha(0.75f);
        chart.getXYPlot().getRenderer().setToolTipGenerator(new XYMultipleYToolTipGenerator());
        return chart;
    }

    public JComponent getMainComponent() {
        return myChartPanel;
    }
}
