/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license, version 2 and without
 * any warranty or technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
//Title:      TASSELMainFrame
//Version:
//Copyright:  Copyright (c) 1997
//Author:     Ed Buckler
//Company:    NCSU
package net.maizegenetics.tassel;

import net.maizegenetics.analysis.imputation.*;
import net.maizegenetics.analysis.popgen.SequenceDiversityPlugin;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.distance.MultiDimensionalScalingPlugin;
import net.maizegenetics.analysis.chart.TableDisplayPlugin;
import net.maizegenetics.analysis.chart.ManhattanDisplayPlugin;
import net.maizegenetics.analysis.chart.QQDisplayPlugin;
import net.maizegenetics.analysis.association.EqtlAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.GenomicSelectionPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.association.WeightedMLMPlugin;
import net.maizegenetics.analysis.popgen.LinkageDiseqDisplayPlugin;
import net.maizegenetics.analysis.popgen.LinkageDisequilibriumPlugin;
import net.maizegenetics.analysis.data.MergeGenotypeTablesPlugin;
import net.maizegenetics.analysis.data.PrincipalComponentsPlugin;
import net.maizegenetics.analysis.data.UnionAlignmentPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.IntersectionAlignmentPlugin;
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.analysis.data.SeparatePlugin;
import net.maizegenetics.analysis.data.SynonymizerPlugin;
import net.maizegenetics.analysis.data.ThinSitesByPositionPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaPropertiesPlugin;
import net.maizegenetics.analysis.filter.FilterSiteNamePlugin;
import net.maizegenetics.analysis.filter.FilterAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterTraitsPlugin;
import net.maizegenetics.analysis.tree.CreateTreePlugin;
import net.maizegenetics.analysis.tree.ArchaeopteryxPlugin;
import net.maizegenetics.analysis.chart.ChartDisplayPlugin;
import net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitterPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.analysis.numericaltransform.TransformDataPlugin;
import net.maizegenetics.gui.PrintHeapAction;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.progress.ProgressPanel;

import org.apache.log4j.Logger;

import javax.swing.*;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.awt.font.TextAttribute;
import java.io.*;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.analysis.data.GenosToABHPlugin;
import net.maizegenetics.analysis.data.GetPositionListPlugin;
import net.maizegenetics.analysis.data.GetTaxaListPlugin;
import net.maizegenetics.analysis.data.HetsToUnknownPlugin;
import net.maizegenetics.analysis.data.SetLowDepthGenosToMissingPlugin;
import net.maizegenetics.analysis.data.SortGenotypeFilePlugin;
import net.maizegenetics.analysis.distance.AMatrixPlugin;
import net.maizegenetics.analysis.distance.DistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.HMatrixPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.gbs.BinaryToTextPlugin;
import net.maizegenetics.analysis.gbs.DiscoverySNPCallerPlugin;
import net.maizegenetics.analysis.gbs.FastqToTagCountPlugin;
import net.maizegenetics.analysis.gbs.MergeMultipleTagCountPlugin;
import net.maizegenetics.analysis.gbs.ModifyTBTHDF5Plugin;
import net.maizegenetics.analysis.gbs.ProductionSNPCallerPlugin;
import net.maizegenetics.analysis.gbs.SAMConverterPlugin;
import net.maizegenetics.analysis.gbs.SeqToTBTHDF5Plugin;
import net.maizegenetics.analysis.gbs.TagCountToFastqPlugin;
import net.maizegenetics.analysis.gbs.UTagCountToTagPairPlugin;
import net.maizegenetics.analysis.gbs.UTagPairToTOPMPlugin;
import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;
import net.maizegenetics.analysis.workflow.WorkflowPlugin;

/**
 * TASSELMainFrame
 *
 */
public class TASSELMainFrame extends JFrame implements ActionListener {

    private static final Logger myLogger = Logger.getLogger(TASSELMainFrame.class);
    public static final String version = "5.2.17";
    public static final String versionDate = "November 12, 2015";
    private DataTreePanel myDataTreePanel;
    //a variable to control when the progress bar was last updated
    private JFileChooser filerSave = new JFileChooser();
    private JFileChooser filerOpen = new JFileChooser();
    private JScrollPane reportPanelScrollPane = new JScrollPane();
    private JTextArea reportPanelTextArea = new JTextArea();
    JScrollPane mainPanelScrollPane = new JScrollPane();
    JPanel mainDisplayPanel = new JPanel();
    private JTextArea mainPanelTextArea = new JTextArea();
    private JTextField myStatusTextField = new JTextField();
    private final ProgressPanel myProgressPanel = ProgressPanel.getInstance();
    private HashMap<JMenuItem, Plugin> myMenuItemHash = new HashMap<JMenuItem, Plugin>();

    public TASSELMainFrame() {
        try {
            loadSettings();
            myDataTreePanel = new DataTreePanel(this);
            myDataTreePanel.setToolTipText("Data Tree Panel");
            addMenuBar();
            initializeMyFrame();

            this.setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage) " + version);

            TasselLogging.basicLoggingInfo();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    //Component initialization
    private void initializeMyFrame() throws Exception {
        getContentPane().setLayout(new BorderLayout());

        int xDim = TasselPrefs.getXDim();
        int yDim = TasselPrefs.getYDim();
        if ((xDim < 50) || (yDim < 50)) {
            Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
            setSize(new Dimension(screenSize.width * 19 / 20, screenSize.height * 19 / 20));
        } else {
            setSize(xDim, yDim);
        }
        setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage)");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                TasselLogging.closeInstance();
                TasselPrefs.putXDim(getWidth());
                TasselPrefs.putYDim(getHeight());
                System.exit(0);
            }
        });

        filerSave.setDialogType(JFileChooser.SAVE_DIALOG);

        JSplitPane dataTreeReportPanelsSplitPanel = new JSplitPane();
        dataTreeReportPanelsSplitPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);

        reportPanelTextArea.setEditable(false);
        reportPanelTextArea.setToolTipText("Report Panel");
        reportPanelTextArea.setLineWrap(true);
        reportPanelTextArea.setWrapStyleWord(true);
        mainPanelTextArea.setDoubleBuffered(true);
        mainPanelTextArea.setEditable(false);
        mainPanelTextArea.setFont(new java.awt.Font("Monospaced", 0, 12));
        mainPanelTextArea.setToolTipText("Main Panel");

        myStatusTextField.setBackground(Color.lightGray);
        myStatusTextField.setBorder(null);

        JSplitPane dataTreeReportMainPanelsSplitPanel = new JSplitPane();

        getContentPane().add(dataTreeReportMainPanelsSplitPanel, BorderLayout.CENTER);

        dataTreeReportMainPanelsSplitPanel.add(dataTreeReportPanelsSplitPanel, JSplitPane.LEFT);

        dataTreeReportPanelsSplitPanel.add(myDataTreePanel, JSplitPane.TOP);

        JSplitPane reportProgressSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        JPanel reportPanel = new JPanel(new BorderLayout());

        reportProgressSplitPane.add(reportPanel, JSplitPane.TOP);

        reportPanel.add(reportPanelScrollPane, BorderLayout.CENTER);

        reportPanelScrollPane.getViewport().add(reportPanelTextArea, null);

        reportProgressSplitPane.add(new JScrollPane(myProgressPanel), JSplitPane.BOTTOM);

        dataTreeReportPanelsSplitPanel.add(reportProgressSplitPane, JSplitPane.BOTTOM);

        dataTreeReportMainPanelsSplitPanel.add(mainDisplayPanel, JSplitPane.RIGHT);

        mainDisplayPanel.setLayout(new BorderLayout());
        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);
        mainDisplayPanel.add(mainPanelScrollPane, BorderLayout.CENTER);

        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);

        getContentPane().add(myStatusTextField, BorderLayout.SOUTH);

        dataTreeReportMainPanelsSplitPanel.setDividerLocation(getSize().width / 4);

        dataTreeReportPanelsSplitPanel.setDividerLocation((int) (getSize().height / 3.5));

        reportProgressSplitPane.setDividerLocation((int) (getSize().height / 3.5));
    }

    private void addMenuBar() {

        JMenuBar jMenuBar = new JMenuBar();

        jMenuBar.add(getFileMenu());
        jMenuBar.add(getDataMenu());
        jMenuBar.add(getImputeMenu());
        jMenuBar.add(getFiltersMenu());
        jMenuBar.add(getAnalysisMenu());
        jMenuBar.add(getResultsMenu());
        jMenuBar.add(getGBSMenu());
        jMenuBar.add(getWorkflowMenu());
        jMenuBar.add(Box.createHorizontalGlue());
        jMenuBar.add(getHelpMenu());

        this.setJMenuBar(jMenuBar);

    }

    //Help | About action performed
    private void helpAbout_actionPerformed(ActionEvent e) {
        AboutBox dlg = new AboutBox(this);
        Dimension dlgSize = dlg.getPreferredSize();
        Dimension frmSize = getSize();
        Point loc = getLocation();
        dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);
        dlg.setModal(true);
        dlg.setVisible(true);
    }

    public void sendMessage(String text) {
        myStatusTextField.setForeground(Color.BLACK);
        myStatusTextField.setText(text);
    }

    public void sendErrorMessage(String text) {
        myStatusTextField.setForeground(Color.RED);
        myStatusTextField.setText(text);
    }

    public void setMainText(String text) {
        mainPanelTextArea.setText(text);
        mainPanelTextArea.setCaretPosition(0);
    }

    public void setNoteText(String text) {
        reportPanelTextArea.setText(text);
        reportPanelTextArea.setCaretPosition(0);
    }

    private void loadSettings() {
        filerOpen.setCurrentDirectory(new File(TasselPrefs.getOpenDir()));
        filerSave.setCurrentDirectory(new File(TasselPrefs.getSaveDir()));
    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {
        myDataTreePanel.addDataSet(theDataSet, defaultNode);
    }

    private void helpButton_actionPerformed(ActionEvent e) {
        HelpDialog theHelpDialog = new HelpDialog(this);
        theHelpDialog.setLocationRelativeTo(this);
        theHelpDialog.setVisible(true);
    }

    public void updateMainDisplayPanel(JPanel panel) {
        mainDisplayPanel.removeAll();
        mainDisplayPanel.add(panel, BorderLayout.CENTER);
        mainDisplayPanel.repaint();
        mainDisplayPanel.validate();
    }

    public DataTreePanel getDataTreePanel() {
        return myDataTreePanel;
    }

    public ProgressPanel getProgressPanel() {
        return myProgressPanel;
    }

    private JMenuItem createMenuItem(Plugin theTP) {
        return createMenuItem(theTP, -1);
    }

    private JMenuItem createMenuItem(Plugin theTP, boolean adjustForIcon) {
        return createMenuItem(theTP, -1, adjustForIcon);
    }

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic) {
        return createMenuItem(theTP, mnemonic, true);
    }

    private static final int ICON_WIDTH_PLUS_GAP = 30;

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic, boolean adjustForIcon) {
        ImageIcon icon = theTP.getIcon();
        JMenuItem menuItem = new JMenuItem(theTP.getButtonName(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        if (adjustForIcon) {
            int pixels = ICON_WIDTH_PLUS_GAP;
            if (icon != null) {
                pixels -= icon.getIconWidth();
                pixels /= 2;
            }
            menuItem.setIconTextGap(pixels);
        }
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setToolTipText(theTP.getToolTipText());
        menuItem.addActionListener(this);
        theTP.addListener(myDataTreePanel);
        myMenuItemHash.put(menuItem, theTP);
        return menuItem;
    }

    private JMenuItem createMenuItem(Action action, int mnemonic) {

        Icon icon = (Icon) action.getValue("SmallIcon");

        JMenuItem menuItem = new JMenuItem(action.getValue("Name").toString(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        int pixels = ICON_WIDTH_PLUS_GAP;
        if (icon != null) {
            pixels -= icon.getIconWidth();
            pixels /= 2;
        }
        menuItem.setIconTextGap(pixels);
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.addActionListener(action);
        return menuItem;

    }

    private JMenuItem createMenuItem(String label, boolean adjustForIcon) {

        JMenuItem menuItem = new JMenuItem(label);
        if (adjustForIcon) {
            menuItem.setIconTextGap(ICON_WIDTH_PLUS_GAP);
        }
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setEnabled(false);
        return menuItem;

    }

    private JMenu getFiltersMenu() {
        JMenu result = new JMenu("Filter");
        result.setMnemonic(KeyEvent.VK_F);
        result.add(createMenuItem(new FilterAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterSiteNamePlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaPropertiesPlugin(this, true)));
        result.add(createMenuItem(new FilterTraitsPlugin(this, true)));
        result.add(createMenuItem(new FilterSiteBuilderPlugin(this, true)));
        return result;
    }

    private JMenu getDataMenu() {

        JMenu result = new JMenu("Data");
        result.setMnemonic(KeyEvent.VK_D);

        result.add(createMenuItem(new FileLoadPlugin(this, true), KeyEvent.VK_L));
        result.add(createMenuItem(new ExportPlugin(this, true)));
        result.add(createMenuItem(new GetTaxaListPlugin(this, true)));
        result.add(createMenuItem(new GetPositionListPlugin(this, true)));
        result.add(createMenuItem(new SortGenotypeFilePlugin(this, true)));
        result.add(createMenuItem(new SynonymizerPlugin(this, true)));
        result.add(createMenuItem(new IntersectionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new UnionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new MergeGenotypeTablesPlugin(this, true)));
        result.add(createMenuItem(new SeparatePlugin(this, true)));
        result.add(createMenuItem(new HetsToUnknownPlugin(this, true)));
        result.add(createMenuItem(new SetLowDepthGenosToMissingPlugin(this, true)));
        result.add(createMenuItem(new TransformDataPlugin(this, true)));
        result.add(createMenuItem(new NumericalGenotypePlugin(this, true)));
        result.add(createMenuItem(new GenosToABHPlugin(this, true)));
        result.add(createMenuItem(new ThinSitesByPositionPlugin(this, true)));
        result.add(createMenuItem(new ClusterGenotypesPlugin(this, true)));
        result.addSeparator();

        JMenuItem delete = new JMenuItem("Delete Dataset");
        delete.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                myDataTreePanel.deleteSelectedNodes();
            }
        });
        delete.setToolTipText("Delete Dataset");
        URL imageURL = TASSELMainFrame.class.getResource("images/trash.gif");
        if (imageURL == null) {
            delete.setIconTextGap(ICON_WIDTH_PLUS_GAP);
        } else {
            int pixels = ICON_WIDTH_PLUS_GAP;
            ImageIcon icon = new ImageIcon(imageURL);
            if (icon != null) {
                pixels -= icon.getIconWidth();
                pixels /= 2;
            }
            delete.setIcon(icon);
            delete.setIconTextGap(pixels);
        }
        delete.setIconTextGap(6);
        result.add(delete);

        return result;
    }
    
    private JMenu getImputeMenu() {

        JMenu result = new JMenu("Impute");
        result.setMnemonic(KeyEvent.VK_I);

        result.add(createMenuItem(new FILLINFindHaplotypesPlugin(this, true), false));
        result.add(createMenuItem(new FILLINImputationPlugin(this, true), false));
        result.add(createMenuItem(new FSFHapImputationPlugin(this, true), false));
        result.add(createMenuItem(new ImputationPlugin(this, true), false));
        result.add(createMenuItem(new RemoveIndelsForBeaglePlugin(this, true), false));
        result.add(createMenuItem(new LDKNNiImputationPlugin(this, true), false));
        result.add(createMenuItem(new ImputationAccuracyPlugin(this, true), false));
        return result;
    }

    private JMenu getAnalysisMenu() {

        JMenu result = new JMenu("Analysis");
        result.setMnemonic(KeyEvent.VK_A);

        result.add(createMenuItem(new SequenceDiversityPlugin(this, true)));
        result.add(createMenuItem(new LinkageDisequilibriumPlugin(this, true)));
        result.add(createMenuItem(new DistanceMatrixPlugin(this, true)));
        result.add(createMenuItem(new CreateTreePlugin(this, true)));
        result.add(createMenuItem(new KinshipPlugin(this, true)));
        result.add(createMenuItem(new AMatrixPlugin(this, true)));
        result.add(createMenuItem(new HMatrixPlugin(this, true)));
        result.add(createMenuItem(new PrincipalComponentsPlugin(this, true)));
        result.add(createMenuItem(new MultiDimensionalScalingPlugin(this, true)));
        result.add(createMenuItem(new FixedEffectLMPlugin(this, true)));
        result.add(createMenuItem(new MLMPlugin(this, true)));
        result.add(createMenuItem(new WeightedMLMPlugin(this, true)));
        result.add(createMenuItem(new GenomicSelectionPlugin(this, true)));
        result.add(createMenuItem(new GenotypeSummaryPlugin(this, true)));
        result.add(createMenuItem(new StepwiseOLSModelFitterPlugin(this, true)));
        result.add(createMenuItem(new EqtlAssociationPlugin(this, true)));
        return result;
    }

    private JMenu getResultsMenu() {

        JMenu result = new JMenu("Results");
        result.setMnemonic(KeyEvent.VK_R);

        result.add(createMenuItem(new TableDisplayPlugin(this, true)));
        result.add(createMenuItem(new ArchaeopteryxPlugin(this, true)));
        result.add(createMenuItem(new LinkageDiseqDisplayPlugin(this, true)));
        result.add(createMenuItem(new ChartDisplayPlugin(this, true)));
        result.add(createMenuItem(new QQDisplayPlugin(this, true)));
        result.add(createMenuItem(new ManhattanDisplayPlugin(this, true)));
        return result;

    }

    private JMenu getFileMenu() {
        JMenu fileMenu = new JMenu();
        fileMenu.setText("File");

        fileMenu.add(createMenuItem(new PreferencesDialog(this, true), false));

        fileMenu.addSeparator();

        JMenuItem exitMenuItem = new JMenuItem("Exit");
        exitMenuItem.addActionListener((ActionEvent e) -> {
            TasselLogging.closeInstance();
            TasselPrefs.putXDim(getWidth());
            TasselPrefs.putYDim(getHeight());
            System.exit(0);
        });
        fileMenu.add(exitMenuItem);

        return fileMenu;
    }

    private JMenu getGBSMenu() {

        JMenu result = new JMenu("GBS");
        result.setMnemonic(KeyEvent.VK_G);

        result.add(createMenuItem(new BinaryToTextPlugin(this, true), false));
        result.add(createMenuItem(new FastqToTagCountPlugin(this, true), false));
        result.add(createMenuItem(new MergeMultipleTagCountPlugin(this, true), false));
        result.addSeparator();
        result.add(getGBSReferenceMenu());
        result.addSeparator();
        result.add(getUNEAKMenu());
        result.addSeparator();
        result.add(createMenuItem(new SeqToTBTHDF5Plugin(this, true), false));
        result.add(createMenuItem(new ModifyTBTHDF5Plugin(this, true), false));
        result.add(createMenuItem(new DiscoverySNPCallerPlugin(this, true), false));
        result.add(createMenuItem(new ProductionSNPCallerPlugin(this, true), false));

        return result;
    }

    private JMenu getGBSReferenceMenu() {

        JMenu result = new JMenu("Reference Genome");

        result.add(createMenuItem(new TagCountToFastqPlugin(this, true), false));
        result.add(createMenuItem("Align to Reference", false));
        result.add(createMenuItem(new SAMConverterPlugin(this, true), false));

        return result;
    }

    private JMenu getUNEAKMenu() {

        JMenu result = new JMenu("UNEAK (No Reference)");

        JMenuItem item = createMenuItem(new UTagCountToTagPairPlugin(this, true));
        Map attributes = item.getFont().getAttributes();
        attributes.put(TextAttribute.STRIKETHROUGH, TextAttribute.STRIKETHROUGH_ON);
        item.setFont(new Font(attributes));
        result.add(item);

        item = createMenuItem(new UTagPairToTOPMPlugin(this, true));
        attributes = item.getFont().getAttributes();
        attributes.put(TextAttribute.STRIKETHROUGH, TextAttribute.STRIKETHROUGH_ON);
        item.setFont(new Font(attributes));
        result.add(item);

        return result;
    }

    private JMenu getWorkflowMenu() {

        JMenu result = new JMenu("Workflow");

        List<WorkflowPlugin> workflows = WorkflowPlugin.getInstances(this);

        for (Plugin current : workflows) {
            result.add(createMenuItem(current));
        }

        return result;
    }

    private JMenu getHelpMenu() {
        JMenu helpMenu = new JMenu();
        helpMenu.setMnemonic(KeyEvent.VK_H);
        helpMenu.setText("Help");

        URL infoImageURL = TASSELMainFrame.class.getResource("images/info.gif");
        Icon infoIcon = null;
        if (infoImageURL != null) {
            infoIcon = new ImageIcon(infoImageURL);
        }
        helpMenu.add(createMenuItem(new AbstractAction("Help Manual", infoIcon) {

            @Override
            public void actionPerformed(ActionEvent e) {
                helpButton_actionPerformed(e);
            }
        }, -1));

        URL aboutImageURL = TASSELMainFrame.class.getResource("images/Tassel_Logo16.png");
        Icon aboutIcon = null;
        if (aboutImageURL != null) {
            aboutIcon = new ImageIcon(aboutImageURL);
        }
        helpMenu.add(createMenuItem(new AbstractAction("About", aboutIcon) {

            @Override
            public void actionPerformed(ActionEvent e) {
                helpAbout_actionPerformed(e);
            }
        }, -1));

        helpMenu.add(createMenuItem(PrintHeapAction.getInstance(this), -1));
        helpMenu.add(createMenuItem(TasselLogging.getInstance(this)));

        return helpMenu;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JMenuItem theMenuItem = (JMenuItem) e.getSource();
        Plugin theTP = this.myMenuItemHash.get(theMenuItem);
        PluginEvent event = new PluginEvent(myDataTreePanel.getSelectedTasselDataSet());
        ProgressPanel progressPanel = getProgressPanel();
        progressPanel.addPlugin(theTP);
        ThreadedPluginListener thread = new ThreadedPluginListener(theTP, event);
        thread.start();
    }
}
