package net.maizegenetics.analysis.association;

import java.util.List;

import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.TableReport;

public interface FixedEffectLM {

	/**
	 * Initializes builders for the site and allele reports. The report names and column headers for a specific analysis are created in this method.
	 */
	public void initializeReportBuilders();

	/**
	 * Solves the linear model
	 */
	public void solve();

	/**
	 * @return		the site report, which contains tests of significance for each site
	 */
	public TableReport siteReport();

	/**
	 * @return		the allele report, which contains estimates and observation numbers for each allele at each site
	 */
	public TableReport alleleReport();

	/**
	 * @return		the site and allele reports as a list of Datum
	 */
	public List<Datum> datumList();
	
	/**
	 * @param permute	if true, conduct an experiment-wise permutation test
	 * @param nperm		the number of permutations to be run
	 */
	public void permutationTest(boolean permute, int nperm);
	
	/**
	 * @param maxP		test results with p > maxP will not be reported.
	 */
	public void maxP(double maxP);
	
	/**
	 * @param savefile	results will be saved to this file instead of memory
	 */
	public void siteReportFilepath(String savefile);
	
	/**
	 * @param savefile	results will be saved to this file instead of memory
	 */
	public void alleleReportFilepath(String savefile);

}