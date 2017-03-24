/*
 *  EndelmanDistanceMatrix
 * 
 *  Created on June 30, 2015
 */
package net.maizegenetics.analysis.distance;

import java.util.Arrays;
import java.util.Optional;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class EndelmanDistanceMatrix {

    private static final Logger myLogger = Logger.getLogger(EndelmanDistanceMatrix.class);

    private static final int DEFAULT_MAX_ALLELES = 2;

    /**
     * Compute Endelman kinship for all pairs of taxa. Missing sites are
     * ignored. http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13
     *
     * @param genotype Genotype Table used to compute kinship
     *
     * @return Endelman Kinship Matrix
     */
    private EndelmanDistanceMatrix() {
        // utility
    }

    /**
     * Compute Endelman Kinship Matrix. Maximum alleles per site to evaluate
     * defaults to 2.
     *
     * @param genotype Genotype Table used to compute kinship
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype) {
        return getInstance(genotype, DEFAULT_MAX_ALLELES, null);
    }

    /**
     * Compute Endelman Kinship Matrix
     *
     * @param genotype Genotype Table used to compute kinship
     * @param maxAlleles maximum alleles per site to evaluate. i.e. Set to 3 to
     * evaluate the three most frequent allele states.
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype, int maxAlleles) {
        return getInstance(genotype, maxAlleles, null);
    }

    /**
     * Compute Endelman Kinship Matrix. Maximum alleles per site to evaluate
     * defaults to 2.
     *
     * @param genotype Genotype Table used to compute kinship
     * @param listener Progress listener
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype, ProgressListener listener) {
        return computeEndelmanDistances(genotype, DEFAULT_MAX_ALLELES, listener);
    }

    /**
     * Compute Endelman Kinship Matrix
     *
     * @param genotype Genotype Table used to compute kinship
     * @param maxAlleles maximum alleles per site to evaluate. i.e. Set to 3 to
     * evaluate the three most frequent allele states.
     * @param listener Progress listener
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype, int maxAlleles, ProgressListener listener) {
        return computeEndelmanDistances(genotype, maxAlleles, listener);
    }

    private static DistanceMatrix computeEndelmanDistances(GenotypeTable genotype, int maxAlleles, ProgressListener listener) {

        if ((maxAlleles < 2) || (maxAlleles > 6)) {
            throw new IllegalArgumentException("EndelmanDistanceMatrix: computeEndelmanDistances: max alleles must be between 2 and 6 inclusive.");
        }

        int numSeqs = genotype.numberOfTaxa();
        long time = System.currentTimeMillis();

        //
        // Sets up parellel stream to divide up sites for processing.
        // Also reduces the distance sums and sum of frequencies into one instance.
        //
        Optional<CountersDistances> optional = stream(genotype, maxAlleles, listener).reduce((CountersDistances t, CountersDistances u) -> {
            t.addAll(u);
            return t;
        });

        if (!optional.isPresent()) {
            return null;
        }
        CountersDistances counters = optional.get();
        double sumpk = counters.mySumPi;
        float[] distances = counters.myDistances;

        //
        // This does the final division of the frequency sum into
        // the distance sums.
        //
        sumpk *= 2.0;
        double[][] result = new double[numSeqs][numSeqs];
        int index = 0;
        for (int t = 0; t < numSeqs; t++) {
            for (int i = 0, n = numSeqs - t; i < n; i++) {
                result[t][t + i] = result[t + i][t] = distances[index] / sumpk;
                index++;
            }
        }

        myLogger.info("EndelmanDistanceMatrix: computeEndelmanDistances time: " + (System.currentTimeMillis() - time) / 1000 + " seconds");
        return new DistanceMatrix(result, genotype.taxa());

    }

    protected static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            if (percent > 100) {
                percent = 100;
            }
            listener.progress(percent, null);
        }

    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate terms of the Endelman equation. These are
    // combined with addAll() to result in one instance at the end.
    //
    private static class CountersDistances {

        private double mySumPi = 0.0;
        private final float[] myDistances;
        private final int myNumTaxa;

        public CountersDistances(int numTaxa) {
            myNumTaxa = numTaxa;
            myDistances = new float[myNumTaxa * (myNumTaxa + 1) / 2];
        }

        public void addAll(CountersDistances counters) {
            float[] otherDistances = counters.myDistances;
            for (int t = 0, n = myDistances.length; t < n; t++) {
                myDistances[t] += otherDistances[t];
            }
            mySumPi += counters.mySumPi;
        }

    }

    //
    // This pre-calculates the number of occurances of the allele
    // for all possible diploid allele values.  Numbers 0 through 7
    // represent A, C, G, T, -, +, N respectively.  First three bits
    // codes the allele.  Remaining six bits codes the diploid
    // allele values. The stored counts are encodings.  Value 7 (bits 111) means
    // it's not a comparable combination because either major allele
    // is unknown or the diploid allele value is unknown.
    // Code 1 (bits 001) is zero count.
    // Code 2 (bits 010) is one count.
    // Code 4 (bits 100) is two count.
    //
    private static final byte[] PRECALCULATED_COUNTS = new byte[512];

    static {
        for (int allele = 0; allele < 8; allele++) {
            for (int a = 0; a < 8; a++) {
                for (int b = 0; b < 8; b++) {
                    int temp = (allele << 6) | (a << 3) | b;
                    if ((allele == 7) | ((a == 7) && (b == 7))) {
                        PRECALCULATED_COUNTS[temp] = 7;
                    } else {
                        if (a == allele) {
                            if (b == allele) {
                                PRECALCULATED_COUNTS[temp] = 4;
                            } else {
                                PRECALCULATED_COUNTS[temp] = 2;
                            }
                        } else if (b == allele) {
                            PRECALCULATED_COUNTS[temp] = 2;
                        } else {
                            PRECALCULATED_COUNTS[temp] = 1;
                        }
                    }
                }
            }
        }
    }

    private static final int NUM_CORES_TO_USE = TasselPrefs.getMaxThreads();

    //
    // Used to report progress.  This is not thread-safe but
    // works well enough for this purpose.
    //
    private static int myNumSitesProcessed = 0;

    //
    // Creates stream from EndelmanSiteSpliterator and Genotype Table
    //
    private static Stream<CountersDistances> stream(GenotypeTable genotypes, int maxAlleles, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new EndelmanSiteSpliterator(genotypes, 0, genotypes.numberOfSites(), maxAlleles, listener), true);
    }

    //
    // Spliterator that splits the sites into halves each time for
    // processing.
    //
    static class EndelmanSiteSpliterator implements Spliterator<CountersDistances> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final int myMaxAlleles;
        private final ProgressListener myProgressListener;
        private int myMinSitesToProcess;

        EndelmanSiteSpliterator(GenotypeTable genotypes, int currentIndex, int fence, int maxAlleles, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myCurrentSite = currentIndex;
            myFence = fence;
            myMaxAlleles = maxAlleles;
            myProgressListener = listener;
            myMinSitesToProcess = Math.max(myNumSites / NUM_CORES_TO_USE, 1000);
        }

        @Override
        public void forEachRemaining(Consumer<? super CountersDistances> action) {

            int numSitesProcessed = myFence - myCurrentSite;
            CountersDistances result = new CountersDistances(myNumTaxa);
            float[] distances = result.myDistances;
            double[] sumpi = new double[1];

            float[] answer1 = new float[32768];
            float[] answer2 = new float[32768];
            float[] answer3 = new float[32768];

            for (; myCurrentSite < myFence;) {

                //
                // This keeps track of number of sites processed.  The blocks
                // of sites may contain entries for minor allele, 2nd minor
                // allele, etc.
                int[] realSites = new int[1];

                //
                // Pre-calculates possible terms and gets counts for
                // three blocks for five (pseudo-)sites.
                //
                Tuple<short[], float[]> firstBlock = getBlockOfSites(myCurrentSite, sumpi, realSites);
                float[] possibleTerms = firstBlock.y;
                short[] alleleCount1 = firstBlock.x;

                Tuple<short[], float[]> secondBlock = getBlockOfSites(myCurrentSite + realSites[0], sumpi, realSites);
                float[] possibleTerms2 = secondBlock.y;
                short[] alleleCount2 = secondBlock.x;

                Tuple<short[], float[]> thirdBlock = getBlockOfSites(myCurrentSite + realSites[0], sumpi, realSites);
                float[] possibleTerms3 = thirdBlock.y;
                short[] alleleCount3 = thirdBlock.x;

                myCurrentSite += realSites[0];

                //
                // Using possible terms, calculates all possible answers
                // for each site block.
                //
                for (int i = 0; i < 32768; i++) {
                    answer1[i] = possibleTerms[(i & 0x7000) >>> 12] + possibleTerms[((i & 0xE00) >>> 9) | 0x8] + possibleTerms[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms[((i & 0x38) >>> 3) | 0x18] + possibleTerms[(i & 0x7) | 0x20];
                    answer2[i] = possibleTerms2[(i & 0x7000) >>> 12] + possibleTerms2[((i & 0xE00) >>> 9) | 0x8] + possibleTerms2[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms2[((i & 0x38) >>> 3) | 0x18] + possibleTerms2[(i & 0x7) | 0x20];
                    answer3[i] = possibleTerms3[(i & 0x7000) >>> 12] + possibleTerms3[((i & 0xE00) >>> 9) | 0x8] + possibleTerms3[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms3[((i & 0x38) >>> 3) | 0x18] + possibleTerms3[(i & 0x7) | 0x20];
                }

                //
                // Iterates through all pair-wise combinations of taxa adding
                // distance comparisons and site counts.
                //
                int index = 0;
                for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                    //
                    // Can skip inter-loop if all fifteen sites for first
                    // taxon is Unknown diploid allele values
                    //
                    if ((alleleCount1[firstTaxa] != 0x7FFF) || (alleleCount2[firstTaxa] != 0x7FFF) || (alleleCount3[firstTaxa] != 0x7FFF)) {
                        for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                            //
                            // Combine first taxon's allele counts with
                            // second taxon's major allele counts to
                            // create index into pre-calculated answers
                            //
                            distances[index] += answer1[alleleCount1[firstTaxa] | alleleCount1[secondTaxa]] + answer2[alleleCount2[firstTaxa] | alleleCount2[secondTaxa]] + answer3[alleleCount3[firstTaxa] | alleleCount3[secondTaxa]];
                            index++;
                        }
                    } else {
                        index += myNumTaxa - firstTaxa;
                    }
                }
            }

            result.mySumPi = sumpi[0];
            action.accept(result);
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);
        }

        private static final int NUM_SITES_PER_BLOCK = 5;

        private Tuple<short[], float[]> getBlockOfSites(int currentSite, double[] sumpi, int[] realSites) {

            int currentSiteNum = 0;

            //
            // This hold possible terms for the Endelman summation given
            // site's allele frequency.  First three bits
            // identifies relative site (0, 1, 2, 3, 4).  Remaining three bits
            // the allele counts encoding.
            //
            float[] possibleTerms = new float[40];

            //
            // This holds count of allele for each taxa.
            // Each short holds count (0, 1, 2, 3) for all four sites
            // at given taxon.  The count encodings are stored in three
            // bits each.
            //
            short[] alleleCount = new short[myNumTaxa];

            //
            // This initializes the counts to 0x7FFF.  That means
            // diploid allele values for the four sites are Unknown.
            //
            Arrays.fill(alleleCount, (short) 0x7FFF);

            while ((currentSiteNum < NUM_SITES_PER_BLOCK) && (currentSite < myFence)) {

                byte[] genotypes = myGenotypes.genotypeAllTaxa(currentSite);
                int[][] alleles = AlleleFreqCache.allelesSortedByFrequencyNucleotide(genotypes);
                int numAlleles = Math.min(alleles[0].length - 1, myMaxAlleles - 1);

                if (numAlleles + currentSiteNum <= NUM_SITES_PER_BLOCK) {

                    //
                    // Calculates total number of haploid alleles that
                    // are not missing.
                    //
                    int totalAlleleCount = 0;
                    for (int i = 0; i < alleles[1].length; i++) {
                        totalAlleleCount += alleles[1][i];
                    }

                    for (int a = 0; a < numAlleles; a++) {

                        byte allele = (byte) alleles[0][a];
                        float alleleFreq = (float) alleles[1][a] / (float) totalAlleleCount;
                        float alleleFreqTimes2 = alleleFreq * 2.0f;
                        sumpi[0] += alleleFreq * (1.0 - alleleFreq);

                        //
                        // Temporarily stores component terms of equation for
                        // individual allele counts (0, 1, 2)
                        //
                        float[] term = new float[3];

                        //
                        // If allele is Unknown, the entire
                        // site is skipped.
                        //
                        if (allele != GenotypeTable.UNKNOWN_ALLELE) {

                            term[0] = 0.0f - alleleFreqTimes2;
                            term[1] = 1.0f - alleleFreqTimes2;
                            term[2] = 2.0f - alleleFreqTimes2;

                            //
                            // Pre-calculates all possible terms of the summation
                            // for this current (pseudo-) site.
                            // Counts (0,0; 0,1; 0,2; 1,1; 1,2; 2,2)
                            //
                            int siteNumIncrement = currentSiteNum * 8;
                            possibleTerms[siteNumIncrement + 1] = term[0] * term[0];
                            possibleTerms[siteNumIncrement + 3] = term[0] * term[1];
                            possibleTerms[siteNumIncrement + 5] = term[0] * term[2];
                            possibleTerms[siteNumIncrement + 2] = term[1] * term[1];
                            possibleTerms[siteNumIncrement + 6] = term[1] * term[2];
                            possibleTerms[siteNumIncrement + 4] = term[2] * term[2];

                            //
                            // Records allele counts for current site in
                            // three bits.
                            //
                            int temp = (allele & 0x7) << 6;
                            int shift = (NUM_SITES_PER_BLOCK - currentSiteNum - 1) * 3;
                            int mask = ~(0x7 << shift) & 0x7FFF;
                            for (int i = 0; i < myNumTaxa; i++) {
                                alleleCount[i] = (short) (alleleCount[i] & (mask | PRECALCULATED_COUNTS[temp | ((genotypes[i] & 0x70) >>> 1) | (genotypes[i] & 0x7)] << shift));
                            }
                        }

                        currentSiteNum++;
                    }
                } else {
                    return new Tuple<>(alleleCount, possibleTerms);
                }

                currentSite++;
                realSites[0]++;
            }

            return new Tuple<>(alleleCount, possibleTerms);

        }

        @Override
        public boolean tryAdvance(Consumer<? super CountersDistances> action) {
            if (myCurrentSite < myFence) {
                forEachRemaining(action);
                return true;
            } else {
                return false;
            }
        }

        @Override
        /**
         * Splits sites
         */
        public Spliterator<CountersDistances> trySplit() {
            int lo = myCurrentSite;
            int mid = lo + myMinSitesToProcess;
            if (mid < myFence) {
                myCurrentSite = mid;
                return new EndelmanSiteSpliterator(myGenotypes, lo, mid, myMaxAlleles, myProgressListener);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return (long) (myFence - myCurrentSite);
        }

        @Override
        public int characteristics() {
            return IMMUTABLE;
        }
    }

}
