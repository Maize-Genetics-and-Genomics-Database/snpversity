package net.maizegenetics.dna.map;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import com.google.common.base.Preconditions;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.*;
import org.apache.log4j.Logger;

import static net.maizegenetics.dna.WHICH_ALLELE.*;

/**
 * A builder for creating immutable PositionList. Can be used for either an in
 * memory or HDF5 list.
 *
 * <p>
 * Example:
 * <pre>   {@code
 *   PositionListBuilder b=new PositionArrayList.Builder();
 *   for (int i = 0; i <size; i++) {
 *       Position ap=new CoreAnnotatedPosition.Builder(chr[chrIndex[i]],pos[i]).refAllele(refSeq[i]).build();
 *       b.add(ap);
 *       }
 *   PositionList instance=b.build();}
 * <p></p>
 * If being built separately from the genotypes, then use validate ordering to make sure sites are added in the
 * intended order.  This list WILL be sorted.
 * <p>Builder instances can be reused - it is safe to call {@link #build()}
 * multiple times to build multiple lists in series. Each new list
 * contains the one created before it.
 *
 * HDF5 Example
 * <p>Example:
 * <pre>   {@code
 *   PositionList instance=new PositionHDF5List.Builder("fileName").build();
 *   }
 *
 * <p>Builder instances can be reused - it is safe to call {@link #build()}
 */
public class PositionListBuilder {

    private static final Logger myLogger = Logger.getLogger(PositionListBuilder.class);

    private ArrayList<Position> myPositions = new ArrayList<>();
    private boolean isHDF5 = false;
    private String genomeVersion = null;
    private IHDF5Reader reader;

    /**
     * Creates a new builder. The returned builder is equivalent to the builder
     * generated by {@link }.
     */
    public PositionListBuilder() {
    }

    /**
     * Creates a new builder with a given number of Positions. This is most
     * useful when the number of sites is known from the beginning and the set
     * method will be used to set positions perhaps out of order. Useful in
     * multithreaded builders.
     */
    public PositionListBuilder(int numberOfPositions) {
        for (int i = 0; i < numberOfPositions; i++) {
            myPositions.add(new GeneralPosition.Builder(Chromosome.UNKNOWN, i).build());
        }
    }

    /**
     * Adds {@code element} to the {@code PositionList}.
     *
     * @param element the element to add
     * @return this {@code Builder} object
     * @throws NullPointerException if {@code element} is null
     */
    public PositionListBuilder add(Position element) {
        if (isHDF5) {
            throw new UnsupportedOperationException("Positions cannot be added to existing HDF5 alignments");
        }
        Preconditions.checkNotNull(element, "element cannot be null");
        myPositions.add(element);
        return this;
    }

    /**
     * Adds each element of {@code elements} to the {@code PositionList}.
     *
     * @param collection collection containing positions to be added to this
     * list
     * @return this {@code Builder} object
     * @throws NullPointerException if {@code elements} is or contains null
     */
    public PositionListBuilder addAll(Collection<? extends Position> collection) {
        if (isHDF5) {
            throw new UnsupportedOperationException("Positions cannot be added to existing HDF5 alignments");
        }
        myPositions.ensureCapacity(myPositions.size() + collection.size());
        for (Position elem : collection) {
            Preconditions.checkNotNull(elem, "elements contains a null");
            myPositions.add(elem);
        }
        return this;
    }

    public PositionListBuilder addAll(PositionListBuilder builder) {
        if (isHDF5) {
            throw new UnsupportedOperationException("Positions cannot be added to existing HDF5 alignments");
        }
        myPositions.ensureCapacity(myPositions.size() + builder.size());
        for (Position elem : builder.myPositions) {
            Preconditions.checkNotNull(elem, "elements contains a null");
            myPositions.add(elem);
        }
        return this;
    }

    /**
     * Replaces the element at the specified position in this list with the
     * specified element.
     *
     * @param index index of the element to replace
     * @param element element to be stored at the specified position
     * @return this {@code Builder} object
     * @throws IndexOutOfBoundsException {@inheritDoc}
     */
    public PositionListBuilder set(int index, Position element) {
        if (isHDF5) {
            throw new UnsupportedOperationException("Positions cannot be edited to existing HDF5 alignments");
        }
        myPositions.set(index, element);
        return this;
    }

    public PositionListBuilder genomeVersion(String genomeVersion) {
        this.genomeVersion = genomeVersion;
        return this;
    }

    /**
     * Returns whether List is already ordered. Important to check this if
     * genotype and sites are separately built, as the PositionArrayList must be
     * sorted, and will be with build.
     */
    public boolean validateOrdering() {
        boolean result = true;
        Position startAP = myPositions.get(0);
        for (Position ap : myPositions) {
            if (ap.compareTo(startAP) < 0) {
                myLogger.error("validateOrdering: " + ap.toString() + " and " + startAP.toString() + " out of order.");
                return false;
            }
            startAP = ap;
        }
        return result;
    }

    /**
     * Returns the size (number of positions) in the current list
     *
     * @return current size
     */
    public int size() {
        return myPositions.size();
    }

    /**
     * Creates a new position list based on an existing HDF5 file.
     */
    public static PositionList getInstance(String hdf5Filename) {
        return new PositionHDF5List(HDF5Factory.openForReading(hdf5Filename));
    }

    /**
     * Creates a new builder based on an existing HDF5 file reader.
     */
    public static PositionList getInstance(IHDF5Reader reader) {
        return new PositionHDF5List(reader);
    }

    /**
     * Generates a generic position list when no position information known
     *
     * @param numSites number of sites
     *
     * @return generic position list
     */
    public static PositionList getInstance(int numSites) {
        PositionListBuilder builder = new PositionListBuilder();
        for (int i = 0; i < numSites; i++) {
            builder.add(new GeneralPosition.Builder(Chromosome.UNKNOWN, i).build());
        }
        return builder.build();
    }

    /**
     * Creates in memory of PositionList from the an array of positions.
     */
    public static PositionList getInstance(List<Position> positions) {
        PositionListBuilder builder = new PositionListBuilder();
        builder.addAll(positions);
        return builder.build();
    }

    /**
     * Creates a positionList in a new HDF5 file.
     */
    public PositionListBuilder(IHDF5Writer h5w, PositionList a) {
        HDF5Utils.createHDF5PositionModule(h5w);
        h5w.setIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES, a.size());
        if (a.hasReference()) {
            h5w.setStringAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_GENOME_VERSION, a.genomeVersion());
            h5w.setBooleanAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_HAS_REFEFERENCE, true);
        }
        String[] lociNames = new String[a.numChromosomes()];
        Map<Chromosome, Integer> locusToIndex = new HashMap<>(10);
        Chromosome[] loci = a.chromosomes();
        for (int i = 0; i < a.numChromosomes(); i++) {
            lociNames[i] = loci[i].getName();
            locusToIndex.put(loci[i], i);
        }
        h5w.createStringVariableLengthArray(Tassel5HDF5Constants.CHROMOSOMES, a.numChromosomes());
        h5w.writeStringVariableLengthArray(Tassel5HDF5Constants.CHROMOSOMES, lociNames);

        int blockSize = 1 << 16;
        h5w.createStringArray(Tassel5HDF5Constants.SNP_IDS, 15, a.numberOfSites(), blockSize, Tassel5HDF5Constants.genDeflation);
        h5w.createIntArray(Tassel5HDF5Constants.CHROMOSOME_INDICES, a.numberOfSites(), Tassel5HDF5Constants.intDeflation);
        h5w.createIntArray(Tassel5HDF5Constants.POSITIONS, a.numberOfSites(), Tassel5HDF5Constants.intDeflation);
        h5w.createIntArray(Tassel5HDF5Constants.REF_ALLELES, a.numberOfSites(), Tassel5HDF5Constants.intDeflation);
        h5w.createIntArray(Tassel5HDF5Constants.ANC_ALLELES, a.numberOfSites(), Tassel5HDF5Constants.intDeflation);

        //This is written in blocks to deal with datasets in the scale for 50M positions
        int blocks = ((a.numberOfSites() - 1) / blockSize) + 1;
        for (int block = 0; block < blocks; block++) {
            int startPos = block * blockSize;
            int length = ((a.numberOfSites() - startPos) > blockSize) ? blockSize : a.numberOfSites() - startPos;
            String[] snpIDs = new String[length];
            int[] locusIndicesArray = new int[length];
            int[] positions = new int[length];
            byte[] refAlleles = new byte[length];
            byte[] ancAlleles = new byte[length];
            for (int i = 0; i < length; i++) {
                Position gp = a.get(i + startPos);
                snpIDs[i] = gp.getSNPID();
                locusIndicesArray[i] = locusToIndex.get(gp.getChromosome());
                positions[i] = gp.getPosition();
                refAlleles[i] = gp.getAllele(Reference);
                ancAlleles[i] = gp.getAllele(Ancestral);
            }
            HDF5Utils.writeHDF5Block(Tassel5HDF5Constants.SNP_IDS, h5w, blockSize, block, snpIDs);
            HDF5Utils.writeHDF5Block(Tassel5HDF5Constants.CHROMOSOME_INDICES, h5w, blockSize, block, locusIndicesArray);
            HDF5Utils.writeHDF5Block(Tassel5HDF5Constants.POSITIONS, h5w, blockSize, block, positions);
            HDF5Utils.writeHDF5Block(Tassel5HDF5Constants.REF_ALLELES, h5w, blockSize, block, refAlleles);
            HDF5Utils.writeHDF5Block(Tassel5HDF5Constants.ANC_ALLELES, h5w, blockSize, block, ancAlleles);
        }
        this.reader = h5w;
        isHDF5 = true;
    }

    /**
     * Returns a newly-created {@code ImmutableList} based on the myPositions of
     * the {@code Builder}.
     */
    public PositionList build() {
        if (isHDF5) {
            return new PositionHDF5List(reader);
        } else {
            Collections.sort(myPositions);
            return new PositionArrayList(myPositions, genomeVersion);
        }
    }

    public PositionList build(GenotypeCallTableBuilder genotypes) {
        sortPositions(genotypes);
        return new PositionArrayList(myPositions, genomeVersion);
    }

    public PositionListBuilder sortPositions(GenotypeCallTableBuilder genotypes) {
        int numPositions = myPositions.size();
        if (numPositions != genotypes.getSiteCount()) {
            throw new IllegalArgumentException("PositionListBuilder: sortPositions: position list size: " + numPositions + " doesn't match genotypes num position: " + genotypes.getSiteCount());
        }
        genotypes.reorderPositions(sort());
        return this;
    }

    public PositionListBuilder sortPositions() {
        sort();
        return this;
    }

    private int[] sort() {

        int numPositions = myPositions.size();

        final int indicesOfSortByPosition[] = new int[numPositions];
        for (int i = 0; i < indicesOfSortByPosition.length; i++) {
            indicesOfSortByPosition[i] = i;
        }

        Swapper swapPosition = new Swapper() {
            @Override
            public void swap(int a, int b) {
                int temp = indicesOfSortByPosition[a];
                indicesOfSortByPosition[a] = indicesOfSortByPosition[b];
                indicesOfSortByPosition[b] = temp;
            }
        };

        IntComparator compPosition = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                return myPositions.get(indicesOfSortByPosition[a]).compareTo(myPositions.get(indicesOfSortByPosition[b]));
            }
        };

        GenericSorting.quickSort(0, indicesOfSortByPosition.length, compPosition, swapPosition);

        ArrayList<Position> temp = new ArrayList<>(numPositions);
        for (int t = 0; t < numPositions; t++) {
            temp.add(myPositions.get(indicesOfSortByPosition[t]));
        }

        myPositions = temp;

        return indicesOfSortByPosition;

    }
}