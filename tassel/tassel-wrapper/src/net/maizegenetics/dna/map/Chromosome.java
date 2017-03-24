package net.maizegenetics.dna.map;

import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Defines the chromosome structure and length. The name and length recorded for
 * each chromosome.
 *
 * @author Terry Casstevens and Ed Buckler
 */
public class Chromosome implements Comparable<Chromosome> {

    private static final long serialVersionUID = -5197800047652332969L;
    public static Chromosome UNKNOWN = new Chromosome("Unknown");
    private final String myName;
    private final int myChromosomeNumber;
    private final int myLength;
    private final GeneralAnnotation myGA;
    private final int hashCode;

    //since there are numerous redundant chromosome, this class use a hash, so that
    //only the pointers are stored.
    private static final ConcurrentMap<Chromosome, Chromosome> CHR_HASH = new ConcurrentHashMap<>(50);

    public static Chromosome getCanonicalChromosome(Chromosome chr) {
        if (CHR_HASH.size() > 1000) {
            CHR_HASH.clear();
        }
        Chromosome canon = CHR_HASH.putIfAbsent(chr, chr);
        return (canon == null) ? chr : canon;
    }

    /**
     *
     * @param name Name of the chromosome
     * @param length Length of chromosome in base pairs
     * @param features Map of features about the chromosome
     */
    public Chromosome(String name, int length, GeneralAnnotation features) {
        myName = parseName(name);
        myLength = length;
        int convChr = Integer.MAX_VALUE;
        try {
            convChr = Integer.parseInt(myName);
        } catch (NumberFormatException ne) {
        }
        myChromosomeNumber = convChr;
        myGA = features;
        hashCode = calcHashCode();
    }

    public Chromosome(String name) {
        this(name, -1, parseAnnotationFromName(name));
    }

    public String getName() {
        return myName;
    }

    /**
     * Returns the integer value of the chromosome (if name is not a number then
     * Integer.MAX_VALUE is returned)
     */
    public int getChromosomeNumber() {
        return myChromosomeNumber;
    }

    public int getLength() {
        return myLength;
    }

    public GeneralAnnotation getAnnotation() {
        return myGA;
    }

    @Override
    public String toString() {
        return getName();
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    private int calcHashCode() {
        int hash = 7;
        hash = 79 * hash + (this.myName != null ? this.myName.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof Chromosome)) {
            return false;
        }
        if (hashCode!= ((Chromosome)obj).hashCode()) return false;
        return (compareTo((Chromosome) obj) == 0);
    }

    @Override
    public int compareTo(Chromosome o) {
        int result = myChromosomeNumber - o.getChromosomeNumber();
        if ((result != 0) || (myChromosomeNumber != Integer.MAX_VALUE)) {
            return result;
        }
        return myName.compareTo(o.getName());
        // compareNonNumericChrom() is too slow.  Use java lexicographic string
        // compare instead. Will re-visit this later
        // return compareNonNumericChrom(myName, o.getName()); 
    }
    /**
     * Takes a string, makes all upper case, removes leading CHROMOSOME/CHR, 
     * returns the resulting string
     * @param chromosome
     * @return the input string minus a leading "chr" or "chromsome" 
     */
    private static String parseName(String name) {
        if (name == null || name == "") return name;
        String parsedName = name.trim();
        parsedName = parsedName.toUpperCase();
        if (parsedName.startsWith("CHROMOSOME")) {
            parsedName = parsedName.replaceFirst("CHROMOSOME","");
        }
        if (parsedName.startsWith("CHR")) {
            parsedName = parsedName.replaceFirst("CHR","");
        }
        int spaceIndex = parsedName.indexOf(" ");
        if (spaceIndex > 0){
            parsedName = parsedName.substring(0,parsedName.indexOf(" "));
        }
        return parsedName;
    }
    
    /**
     * Takes a chromosome name, looks for the first space, returns
     * the data beyond as an annotation.  This takes care of lines in
     * a fasta file that look like this:
     * >3 This is a description
     * @param name - the string chromosome passed in
     * @return Annotations built from the string beyond the name
     */
    private static GeneralAnnotation parseAnnotationFromName(String name) {
        if (name == null || name == "") return null;
        GeneralAnnotation myAnnotations = null;
        String currChrDesc = null;
        int spaceIndex = name.indexOf(" ");
        if (spaceIndex > 0) {                   
            currChrDesc = name.substring(name.indexOf(" ") + 1);
            myAnnotations = GeneralAnnotationStorage.getBuilder().addAnnotation("Description", currChrDesc).build();
        } 

        return myAnnotations;
    }
    /**
     * Compares non-numeric chromsomes.  "chr"/"chromsome" was already parsed off name.
     * Strictly numeric chromosomes were already handles.
     * This code compares chromsomes where one or both contain non-numeric 
     * characters.  A strict name comparison does not work as "10" is 
     * considered less than "9" by java string compare, but not by scientists.
     * Same problem with 9A and 10D - 10D is considered smaller by java.
     * I'm doing this for the wheat folks.
     * @param chromA
     * @param chromB
     * @return
     */
    private static int compareNonNumericChrom(String chromA, String chromB) {
        int chromAInt = Integer.MAX_VALUE;
        int chromBInt = Integer.MAX_VALUE;
        String chromAIntStr = null;
        String chromBIntStr = null;
        
        Matcher matcher1 = Pattern.compile("\\D+").matcher(chromA);
        Matcher matcher2 = Pattern.compile("\\D+").matcher(chromB);
        int startAlphaChromA = -1;
        if (matcher1.find()) { // find first non-numeric in chromA
            startAlphaChromA = matcher1.start();
        }
        int startAlphaChromB = -1;
        if (matcher2.find()){ // first non-numeric in chromB
            startAlphaChromB = matcher2.start();
        }
 
        if (startAlphaChromA != 0) { // -1 means NO non-numeric chars, > 0 is start of non-numeric
            int endString = startAlphaChromA > 0 ? startAlphaChromA:chromA.length();
            chromAIntStr = chromA.substring(0,endString);           
        }
        if (startAlphaChromB != 0) {
            int endString = startAlphaChromB > 0 ? startAlphaChromB:chromB.length();
            chromBIntStr = chromB.substring(0,endString);           
        } 
 
        //Compare the integer portions
        if (chromBIntStr != null && chromAIntStr != null) {
            int convChr = Integer.MAX_VALUE;
            try {
                convChr = Integer.parseInt(chromAIntStr);
            } catch (NumberFormatException ne) {
            }
            chromAInt = convChr;
            try {
                convChr = Integer.parseInt(chromBIntStr);
            } catch (NumberFormatException ne) {
            }
            chromBInt = convChr;
            int result = chromAInt - chromBInt;
            if ((result != 0) || ( result == 0 && chromAInt == Integer.MAX_VALUE )) {
                return result;
            }
        }
        // Either one or both chrom strings do not start with a number, or they do
        // but it is the same number.  Comparison must be based on the alpha part
        String chromAAlpha = null;
        String chromBAlpha = null;
 
        if (startAlphaChromA < 0 || startAlphaChromB < 0 ) return chromA.compareTo(chromB);
        chromAAlpha = chromA.substring(startAlphaChromA);
        chromBAlpha = chromB.substring(startAlphaChromB);
        return chromAAlpha.compareTo(chromBAlpha);
    }
}
