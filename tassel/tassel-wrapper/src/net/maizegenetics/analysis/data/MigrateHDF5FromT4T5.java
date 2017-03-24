package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.HDF5ObjectType;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.dna.snp.HapMapHDF5Constants;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.List;

/**
 * Provides a migration tool from TASSEL4 HDF5 to TASSEL5 HDF5
 *
 * @author Ed Buckler
 */
public class MigrateHDF5FromT4T5 {
    public static void copyGenotypes(String t4File, String newT5File) {
        IHDF5Reader reader=HDF5Factory.openForReading(t4File);
        IHDF5Writer writer=HDF5Factory.open(newT5File);

        writer.createGroup(Tassel5HDF5Constants.GENOTYPES_MODULE);
        HDF5Utils.unlockHDF5GenotypeModule(writer);
        HDF5Utils.createHDF5TaxaModule(writer);
        HDF5Utils.unlockHDF5TaxaModule(writer);
        int numTaxa = 0;
        HDF5Utils.writeHDF5GenotypesAlleleStates(writer,NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        HDF5Utils.writeHDF5GenotypesMaxNumAlleles(writer,NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        HDF5Utils.writeHDF5GenotypesRetainRareAlleles(writer,false);
        List<HDF5LinkInformation> fields = reader.getAllGroupMemberInformation(HapMapHDF5Constants.GENOTYPES, true);
        for (HDF5LinkInformation is : fields) {
            if (is.isDataSet() == false) continue;
            String taxonName=is.getName();
            System.out.println(taxonName);
            //This is two step copy & then rename.  I couldn't get it to work with one step - it should.
            reader.copy(HapMapHDF5Constants.GENOTYPES+"/"+taxonName, writer,
                    Tassel5HDF5Constants.GENOTYPES_MODULE+"/"+taxonName+"/");
            writer.move(Tassel5HDF5Constants.GENOTYPES_MODULE+"/"+taxonName+"/"+taxonName,
                    Tassel5HDF5Constants.getGenotypesCallsPath(taxonName));
            //copy depth if it exists
            if(reader.exists(HapMapHDF5Constants.DEPTH+"/"+taxonName)) {
                reader.copy(HapMapHDF5Constants.DEPTH+"/"+taxonName, writer,
                        Tassel5HDF5Constants.GENOTYPES_MODULE+"/"+taxonName+"/");
                writer.move(Tassel5HDF5Constants.GENOTYPES_MODULE+"/"+taxonName+"/"+taxonName,
                        Tassel5HDF5Constants.getGenotypesDepthPath(taxonName));
            }
            HDF5Utils.addTaxon(writer,new Taxon(taxonName));
            numTaxa++;
        }
        HDF5Utils.writeHDF5GenotypesNumTaxa(writer,numTaxa);
        HDF5Utils.writeHDF5TaxaNumTaxa(writer,numTaxa);

        //Position module
        writer.createGroup(Tassel5HDF5Constants.POSITION_MODULE);
        int numSites = reader.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SITES);
        HDF5Utils.writeHDF5PositionNumSite(writer,numSites);
        System.out.println(reader.exists(HapMapHDF5Constants.POSITIONS));
        reader.copy(HapMapHDF5Constants.POSITIONS, writer, Tassel5HDF5Constants.POSITIONS);
        reader.copy(HapMapHDF5Constants.LOCI, writer, Tassel5HDF5Constants.CHROMOSOMES);
        reader.copy(HapMapHDF5Constants.LOCUS_INDICES, writer, Tassel5HDF5Constants.CHROMOSOME_INDICES);
        reader.copy(HapMapHDF5Constants.SNP_IDS, writer, Tassel5HDF5Constants.SNP_IDS);

        //Precalculated Stats
        writer.createGroup(Tassel5HDF5Constants.GENO_DESC);
        reader.copy(HapMapHDF5Constants.ALLELE_CNT, writer, Tassel5HDF5Constants.ALLELE_CNT);
        reader.copy(HapMapHDF5Constants.MAF, writer, Tassel5HDF5Constants.MAF);
        reader.copy(HapMapHDF5Constants.SITECOV, writer, Tassel5HDF5Constants.SITECOV);
        reader.copy(HapMapHDF5Constants.ALLELE_FREQ_ORD, writer, Tassel5HDF5Constants.ALLELE_FREQ_ORD);
        reader.copy(HapMapHDF5Constants.TAXACOV, writer, Tassel5HDF5Constants.TAXACOV);
        reader.copy(HapMapHDF5Constants.TAXAHET, writer, Tassel5HDF5Constants.TAXAHET);
        
        HDF5Utils.lockHDF5GenotypeModule(writer);
        HDF5Utils.lockHDF5TaxaModule(writer);
        reader.close();
        writer.close();
    }

}
