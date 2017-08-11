package findhybrids;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.UnknownHostException;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.jtr.transliterate.CharacterParseException;

import MELT.MELTIllumina.checkStats.CheckAnnotation;
import findhybrids.ExonInfo.ExonType;
import findhybrids.HybridCheckMethods.DistributionException;
import findhybrids.HybridCheckMethods.NoDiscordantPairs;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import jaligner.matrix.MatrixLoaderException;
import mergehybrids.MergeHybrids;

public class FindHybrids {

	private static BufferedWriter outFile;
	
	public static void main(String[] args) throws IOException, MatrixLoaderException, NumberFormatException, DistributionException, NoDiscordantPairs, InterruptedException, URISyntaxException, CharacterParseException {

		//Take as first input the sam files of each read aligned separately to SINE and bam file of properly aligned reads
		/*
		 * 0 = Mate 1 -> SINE (.sam - bowtie2) -- Already filtered based on flag 4 (Only MAPPED mates)
		 * 1 = Mate 2 -> SINE (.sam - bowtie2) -- Already filtered based on flag 4 (Only MAPPED mates)
		 * 2 = Reads aligned pair-wise to Ref (.bam - bwa) -- NAME Sorted
		 * 3 = Gene annotation
		 * 4 = output filtered BAM
		 * 5 = DOG fasta
		 * 6 = MEI fasta
		 * 7 = mei Name
		 * 8 = Working DIR
		 * 
		 */
		
		File firstMatesFile = new File(args[0]);
		File secondMatesFile = new File(args[1]);
		File refAlignedFile = new File(args[2]);
		Map<String, IntervalTree<ExonInfo>> annotation = buildAnnotation(new BufferedReader(new FileReader(new File(args[3]))));
		File outputBam = new File(args[4]);
		File fasta = new File(args[5]);
		File fastaIndex = new File(args[5] + ".fai");
		File meiFasta = new File(args[6]);
		File meiFastaIndex = new File(args[6] + ".fai");
		String meiName = args[7];
		
		outFile = new BufferedWriter(new FileWriter(new File(args[8] + "/putative_hits.txt")));
		
		Map<String, Integer> validChrs = processFAI(fastaIndex);
		IndexedFastaSequenceFile dogIndexedFasta = new IndexedFastaSequenceFile(fasta, new FastaSequenceIndex(fastaIndex));
		IndexedFastaSequenceFile meiIndexedFasta = new IndexedFastaSequenceFile(meiFasta, new FastaSequenceIndex(meiFastaIndex));
		
		//Grab traces that have only 1 mate aligned to Mobile Element
		BuildSINEAligned buildFinal = new BuildSINEAligned(firstMatesFile, secondMatesFile);
		Map<String, SAMRecord> finalReads = buildFinal.getFinalReads();
		//Grab reads from the name sorted bam with mates that map to SINE
		buildFilteredBam(finalReads, refAlignedFile, outputBam);
		//Grab putative locations of SINE exonization that intersect with known genes
		List<PutativeHit> putativeHits = GrabPutativeLocations.GrabLocations(outputBam, validChrs);
		
		printSites(putativeHits, dogIndexedFasta, meiIndexedFasta, annotation, new File(outputBam.getAbsolutePath() + "." + meiName + "Assembly/"), meiFastaIndex, meiName);
		outFile.close();
		
	}
	
	public static Map<String, Integer> processFAI(File refSeqIndex) throws NumberFormatException, IOException {

		Map<String,Integer> validChrs = new LinkedHashMap<String, Integer>();
		if (refSeqIndex.isFile()) {
			
			BufferedReader faiReader = new BufferedReader(new FileReader(refSeqIndex));
			String line;
			String data[];
			
			while ((line = faiReader.readLine()) != null) {
				
				data = line.split("\t");
				int length = Integer.parseInt(data[1]); 
				if (length > 1000000) {
					validChrs.put(data[0], length);
				}
			}
			faiReader.close();
		}
		return validChrs;
	}
	private static void buildFilteredBam (Map<String, SAMRecord> finalReads, File refAlignedFile, File outputBam) throws IOException {
		
		SamReader refAlignedReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(refAlignedFile);
		SAMFileHeader newHeader = refAlignedReader.getFileHeader();
		SAMFileWriterFactory samFactory = new SAMFileWriterFactory().setCreateIndex(true);
		newHeader.setSortOrder(SortOrder.coordinate);
		SAMFileWriter pairWriter = samFactory.makeBAMWriter(newHeader, false, outputBam);
				
		SAMRecordIterator samItr = refAlignedReader.iterator();

		while (samItr.hasNext()) {

			SAMRecord current = samItr.next();
			
			String name = current.getReadName();
			boolean isFirstMate = current.getFirstOfPairFlag();
			boolean isSuppAlign = current.getSupplementaryAlignmentFlag();
			boolean notPrimaryAlign = current.getNotPrimaryAlignmentFlag();
			boolean isPairedRead = current.getReadPairedFlag();
			
			if (finalReads.containsKey(name) && isSuppAlign == false && notPrimaryAlign == false && isPairedRead == true) {
				
				SAMRecord MEAligned = finalReads.get(name);
				
				boolean isFirstMateME = MEAligned.getFirstOfPairFlag();
				
				if ((isFirstMate == true && isFirstMateME == false) || (isFirstMate == false && isFirstMateME == true)) {
					
					current.setAttribute("OS", MEAligned.getAlignmentStart());
					current.setAttribute("OE", MEAligned.getAlignmentEnd());
					current.setAttribute("OC", MEAligned.getCigarString());
					current.setAttribute("R2", MEAligned.getReadString());
					current.setAttribute("Q2", MEAligned.getBaseQualityString());
					current.setAttribute("OF", MEAligned.getFlags());
					current.setAttribute("OQ", MEAligned.getMappingQuality());
					current.setAttribute("US", "TRUE");
					
					pairWriter.addAlignment(current);
					
				} else if ((isFirstMate == true && isFirstMateME == true) || (isFirstMate == true && isFirstMateME == true)) {
					
					current.setAttribute("US", "FALSE");
					pairWriter.addAlignment(current);
					
				}
				
				
			}
						
		}
		
		refAlignedReader.close();
		samItr.close();
		pairWriter.close();
		
	}
	
	public static Map<String,IntervalTree<ExonInfo>> buildAnnotation(BufferedReader humRef) throws IOException {
		
		Hashtable<String,IntervalTree<ExonInfo>> geneData = new Hashtable<String,IntervalTree<ExonInfo>>();
		String data[];
		String lengthsArray[];
		String startsArray[];
		String myLine;
		
		while ((myLine = humRef.readLine()) != null) {
			
			data = myLine.split("\t");
			lengthsArray = data[10].split(",");
			startsArray = data[11].split(",");
			int geneStart = Integer.parseInt(data[1]);
			int geneStop = Integer.parseInt(data[2]);
			int codingStart = Integer.parseInt(data[6]);
			int codingStop = Integer.parseInt(data[7]);	
			String geneName = data[3];
			String chr = data[0];
			ExonInfo currentInfo;
			int storedExonNum = 0;
			
			//Determine strandedness
			boolean isForward;
			if (data[5].equals("+")) {
				isForward = true;
			} else {
				isForward = false;
			}
			
			if (codingStart == codingStop) {continue;}
			
			//Add a fake UTR to each gene:
			
			if (isForward) {
				currentInfo = new ExonInfo(geneName, isForward, ExonType.FAKE_FUTR, 1);
				addExonInfo(geneData, chr, currentInfo, (geneStart - 2001), (geneStart - 1));
				currentInfo = new ExonInfo(geneName, isForward, ExonType.FAKE_TUTR, startsArray.length);
				addExonInfo(geneData, chr, currentInfo, (geneStop + 1), (geneStop + 2001));
			} else {
				currentInfo = new ExonInfo(geneName, isForward, ExonType.FAKE_TUTR, startsArray.length);
				addExonInfo(geneData, chr, currentInfo, (geneStart - 2001), (geneStart - 1));
				currentInfo = new ExonInfo(geneName, isForward, ExonType.FAKE_FUTR, 1);
				addExonInfo(geneData, chr, currentInfo, (geneStop + 1), (geneStop + 2001));
			}
			
			for (int exons = 0; exons <= (startsArray.length - 1); exons++) {
				
				int currentStart = geneStart + Integer.parseInt(startsArray[exons]);
				int currentStop = currentStart + Integer.parseInt(lengthsArray[exons]);
				currentStart++;
				
				if (isForward) {
					storedExonNum = exons + 1;
				} else {
					storedExonNum = startsArray.length - (exons);
				}
				
				if (isForward) {
					
					//Add or subtract 2kbp to 3' and 5' UTR due to poor annotation of dog UTRs
					
					if (currentStart < codingStart && currentStop > codingStart) { //EXON with mixed coding & non-coding (5')
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.FUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, codingStart);
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.EXON, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, codingStart + 1, currentStop);
						
					} else if (currentStop <= codingStart) {
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.FUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, currentStop);
	
						
					} else if (currentStop > codingStop && currentStart < codingStop) { //EXON with mixed coding & non-coding (3')
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.EXON, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, codingStop);
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.TUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, codingStop + 1, currentStop);
						
					} else if (currentStart >= codingStop) {
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.TUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, currentStop);

					} else {
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.EXON, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, currentStop);
						
					}
					
				} else {
					
					if (currentStart < codingStart && currentStop > codingStart) { //EXON with mixed coding & non-coding (3')
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.TUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, codingStart);					
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.EXON, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, codingStart + 1, currentStop);
						
					} else if (currentStop <= codingStart) { //JUST UTR
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.TUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, currentStop);
						
					} else if (currentStop > codingStop && currentStart < codingStop) { //EXON with mixed coding & non-coding (5')
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.EXON, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, codingStop);
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.FUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, codingStop + 1, currentStop);
						
					} else if (currentStart >= codingStop) { //JUST UTR
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.FUTR, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, currentStop);
												
					} else {
						
						currentInfo = new ExonInfo(geneName, isForward, ExonType.EXON, storedExonNum);
						addExonInfo(geneData, chr, currentInfo, currentStart, currentStop);
						
					}
										
				}
						
			}
						
		}
				
		humRef.close();
		return geneData;
		
	}
	private static void addExonInfo(Map<String,IntervalTree<ExonInfo>> geneData, String chr, ExonInfo currentInfo, int start, int stop) {
		if (geneData.containsKey(chr)) {
			IntervalTree<ExonInfo> positionInfo = geneData.get(chr);
			positionInfo.put(start, stop, currentInfo);
			geneData.put(chr, positionInfo);
		} else {
			IntervalTree<ExonInfo> positionInfo = new IntervalTree<ExonInfo>();
			positionInfo.put(start, stop, currentInfo);
			geneData.put(chr, positionInfo);
		}
	}
	
	private static void printSites (List<PutativeHit> hits, IndexedFastaSequenceFile dogIndexedFasta, IndexedFastaSequenceFile meiIndexedFasta, Map<String, IntervalTree<ExonInfo>> annotation, File SINEAssemblies, File firstMatesFile, String meiName) throws UnknownHostException, IOException, MatrixLoaderException, NumberFormatException, DistributionException, NoDiscordantPairs, InterruptedException, URISyntaxException, CharacterParseException {
		
		DecimalFormat format = new DecimalFormat("#.###");
		for (PutativeHit hit : hits) {
			
			HybridCheckMethods methods = new HybridCheckMethods(hit, dogIndexedFasta, meiIndexedFasta, meiName, annotation, SINEAssemblies, firstMatesFile);
			CheckAnnotation finalHit = methods.call();
			if (methods.getGeneInfo() != null) {
				if (methods.getGeneInfo().getType() == ExonType.EXON) {
					
					int start = determineStart(finalHit.getLeftposition(), finalHit.getRightposition(), finalHit.getQual(), finalHit.getTsd());					
					outFile.write(hit.getChr() + "\t" + start + "\t" + start + "\t" + hit.getTotalReads() + "\t" + format.format(hit.getPercDisc()) + "\t" + methods.getTotalSplit() + "\t" + methods.getTotalPassSplit() + "\t" + methods.getTotalDiscordant() + "\t" + finalHit.getQual() + "\t" + (int) finalHit.getInsertSizeStart() + "\t" + (int) finalHit.getInsertSizeStop() + "\t" + methods.getGeneInfo().getGeneName() + "\t" + methods.getGeneInfo().getType() + "\t" + methods.getGeneInfo().getExonNum() + "\t" + MergeHybrids.getDirection(finalHit.getTotalLeft(), finalHit.getTotalRight()) + "\n");
					outFile.flush();

				}
			}
		}
	}
	private static int determineStart (double left, double right, int qual, String tsd) {
		
		int start;
		
		if (qual == 5) {
			if (tsd.matches("d[ATGC]+") || tsd.matches("0tsd")) {
				start = (int) left;
			} else {
				start = (int) right;
			}
		} else if (qual == 3) {
			start = (int) left;
		} else if (qual == 2) {
			if (right == 0) {
				start = (int) left;
			} else {
				start = (int) right;
			}
		} else {
			start = (int) right;
		}
				
		return start;
	}

}
