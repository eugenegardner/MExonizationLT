package findhybrids;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.DNASequence;
import org.jtr.transliterate.CharacterParseException;

import MELT.MELTIllumina.checkStats.CheckAnnotation;
import MELT.MELTIllumina.checkStats.CheckAnnotation.Strand;
import MELT.MELTIllumina.checkStats.CheckClasses;
import MELT.MELTIllumina.classification.generic.GenericResults;
import MELT.MELTIllumina.classification.generic.GenericU;
import MELT.utilities.SAMTools;
import MELT.utilities.SAMTools.AlignmentDecision;
import MELT.utilities.SAMTools.polADecision;
import findhybrids.ExonInfo.ExonType;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import jaligner.Sequence;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import jaligner.ui.filechooser.NamedInputStream;

public class HybridCheckMethods {
	
	// Information about the insertion:
	private boolean posForSplit;
	private int totalPassSplit;
	private int totalSplit;
	private int totalDiscordant;
	private PutativeHit exonization;
	private Map<String, IntervalTree<ExonInfo>> annotation;
	private File SINEDumpFolder;
	
	// SAMTools iterator to use for determine pair status:
	private Iterator<SAMRecord> verticalItr;
	private Hashtable<Character,Integer> qualityHash;
	
	// SAMTools handle for running specific tools:
	private SAMTools tools;
	private IndexedFastaSequenceFile reference;
	
	// Handles for supplemental classes:
	private CheckClasses holder;
	private CheckAnnotation annote;
	
	// Handles for subclasses of CheckClasses:
	private CheckClasses.lengthInfo myLength;
	private CheckClasses.supInfo mySupInfo;
	private ArrayList<Integer> starts;
	private ArrayList<Integer> stops;
	
	// Lists to hold information about site:
	private List<SAMRecord> lefts;
	private List<Integer> leftPos;
	private List<SAMRecord> rights;
	private List<Integer> rightPos;
	private ExonInfo geneInfo;
	
	// Loaders for mobile element determination:
	private Sequence meiLoaded;
	private Matrix matrix;
	@SuppressWarnings("unused")
	private int total2X;
	
	// Patterns used to ID specific split read information:
	private Pattern cigarCaptureRight = Pattern.compile("(\\S+[A-Z]{1})(\\d+)S$");
	private Pattern cigarCaptureLeft = Pattern.compile("^(\\d+)S(\\S+)");
	
	// HashMap/Tables that hold mate information:
	private Map<String,SAMRecord> discordantPairs;
	
	// HashMap for calling InDels in Species Finder:
	private HashMap<Integer,CheckClasses.InDelInfo> deletions;
	private HashMap<Integer,CheckClasses.InDelInfo> insertions;
	private List<CheckClasses.cutter> toCut;
	
	//Doubles for determining SR status of DPs
	private double splitDP;
	private double totalDP;
	
	//Handles for making SINE assemblies:
	private SINERecordFactory SINEfactory;
	private File firstMatesFile;
	private String species;
	
	// Constructor method, a little different than my stuff before, just sets up all handles needed for the entire class
	public HybridCheckMethods(PutativeHit exonization, IndexedFastaSequenceFile reference, IndexedFastaSequenceFile mei, String species, Map<String, IntervalTree<ExonInfo>> annotation, File SINEDumpFolder, File firstMatesFile) throws MatrixLoaderException, IOException {
		
		// Declare Variables:
		this.reference = reference;
		this.exonization = exonization;
		this.annote = new CheckAnnotation();
		this.annotation = annotation;
		this.SINEDumpFolder = SINEDumpFolder;
		this.firstMatesFile = firstMatesFile;
		this.species = species;
		matrix = MatrixLoader.load(new NamedInputStream("NAMatrix",this.getClass().getResourceAsStream("/MELT/MELTIllumina/importFiles/nuc-4_4.txt")));
		totalPassSplit = 0;
		totalSplit = 0;
		totalDiscordant = 0;
		holder = new CheckClasses();
		myLength = holder.new lengthInfo();
		mySupInfo = holder.new supInfo();
		tools = new SAMTools();
		// Get hashtable for insert sequences
		lefts = new ArrayList<SAMRecord>();
		leftPos = new ArrayList<Integer>();
		rights = new ArrayList<SAMRecord>();
		rightPos = new ArrayList<Integer>();
		starts = new ArrayList<Integer>();
		stops = new ArrayList<Integer>();
		discordantPairs = new HashMap<String,SAMRecord>();
		splitDP = 0;
		totalDP = 0;
		qualityHash = tools.loadQualities(false);
				
		String transSeq = new String(mei.getSequence(species).getBases()).toUpperCase();
		meiLoaded = new Sequence(transSeq);
		meiLoaded.setId(species);
		
		deletions = new HashMap<Integer,CheckClasses.InDelInfo>();
		insertions = new HashMap<Integer,CheckClasses.InDelInfo>();
		toCut = new ArrayList<CheckClasses.cutter>();
		
	}
	
	public CheckAnnotation call() throws IOException, findhybrids.HybridCheckMethods.DistributionException, findhybrids.HybridCheckMethods.NoDiscordantPairs, NumberFormatException, InterruptedException, MatrixLoaderException, URISyntaxException, CharacterParseException {
		
		// Do actual computation
		verticalItr = exonization.getSupportingReads().iterator();
		posForSplit = useSplitReads();
		catalogDisc(posForSplit);
		
		//Get length:
		myLength = getLength();
		annote.setInsertSizeStart(myLength.getInsertSizeStart());
		annote.setInsertSizeStop(myLength.getInsertSizeStop());
	
		// Do all the statistic finding required by the class associated with this class:		
		if (posForSplit == true) {
			
			splitReads polyAInfo = assessSplit();

			if (totalPassSplit > 0) {

				// Check to see if both pairs are defined:
				annote.setStrand(determineStrandednessWithSplit(polyAInfo));
				
				// If insert has split reads, use these rather than discordant pairs to determine insert length.
				if (annote.getStrand() == Strand.POS) {
					
				} else if (annote.getStrand() == Strand.NEG) {

				} else {
					if (getStrandednessWithDiscordant().equals("+")) {
						annote.setStrand(Strand.POS);
					} else {
						annote.setStrand(Strand.NEG);
					}
				}
				
				determineTSD();
				
			//This is if, after checking for GOOD SRs, we don't actually find any
			} else {
				recoverTSD();
				if (getStrandednessWithDiscordant().equals("+")) {
					annote.setStrand(Strand.POS);
				} else {
					annote.setStrand(Strand.NEG);
				}
				annote.setTwinPrime(-1);
				// Check Insert Length Info
			}
				
		} else {
			
			if (getStrandednessWithDiscordant().equals("+")) {
				annote.setStrand(Strand.POS);
			} else {
				annote.setStrand(Strand.NEG);
			}
			annote.setTwinPrime(-1);
			
		}
		
		// Determine mutations
		int pos = determineStart(annote.getLeftposition(), annote.getRightposition(), annote.getQual(), annote.getTsd());
		SINEfactory = new SINERecordFactory(SINEDumpFolder, exonization.getChr() + "_" + pos, firstMatesFile, species);
		String consensus = getSpecies().replaceAll("N", "");
		annote.setConsensus(consensus);
		GenericResults result = getGenericSpecies(consensus);
		annote.setInsertSizeStart(result.getStart());
		annote.setInsertSizeStop(result.getStop());
		
		SINEfactory.closeWriter();
		
		//Get gene location
		geneInfo = checkGene();
				
		//Set Discordant Pair stats:
		annote.setDpInfo((splitDP/totalDP) * 100);
		return annote;
	
	}
	
	//Get number of SRs:
	public int getTotalSplit () {
		return totalSplit;
	}
	public int getTotalPassSplit() {
		return totalPassSplit;
	}
	public int getTotalDiscordant() {
		return totalDiscordant;
	}
	public ExonInfo getGeneInfo() {
		return geneInfo;
	}
	
	// Catalog split and discordant reads (place them into different groups):
	private boolean useSplitReads() {
		
		int counter = 0;
		Pattern cigarCaptureRight = Pattern.compile("\\S+[A-Z]{1}(\\d+)S$");
		Pattern cigarCaptureLeft = Pattern.compile("^(\\d+)S\\S+");
		
		while(verticalItr.hasNext() == true) {
			
			SAMRecord currentRecord = verticalItr.next();
			if (currentRecord.getReadUnmappedFlag() == true || currentRecord.getMappingQuality() == 0) {
				continue;
			}
			String cigar = currentRecord.getCigarString();
			
			//I believe all reads should be treated as discordant. Need to check!
			if (currentRecord.getMateReferenceName() != currentRecord.getReferenceName()) {
				totalDiscordant++;
			}
			Matcher right = cigarCaptureRight.matcher(cigar);
			Matcher left = cigarCaptureLeft.matcher(cigar);
			if (right.matches() && left.matches()) {
				//continue;
			} else if (right.matches() || left.matches()) {
				splitDP++;
				totalDP++;
			} else {
				totalDP++;
			}
			discordantPairs.put(currentRecord.getReadName(), currentRecord);
			
			counter += catalogReads(cigar, currentRecord);
			
//			}
			
		}
		
		if (counter == 0) {
		
			return false;
		
		} else {
			
			return true;
			
		}
		
	}
	private int catalogReads (String cigar, SAMRecord currentRecord) {
		
		// Store Match for Left and Right
		Matcher cigarRight = cigarCaptureRight.matcher(cigar);
		Matcher cigarLeft = cigarCaptureLeft.matcher(cigar);
		
		if (cigarRight.matches() && cigarLeft.matches()) {
			
			return 0;
			
		} else if (cigarRight.matches()) {
			
			String qualityString = currentRecord.getBaseQualityString().substring(currentRecord.getReadLength() - Integer.parseInt(cigarRight.group(2)), currentRecord.getReadLength());
			double quality = tools.qualityString(qualityString, qualityHash);
			if (quality > 5) {
				totalSplit++;
				rights.add(currentRecord);
			} else {
				return 0;
			}
					
		} else if (cigarLeft.matches()) {
				
			String qualityString = currentRecord.getBaseQualityString().substring(0, Integer.parseInt(cigarLeft.group(1)) - 1);
			double quality = tools.qualityString(qualityString, qualityHash);
			if (quality > 5) {
				totalSplit++;
				lefts.add(currentRecord);
			} else {
				return 0;
			}
		} else {
			return 0;
		}

		return 1;
		
	}
	private void recoverTSD () throws NoDiscordantPairs {
		
		int counter = 0;
		ArrayList<Integer> rightStart = new ArrayList<Integer>();
		ArrayList<Integer> leftStop = new ArrayList<Integer>();
		
		for(Map.Entry<String, SAMRecord> entry : discordantPairs.entrySet()) {
			
			counter++;
			SAMRecord currentRecord = entry.getValue();
			boolean readStrand = currentRecord.getReadNegativeStrandFlag();
			
			if (readStrand == true) {
				rightStart.add(currentRecord.getAlignmentStart());
			} else {
				leftStop.add(currentRecord.getAlignmentEnd());
			}
			
		}
		
		if (counter > 0) {
			Collections.sort(rightStart);
			Collections.sort(leftStop);
			ReferenceSequence refContig;
			if (rightStart.size() == 0 || leftStop.size() == 0) {
				if (rightStart.size() == 0) {
					annote.setRightposition(0);
					annote.setLeftposition(leftStop.get(leftStop.size() - 1) - 1);
					annote.setQual(2);
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long) annote.getLeftposition(), (long) annote.getLeftposition());
					annote.setAlt(new String(refContig.getBases()));
				} else if (leftStop.size() == 0) {
					annote.setLeftposition(0);
					annote.setRightposition(rightStart.get(0));
					annote.setQual(2);
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long) annote.getRightposition(), (long) annote.getRightposition());
					annote.setAlt(new String(refContig.getBases()));
				}
			} else {
				long rightposition = leftStop.get(leftStop.size() - 1);
				long leftposition = rightStart.get(0);
				if (leftposition < rightposition) {
					long distance = rightposition - leftposition;
					if (distance <= 30) {
						annote.setRightposition(rightposition + 1);
						annote.setLeftposition(leftposition);
						refContig = reference.getSubsequenceAt(exonization.getChr(), rightposition+1, rightposition+1);
						annote.setAlt(new String(refContig.getBases()));
						refContig = reference.getSubsequenceAt(exonization.getChr(), leftposition, rightposition);
						annote.setTsd(new String(refContig.getBases()));
						annote.setQual(2);
					} else {
						annote.setRightposition(Math.round((rightposition + leftposition) / 2));
						annote.setLeftposition(0);
						refContig = reference.getSubsequenceAt(exonization.getChr(), rightposition, rightposition);
						annote.setAlt(new String(refContig.getBases()));
						annote.setQual(1);
					}
					
				} else {
					annote.setRightposition(Math.round((rightposition+leftposition) / 2));
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long)annote.getRightposition(), (long)annote.getRightposition());
					annote.setAlt(new String(refContig.getBases()));
					annote.setLeftposition(0);
					annote.setQual(0);
				}
				
			}
			
		} else {
			throw new NoDiscordantPairs();
		}
		
	}
	
	private boolean catalogDisc(boolean determineTSD) throws DistributionException {
		
		ArrayList<SAMRecord> speciesCheck = new ArrayList<SAMRecord>();
		ArrayList<Integer> lengthCheck = new ArrayList<Integer>();
		Hashtable<Integer,Boolean> dirCheck = new Hashtable<Integer,Boolean>();
		ArrayList<Integer> rightStart = new ArrayList<Integer>();
		ArrayList<Integer> leftStop = new ArrayList<Integer>();
		int counter = 0;
		int totalLeft = 0;
		int totalRight = 0;
		DescriptiveStatistics avgReadLength = new DescriptiveStatistics();
		DescriptiveStatistics avgStart = new DescriptiveStatistics();
		DescriptiveStatistics avgEnd = new DescriptiveStatistics();
		
		for (Map.Entry<String, SAMRecord> entry : discordantPairs.entrySet()) {
			
			counter++;
			SAMRecord currentRecord = entry.getValue();
			avgReadLength.addValue(currentRecord.getReadLength());
			boolean readStrand = currentRecord.getReadNegativeStrandFlag();
			if (readStrand == false) {
				totalLeft++;
			} else {
				totalRight++;
			}
			//Get width of the DP region
			avgStart.addValue(currentRecord.getStart());
			avgEnd.addValue(currentRecord.getEnd());
			// add insertion based tags to the relevant supInfo (OS, OE, ST, SE, QA, OC):

			boolean isNegative = currentRecord.getReadNegativeStrandFlag();
			// Add info to the extra fields for later!
			dirCheck.put(currentRecord.getIntegerAttribute("OE"), isNegative);
			// Add in a predictable pattern to the array, start than stop!
			
			speciesCheck.add(currentRecord);
			
			int meiStart = currentRecord.getIntegerAttribute("OS");
			int meiStop = currentRecord.getIntegerAttribute("OE");
			String cigar = currentRecord.getStringAttribute("OC");
			
			Matcher right = cigarCaptureRight.matcher(cigar);
			Matcher left = cigarCaptureLeft.matcher(cigar);
//			private Pattern cigarCaptureRight = Pattern.compile("(\\S+[A-Z]{1})(\\d+)S$");
//			private Pattern cigarCaptureLeft = Pattern.compile("^(\\d+)S(\\S+)");
			//2 right 1 left
			if (right.matches() && left.matches()) {
				if (Integer.parseInt(right.group(2)) >= 5 && Integer.parseInt(left.group(1)) >= 5) {
					starts.add(meiStart);
					stops.add(meiStop);
				}
			} else if (right.matches() && !left.matches()) {
				if (Integer.parseInt(right.group(2)) >= 5) {
					stops.add(meiStop);
				}
			} else if (!right.matches() && left.matches()) {
				if (Integer.parseInt(left.group(1)) >= 5) {
					starts.add(meiStart);
			
				}
			}
			if (meiStop == meiLoaded.length()) {
				stops.add(meiStop);
			}
			if (meiStart == 1) {
				starts.add(meiStart);
			}
			
				
			if (determineTSD == false) {
				if (readStrand == true) {
					rightStart.add(currentRecord.getAlignmentStart());
				} else {
					leftStop.add(currentRecord.getAlignmentEnd());
				}
			}
		}
		
		if (totalLeft == 0) {
			annote.setTotalLeft(0);
			annote.setTotalRight(totalRight);
		} else if (totalRight == 0) {
			annote.setTotalLeft(totalLeft);
			annote.setTotalRight(0);
		} else {
			annote.setTotalLeft(totalLeft);
			annote.setTotalRight(totalRight);
		}
		
		// Check to see if the standard deviation of the trimmed mean of both start and stop is below our threshold (current 35.0 bp)
		double[] startSort = avgStart.getSortedValues();
		double[] endSort = avgEnd.getSortedValues();
		avgStart.clear();
		avgEnd.clear();
		double toTrim = counter * 0.1;
		int start = (int) Math.round(toTrim);
		for (int x = 0 + start; x < (counter - start); x++) {
			avgStart.addValue(startSort[x]);
			avgEnd.addValue(endSort[x]);
		}
		
		if (determineTSD == false) {
			Collections.sort(rightStart);
			Collections.sort(leftStop);
			ReferenceSequence refContig;
			if (rightStart.size() == 0 || leftStop.size() == 0) {
				if (rightStart.size() == 0) {
					annote.setRightposition(0);
					annote.setLeftposition(leftStop.get(leftStop.size() - 1) - 1);
					annote.setQual(2);
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long) annote.getLeftposition(), (long) annote.getLeftposition());
					annote.setAlt(new String(refContig.getBases()));
				} else if (leftStop.size() == 0) {
					annote.setLeftposition(0);
					annote.setRightposition(rightStart.get(0));
					annote.setQual(2);
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long) annote.getRightposition(), (long) annote.getRightposition());
					annote.setAlt(new String(refContig.getBases()));
				}
			} else {
				long rightposition = leftStop.get(leftStop.size() - 1);
				long leftposition = rightStart.get(0);
				if (leftposition < rightposition) {
					long distance = rightposition - leftposition;
					if (distance <= 30) {
						annote.setRightposition(rightposition + 1);
						annote.setLeftposition(leftposition);
						refContig = reference.getSubsequenceAt(exonization.getChr(), rightposition+1, rightposition+1);
						annote.setAlt(new String(refContig.getBases()));
						refContig = reference.getSubsequenceAt(exonization.getChr(), leftposition, rightposition);
						annote.setTsd(new String(refContig.getBases()));
						annote.setQual(2);
					} else {
						annote.setRightposition(Math.round((rightposition + leftposition) / 2));
						annote.setLeftposition(0);
						refContig = reference.getSubsequenceAt(exonization.getChr(), rightposition, rightposition);
						annote.setAlt(new String(refContig.getBases()));
						annote.setQual(1);
					}
					
				} else {
					annote.setRightposition(Math.round((rightposition+leftposition) / 2));
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long)annote.getRightposition(), (long)annote.getRightposition());
					annote.setAlt(new String(refContig.getBases()));
					annote.setLeftposition(0);
					annote.setQual(0);
				}
				
			}
		}
		
		mySupInfo.setDirCheck(dirCheck);
		mySupInfo.setLengthCheck(lengthCheck);
		mySupInfo.setSpeciesCheck(speciesCheck);
		return false;
		
	}
	
	//Determine where split reads go and help figure out Strandedness
	private splitReads assessSplit() throws CompoundNotFoundError {
		
		splitReads info = new splitReads();

		polADecision polyForDecision;
		
		for (SAMRecord left : lefts) {
			
			Matcher matched = cigarCaptureLeft.matcher(left.getCigarString());
			if (matched.matches()) {
				
				Sequence pulledCigar = new Sequence(left.getReadString().substring(0, Integer.parseInt(matched.group(1)) - 1));
				Sequence pulledCigarRev = new Sequence(new DNASequence(pulledCigar.getSequence()).getReverseComplement().getSequenceAsString());
				
				polyForDecision = tools.callPolyA(pulledCigarRev);
				if (polyForDecision.getPolDec() == true) {
					info.incLeftPolyTrue();
					leftPos.add(left.getAlignmentStart());
					totalPassSplit++;
				} else {
					if (testSplit(pulledCigarRev)) {
						info.incLeftPolyFalse();
						leftPos.add(left.getAlignmentStart());
						totalPassSplit++;
					}
				}
			}
		}
				
		for (SAMRecord right : rights) {

			Matcher matched = cigarCaptureRight.matcher(right.getCigarString());
			if (matched.matches()) {
				Sequence pulledCigar = new Sequence(right.getReadString().substring(right.getReadLength() - Integer.parseInt(matched.group(2)), right.getReadLength()));
				polyForDecision = tools.callPolyA(pulledCigar);
				if (polyForDecision.getPolDec() == true) {
					rightPos.add(right.getAlignmentEnd());
					info.incRightPolyTrue();
					totalPassSplit++;
				} else {
					if (testSplit(pulledCigar)) {
						info.incRightPolyFalse();
						rightPos.add(right.getAlignmentEnd());
						totalPassSplit++;
					}
				}
			}
		}
		
		
		return info;
		
	}
	private class splitReads {
		
		int leftPolyTrue;
		int leftPolyFalse;
		int rightPolyTrue;
		int rightPolyFalse;
		
		splitReads () {
			this.leftPolyTrue = 0;
			this.leftPolyFalse = 0;
			this.rightPolyTrue = 0;
			this.rightPolyFalse = 0;
		}

		public void incLeftPolyTrue() {
			leftPolyTrue++;
		}
		public void incLeftPolyFalse() {
			leftPolyFalse++;
		}
		public void incRightPolyTrue() {
			rightPolyTrue++;
		}
		public void incRightPolyFalse() {
			rightPolyFalse++;
		}
		public int getLeftPolyTrue() {
			return leftPolyTrue;
		}
		public int getLeftPolyFalse() {
			return leftPolyFalse;
		}
		public int getRightPolyTrue() {
			return rightPolyTrue;
		}
		public int getRightPolyFalse() {
			return rightPolyFalse;
		}

	}

	private boolean testSplit(Sequence toTest) throws CompoundNotFoundError {
		
		boolean pass = false;
		
		int forLength = toTest.length();
		if (forLength > 5) {
			Sequence revLoaded = new Sequence(new DNASequence(toTest.getSequence()).getReverseComplement().getSequenceAsString());
			AlignmentDecision insDecision = tools.alignmentMatrix(toTest,revLoaded,meiLoaded,forLength,90, matrix);
			if (insDecision == AlignmentDecision.FOR || insDecision == AlignmentDecision.REV) {
				pass = true;
			}
		}
				
		return pass;
		
	}
	
	// Get Strand Methods:
	// Additional Method for getting strandedness if no split reads:
	private String getStrandednessWithDiscordant () {
		
		String myStrandedness = "";
		Hashtable<Integer,Boolean> dirCheck = mySupInfo.getDirCheck();
		List<Integer> keys = new LinkedList<Integer>(dirCheck.keySet());
		Collections.sort(keys);
		
		for(Integer key: keys) {
			if (dirCheck.get(key) == true) {
				myStrandedness = "+";
				break;
			} else {
				myStrandedness = "-";
				break;
			}
		}
						
		return myStrandedness;
	}
	//Determine strand with SRs:
	public Strand determineStrandednessWithSplit(splitReads polAInfo) {
		Strand strand;
		//Check to see if both are set:
		if ((polAInfo.getLeftPolyTrue() > 0 || polAInfo.getLeftPolyFalse() > 0) && (polAInfo.getRightPolyTrue() > 0 || polAInfo.getRightPolyFalse() > 0)) {
			if (polAInfo.getLeftPolyTrue() < polAInfo.getLeftPolyFalse() && polAInfo.getRightPolyTrue() > polAInfo.getRightPolyFalse()) {
				strand = Strand.NEG;
			} else if (polAInfo.getLeftPolyTrue() > polAInfo.getLeftPolyFalse() && polAInfo.getRightPolyTrue() < polAInfo.getRightPolyFalse()) {
				strand = Strand.POS;
			} else {
				strand = Strand.NULL;
			}
		// If not, find out which one is (left first):	
		} else if ((polAInfo.getLeftPolyTrue() > 0 || polAInfo.getLeftPolyFalse() > 0) && (polAInfo.getRightPolyTrue() == 0 && polAInfo.getRightPolyFalse() == 0)) {
			// Now ask which one is present:
			if (polAInfo.getLeftPolyTrue() > polAInfo.getLeftPolyFalse()) {
				strand = Strand.POS;
			} else {
				strand = Strand.NEG;
			}
		// Now Right:
		} else if ((polAInfo.getLeftPolyTrue() == 0 && polAInfo.getLeftPolyFalse() == 0) && (polAInfo.getRightPolyTrue() > 0 || polAInfo.getRightPolyFalse() > 0)) {
			if (polAInfo.getRightPolyTrue() > polAInfo.getRightPolyFalse()) {
				strand = Strand.NEG;
			} else {
				strand = Strand.POS;
			}
		} else {
			strand = Strand.NULL;
		}
		return strand;
	}
	
	//Determine where the TSD is using SRs
	private void determineTSD() {
		
		Map<Integer,Integer> leftMap;
		Map<Integer,Integer> rightMap;
		List<Integer> leftModes = new ArrayList<Integer>();
		List<Integer> rightModes = new ArrayList<Integer>();
		double leftposition;
		double rightposition;
		ReferenceSequence refContig;
		int rightToUse = 0;
		int leftToUse = 0;
		
		rightMap = getMode(rightPos);
		rightModes.addAll(rightMap.keySet());
		leftMap = getMode(leftPos);
		leftModes.addAll(leftMap.keySet());
		if (leftPos.size() == 0 && rightPos.size() > 0) {
			annote.setQual(4);
		} else if (leftPos.size() > 0 && rightPos.size() == 0) {
			annote.setQual(3);
		} else {
			if (rightModes.size() == 1 && leftModes.size() == 1) {
				if (rightModes.get(0).equals(leftModes.get(0))) {
					annote.setQual(5);
				} else if ((leftModes.get(0) < (rightModes.get(0) - 30)) || (rightModes.get(0) < (leftModes.get(0) - 30))) {
					if (rightMap.get(rightModes.get(0)) > leftMap.get(leftModes.get(0))) {
						annote.setQual(4);
					} else if (rightMap.get(rightModes.get(0)) < leftMap.get(leftModes.get(0))) {
						annote.setQual(3);
					} else {
						annote.setQual(1);
					}
				} else {
					rightToUse = 0;
					leftToUse = 0;
					annote.setQual(5);
				}
			} else {
				if ((rightModes.size() > 1 && rightModes.size() < 4) && leftModes.size() == 1) {
					leftposition = leftModes.get(0) + 1;
					leftToUse = 0;
					for (int x = 0; x <= rightModes.size(); x++) {
						if (x < rightModes.size()) {
							if ((rightModes.get(x) < (leftposition + 30)) && (rightModes.get(x) > leftposition)) {
								rightToUse = x;
								annote.setQual(5);
								break;
							}
						} else {
							annote.setQual(3);
						}
					}
				} else if (rightModes.size() == 1 && (leftModes.size() > 1 && leftModes.size() < 4)) {
					rightposition = rightModes.get(0);
					rightToUse = 0;
					for (int x = 0; x <= leftModes.size(); x++) {
						if (x < leftModes.size()) {
							if ((leftModes.get(x) > (rightposition - 30)) && (leftModes.get(x) < rightposition)) {
								leftToUse = x;
								annote.setQual(5);
								break;
							}
						} else {
							annote.setQual(4);						
						}
					}
				} else if (rightModes.size() >= 4 && leftModes.size() == 1) {
					annote.setQual(3);
				} else if (leftModes.size() >= 4 && rightModes.size() == 1) {
					annote.setQual(4);
				} else {
					annote.setQual(1);
				}
			}
		}
		
		switch (annote.getQual()) {
		case 1:
			annote.setRightposition((leftModes.get(0) + rightModes.get(0)) / 2);
			refContig = reference.getSubsequenceAt(exonization.getChr(), (long)annote.getRightposition(), (long)annote.getRightposition());
			annote.setAlt(new String(refContig.getBases()).toUpperCase());
			break;
		case 3:
			leftposition = leftModes.get(0) - 1;
			annote.setLeftposition(leftposition);
			annote.setRightposition(0);
			refContig = reference.getSubsequenceAt(exonization.getChr(), (long)annote.getLeftposition(), (long)annote.getLeftposition());
			annote.setAlt(new String(refContig.getBases()).toUpperCase());
			break;
		case 4:
			rightposition = rightModes.get(0);
			annote.setRightposition(rightposition);
			refContig = reference.getSubsequenceAt(exonization.getChr(), (long)rightposition, (long)rightposition);
			annote.setAlt(new String(refContig.getBases()).toUpperCase());
			annote.setLeftposition(0);
			break;
		case 5:
			rightposition = rightModes.get(rightToUse);
			leftposition = leftModes.get(leftToUse);
			if (rightposition == leftposition) {
				annote.setRightposition(rightposition);
				annote.setLeftposition(leftposition);
				refContig = reference.getSubsequenceAt(exonization.getChr(), (long)rightposition, (long)rightposition);
				annote.setTsd(new String(refContig.getBases()).toUpperCase());
				refContig = reference.getSubsequenceAt(exonization.getChr(), (long)rightposition, (long)rightposition);
				annote.setAlt(new String(refContig.getBases()).toUpperCase());
			} else {
				if (rightposition < leftposition) {
					//This is a TS Deletion
					annote.setRightposition(leftposition);
					annote.setLeftposition(rightposition);
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long)rightposition + 1, (long)leftposition - 1);
					if (new String(refContig.getBases()).length() > 0) {
						annote.setTsd("d" + new String(refContig.getBases()).toUpperCase());
					} else {
						annote.setTsd("0tsd");
					}
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long)leftposition, (long)leftposition);
					annote.setAlt(new String(refContig.getBases()).toUpperCase());
				} else {
					//This is a TS Duplication
					annote.setRightposition(rightposition);
					annote.setLeftposition(leftposition);
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long)leftposition, (long)rightposition);
					if (new String(refContig.getBases()).length() > 0) {
						annote.setTsd(new String(refContig.getBases()).toUpperCase());
					} else {
						annote.setTsd("0tsd");
					}
					refContig = reference.getSubsequenceAt(exonization.getChr(), (long)rightposition, (long)rightposition);
					annote.setAlt(new String(refContig.getBases()).toUpperCase());
				}
				
			}
			
			break;
		default:
			System.err.println("Fatal Error involving TSD!");
			System.exit(1);
			break;
		}
		
	}
	
	// Determine length of insertion using DPs...
	private CheckClasses.lengthInfo getLength() {
		
		CheckClasses.lengthInfo myLength = holder.new lengthInfo();
				
		if (starts.size() > 0) {
			
			Collections.sort(starts);
			Map<Integer, Integer> startMap = getMode(starts);
			List<Integer> startModes = new ArrayList<Integer>(); 
			startModes.addAll(startMap.keySet());
					
			if (startModes.size() >= 1) {
				myLength.setInsertSizeStart(startModes.get(0));
			} else {
				myLength.setInsertSizeStart(0);
			}
				
		} else {
			myLength.setInsertSizeStart(0);
		}
			
		if (stops.size() > 0) {
			
			Collections.sort(stops);
			Map<Integer, Integer> stopMap = getMode(stops);
			List<Integer> stopModes = new ArrayList<Integer>(); 
			stopModes.addAll(stopMap.keySet());
			
			if (stopModes.size() >= 1) {
				myLength.setInsertSizeStop(stopModes.get(0));
			} else {
				myLength.setInsertSizeStop(meiLoaded.length() - 1);
			}
		} else {
			myLength.setInsertSizeStop(meiLoaded.length() - 1);
		}
		
		return myLength;
		
	}
	
	// Get species information
	private String getSpecies() throws IOException, InterruptedException, MatrixLoaderException, URISyntaxException, NumberFormatException {
		
		Hashtable<Character,Integer> qualityHash = tools.loadQualities(false);
		HashMap<Integer,String> mutations = new LinkedHashMap<Integer,String>();
		
		//Populate the map, so that I can add stuff to it
		for(int parser = 1; parser <= meiLoaded.length(); parser++) {
			mutations.put(parser, "");
		}
		
		for (SAMRecord thisOne : mySupInfo.getSpeciesCheck()) {
			
			char[] seqSubArray = thisOne.getStringAttribute("R2").toCharArray();
			char[] seqQualArray = thisOne.getStringAttribute("Q2").toCharArray();
			
			SINEfactory.AddSINEAlignment(thisOne);
			
			List<CheckClasses.AlignmentBlob> blobs = getReadBlocks(buildCigar(thisOne.getStringAttribute("OC")), thisOne.getIntegerAttribute("OS"), thisOne.getReadName());			
						
			for (CheckClasses.AlignmentBlob blo : blobs) {
				
				for(int refParser = blo.getReferenceStart(), seqParser = blo.getReadStart() - 1; refParser <= ((blo.getReferenceStart() + blo.getLength()) - 1); refParser++, seqParser++) {
					if (tools.qualityString(String.valueOf(seqQualArray[seqParser]), qualityHash) > 10.0) {
						String newString = mutations.get(refParser) + seqSubArray[seqParser];
						newString = newString.toUpperCase();
						mutations.put(refParser, newString);
					}

				}
				
			}
					
		}
		
		// Generate a pileup of all the read information at this site, and decide consensus (all the while caring about ins/del status):
		String consensus = getPileup(mutations);
		
		// Check for deletions in the consensus:
		checkDel(consensus);
			
		// Check for insertions in the consensus:
		checkIns(consensus);
		
		// Now apply InDels to the consensus:
		// Have to apply the change to the consensus, and the subtract or add to each element in Array
		
		String meiCut = meiLoaded.getSequence();
		for (int x = 0; x < toCut.size(); x++) {
			CheckClasses.cutter current = toCut.get(x);
			if (current.getType().equals("insertion")) {
				String before = consensus.substring(0, current.getLocation());
				String after = consensus.substring(current.getLocation());
				int length = current.getSequence().length();
				consensus = before + current.getSequence() + after;
				before = meiCut.substring(0, current.getLocation());
				after = meiCut.substring(current.getLocation());
				length = current.getSequence().length();
				meiCut = before + current.getSequence() + after;
				// Now go through anything after and add length to every position
				for (int y = x+1; y < toCut.size(); y++) {
					CheckClasses.cutter fixed = toCut.get(y);
					fixed.setLocation(fixed.getLocation() + length);
				}
			} else {
				String before = consensus.substring(0, current.getLocation());
				String after = consensus.substring(current.getLocation() + current.getSequence().length() - 1, meiCut.length());
				int length = current.getSequence().length();
				consensus = before + after;
				before = meiCut.substring(0, current.getLocation());
				after = meiCut.substring(current.getLocation() + current.getSequence().length(), meiCut.length());
				length = current.getSequence().length();
				meiCut = before + after;
				// Now go through anything after and subtract length from every position
				for (int y = x+1; y < toCut.size(); y++) {
					CheckClasses.cutter fixed = toCut.get(y);
					fixed.setLocation(fixed.getLocation() - length);
				}
			}
		}
		
		return consensus;
		
	}
	private GenericResults getGenericSpecies(String consensus) throws IOException, InterruptedException, MatrixLoaderException {
				
		GenericU getType = new GenericU(consensus, meiLoaded);
		GenericResults results = getType.getPosResults();
		String diffs = results.getDifferences();
		
		Pattern dBoth = Pattern.compile("d(\\S+\\|)d(\\d+\\-\\d+)");
		Pattern dLeft = Pattern.compile("d(\\S+)");
		Pattern dRight = Pattern.compile("(\\S+\\|)d(\\d+\\-\\d+)");
		Pattern dRightAlt = Pattern.compile("(\\S+\\|)d(\\d+)");
		Pattern dBothAlt = Pattern.compile("d(\\S+\\|)d(\\d+)");
		Matcher matchBoth = dBoth.matcher(diffs);
		Matcher matchLeft = dLeft.matcher(diffs);
		Matcher matchRight = dRight.matcher(diffs);
		Matcher matchRightAlt = dRightAlt.matcher(diffs);
		Matcher matchBothAlt = dBothAlt.matcher(diffs);
		
		if (matchBoth.matches()) {
			results.setDifferences("n" + matchBoth.group(1) + "n" + matchBoth.group(2));
		} else if (matchBothAlt.matches()) {
			results.setDifferences("n" + matchBothAlt.group(1) + "n" + matchBothAlt.group(2));
		} else if (matchRightAlt.matches()) {
			results.setDifferences(matchRightAlt.group(1) + "n" + matchRightAlt.group(2));
		} else if (matchLeft.matches()) {
			results.setDifferences("n" + matchLeft.group(1));
		} else if (matchRight.matches()) {
			results.setDifferences(matchRight.group(1) + "n" + matchRight.group(2));
		}
		
		return results;
		
	}
	
	// private methods for getting consensus sequence for the species finder:
	private List<CheckClasses.AlignmentBlob> getReadBlocks(Cigar cigar, int readStart, String readName) {
		
		List<CheckClasses.AlignmentBlob> blobs = new ArrayList<CheckClasses.AlignmentBlob>();
				
		int currentRefPos = readStart;
		int currentReadPos = 1;
		
		for (CigarElement ele : cigar.getCigarElements()) {
			
			if (ele.getOperator().toString().equals("S")) {
				currentReadPos += ele.getLength();
				continue;
			} else if (ele.getOperator().toString().equals("M")) {
				CheckClasses.AlignmentBlob blob = holder.new AlignmentBlob();
				blob.setLength(ele.getLength());
				blob.setReadStart(currentReadPos);
				blob.setReferenceStart(currentRefPos);
				blobs.add(blob);
				currentReadPos+= ele.getLength();
				currentRefPos+= ele.getLength();
			} else if (ele.getOperator().toString().equals("D")) {
				if (deletions.containsKey(currentRefPos)) {
					CheckClasses.InDelInfo delInfo = deletions.get(currentRefPos);
					delInfo.setTotal(delInfo.getTotal() + 1);
					deletions.put(currentRefPos, delInfo);
				} else {
					CheckClasses.InDelInfo delInfo = holder.new InDelInfo();
					delInfo.setTotal(1);
					delInfo.setLength(ele.getLength());
					deletions.put(currentRefPos, delInfo);
				}
				currentRefPos+= ele.getLength();
				continue;
			} else if (ele.getOperator().toString().equals("I")) {
				if (insertions.containsKey(currentRefPos)) {
					CheckClasses.InDelInfo insInfo = insertions.get(currentRefPos);
					insInfo.setTotal(insInfo.getTotal() + 1);
					CheckClasses.InDelInfo.Info info = insInfo.new Info();
					info.setInsLen(ele.getLength());
					info.setPosition(currentReadPos);
					insInfo.addToInsPos(readName, info);
					insertions.put(currentRefPos, insInfo);
				} else {
					CheckClasses.InDelInfo insInfo = holder.new InDelInfo();
					insInfo.setTotal(1);
					CheckClasses.InDelInfo.Info info = insInfo.new Info();
					info.setInsLen(ele.getLength());
					info.setPosition(currentReadPos);
					insInfo.addToInsPos(readName, info);
					insertions.put(currentRefPos, insInfo);
				}
				currentReadPos+=ele.getLength();
				continue;
			} else {
				continue;
			}
		}
		
		return blobs;
		
	}
	private String getPileup (HashMap<Integer,String> mutations) {
		
		String consensus = "";

		for (int x = 1; x <= meiLoaded.length(); x++) {
			
			if (mutations.get(x).length() > 0 && x >= myLength.getInsertSizeStart() && x <= myLength.getInsertSizeStop()) {
				
				if (mutations.get(x).length() >= 2) {
					total2X++;
				}
				
				double aTotal = 0;
				double tTotal = 0;
				double cTotal = 0;
				double gTotal = 0;
				char bases[] =  mutations.get(x).toCharArray();
				double numBases = bases.length;
				
				for (char base : bases) {
						
					if (base == 'A') {
						aTotal++;
					} else if (base == 'T') {
						tTotal++;
					} else if (base == 'G') {
						gTotal++;
					} else {
						cTotal++;
					}
						
				}
				aTotal /= numBases;
				tTotal /= numBases;
				cTotal /= numBases;
				gTotal /= numBases;
					
				if (aTotal >= tTotal && aTotal >= cTotal && aTotal >= gTotal) {
					consensus += "A";
				} else if (tTotal >= aTotal && tTotal >= cTotal && tTotal >= gTotal) {
					consensus += "T";
				} else if (cTotal >= aTotal && cTotal >= tTotal && cTotal >= gTotal) {
					consensus += "C";
				} else {
					consensus += "G";
				}
				
			} else {
				
				//consensus += meiLoaded.getSequence().substring(x-1, x);
				consensus += "N";	
				
			}
			
		}
		
		return consensus;
		
	}
	private void checkDel(String consensus) {
		
		for (Map.Entry<Integer, CheckClasses.InDelInfo> entry : deletions.entrySet()) {
			
			double totalHits = entry.getValue().getTotal();
			int delLoc = entry.getKey();
			double totalReads = 0;
			
			for (SAMRecord thisOne : mySupInfo.getSpeciesCheck()) {
				
				if (thisOne.getIntegerAttribute("OS") < delLoc && thisOne.getIntegerAttribute("OE") > delLoc) {
					
					totalReads++;
					
				}
				
			}
			
			double calc = totalHits / totalReads * 100;
			
			if ((calc > 50) && totalReads > 2) {
				
				String del = consensus.substring(delLoc, delLoc + entry.getValue().getLength());
				CheckClasses.cutter cut = holder.new cutter();
				cut.setLocation(delLoc);
				cut.setSequence(del);
				cut.setType("deletion");
				// add items in order to the array
				if (toCut.isEmpty()) {
					toCut.add(cut);	
				} else {
					int size = toCut.size();
					if (size == 1) {
						if (toCut.get(0).getLocation() > delLoc) {
							toCut.add(0, cut);
						} else {
							toCut.add(cut);
						}
					} else {
						// Walk through each array element to see where the appropriate place is
						int cutSize = toCut.size();
						for (int x = 0; x < cutSize; x++) {
							if (toCut.get(x).getLocation() > delLoc) {
								toCut.add(x, cut);
								break;
							}
							if (x == toCut.size() - 1) {
								toCut.add(cut);
							}
						}
					}
				}	
			}
		}
	}
	private void checkIns(String consensus) {
		
		for (Map.Entry<Integer, CheckClasses.InDelInfo> entry : insertions.entrySet()) {
			
			double totalHits = entry.getValue().getTotal();
			int insLoc = entry.getKey();
			double totalReads = 0;
			Map<String, Integer> insConsensus = new HashMap<String, Integer>();
			
			for (SAMRecord thisOne : mySupInfo.getSpeciesCheck()) {
				
				if (thisOne.getIntegerAttribute("OS") < insLoc && thisOne.getIntegerAttribute("OE") > insLoc) {
					
					totalReads++;
					if (entry.getValue().getInsPos().containsKey(thisOne.getReadName())) {
						int insReadLoc = entry.getValue().getInsPos().get(thisOne.getReadName()).getPosition();
						int insLength = entry.getValue().getInsPos().get(thisOne.getReadName()).getInsLen();
						String insString = thisOne.getStringAttribute("R2").substring(insReadLoc, insReadLoc + insLength);
						if (insConsensus.containsKey(insString)) {
							int adjust = insConsensus.get(insString);
							adjust++;
							insConsensus.put(insString, adjust);
						} else {
							insConsensus.put(insString, 1);
						}
					}
					
				}
				
			}
			
			double calc = totalHits / totalReads * 100;
			
			if ((calc >= 50) && totalReads > 2) {
			
				//  Decide the proper insertion to place into the consensus:
				String highest = "";
				int currentScore = 0;
				
				for (Map.Entry<String, Integer> comp : insConsensus.entrySet()) {
					
					if (comp.getValue() > currentScore) {
						
						currentScore = comp.getValue();
						highest = comp.getKey();
						
					}
					
				}
				
				CheckClasses.cutter cut = holder.new cutter();
				cut.setLocation(insLoc);
				cut.setSequence(highest);
				cut.setType("insertion");
				if (toCut.isEmpty()) {
					toCut.add(cut);	
				} else {
					int size = toCut.size();
					if (size == 1) {
						if (toCut.get(0).getLocation() > insLoc) {
							toCut.add(0, cut);
						} else {
							toCut.add(cut);
						}
					} else {
						// Walk through each array element to see where the appropriate place is
						int cutSize = toCut.size();
						for (int x = 0; x < cutSize; x++) {
							
							if (toCut.get(x).getLocation() > insLoc) {
								toCut.add(x, cut);
								break;
							}
							if (x == (toCut.size() - 1)) {
								toCut.add(cut);
							}
						}
					}
				}	
			}
		}		
	}
	
	//Get Gene info specific to RNA-Transposon Hybrids:
	private ExonInfo checkGene() throws UnknownHostException, IOException {
		
		IntervalTree<ExonInfo> currentGeneInfo = annotation.get(exonization.getChr());
		ExonInfo hit = null;
		int searchPos;

		if (annote.getQual() == 5) {
			if (annote.getTsd().matches("d[ATGC]+") || annote.getTsd().matches("0tsd")) {
				searchPos = (int) annote.getLeftposition();
			} else {
				searchPos = (int) annote.getRightposition();
			}
		} else if (annote.getQual() == 3) {
			searchPos = (int) annote.getLeftposition();
		} else if (annote.getQual() == 2) {
			if (annote.getRightposition() == 0) {
				searchPos = (int) annote.getLeftposition();
			} else {
				searchPos = (int) annote.getRightposition();
			}
		} else {
			searchPos = (int) annote.getRightposition();
		}
		
		Iterator<Node<ExonInfo>> geneOverlap = currentGeneInfo.overlappers(searchPos - 5, searchPos + 5);
						
		while (geneOverlap.hasNext()) {
			
			ExonInfo currentHit = geneOverlap.next().getValue();
			if (hit == null) {
				hit = currentHit;
			} else if (hit.getType() == ExonType.EXON) {
				continue;
			} else if ((hit.getType() == ExonType.TUTR || hit.getType() == ExonType.FUTR) && currentHit.getType() == ExonType.EXON) {
				hit = currentHit;						
			} else if (((hit.getType() == ExonType.FAKE_TUTR || hit.getType() == ExonType.FAKE_FUTR) && (currentHit.getType() == ExonType.EXON || (currentHit.getType() == ExonType.FUTR || currentHit.getType() == ExonType.TUTR)))) {
				hit = currentHit;
			}
			
		}
		
		return hit;
		
	}
	
	// Private mode method
	public static Map<Integer,Integer> getMode (List<Integer> toCheck) {
		
		int maxCount=0;
		Map<Integer,Integer> modes = new HashMap<Integer,Integer>();
		Map<Integer,Integer> searched = new HashMap<Integer,Integer>();
		
		for (int check : toCheck) {
			
			//Want to make sure that due to sliding of Breakpoint, that SRs within 1/2 bp of each other get counted together
			
			if (searched.containsKey(check)) {
				
				int replace = searched.get(check);
				replace++;
				searched.put(check, replace);
				
			} else {
				
				searched.put(check, 1);
				
			}
			
		}
		
		// get modes that are within 1 bp of each other, and use the highest one.
		
		Map<Integer,Integer> decider = new HashMap<Integer,Integer>();
		List<Integer> sites = new ArrayList<Integer>(searched.keySet());
		Collections.sort(sites);
		int last = -9999;
		for (int merger : sites) {
			
			int min_one = merger - 1;
			
			if (min_one == last) {
				
				int test = searched.get(merger);
				int repl = searched.get(min_one);
				
				decider.put(merger, test);
				decider.put(min_one, repl);
				
				if (searched.containsKey(merger)) {
					searched.remove(merger);
				}
				if (searched.containsKey(min_one)) {
					searched.remove(min_one);
				}
				
			} else {
				
				int held = 0;
				int highest = 0;
				int total = 0;
				for(Map.Entry<Integer, Integer> entry : decider.entrySet()) {
					if (entry.getValue() > held) {
						highest = entry.getKey();
						held = entry.getValue();
					}
					total+=entry.getValue();
				}
				searched.put(highest, total);
				decider.clear();
			}
			
		}
		
		for (Map.Entry<Integer, Integer> entry : searched.entrySet()) {
			
			if (entry.getValue() > maxCount) {
				
				modes.clear();
				modes.put(entry.getKey(), entry.getValue());
				maxCount = entry.getValue();
				
			} else if (entry.getValue() == maxCount) {
				
				modes.put(entry.getKey(), entry.getValue());
				
			}
			
		}
		
		return modes;
		
	}

	//Private cigar builder:
	Cigar buildCigar(String cigarString) {
		List<CigarElement> cigarElements = new ArrayList<CigarElement>();
		Pattern cigarPatt = Pattern.compile("(\\d+)([MIDNSHPX=])(\\S*)");
		while (cigarString.length() > 0) {
			Matcher cigarMatch = cigarPatt.matcher(cigarString);
			if (cigarMatch.matches()) {
				cigarString = cigarMatch.group(3);
				CigarElement ele;
				if (cigarMatch.group(2).equals("=")) {
					ele = new CigarElement(Integer.parseInt(cigarMatch.group(1)), CigarOperator.EQ);
				} else {
					ele = new CigarElement(Integer.parseInt(cigarMatch.group(1)), CigarOperator.characterToEnum(cigarMatch.group(2).toCharArray()[0]));
				}
				cigarElements.add(ele);
			}
		}
		return new Cigar(cigarElements);
	}
	
	// Exception for length issue:
	public class DistributionException extends Exception {
		
		private static final long serialVersionUID = 1L;
		
		public DistributionException() {
			super();
		}
		
	}
	public class NoDiscordantPairs extends Exception {
		
		private static final long serialVersionUID = 1L;
		
		public NoDiscordantPairs() {
			super();
		}
		
	}

	private int determineStart (double left, double right, int qual, String tsd) {
		
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
