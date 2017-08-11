package mergehybrids;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import MELT.MELTIllumina.checkStats.CheckClasses;
import MELT.MELTIllumina.classification.generic.GenericResults;
import MELT.MELTIllumina.classification.generic.GenericU;
import MELT.utilities.SAMTools;
import findhybrids.HybridCheckMethods;
import findhybrids.SINERecordFactory;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jaligner.Sequence;
import jaligner.matrix.MatrixLoaderException;

public class MergeHybridsMethods {

	private List<SAMRecord> reads;
	private Sequence meiLoaded;

	// HashMap for calling InDels in Species Finder:
	private HashMap<Integer,CheckClasses.InDelInfo> deletions;
	private HashMap<Integer,CheckClasses.InDelInfo> insertions;
	private List<CheckClasses.cutter> toCut;

	private Pattern cigarCaptureRight = Pattern.compile("(\\S+[A-Z]{1})(\\d+)S$");
	private Pattern cigarCaptureLeft = Pattern.compile("^(\\d+)S(\\S+)");
	private ArrayList<Integer> starts;
	private ArrayList<Integer> stops;
	
	private SAMTools tools;
	private CheckClasses holder;
	
	private int insertStart;
	private int insertStop;
	private int totalLeft;
	private int totalRight;
	
	private int meiStart;
	private int meiStop;
	private String consensus;
	
	public MergeHybridsMethods(List<SAMRecord> reads, IndexedFastaSequenceFile mei, String species, File SINEDumpFolder, FinalHit hit) throws IOException, InterruptedException, MatrixLoaderException, NumberFormatException, URISyntaxException {
		
		this.reads = reads;
		String transSeq = new String(mei.getSequence(species).getBases()).toUpperCase();
		meiLoaded = new Sequence(transSeq);
		meiLoaded.setId(species);

		starts = new ArrayList<Integer>();
		stops = new ArrayList<Integer>();
		deletions = new HashMap<Integer,CheckClasses.InDelInfo>();
		insertions = new HashMap<Integer,CheckClasses.InDelInfo>();
		toCut = new ArrayList<CheckClasses.cutter>();
		tools = new SAMTools();

		holder = new CheckClasses();
		
		consensus = grabInformation();
		GenericResults result = getGenericSpecies(consensus);
		if (totalLeft > 0) {
			meiStart = result.getStart();
		} else {
			meiStart = -1;
		}
		if (totalRight > 0) {
			meiStop = result.getStop();
		} else {
			meiStop = -1;
		}
		
	}

	public String getConsensus() {
		return consensus;
	}
	public int getMEIStart() {
		return meiStart;
	}
	public int getMEIStop() {
		return meiStop;
	}
	
	public String grabInformation () throws IOException, NumberFormatException, InterruptedException, MatrixLoaderException, URISyntaxException {
		
		Hashtable<Character,Integer> qualityHash = tools.loadQualities(false);
		HashMap<Integer,String> mutations = new LinkedHashMap<Integer,String>();
		
		//Populate the map, so that I can add stuff to it
		for(int parser = 1; parser <= meiLoaded.length(); parser++) {
			mutations.put(parser, "");
		}
		
		for (SAMRecord record : reads) {
		
			//Add start/stop information:
			int meiStart = record.getStart();
			int meiStop = record.getEnd();
			String cigar = record.getCigarString();
			
			//Do first species information:
			char[] seqSubArray = record.getReadString().toCharArray();
			char[] seqQualArray = record.getBaseQualityString().toCharArray();
			
			List<CheckClasses.AlignmentBlob> blobs = getReadBlocks(buildCigar(cigar), meiStart, record.getReadName());			
						
			for (CheckClasses.AlignmentBlob blo : blobs) {
				
				for(int refParser = blo.getReferenceStart(), seqParser = blo.getReadStart() - 1; refParser <= ((blo.getReferenceStart() + blo.getLength()) - 1); refParser++, seqParser++) {
					if (tools.qualityString(String.valueOf(seqQualArray[seqParser]), qualityHash) > 10.0) {
						String newString = mutations.get(refParser) + seqSubArray[seqParser];
						newString = newString.toUpperCase();
						mutations.put(refParser, newString);
					}

				}
				
			}
			
			//Now build length info stuff:
			Matcher right = cigarCaptureRight.matcher(cigar);
			Matcher left = cigarCaptureLeft.matcher(cigar);
//			private Pattern cigarCaptureRight = Pattern.compile("(\\S+[A-Z]{1})(\\d+)S$");
//			private Pattern cigarCaptureLeft = Pattern.compile("^(\\d+)S(\\S+)");
			//2 right 1 left
			if (right.matches() && left.matches()) {
				if (Integer.parseInt(right.group(2)) >= 5 && Integer.parseInt(left.group(1)) >= 5) {
					starts.add(meiStart);
					stops.add(meiStop);
					totalLeft++;
					totalRight++;
				}
			} else if (right.matches() && !left.matches()) {
				if (Integer.parseInt(right.group(2)) >= 5) {
					stops.add(meiStop);
					totalRight++;
				}
			} else if (!right.matches() && left.matches()) {
				if (Integer.parseInt(left.group(1)) >= 5) {
					starts.add(meiStart);
					totalLeft++;
				}
			}
			if (meiStop == meiLoaded.length()) {
				stops.add(meiStop);
				totalRight++;
			}
			if (meiStart == 1) {
				starts.add(meiStart);
				totalLeft++;
			}

		}
		
		//get length
		getLength();
		//get the consensus sequence
		return getSpecies(mutations).replaceAll("N", "");
		
	}
	
	private void getLength() {
						
		if (starts.size() > 0) {
			
			Collections.sort(starts);
			Map<Integer, Integer> startMap = HybridCheckMethods.getMode(starts);
			List<Integer> startModes = new ArrayList<Integer>(); 
			startModes.addAll(startMap.keySet());
					
			if (startModes.size() >= 1) {
				insertStart = startModes.get(0);
			} else {
				insertStart = 0;
			}
			
		} else {
			insertStart = 0;
		}
			
		if (stops.size() > 0) {
			
			Collections.sort(stops);
			Map<Integer, Integer> stopMap = HybridCheckMethods.getMode(stops);
			List<Integer> stopModes = new ArrayList<Integer>(); 
			stopModes.addAll(stopMap.keySet());
			
			if (stopModes.size() >= 1) {
				insertStop = stopModes.get(0);
			} else {
				insertStop = meiLoaded.length() - 1;
			}
		} else {
			insertStop = meiLoaded.length() - 1;
		}
		
	}	
	
	private String getSpecies(HashMap<Integer,String> mutations) throws IOException, InterruptedException, MatrixLoaderException, URISyntaxException, NumberFormatException {
		
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
			
			if (mutations.get(x).length() > 0 && x >= insertStart && x <= insertStop) {
				
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
			
			for (SAMRecord thisOne : reads) {
				
				if (thisOne.getStart() < delLoc && thisOne.getEnd() > delLoc) {
					
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
			
			for (SAMRecord thisOne : reads) {
				
				if (thisOne.getStart() < insLoc && thisOne.getEnd() > insLoc) {
					
					totalReads++;
					if (entry.getValue().getInsPos().containsKey(thisOne.getReadName())) {
						int insReadLoc = entry.getValue().getInsPos().get(thisOne.getReadName()).getPosition();
						int insLength = entry.getValue().getInsPos().get(thisOne.getReadName()).getInsLen();
						String insString = thisOne.getReadString().substring(insReadLoc, insReadLoc + insLength);
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
	
}
