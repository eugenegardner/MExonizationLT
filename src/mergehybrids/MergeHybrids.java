package mergehybrids;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import MELT.MELTIllumina.checkStats.CheckAnnotation.Strand;
import MELT.utilities.Combine;
import findhybrids.ExonInfo;
import findhybrids.FindHybrids;
import findhybrids.IGVPrinter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import jaligner.matrix.MatrixLoaderException;
import mergehybrids.CheckAssembly.MEIType;
import mergehybrids.CheckAssembly.MatchedMEI;
import mergehybrids.FinalHit.Sample;

public class MergeHybrids {

	private static List<String> allSamples;
	
	public static void main(String[] args) throws IOException, MatrixLoaderException, NumberFormatException, InterruptedException, URISyntaxException {

		File putativeList = new File(args[0]);
		Map<String, IntervalTree<String>> sineMap = loadSINES(new File (args[1]));
		Map<String, IntervalTree<String>> transRef = buildTransRef(new BufferedReader(new FileReader(new File(args[2]))));
		File fasta = new File(args[3]);
		File fastaIndex = new File(args[3] + ".fai");
		File meiFasta = new File(args[4]);
		File meiFastaIndex = new File(args[4] + ".fai");
		Map<String, IntervalTree<ExonInfo>> annotation = FindHybrids.buildAnnotation(new BufferedReader(new FileReader(new File(args[5]))));
		List<Integer> sites = loadSampleHash(args[6]);
		File assembliesDir = new File("SINEAssemblies/");
		if (!assembliesDir.isDirectory()) {
			assembliesDir.mkdir();
		}
		
		IndexedFastaSequenceFile dogIndexedFasta = new IndexedFastaSequenceFile(fasta, new FastaSequenceIndex(fastaIndex));
		IndexedFastaSequenceFile meiIndexedFasta = new IndexedFastaSequenceFile(meiFasta, new FastaSequenceIndex(meiFastaIndex));
		
		CheckAssembly checkAssembly = new CheckAssembly(dogIndexedFasta, annotation, sineMap, transRef);
		
		Map<String, FinalHit> merged = new HashMap<String, FinalHit>();
		allSamples = new ArrayList<String>();
		
		
		//Go through and collect everything first:
		BufferedReader listReader = new BufferedReader(new FileReader(putativeList));
		String line;
		
		while ((line = listReader.readLine()) != null) {
			
			File putativeFile = new File(line);
			merged = assessHits(putativeFile, merged);
			
		}
		
		//Try to connect to IGV
		IGVPrinter igvPrinterDog;
		IGVPrinter igvPrinterSINE;
		try {
			igvPrinterDog = new IGVPrinter(60151, dogIndexedFasta);
			igvPrinterSINE = new IGVPrinter(60152, dogIndexedFasta);
			System.out.println("Connected to IGV.");
		} catch (Exception e) {
			igvPrinterDog = null;
			igvPrinterSINE = null;
		}
		
		listReader.close();
		
		//Now merge adjacent hits:
		List<FinalHit> finalHits = mergeHits(merged.values());
		
		//Build the header for the final file:
		Collections.sort(allSamples);
		System.out.print("chr\tstart\tstop\tgene\texon\tmeiBreakLeft\tmeiBreakRight\tmeiStart\tmeiStop\tidentity\talignmentScore\tmeiType\t");
		String joinedHeader = Combine.combineList(allSamples, "\t");
		System.out.println(joinedHeader);
		
		SamReaderFactory samFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SAMFileWriterFactory SINEWriterFactory = new SAMFileWriterFactory().setCreateIndex(true);
		
		//Go through each hit and find the SINE element it closest matches to.
		for (FinalHit hit : finalHits) {
			if (sites.contains(hit.getStart())) {
			Map<String, Integer> include = setSampleTree();
			//Check samples this hit is in:
			for (Sample sample : hit.getSamples()) {
				String sampleName = sample.getSample();
				int total = include.get(sampleName);
				total += sample.getTotalReads();
				include.put(sampleName, total);
			}
			List<String> sampleList = new ArrayList<String>();
			for (Map.Entry<String, Integer> entry : include.entrySet()) {
				sampleList.add(entry.getValue().toString());
			}
			String samples = Combine.combineList(sampleList, "\t");
			
			//Make a list of all reads contained in this 'hit'
			List<SAMRecord> reads = new ArrayList<SAMRecord>();
			SAMFileHeader useHeader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(hit.getSamples().get(0).getSINEreads()).getFileHeader();
			SAMFileWriter SINEwriter = SINEWriterFactory.makeBAMWriter(useHeader, false, new File(assembliesDir + "/" + hit.getChr() + "_" + hit.getStart() + ".sorted.bam"));
			for (Sample sample : hit.getSamples()) {
				SAMRecordIterator itr = samFactory.open(SamInputResource.of(sample.getSINEreads())).iterator();
				while (itr.hasNext()) {
					SAMRecord current = itr.next();
					reads.add(current);
					SINEwriter.addAlignment(current);
				}
			}
			SINEwriter.close();
			
			MergeHybridsMethods methods = new MergeHybridsMethods(reads, meiIndexedFasta, "SINEC", assembliesDir, hit);
			
			MatchedMEI highestMatch;
			ExonInfo currentInfo;
			if (hit.getStrand() == Strand.POS) {
				if (hit.getExon().size() > 1) {
					currentInfo = new ExonInfo(hit.getGene(), false, null, hit.getExon().get(0));
					highestMatch = checkAssembly.checkAlignments(hit.getChr(), hit.getStart(), currentInfo, methods.getConsensus(), hit.getStrand());
				} else {
					currentInfo = new ExonInfo(hit.getGene(), false, null, hit.getExon().get(0));
					highestMatch = checkAssembly.checkAlignments(hit.getChr(), hit.getStop(), currentInfo, methods.getConsensus(), hit.getStrand());
				}
			} else {
				if (hit.getExon().size() > 1) {
					currentInfo = new ExonInfo(hit.getGene(), false, null, hit.getExon().get(hit.getExon().size() - 1));
					highestMatch = checkAssembly.checkAlignments(hit.getChr(), hit.getStop(), currentInfo, methods.getConsensus(), hit.getStrand());
				} else {
					currentInfo = new ExonInfo(hit.getGene(), false, null, hit.getExon().get(0));
					highestMatch = checkAssembly.checkAlignments(hit.getChr(), hit.getStop(), currentInfo, methods.getConsensus(), hit.getStrand());
				}
			}
			System.out.println(methods.getConsensus());
			if (highestMatch == null) {
				System.out.println(hit.getChr() + "\t" + hit.getStart() + "\t" + hit.getStop() + "\t" + hit.getGene() + "\t" + hit.getExon() + "\t" + methods.getMEIStart() + "\t" + methods.getMEIStop() + "\tnull\tnull\tnull\tnull\tnull\t" + samples);

				if (igvPrinterDog != null && igvPrinterSINE != null) {
					igvPrinterSINE.viewConsensus(methods.getConsensus(), new File("/local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/analysis/rna-seq/SINEAssemblies/"), hit.getChr(), hit.getStart());
					igvPrinterDog.viewSite(hit, include, null);
				}
			} else {
				if (highestMatch.getType() != MEIType.GENE) {
					if (highestMatch.getType() == MEIType.NONREFUTR && highestMatch.getPercIdentity() >= 100 && highestMatch.getAlignmentScore() > 300) {
						System.out.println(hit.getChr() + "\t" + hit.getStart() + "\t" + hit.getStop() + "\t" + hit.getGene() + "\t" + hit.getExon() + "\t" + methods.getMEIStart() + "\t" + methods.getMEIStop() + "\t" + highestMatch.getStart() + "\t" + highestMatch.getStop() + "\t" + highestMatch.getPercIdentity() + "\t" + highestMatch.getAlignmentScore() + "\t" + highestMatch.getType() + "\t" + samples);
					} else if (highestMatch.getType() == MEIType.REFUTR && highestMatch.getPercIdentity() >= 98.5 && highestMatch.getAlignmentScore() > 250) {
						System.out.println(hit.getChr() + "\t" + hit.getStart() + "\t" + hit.getStop() + "\t" + hit.getGene() + "\t" + hit.getExon() + "\t" + methods.getMEIStart() + "\t" + methods.getMEIStop() + "\t" + highestMatch.getStart() + "\t" + highestMatch.getStop() + "\t" + highestMatch.getPercIdentity() + "\t" + highestMatch.getAlignmentScore() + "\t" + highestMatch.getType() + "\t" + samples);
					} else if (highestMatch.getType() == MEIType.REF || highestMatch.getType() == MEIType.NONREF) {
						System.out.println(hit.getChr() + "\t" + hit.getStart() + "\t" + hit.getStop() + "\t" + hit.getGene() + "\t" + hit.getExon() + "\t" + methods.getMEIStart() + "\t" + methods.getMEIStop() + "\t" + highestMatch.getStart() + "\t" + highestMatch.getStop() + "\t" + highestMatch.getPercIdentity() + "\t" + highestMatch.getAlignmentScore() + "\t" + highestMatch.getType() + "\t" + samples);
					} else {
						System.out.println(hit.getChr() + "\t" + hit.getStart() + "\t" + hit.getStop() + "\t" + hit.getGene() + "\t" + hit.getExon() + "\t" + methods.getMEIStart() + "\t" + methods.getMEIStop() + "\tnull\tnull\tnull\tnull\tnull\t" + samples);
					}
					
					if (igvPrinterDog != null && igvPrinterSINE != null) {
						igvPrinterSINE.viewConsensus(methods.getConsensus(), new File("/local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/analysis/rna-seq/SINEAssemblies/"), hit.getChr(), hit.getStart());
						if (highestMatch != null) {
							if (highestMatch.getType() == MEIType.REF || highestMatch.getType() == MEIType.REFUTR) {
								igvPrinterSINE.viewRefSine(highestMatch.getChr(), highestMatch.getStart(), highestMatch.getStop(), new File("/local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/analysis/rna-seq/SINEAssemblies/"));
							} else if (highestMatch.getType() == MEIType.NONREF || highestMatch.getType() == MEIType.NONREFUTR) {
								igvPrinterSINE.viewSINE(hit.getChr(), highestMatch.getStart(), sineMap, new File("/local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/analysis/rna-seq/SINEAssemblies/"));
							}
						}
						igvPrinterDog.viewSite(hit, include, highestMatch);
					}
				}
			}
			}
		}
	}
	
	private static Map<String, FinalHit> assessHits (File putativeFile, Map<String, FinalHit> merged) throws IOException {
		
		BufferedReader putReader = new BufferedReader(new FileReader(putativeFile));
		String line;
		String data[];
		Matcher nameMatch = Pattern.compile("\\S*(SRR\\d+)\\/putative_hits\\.txt").matcher(putativeFile.getAbsolutePath());
		String sampleName = null;
		if (nameMatch.matches()) {
			sampleName = nameMatch.group(1);
			if (!allSamples.contains(sampleName)) {
				allSamples.add(sampleName);
			}
		}
		File sampleDir = putativeFile.getParentFile();
		
		while ((line = putReader.readLine()) != null) {
			
			data = line.split("\\t");
			
			String chr = data[0];
			int start = Integer.parseInt(data[1]);
			int stop = Integer.parseInt(data[2]);
			String gene = data[11];
			int exonNum = Integer.parseInt(data[13]);
			Strand strand = Strand.valueOf(data[14]);
			int totalReads = (int) Double.parseDouble(data[3]);
			
			String key = gene + "_" + exonNum;
			
			if (merged.containsKey(key)) {
				merged.put(key, assessHit(start, stop, exonNum, merged.get(key), sampleName, sampleDir, strand, totalReads));
			} else {
				merged.put(key, createFinalHit(chr, start, stop, gene, exonNum, sampleName, sampleDir, strand, totalReads));
			}
									
		}
		
		putReader.close();
		return merged;
		
	}
	
	public static FinalHit createFinalHit (String chr, int start, int stop, String gene, int exonNum, String sampleName, File parentDir, Strand strand, int totalReads) throws IOException {
		
		FinalHit hit = new FinalHit(chr, start, stop, gene, exonNum, sampleName, start, totalReads, parentDir, strand);
		return hit;
		
	}
	public static FinalHit assessHit (int start, int stop, int exonNum, FinalHit toAssess, String sampleName, File parentDir, Strand strand, int totalReads) {
		
		toAssess.setStart(start);
		toAssess.setStop(stop);
		toAssess.addExon(exonNum);
		toAssess.addSample(sampleName, start, totalReads, parentDir);
		toAssess.addStrand(strand);
		return toAssess;
		
	}
	public static FinalHit assessHit (int start, int stop, int exonNum, FinalHit toAssess, Sample sample, Strand strand) {
		
		toAssess.setStart(start);
		toAssess.setStop(stop);
		toAssess.addExon(exonNum);
		toAssess.addSample(sample);
		toAssess.addStrand(strand);
		return toAssess;
		
	}
	
	private static List<FinalHit> mergeHits (Collection<FinalHit> mergedValues) {
		
		List<FinalHit> finalHits = new ArrayList<FinalHit>();
		List<FinalHit> sortedHits = new ArrayList<FinalHit>();
		sortedHits.addAll(mergedValues);
		Collections.sort(sortedHits, new MergedRecordComparator());
				
		for (int x = 0; x < sortedHits.size(); x++) {
		
			FinalHit currentHit = sortedHits.get(x);
			String currentGene = currentHit.getGene();
			int currentExon = currentHit.getExon().get(0);
			
			boolean keepIterating = true;
			if ((x + 1) < sortedHits.size()) {
				for (int y = x + 1; keepIterating; y++) {
					
					FinalHit checkHit = sortedHits.get(y);
					String checkGene = checkHit.getGene();
					int checkExon = checkHit.getExon().get(0);
					
					if (currentGene.equals(checkGene) && checkExon == (currentExon + 1)) {
						
						for (Sample sample : checkHit.getSamples()) {
							currentHit = assessHit(checkHit.getStart(), checkHit.getStart(), checkExon, currentHit, sample, checkHit.getStrand());
						}
						currentExon++;
						x++;
					} else if (currentGene.equals(checkGene) && checkExon == (currentExon - 1)) {
						
						for (Sample sample : checkHit.getSamples()) {
							currentHit = assessHit(checkHit.getStart(), checkHit.getStart(), checkExon, currentHit, sample, checkHit.getStrand());
						}
						currentExon--;
						x++;
					} else {
						keepIterating = false;
					}
					
				}
			}
			
			finalHits.add(currentHit);			
			
		}
		
		return finalHits;
		
	}
		
	private static Map<String, IntervalTree<String>> loadSINES (File SINEs) throws IOException {
		
		Map<String, IntervalTree<String>> sineMap = new HashMap<String, IntervalTree<String>>();
		
		BufferedReader SINEreader = new BufferedReader(new FileReader(SINEs));
		String line;
		String data[];
		
		while ((line = SINEreader.readLine()) != null) {
			
			data = line.split("\\t");
			String chr = data[0];
			int start = Integer.parseInt(data[1]);
			int stop = Integer.parseInt(data[1]);
			
			if (sineMap.containsKey(chr)) {
				IntervalTree<String> positionInfo = sineMap.get(chr);
				positionInfo.put(start, stop, data[2]);
				sineMap.put(chr, positionInfo);
			} else {
				IntervalTree<String> positionInfo = new IntervalTree<String>();
				positionInfo.put(start, stop, data[2]);
				sineMap.put(chr, positionInfo);
			}
			
		}
		
		SINEreader.close();
		
		return sineMap;
		
	}
	private static Map<String, IntervalTree<String>> buildTransRef (BufferedReader insMasker) throws IOException {
		
		Map<String,IntervalTree<String>> insTree = new HashMap<String, IntervalTree<String>>();
		String myLine;
		String data[];
		
		while((myLine = insMasker.readLine()) != null) {
			data = myLine.split("\t");
			if (insTree.containsKey(data[0])) {
				insTree.get(data[0]).put(Integer.parseInt(data[1]), Integer.parseInt(data[2]), data[5]);
			} else {
				IntervalTree<String> inter = new IntervalTree<String>();
				insTree.put(data[0], inter);
				insTree.get(data[0]).put(Integer.parseInt(data[1]), Integer.parseInt(data[2]), data[5]);
			}
			
		}
		
		return insTree;
			
	}
	public static Strand getDirection(int totalLeft, int totalRight) {
		
		Strand strand;
		//totalLeft is # of reads FACING left;
		//totalRight is # of reads FACING right;
		
		if (totalLeft > totalRight) {
			strand = Strand.POS;
		} else {
			strand = Strand.NEG;
		}
		
		return strand;
		
	}
	
	private static Map<String, Integer> setSampleTree () {
		
		Map<String, Integer> sampleHash = new LinkedHashMap<String, Integer>();
		for (String sample : allSamples) {
			
			sampleHash.put(sample, 0);
			
		}
		
		return sampleHash;
		
	}

	private static List<Integer> loadSampleHash(String file) throws IOException {
		
		List<Integer> sites = new ArrayList<Integer>();
		BufferedReader siteReader = new BufferedReader(new FileReader(new File(file)));
		String line;
		while ((line = siteReader.readLine()) != null) {
			int site = Integer.parseInt(line);
			sites.add(site);
		}
		siteReader.close();
		return sites;
		
	}
	
}
