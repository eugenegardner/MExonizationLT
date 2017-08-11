package mergehybrids;

import java.util.Iterator;
import java.util.Map;

import org.biojava3.core.sequence.DNASequence;

import MELT.MELTIllumina.checkStats.CheckAnnotation.Strand;
import findhybrids.ExonInfo;
import findhybrids.ExonInfo.ExonType;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import jaligner.ui.filechooser.NamedInputStream;

public class CheckAssembly {

	private static final float GAP_PENALTY = 50;
	private static final float EXTEND_PENALTY = 5;
	
	private IndexedFastaSequenceFile humRef;
	private Map<String, IntervalTree<ExonInfo>> annotation;
	private Map<String, IntervalTree<String >>sineMap;
	private Map<String, IntervalTree<String>> transRef;
	private Matrix matrix;
	
	public CheckAssembly (IndexedFastaSequenceFile humRef, Map<String, IntervalTree<ExonInfo>> annotation, Map<String, IntervalTree<String>> sineMap, Map<String, IntervalTree<String>> transRef) throws MatrixLoaderException {
		
		this.humRef = humRef;
		this.annotation = annotation;
		this.sineMap = sineMap;
		this.transRef = transRef;
		
		matrix = MatrixLoader.load(new NamedInputStream("NAMatrix",this.getClass().getResourceAsStream("/importFiles/nuc-4_4.txt")));
		
	}

	public Map<String, IntervalTree<String>> getSineMap () {
		return sineMap;
	}
	
	public MatchedMEI checkAlignments (String chr, int position, ExonInfo geneInfo, String assembly, Strand strand) {
		
		int lastExon = 0;
		int nextExon = 0;
		ExonInfo recursiveInfo = null;

		boolean caught = false;
		if (strand == Strand.POS) {
			lastExon = position;
			Iterator<Node<ExonInfo>> forwardItr = annotation.get(chr).iterator(position, position);
			while (forwardItr.hasNext() && caught == false) {
				Node<ExonInfo> test = forwardItr.next();
				if (test.getValue().getGeneName().equals(geneInfo.getGeneName()) && test.getValue().getExonNum() != geneInfo.getExonNum()) {
					caught = true;
					nextExon = test.getStart();
					recursiveInfo = test.getValue();
				}
			}
			
		} else {
			Iterator<Node<ExonInfo>> reverseItr = annotation.get(chr).reverseIterator(position, position);
			while (reverseItr.hasNext() && caught == false) {
				Node<ExonInfo> test = reverseItr.next();
				if (test.getValue().getGeneName().equals(geneInfo.getGeneName()) && test.getValue().getExonNum() != geneInfo.getExonNum()) {
					caught = true;
					lastExon = test.getEnd();
					recursiveInfo = test.getValue();
				}				
			}
			nextExon = position;
		}
		
		if (caught == false) {
			if (strand == Strand.POS) {
				nextExon = position + 30000;
			} else {
				lastExon = position - 30000;
			}
		}
		
		IntervalTree<String> currentNonrefIntervalTree = sineMap.get(chr);
		IntervalTree<String> currentRefIntervalTree = transRef.get(chr);		
		
		Iterator<Node<String>> nonRef = currentNonrefIntervalTree.overlappers(lastExon, nextExon);
		Iterator<Node<String>> ref = currentRefIntervalTree.overlappers(lastExon, nextExon);
		MatchedMEI maxNonRef = checkForMatches(nonRef, assembly, chr, MEIType.NONREF);
		MatchedMEI maxRef = checkForMatches(ref, assembly, chr, MEIType.REF);
		
		//If null, try to recover with recursion -- will go until the end of the gene, at which point will throw a NoNextExon error!
		double maxHit;
		if (maxNonRef == null || maxRef == null) {
			maxHit = 0;
		} else if (maxNonRef.getPercIdentity() >= 98.5 && maxRef.getPercIdentity() >= 95) {
			maxHit = Math.max(maxNonRef.getPercIdentity(), maxRef.getPercIdentity());
		} else if (maxNonRef.getPercIdentity() >= 98.5 && maxRef.getPercIdentity() < 95) {
			maxHit = maxNonRef.getPercIdentity();
			maxRef = null;
		} else if (maxNonRef.getPercIdentity() < 98.5 && maxRef.getPercIdentity() >= 95) {
			maxHit = maxRef.getPercIdentity();
			maxNonRef = null;
		} else {
			maxHit = 0;
		}
				
		if (maxNonRef == null && maxRef == null) {
			if (caught == true) {
				if (strand == Strand.POS) {
					return checkAlignments(chr, nextExon, recursiveInfo, assembly, strand);
				} else {
					return checkAlignments(chr, lastExon, recursiveInfo, assembly, strand);
				}
			} else {
				return null;
			}
		} else if (maxNonRef == null && maxRef != null) {
			if (maxRef.getPercIdentity() >= 95) { //Don't want bad hits.
				return checkGene(maxRef, caught);
			} else {
				if (caught == true) {
					if (strand == Strand.POS) {
						return checkAlignments(chr, nextExon, recursiveInfo, assembly, strand);
					} else {
						return checkAlignments(chr, lastExon, recursiveInfo, assembly, strand);
					}
				} else {
					return null;
				}
			}
		} else if (maxNonRef != null && maxRef == null) {
			if (maxNonRef.getPercIdentity() >= 98.5) { //Don't want bad hits.
				return checkGene(maxNonRef, caught);
			} else {
				if (caught == true) {
					if (strand == Strand.POS) {
						return checkAlignments(chr, nextExon, recursiveInfo, assembly, strand);
					} else {
						return checkAlignments(chr, lastExon, recursiveInfo, assembly, strand);
					}
				} else {
					return null;
				}
			}
		} else { //This should be if both are non null
			
			if (maxHit >= 95) { //Don't want crap
				if (maxNonRef.getPercIdentity() > maxRef.getPercIdentity()) {
					return checkGene(maxNonRef, caught);
				} else if (maxNonRef.getPercIdentity() < maxRef.getPercIdentity()) {
					return checkGene(maxRef, caught);
				} else { //This should default to choosing nonRef site if it is higher!
					return checkGene(maxNonRef, caught);
				}
			} else {
				if (caught == true) {
					if (strand == Strand.POS) {
						return checkAlignments(chr, nextExon, recursiveInfo, assembly, strand);
					} else {
						return checkAlignments(chr, lastExon, recursiveInfo, assembly, strand);
					}
				} else {
					return null;
				}
			}
		}
		
	}
		
	private MatchedMEI checkForMatches (Iterator<Node<String>> iterator, String assembly, String chr, MEIType type) {
		
		float maxHit = 0;
		MatchedMEI currentMatch = null;
		
		while (iterator.hasNext()) {
			
			Node<String> currentMEI = iterator.next();
			String reference;
			if (type == MEIType.NONREF) {
				reference = currentMEI.getValue();
			} else {
				if (currentMEI.getValue().matches("\\+")) {
					reference = humRef.getSubsequenceAt(chr, currentMEI.getStart(), currentMEI.getEnd()).getBaseString();
				} else {
					reference = new DNASequence(reference = humRef.getSubsequenceAt(chr, currentMEI.getStart(), currentMEI.getEnd()).getBaseString()).getReverseComplement().getSequenceAsString();
				}
			}
			reference = reference.toUpperCase();
			Alignment result = alignMEIs(assembly, reference);
			double identity = calculateIdentity(result.getMarkupLine(), assembly.length());
//			System.out.println(chr + "\t" + currentMEI.getStart() + "\t" + currentMEI.getEnd() + "\t" + identity + "\t" + result.getScore());
			if (result.getScore() >= 190 && result.getScore() > maxHit) {
				currentMatch = new MatchedMEI(chr, currentMEI.getStart(), currentMEI.getEnd(), identity, result.getScore(), type);
				maxHit = result.getScore();
			}
			
		}
		
		return currentMatch;
		
	}
	
	private Alignment alignMEIs (String query, String reference) {
		
		Sequence querySeq = new Sequence(query);
		Sequence refSeq = new Sequence(reference);
		Alignment forwardAlignment = SmithWatermanGotoh.align(querySeq, refSeq, matrix, GAP_PENALTY, EXTEND_PENALTY);
//		System.out.println(query);
//		System.out.print(String.valueOf(forwardAlignment.getSequence1()) + "\n" + String.valueOf(forwardAlignment.getMarkupLine()) + "\n" + String.valueOf(forwardAlignment.getSequence2()) + "\n");
//		System.out.println(reference);
		return forwardAlignment;
		
	}
	private double calculateIdentity(char[] markup, double expectedLength) {
		
		double alignmentLength = markup.length;
		double factor;
		if (alignmentLength > expectedLength) {
			factor = expectedLength / alignmentLength;
		} else {
			factor = alignmentLength / expectedLength;
		}
		double identity = (getSimilarity(markup) * factor) * 100;
		return identity;
		
	}
	private double getSimilarity (char[] markup) {
		
		double matches = 0;
		double length = markup.length;
		
		for (char c : markup) {
			if (c == '|') {
				matches++;
			}
		}
		
		return matches / length;
		
	}
	
	private MatchedMEI checkGene (MatchedMEI matched, boolean caught) {
		
		IntervalTree<ExonInfo> currentGeneInfo = annotation.get(matched.getChr());
		ExonInfo hit = null;
		
		Iterator<Node<ExonInfo>> geneOverlap = currentGeneInfo.overlappers(matched.getStart(), matched.getStop());
		
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
		
		if (hit != null) {
			if (hit.getType() == ExonType.EXON || hit.getType() == ExonType.FUTR || hit.getType() == ExonType.TUTR) {
				matched.setMEIType(MEIType.GENE);
			}
		}
		if (caught == false) {
			if (matched.getType() == MEIType.REF) {
				matched.setMEIType(MEIType.REFUTR);
			} else if (matched.getType() == MEIType.NONREF) {
				matched.setMEIType(MEIType.NONREFUTR);
			}
		}
		
		return matched;
		
	}
	
	public class MatchedMEI {
		
		private String chr;
		private int start;
		private int stop;
		private double percIdentity;
		private float alignmentScore;
		private MEIType type;
		
		public MatchedMEI(String chr, int start, int stop, double percIdentity, float alignmentScore, MEIType type) {
			
			this.chr = chr;
			this.start = start;
			this.stop = stop;
			this.percIdentity = percIdentity;
			this.alignmentScore = alignmentScore;
			this.type = type;
			
		}

		public void setMEIType(MEIType type) {
			this.type = type;
		}
 		public String getChr() {
			return chr;
		}
		public int getStart() {
			return start;
		}
		public int getStop() {
			return stop;
		}
		public double getPercIdentity() {
			return percIdentity;
		}
		public float getAlignmentScore() {
			return alignmentScore;
		}
		public MEIType getType() {
			return type;
		}
		
	}
	
	public enum MEIType {
		REF,NONREF,GENE, REFUTR, NONREFUTR;
	}

}
