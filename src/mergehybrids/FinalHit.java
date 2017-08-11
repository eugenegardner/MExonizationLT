package mergehybrids;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import MELT.MELTIllumina.checkStats.CheckAnnotation.Strand;

public class FinalHit {

	private String chr;
	private int start;
	private int stop;
	private String gene;
	private List<Integer> exon;
	private List<Sample> samples;
	private int totalPos;
	private int totalNeg;
	
	public FinalHit (String chr, int start, int stop, String gene, int exonNum, String sample, int position, int totalReads, File parentDir, Strand strand) {
		
		this.chr = chr;
		this.start = start;
		this.stop = stop;
		this.gene = gene;
		this.exon = new ArrayList<Integer>();
		this.exon.add(exonNum);
		this.samples = new ArrayList<Sample>();
		samples.add(new Sample(sample, chr, position, totalReads, parentDir));
		if (strand == Strand.POS) {
			totalPos = 1;
			totalNeg = 0;
		} else {
			totalPos = 0;
			totalNeg = 1;
		}
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
	public String getGene() {
		return gene;
	}
	public List<Integer> getExon() {
		return exon;
	}
	public List<Sample> getSamples() {
		return samples;
	}
	public Strand getStrand() {
		if (totalPos > totalNeg) {
			return Strand.POS;
		} else if (totalPos < totalNeg) {
			return Strand.NEG;
		} else {
			return null;
		}
	}
	
	public void setStart (int start) {
		if (start < this.start) {
			this.start = start;
		}
	}
	public void setStop (int stop) {
		if (stop > this.stop) {
			this.stop = stop;
		}
	}
	public void addSample(String sample, int position, int totalReads, File parentDir) {
		samples.add(new Sample(sample, chr, position, totalReads, parentDir));
	}
	public void addSample(Sample sample) {
		samples.add(sample);
	}
	
	public void addExon(int exonNum) {
	
		if (!exon.contains(exonNum)) {
			exon.add(exonNum);
		}
		
	}
	public void addStrand(Strand strand) {
		if (strand == Strand.POS) {
			totalPos++;
		} else {
			totalNeg++;
		}
	}
	
 	public static class Sample {
		
		private String sample;
		private int position;
		private int totalReads;
		private File SINEreads;
		private File parentDir;
		
		public Sample (String sample, String chr, int position, int totalReads, File parentDir) {
			
			this.sample = sample;
			this.position = position;
			SINEreads = new File(parentDir.getAbsolutePath() + "/hybrids.sorted.bam.SINEAssembly/" + chr + "_" + position + ".sorted.bam");
			this.parentDir = parentDir;
			this.totalReads = totalReads;
			
		}

		public String getSample() {
			return sample;
		}

		public int getPosition() {
			return position;
		}

		public File getSINEreads() {
			return SINEreads;
		}
		
		public int getTotalReads() {
			return totalReads;
		}
		
		public File getParentDir() {
			return parentDir;
		}
				
	}

}
