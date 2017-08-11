package findhybrids;

import java.util.List;

import htsjdk.samtools.SAMRecord;

public class PutativeHit {
		
	String chr;
	int start;
	int stop;
	double totalReads;
	double percDisc;
	ExonInfo exonInfo;
	List<SAMRecord> supportingReads;
		
	public PutativeHit(String chr, int start, int stop, double totalReads, double discReads, List<SAMRecord> supportingReads) {
		
		this.chr = chr;
		this.start = start;
		this.stop = stop;
		this.totalReads = totalReads;
		percDisc = discReads/totalReads;
		this.supportingReads = supportingReads;
		
	}

	public List<SAMRecord> getSupportingReads() {
		return supportingReads;
	}
	public void setExonInfo(ExonInfo exonInfo) {
		this.exonInfo = exonInfo;
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
	public double getPercDisc() {
		return percDisc;
	}
	public ExonInfo getExonInfo() {
		return exonInfo;
	}
	public double getTotalReads() {
		return totalReads;
	}
	
}
