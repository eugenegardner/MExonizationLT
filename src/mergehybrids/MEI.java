package mergehybrids;

import mergehybrids.CheckAssembly.MEIType;

public class MEI {

	private int start;
	private int stop;
	private double identity;
	private int alignmentScore;
	private MEIType type;
	
	public MEI (int start, int stop, double identity, int alignmentScore, MEIType type) {
		
		this.start = start;
		this.stop = stop;
		this.identity = identity;
		this.alignmentScore = alignmentScore;
		this.type = type;
					
	}

	public int getStart() {
		return start;
	}
	public int getStop() {
		return stop;
	}
	public double getIdentity() {
		return identity;
	}
	public int getAlignmentScore() {
		return alignmentScore;
	}
	public MEIType getType() {
		return type;
	}
	
	
}
