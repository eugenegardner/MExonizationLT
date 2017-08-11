package findhybrids;

public class ExonInfo {

	private String geneName;
	private boolean isForward;
	private ExonType type;
	private int exonNum;
	
	public ExonInfo (String geneName, boolean isForward, ExonType type, int exonNum) {
		
		this.geneName = geneName;
		this.isForward = isForward;
		this.type = type;
		this.exonNum = exonNum;
		
	}
	
	public String getGeneName() {
		return geneName;
	}
	public boolean isForward() {
		return isForward;
	}
	public ExonType getType() {
		return type;
	}
	public int getExonNum() {
		return exonNum;
	}
	
	public enum ExonType {
		
		EXON,FUTR,TUTR,FAKE_FUTR,FAKE_TUTR;
		
	}
	
}
