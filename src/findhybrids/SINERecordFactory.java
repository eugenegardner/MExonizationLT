package findhybrids;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import MELT.MELTIllumina.mergeSAM.MergeSAMMain;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

public class SINERecordFactory {

	private SAMFileWriter pairWriter;
	private SAMFileHeader useHeader;
	private String meiName;
	
	public SINERecordFactory(File folder, String name, File meiIndex, String meiName) throws IOException {
		
		if (!folder.isDirectory()) {
			folder.mkdir();
		}
		this.meiName = meiName;
		SAMSequenceDictionary seqDic = MergeSAMMain.buildSeqDictionary(meiIndex);
		SAMFileHeader useHeader = new SAMFileHeader();
		useHeader.setSequenceDictionary(seqDic);
		useHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		useHeader.setSortOrder(SortOrder.coordinate);
		SAMFileWriterFactory samFactory = new SAMFileWriterFactory().setCreateIndex(true);
		pairWriter = samFactory.makeBAMWriter(useHeader, false, new File(folder.getAbsolutePath() + "/" + name + ".sorted.bam"));
		
	}
	
	//Method to dump reads
	public void AddSINEAlignment(SAMRecord matchedNormal) {
			
		SAMRecord toDump = new SAMRecord(useHeader);
		Cigar cigar = buildCigar(matchedNormal.getStringAttribute("OC"));
		String sequence = matchedNormal.getStringAttribute("R2");
		String qual = matchedNormal.getStringAttribute("Q2");
		String name = matchedNormal.getReadName();
		int start = matchedNormal.getIntegerAttribute("OS");
		int flags = matchedNormal.getIntegerAttribute("OF");
		int mapQ = matchedNormal.getIntegerAttribute("OQ");
		
		toDump.setReadName(name);
		toDump.setReferenceName(meiName);
		toDump.setAlignmentStart(start);
		toDump.setFlags(flags);
		toDump.setReadString(sequence);
		toDump.setBaseQualityString(qual);
		toDump.setCigar(cigar);
		toDump.setMateAlignmentStart(0);
		toDump.setMateReferenceIndex(-1);
		toDump.setMappingQuality(mapQ);
							
		pairWriter.addAlignment(toDump);
								
	}
	
	public void closeWriter() {
		pairWriter.close();
	}
	
	private Cigar buildCigar(String cigarString) {
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
