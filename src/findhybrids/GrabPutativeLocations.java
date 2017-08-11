package findhybrids;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class GrabPutativeLocations {

	//Walk across the genome and grab signiciant hits (genomeWalker)
	public static List<PutativeHit> GrabLocations (File outputBam, Map<String, Integer> chrs) throws IOException {
		
		return genomeWalker(outputBam, chrs);
		
	}	
	
	private static List<PutativeHit> genomeWalker (File outputBam, Map<String, Integer> chrs) throws IOException {
		
		List<PutativeHit> regions = new ArrayList<PutativeHit>();
		int insertCounter = 0;		
		
		for (Map.Entry<String, Integer> chrEntry : chrs.entrySet()) {
						
			SamReader refAlignedReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(outputBam).index(new File(outputBam.getAbsolutePath() + ".bai")));
						
			SAMRecordIterator itr = refAlignedReader.queryContained(chrEntry.getKey(), 0, chrEntry.getValue());
			
			int currentStart;
			int currentStop;
			int previousStart = 0;
			int previousStop = 0;
			int finalStart = 0;
			int finalStop = 0;
			double numReads = 0;
			double numDisc = 0;
			List<SAMRecord> supportingReads = new ArrayList<SAMRecord>();
			
			while(itr.hasNext() == true) {
				
				SAMRecord currentRead = itr.next();
				if (currentRead.getStringAttribute("US").matches("FALSE")) {continue;}
				
				if (currentRead.getMappingQuality() != 0) {
					currentStart = currentRead.getAlignmentStart();
					currentStop = currentRead.getAlignmentEnd();
					if ((currentStart >= (previousStart - 500) && currentStart <= (previousStop + 500)) && insertCounter != 0) {
						
						finalStop = currentStop;
						previousStart = currentStart;
						previousStop = currentStop;
						numReads++;
						supportingReads.add(currentRead);
						if (currentRead.getInferredInsertSize() >= 500 || currentRead.getInferredInsertSize() <= -500 || (currentRead.getMateReferenceName() != currentRead.getReferenceName())) {
							numDisc++;
						}
						
					} else {
						if (insertCounter != 0 && numReads > 4 && ((numDisc / numReads) >= .5)) {
							regions.add(new PutativeHit(chrEntry.getKey(), finalStart, finalStop, numReads, numDisc, supportingReads));
						}
						insertCounter++;
						finalStart = currentStart;
						finalStop = currentStop;
						previousStart = currentStart;
						previousStop = currentStop;
						numReads = 1;
						supportingReads = new ArrayList<SAMRecord>();
						supportingReads.add(currentRead);
						
						if (currentRead.getInferredInsertSize() >= 500 || currentRead.getInferredInsertSize() <= -500 || (currentRead.getMateReferenceName() != currentRead.getReferenceName())) {
							numDisc = 1;
						} else {
							numDisc = 0;
						}
						
					}			
				}
			}
			
			if (numReads > 3 && ((numDisc / numReads) >= .2)) {
				regions.add(new PutativeHit(chrEntry.getKey(), finalStart, finalStop, numReads, numDisc, supportingReads));
			}
			refAlignedReader.close();
			itr.close();
			
		}

		return regions;
		
	}
	


}
