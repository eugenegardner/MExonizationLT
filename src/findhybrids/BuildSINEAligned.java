package findhybrids;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BuildSINEAligned {

	private Map<String, SAMRecord> finalReads;
	
	public BuildSINEAligned (File firstMates, File secondMates) {
		
		SamReaderFactory samFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		
		SamReader mateOneReader = samFactory.open(firstMates);
		SamReader mateTwoReader = samFactory.open(secondMates);
		Map<String, SAMRecord> firstMateHash = buildMaps(mateOneReader, Mate.FIRST);
		Map<String, SAMRecord> secondMateHash = buildMaps(mateTwoReader, Mate.SECOND);
		finalReads = filterMates(firstMateHash, secondMateHash);
		
				
	}
	
	public Map<String, SAMRecord> getFinalReads() {
		return finalReads;
	}
	
	private Map<String, SAMRecord> buildMaps (SamReader currentReader, Mate mate) {
		
		SAMRecordIterator itr = currentReader.iterator();
		Map<String, SAMRecord> reads = new HashMap<String, SAMRecord>();
		
		while (itr.hasNext()) {
			
			SAMRecord read = itr.next();
			String name = read.getReadName();
			read.setReadPairedFlag(true);
			read.setFirstOfPairFlag(mate.isFirst());
			reads.put(name, read);
			
		}
		
		return reads;
		
	}
	private Map<String, SAMRecord> filterMates (Map<String, SAMRecord> firstMates, Map<String, SAMRecord> secondMates) {
		
		Map<String, SAMRecord> reads = new HashMap<String, SAMRecord>();
		
		for (Map.Entry<String, SAMRecord> entry : firstMates.entrySet()) {
			
			if (!secondMates.containsKey(entry.getKey())) {
				reads.put(entry.getKey(), entry.getValue());
			}
			
		}
		for (Map.Entry<String, SAMRecord> entry : secondMates.entrySet()) {
			
			if (!firstMates.containsKey(entry.getKey())) {
				reads.put(entry.getKey(), entry.getValue());
			}
			
		}
		
		return reads;
		
	}
	
	enum Mate {
		FIRST(true),
		SECOND(false);
		
		boolean isFirst;
		Mate(boolean isFirst) {
			this.isFirst = isFirst;
		}
		private boolean isFirst () {
			return isFirst;
		}
	}
	
}
