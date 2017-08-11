package quantRegions;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import MELT.utilities.EvaluateIndex;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Quant {

	public static void main(String[] args) {

		File bamFile = new File(args[0]);
		File bamIndex = EvaluateIndex.returnIndex(bamFile);
		Matcher posMatcher = Pattern.compile("(chr[0-9XY]{1,2}):(\\d+)\\-(\\d+)").matcher(args[1]);
		
		if (posMatcher.matches()) {
			
			String chr = posMatcher.group(1);
			int start = Integer.parseInt(posMatcher.group(2));
			int end = Integer.parseInt(posMatcher.group(3));
			
			int totalFragments = 0;
			
			SAMRecordIterator recordItr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(bamFile).index(bamIndex)).queryContained(chr, start, end);
			
			List<String> alreadySeen = new ArrayList<String>();
			
			while (recordItr.hasNext()) {
				
				SAMRecord current = recordItr.next();
				if (alreadySeen.contains(current.getReadName()) == false && current.getMappingQuality() >= 10 && current.getDuplicateReadFlag() == false) {
					totalFragments++;
					alreadySeen.add(current.getReadName());
				}
				
			}
			
			recordItr.close();
			
			System.out.println(bamFile.getName() + "\t" + totalFragments);
			
		} else {
			
			System.err.println("Coordinates does not conform to standard chrZZ:12345-67890!\nExiting!!!\n");
			System.exit(1);
			
		}
		
		
	}

}
