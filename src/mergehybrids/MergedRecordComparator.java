package mergehybrids;

import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import MELT.utilities.VCFRecord;

public class MergedRecordComparator implements Comparator<FinalHit> {

	@Override
	public int compare(FinalHit record1, FinalHit record2) {
		
		String chrOne = record1.getChr();
		String chrTwo = record2.getChr();
		
		//First remove any 'chr' at the beginning of the chr
		Pattern chrPatt = Pattern.compile("chr(\\S+)");
		Matcher chrMatcher = chrPatt.matcher(chrOne);
		if (chrMatcher.matches()) {
			chrOne = chrMatcher.group(1);
		}
		chrMatcher = chrPatt.matcher(chrTwo);
		if (chrMatcher.matches()) {
			chrTwo = chrMatcher.group(1);
		}
		
		//Test chr first, if one is int, and other is not
		boolean chrOneIsInt = true;
		boolean chrTwoIsInt = true;
		
		try {
			Integer.parseInt(chrOne);
		} catch (NumberFormatException e) {
			chrOneIsInt = false;
		}
		try {
			Integer.parseInt(chrTwo);
		} catch (NumberFormatException e) {
			chrTwoIsInt = false;
		}
		if (chrOneIsInt == true && chrTwoIsInt == false) {
			return -1;
		} else if (chrOneIsInt == false && chrTwoIsInt == true) {
			return 1;
		} else if (chrOneIsInt == true && chrTwoIsInt == true) {
			int cmp = Integer.parseInt(chrOne) - Integer.parseInt(chrTwo);
			if (cmp != 0) {
				return cmp;
			} else {
				return record1.getStart() - record2.getStart();
			}
		} else {
			//Now compare the two chromosomes:
			int cmp = chrOne.compareTo(chrTwo);
			if (cmp != 0) {
				return cmp;
			} else {
				// Now compare positions:
				return record1.getStart() - record2.getStart();
			}		
		}
	}
	
}
