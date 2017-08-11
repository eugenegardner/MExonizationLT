package ManualView;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class ManualView {

	public static void main (String args[]) throws IOException {
		
		File bams = new File(args[0]);
		File sites = new File(args[1]);
		
		Map<String, IntervalTree<File>> bamInt = loadBams(bams);
				
		BufferedReader siteReader = new BufferedReader(new FileReader(sites));
		
		String line;
		String data[];
		
		Socket socket = new Socket("devinelab-lx.igs.umaryland.edu", 60152);
		PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		
		while ((line = siteReader.readLine()) != null) {
			
			data = line.split("\\t");
			if (data[0].matches("chr")) {continue;}
			
			String chr = data[0];
			int start = Integer.parseInt(data[2]);
			int end = Integer.parseInt(data[3]);
			
			out.println("new");
			out.println("setSleepInterval 50");
			
			Iterator<Node<File>> itr = bamInt.get(chr).overlappers(start, end);
			
			System.out.println("Interval " + chr + ":" + start + "-" + end);
			System.out.println("meiStart: " + data[7]);
			System.out.println("meiEnd  : " + data[8]);
			
			while (itr.hasNext()) {
				
				Node<File> node = itr.next();
				File current = node.getValue();
				System.out.println("loading " + current.getName());
				out.println("load " + current.getAbsolutePath());
				
			}
			
			String s = "";
			while (s.length() == 0){ //Hit a character to go to next site
				s = br.readLine();   //Hit a character to go to next site
			}     					 //Hit a character to go to next site
			s = "";					 //Hit a character to go to next site
			
		}
		
		socket.close();
		siteReader.close();
		
	}
	
	private static Map<String, IntervalTree<File>> loadBams(File bams) throws IOException {
		
		Map<String, IntervalTree<File>> bamInt = new HashMap<String, IntervalTree<File>>();
		BufferedReader fileReader = new BufferedReader(new FileReader(bams));
		
		String line;
		Pattern grab = Pattern.compile("\\S+(chr[0-9XY]{1,2})_(\\d+)\\.sorted.bam");
		
		while((line = fileReader.readLine()) != null) {

			Matcher bamMatch = grab.matcher(line);
			File currentFile = new File(line);
			
			if (bamMatch.matches()) {
				String chr = bamMatch.group(1);
				int pos = Integer.parseInt(bamMatch.group(2));
				
				if (bamInt.containsKey(chr)) {
					bamInt.get(chr).put(pos, pos, currentFile);
				} else {
					IntervalTree<File> inter = new IntervalTree<File>();
					bamInt.put(chr, inter);
					bamInt.get(chr).put(pos, pos, currentFile);
				}
			}
						
		}
		
		fileReader.close();
		
		return bamInt;		
		
	}
	
}
