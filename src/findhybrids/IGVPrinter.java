package findhybrids;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.ExecuteException;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import mergehybrids.CheckAssembly.MatchedMEI;
import mergehybrids.FinalHit;

public class IGVPrinter {

	private Socket socket;
	private PrintWriter out;
	private BufferedReader br;
	private IndexedFastaSequenceFile reference;
	
	public IGVPrinter (int port, IndexedFastaSequenceFile reference) throws UnknownHostException, IOException {

		socket = new Socket("devinelab-lx.igs.umaryland.edu", port); 
		out = new PrintWriter(socket.getOutputStream(), true);
		br= new BufferedReader(new InputStreamReader(System.in));
		this.reference = reference;
		
	}
	
	public void viewSite (FinalHit hit) throws IOException {
		
		String s = "";
		out.println("goto " + hit.getChr() + ":" + hit.getStart() + "-" + hit.getStop()); //socket printer
		while (s.length() == 0){ //Hit a character to go to next site
			s = br.readLine();   //Hit a character to go to next site
		}     					 //Hit a character to go to next site
		s = "";					 //Hit a character to go to next site
	}
	public void viewSite(FinalHit hit, Map<String, Integer> include, MatchedMEI mei) throws IOException {
		
		String s = "";
		out.println("new");
		out.println("setSleepInterval 50");
		if (mei == null) {
			out.println("goto " + hit.getChr() + ":" + hit.getStart() + "-" + hit.getStop()); //socket printer
		} else if (mei.getStart() < hit.getStart()) {
			out.println("goto " + hit.getChr() + ":" + mei.getStart() + "-" + hit.getStop()); //socket printer
		} else {
			out.println("goto " + hit.getChr() + ":" + hit.getStart() + "-" + mei.getStop()); //socket printer
		}
		out.println("load /local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/VCFFiles/indexedVCF/SINEC.pass.vcf.gz");
		out.println("load /home/eugene.gardner/Old_MELT_Versions/MELTv1.2_Internal/me_refs/CanFam3.1/test/SINEC.bed");
		out.println("load /home/eugene.gardner/Old_MELT_Versions/MELTv1.2_Internal/add_bed_files/CanFam3.1/canFam3.genes.bed");
		for (Map.Entry<String, Integer> entry : include.entrySet()) {
			if (entry.getValue() > 0) {
				out.println("load /local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/analysis/rna-seq/" + entry.getKey() + "/hybrids.sorted.bam");
				out.println("load /local/aberdeen2rw/DEVINE/DOG/AnalysisApr2016/FinalVCF/FinalVCF_Aug05/analysis/rna-seq/" + entry.getKey() + "/" + entry.getKey() + ".sorted.bam");
				out.println("viewaspairs hybrids.sorted.bam");
				System.out.println(entry.getKey() + " : " + entry.getValue());
			}
		}
		while (s.length() == 0){ //Hit a character to go to next site
			s = br.readLine();   //Hit a character to go to next site
		}     					 //Hit a character to go to next site
		s = "";					 //Hit a character to go to next site
		
	}
	
	public void viewSINE (String chr, int start, Map<String, IntervalTree<String>> sines, File assemblyPath) throws ExecuteException, IOException {
		
		if (!new File(assemblyPath.getAbsolutePath() + ".nonref/").isDirectory()) {
			new File(assemblyPath.getAbsolutePath() + ".nonref/").mkdir();
		}
		
		DefaultExecutor executor = new DefaultExecutor();
		
		String faPath = assemblyPath.getAbsolutePath() + ".nonref/" + chr + "_" + start + ".ref";
		
		BufferedWriter refWriter = new BufferedWriter(new FileWriter(new File(faPath + ".fa")));
		
		String toWrite = sines.get(chr).find(start, start).getValue();
		refWriter.write(">" + chr + "_" + start + "\n");
		refWriter.write(toWrite);
		refWriter.newLine();
		refWriter.flush();
		refWriter.close();		
		
		executor.execute(CommandLine.parse("bowtie2 -f --local --quiet -x /home/eugene.gardner/Old_MELT_Versions/MELTv1.2_Internal/me_refs/CanFam3.1/test/SINEC -U " + faPath + ".fa -S " + faPath + ".sam"));
		executor.execute(CommandLine.parse("samtools view -Sbo " + faPath + ".bam " + faPath + ".sam"));
		executor.execute(CommandLine.parse("samtools index " + faPath + ".bam"));
		
		System.out.println("\tLoading maximum nonref SINE at " + start);
		out.println("load " + faPath + ".bam");
				
	}
	
	public void viewRefSine (String chr, int start, int stop, File assemblyPath) throws IOException {
		
		if (!new File(assemblyPath.getAbsolutePath() + ".ref/").isDirectory()) {
			new File(assemblyPath.getAbsolutePath() + ".ref/").mkdir();
		}
		
		DefaultExecutor executor = new DefaultExecutor();
		
		String faPath = assemblyPath.getAbsolutePath() + ".ref/" + chr + "_" + start + ".ref";
		
		BufferedWriter refWriter = new BufferedWriter(new FileWriter(new File(faPath + ".fa")));
		String toWrite = reference.getSubsequenceAt(chr, start, stop).getBaseString().toUpperCase();
		refWriter.write(">" + chr + "_" + start + "-" + stop + "\n");
		refWriter.write(toWrite);
		refWriter.newLine();
		refWriter.flush();
		refWriter.close();		
		
		executor.execute(CommandLine.parse("bowtie2 -f --local --quiet -x /home/eugene.gardner/Old_MELT_Versions/MELTv1.2_Internal/me_refs/CanFam3.1/test/SINEC -U " + faPath + ".fa -S " + faPath + ".sam"));
		executor.execute(CommandLine.parse("samtools view -Sbo " + faPath + ".bam " + faPath + ".sam"));
		executor.execute(CommandLine.parse("samtools index " + faPath + ".bam"));
		
		System.out.println("\tLoading maximum ref SINE at " + start + "-" + stop);
		out.println("load " + faPath + ".bam");
		
	}
	
	public void viewConsensus (String consensus, File assemblyPath, String chr, int pos) throws ExecuteException, IOException {
		
		out.println("new");
		out.println("setSleepInterval 50");
		out.println("load " + assemblyPath.getAbsolutePath() + "/" + chr + "_" + pos + ".sorted.bam");
		
		if (!new File(assemblyPath.getAbsolutePath() + ".cons/").isDirectory()) {
			new File(assemblyPath.getAbsolutePath() + ".cons/").mkdir();
		}
		
		DefaultExecutor executor = new DefaultExecutor();
		
		String faPath = assemblyPath.getAbsolutePath() + ".cons/" + chr + "_" + pos + ".cons";
		
		BufferedWriter refWriter = new BufferedWriter(new FileWriter(new File(faPath + ".fa")));
		refWriter.write(">consensus\n");
		refWriter.write(consensus);
		refWriter.newLine();
		refWriter.flush();
		refWriter.close();		
		
		executor.execute(CommandLine.parse("bowtie2 -f --local --quiet -x /home/eugene.gardner/Old_MELT_Versions/MELTv1.2_Internal/me_refs/CanFam3.1/test/SINEC -U " + faPath + ".fa -S " + faPath + ".sam"));
		executor.execute(CommandLine.parse("samtools view -Sbo " + faPath + ".bam " + faPath + ".sam"));
		executor.execute(CommandLine.parse("samtools index " + faPath + ".bam"));
		
		out.println("load " + faPath + ".bam");
		
	}
	
	public void closeHandles () throws IOException {
		out.close();
		socket.close();
	}

	
	
	
	
}
