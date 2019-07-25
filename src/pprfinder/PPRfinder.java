/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 *
 * @author ian
 */
public class PPRfinder {
    
    public static HashMap<String,PPR> PPRs = new HashMap();
    public static HashMap<String,Boolean> Gapless_Motif_Combinations = new HashMap();
    public static HashMap<String,Boolean> Motif_Combinations = new HashMap();
    
    /**
     * @param args the command line arguments
     * args[0] protein sequence file (fasta)
     * args[1] --domtout table from HMMER3.2
     */
    
    public static void main(String[] args) {
        
        //check arguments
        if (args.length != 2 || !(new File(args[0]).exists()) || !(new File(args[1]).exists())){
            System.out.println("Program: PPRfinder (for characterising PPR motifs in proteins)");
            System.out.println("Version: 1.0");
            System.out.println("Usage: java -jar PPRfinder.jar <orfs.fasta> <orfs.domt>");
            System.out.println("<orfs.fasta>: fasta file of protein sequences");
            System.out.println("<orfs.domt>: domain table output from running hmmsearch on <orfs.fasta> with PPR HMMs");
            System.exit(0);
        }
        
        try {
            //read --domtout table
            System.out.println("reading " + args[1]);
            BufferedReader reader = new BufferedReader(new FileReader(args[1]));
            String s;
            String pprid, orfid;
            PPR p;    
            ORF o;
            Motif m;
            String[] fields;
             
            while ( (s=reader.readLine()) != null) {
                if (s.startsWith("#")) {
                    continue;
                }
                fields = s.split(" +");
                
                //ignore coordinates in ORF name, to collect all motifs on same strand, whatever frame they're in
                pprid = fields[0].substring(0,fields[0].lastIndexOf(":"));
                orfid = fields[0].substring(fields[0].lastIndexOf(":")+1);
                p = PPRs.get(pprid);
                if (p == null) {
                    p = new PPR(pprid);
                    PPRs.put(pprid,p);
                }
                o = p.addORF(orfid);
                //System.out.println("adding ORF " + o.range + " to " + pprid);
                m = new Motif(fields);
                o.addMotif(m);
                //System.out.println(m.toString());
            }
            reader.close();
            System.out.println("read " + PPRs.size() + " PPRs");
            
            //read fasta sequence file
            System.out.println("reading sequences from " + args[0]);
            reader = new BufferedReader(new FileReader(args[0]));
            StringBuilder sb = new StringBuilder(3000);
            p = null;
            fields = null;
            String range = null;
            while ( (s=reader.readLine()) != null) {
                if (s.startsWith(">")) {
                    //finish previous ORF if there is one
                    if ( p != null ) {
                        p.addORFsequence(range,sb.toString());
                    }
                    pprid = s.substring(1,s.lastIndexOf(":"));
                    p = PPRs.get(pprid);
                    range = s.substring(s.lastIndexOf(":")+1);
                    sb.setLength(0);
                }
                else if (p != null) {
                    sb.append(s);
                }
            }
            if ( p != null ) {
                p.setSequence(sb.toString());
                sb.setLength(0);
            }
            reader.close();
            
            //read valid motif combinations
            InputStream is = PPRfinder.class.getClassLoader().getResourceAsStream("resources/motif_combinations.txt");
            reader = new BufferedReader(new InputStreamReader(is));
            while ( (s=reader.readLine()) != null) {
                fields = s.split("\t");
                Motif_Combinations.put(fields[0], fields[1].equals("1") ? Boolean.TRUE:Boolean.FALSE);
            }
            reader.close();
                        
            is = PPRfinder.class.getClassLoader().getResourceAsStream("resources/gapless_motif_combinations.txt");
            reader = new BufferedReader(new InputStreamReader(is));
            while ( (s=reader.readLine()) != null) {
                fields = s.split("\t");
                Gapless_Motif_Combinations.put(fields[0], fields[1].equals("1") ? Boolean.TRUE:Boolean.FALSE);
            }
            reader.close();
            
            //System.out.println("inferring motifs... ");
            BufferedWriter motif_writer = new BufferedWriter(new FileWriter(args[1]+"_motifs.txt"));
            BufferedWriter orf_writer = new BufferedWriter(new FileWriter(args[1]+"_pprs.fa"));
            BufferedWriter beads_writer = new BufferedWriter(new FileWriter(args[1]+"_beads.txt"));
            
            String beads;
            double score;
            
            String sequence;
            
            for (PPR ppr : PPRs.values()) {
                //System.out.println(ppr.getID());
                //System.out.println("merging motifs...");
                ppr.mergeMotifs();
                /*for (Motif motif : ppr.motifs) {
                    System.out.println(motif.toString());
                }*/
                if (ppr.motifs.isEmpty()) continue;
                //System.out.println("constructing beads...");
                beads = ppr.getBestString();
                //System.out.println("calculating score...");
                score = ppr.best_sob.getScore();
                //System.out.println("generating sequence ...");
                sequence = ppr.getSequence();
                if (score >= 40 && (ppr.best_sob.hasAdjacentMotifs() || ppr.best_sob.hasDYW())) {
                    beads_writer.write(ppr.id + "\t" + sequence.length()  + "\t" + beads + "\t" + ppr.guessSubclass() + "\t" + score);
                    beads_writer.newLine();  
                    beads_writer.flush(); 
                    orf_writer.write(">" + ppr.id);
                    orf_writer.newLine();
                    orf_writer.write(sequence);
                    orf_writer.newLine();
                    orf_writer.flush(); 
                    ppr.writeMotifs(motif_writer);
                    motif_writer.flush();
                }
            }
                        
            orf_writer.close();
            motif_writer.close();
            beads_writer.close();
            
        }
        
        catch(Exception ex) {
            ex.printStackTrace();
            System.out.println(ex.toString());
        }
    }
}
