/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

class StartComparator implements Comparator<Motif> {
    @Override
    public int compare(Motif m1, Motif m2) {
        if (m1.alignment_start == m2.alignment_start) {
            return Double.compare(m2.score,m1.score);
        }
        return m1.alignment_start - m2.alignment_start;
    }
}

class TypeComparator implements Comparator<Motif> {
    @Override
    public int compare(Motif m1, Motif m2) {
        if (m1.type.equals(m2.type)) {
            return Integer.compare(m1.alignment_start,m2.alignment_start);
        }
        else return m1.type.compareTo(m2.type);
    }
}

class ScoreComparator implements Comparator<String_of_Beads> {
    @Override
    public int compare(String_of_Beads sob1, String_of_Beads sob2) {
        return Double.compare(sob1.getScore(),sob2.getScore());
    }
}
/**
 *
 * @author ian
 */
public class PPR {
    String id;
    String sequence;
    ArrayList<ORF> orfs;
    ArrayList<Motif> motifs;
    String_of_Beads best_sob;
    String_of_Beads best_valid_sob;
       
    public PPR(String name) {
        id = name;
        orfs = new ArrayList();
        motifs = new ArrayList();
    }
    
    public String getID() {
        return id;
    }
    
    public void setSequence(String data) {
        if (data.endsWith("*")) {
            sequence = data.substring(0,data.length()-1);  // removes trailing *
        }
        else {
            sequence = data.substring(0,data.length());
        }
    }
    
    public String getSequence() {
        orfs.sort(Comparator.comparingInt(orf -> orf.start));
        int seq_start = orfs.get(0).start;
        orfs.sort(Comparator.comparingInt(orf -> orf.end));
        int seq_end = orfs.get(orfs.size()-1).end;
        
        //System.out.println("expected sequence length of " + id + ": " + String.valueOf(seq_end-seq_start+1));
        //System.out.println("num ORFs " + String.valueOf(orfs.size()));
        //System.out.println("num motifs " + String.valueOf(best_sob.motifs.size()));
        
        /*for (Motif m : best_sob.motifs) {
            System.out.println(m.toString());
        }*/
        
        StringBuilder seq = new StringBuilder();
        char aa;
        boolean appended;
        ORF orf;
        Motif m;
        for (int pos = seq_start; pos <= seq_end; pos++ ) {
            appended = false;
            //check if this position is in a motif
            for (int i = 0; i < best_sob.motifs.size(); i++) {
                m = best_sob.motifs.get(i);
                if (m.contains(pos)) {
                    seq.append(m.getaa(this,pos));
                    appended = true;
                    break;
                }
            }
            if (!appended) {
                for (int o = 0; o < orfs.size(); o++) {
                    if (orfs.get(o).contains(pos)) {
                        seq.append(orfs.get(o).getaa(pos));
                        appended = true;
                        break;
                    }
                }
            }
            if (!appended) seq.append('X');
        }
        sequence = seq.toString();
        return sequence;
    }
    
    ORF addORF(String orfid){
        ORF o = getORF(orfid);
        if (o == null) {
            o = new ORF(orfid);
            orfs.add(o);
        }
        return o;
    }
    
    ORF getORF(String name) {
        for (ORF o: orfs){
            if (o.range.equals(name)) {
                return o;
            }
        }
        //System.out.println("couldn't find " + name + " in ORFs in " + id);
        return null;
    }
    
    public void addORFsequence(String range, String seq){
        ORF o = getORF(range);
        if (o != null) o.addSequence(seq);
    }
    
    public int end(){
        orfs.sort(Comparator.comparingInt(orf -> orf.end));
        return orfs.get(orfs.size()-1).end;
    }
    
    public int length() {
        
        orfs.sort(Comparator.comparingInt(orf -> orf.start));
        int seq_start = orfs.get(0).start;
        orfs.sort(Comparator.comparingInt(orf -> orf.end));
        int seq_end = orfs.get(orfs.size()-1).end;
        
        return seq_end - seq_start + 1;
    }
    
    public boolean hasValidSob() {
        if (best_valid_sob == null) {
            return false;
        }
        return true;
    }
      
    public void mergeMotifs(){
        
        for (ORF o : orfs){
            motifs.addAll(o.motifs);
        }
        
        if (motifs.size() < 2) return;
        
        Collections.sort(motifs, new TypeComparator());
        
        List<Motif> dead_motifs = new ArrayList();
        List<Motif> live_motifs = new ArrayList();
        
        
        Motif m1, m2;
        double hmm_coverage, alignment_coverage;
        for (int i=1; i<motifs.size(); i++){
            m1 = motifs.get(i-1);
            m2 = motifs.get(i);
            
            //check if merge is possible
            if (!m1.type.equals(m2.type)) continue;
            hmm_coverage = m1.hmm_end-m1.hmm_start+m2.hmm_end-m2.hmm_start;
            alignment_coverage = m2.alignment_end-m1.alignment_start;
            if ((hmm_coverage/Motif.getCanonicalLength(m2) > 1.5) || (alignment_coverage/Motif.getCanonicalLength(m2) > 1.5)) continue;
            
            //merge
            //System.out.println("merging " + id + " " + m2.type);
            Motif new_motif = new Motif(m2.type);
            new_motif.hmm_start = m1.hmm_start;
            new_motif.hmm_end = m2.hmm_end;
            new_motif.alignment_start = m1.alignment_start;
            new_motif.alignment_end = m2.alignment_end;
            new_motif.score = m1.score + m2.score;
            new_motif.segments = new ArrayList();
            new_motif.segments.addAll(m1.segments);
            new_motif.segments.addAll(m2.segments);
                    
            live_motifs.add(new_motif);
            dead_motifs.add(m1);
            dead_motifs.add(m2);
        }
        motifs.removeAll(dead_motifs);
        motifs.addAll(live_motifs);
        
    }
    
    public String getBestString() {
        
        Collections.sort(motifs, new StartComparator());
        ArrayList<String_of_Beads> sobs = new ArrayList();
        ArrayList<String_of_Beads> new_sobs = new ArrayList();
        ArrayList<String_of_Beads> prunings = new ArrayList();
        
        String_of_Beads new_sob;
        String_of_Beads extended_sob;
                
        for (Motif m : motifs) {
            new_sob = new String_of_Beads(m);
            new_sob.inferMotifs(this);
            if (new_sob.isValid()) {
                new_sobs.add(new_sob);
            }
            for (String_of_Beads sob : sobs) {
                if (sob.canAddMotif(m)) {
                    extended_sob = new String_of_Beads(sob.getMotifs(),m);
                    extended_sob.inferMotifs(this);
                    if (extended_sob.isValid()) {
                        new_sobs.add(extended_sob);
                    }
                }
            }
            sobs.addAll(new_sobs);
            new_sobs.clear();
            
            Collections.sort(sobs, new ScoreComparator());
            
            // prune sobs to avoid running out of memory...
            if (sobs.size() > 31) {
                String_of_Beads s;
                for (int n = 0; n < sobs.size()-1; n++) {
                    s = sobs.get(n);
                    if (n < sobs.size()/4) {
                        prunings.add(sobs.get(n));
                    }
                }
                sobs.removeAll(prunings);
                prunings.clear();
            }
            //System.out.println(id + " sobs contains " + sobs.size() + " members ");
        }
            
        best_sob = sobs.get(sobs.size()-1);
        best_sob.inferMotifs(this);
        
        return best_sob.beads(this);
    }
    
    public int get_max_motif_length(String motif_type) {
        
        int max_motif_length = 0;
        
        for (Motif m : motifs) {
            if (m.hmmer_alignment != null && m.inferred_start <= m.alignment_start && m.getType().equals(motif_type)) {
                max_motif_length = Math.max(max_motif_length,m.get_motif_length());
            }
        }
        
        return max_motif_length;
    }
    
    public String guessSubclass() {
        
        if (best_sob.hasMotifType("P")) {
            return "P";
        }
        
        return "PLS";
    }
    
    /**
     *
     * @param motif_writer
     * @throws IOException
     */
    public void writeMotifs(BufferedWriter motif_writer) throws IOException {
        
        String_of_Beads sob = best_sob;
        
        orfs.sort(Comparator.comparingInt(orf -> orf.start));
        int seq_start = orfs.get(0).start;
        
        String motif_sequence;
        
        try {
            
            for (Motif m : sob.getMotifs()) {
                motif_writer.write(id);
                motif_writer.write("\t");
                motif_writer.write(String.valueOf(m.inferred_start - seq_start + 1));
                motif_writer.write("\t");
                motif_writer.write(String.valueOf(m.inferred_end - seq_start + 1));
                motif_writer.write("\t");
                motif_writer.write(String.valueOf(m.score));
                motif_writer.write("\t");
                motif_sequence = sequence.substring(m.inferred_start - seq_start,m.inferred_end - seq_start + 1);
                motif_writer.write(motif_sequence);
                motif_writer.write("\t");
                motif_writer.write(motif_sequence.charAt(1));
                motif_writer.write("\t");
                motif_writer.write(motif_sequence.charAt(4));
                motif_writer.write("\t");
                motif_writer.write(motif_sequence.charAt(motif_sequence.length()-1));
                motif_writer.write("\t");
                motif_writer.write(m.getType());
                motif_writer.newLine();    
            }
        
        } catch (StringIndexOutOfBoundsException ex) {
               System.out.println("StringIndexOutOfBoundsException for " + id);
        }
    
    }
    
}

