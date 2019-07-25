/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.util.ArrayList;

/**
 *
 * @author ian
 */
public class ORF {
    
    String range;
    String sequence;
    
    int start;
    int end;
    
    ArrayList<Motif> motifs;
    
    public ORF(String name){
        range = name;
        String[] coords = range.split("-");
        start = Integer.parseInt(coords[0]);
        end = Integer.parseInt(coords[1]);
        motifs = new ArrayList();
    }
    
    public void addMotif(Motif m) {
        
        // filter to ignore DYW variants
        
        //if (m.type.equals("K123")) return;
        //if (m.type.equals("WW")) return;
        
        //bitscore filter
        
        if (m.type.equals("DYW") && m.score < 30) m.hidden = true;
        else if (m.type.equals("KPAKA") && m.score < 30) m.hidden = true;
        //else if (m.type.equals("DYW:PGWW") && m.score < 90) m.hidden = true;
        else if (m.type.equals("SS") && m.score < 10) m.hidden = true;
        else if (m.score < 0) m.hidden = true;
        
        if (!m.hidden) motifs.add(m);
    }
    
    public void addSequence(String s){
        sequence = s;
    }
    
    public boolean contains(int pos){
        if (pos >= start && pos <= end){
            return true;
        }
        return false;
    }
    
    public char getaa(int pos){
        return sequence.charAt(pos-start);
    }
    
    @Override
    public String toString(){
        return range;
    }
    
}
