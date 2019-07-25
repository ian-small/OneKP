/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

/**
 *
 * @author ian
 */
public class MotifSegment {
    
    String range;
    int orf_start;
    int orf_end;
    int motif_segment_start;
    int motif_segment_end;
    
    public MotifSegment(String data){
        //get coordinates from ORF name
        range = data.substring(data.lastIndexOf(":")+1);
        String[] coords = range.split("-");
        orf_start = Integer.parseInt(coords[0]);
        orf_end = Integer.parseInt(coords[1]);
    }
    
    public boolean contains(int pos){
        if (pos >= motif_segment_start && pos <= motif_segment_end) return true;
        return false;
    }
    
}