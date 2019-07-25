/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pprfinder;

import java.util.HashMap;
import java.util.Random;

/**
 *
 * @author ian
 */
public class Blosum62 {
    
    String[] aas = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};
    double [] proportions = {0.0791, 0.0512, 0.0441, 0.0519, 0.0267, 0.0319, 0.0695, 0.0662, 0.0246, 0.0552, 0.1050, 0.0639, 0.0385, 0.0394, 0.0294, 0.0670, 0.0422, 0.0121, 0.0349, 0.0672};
    
    int[][] blosum = new int [][]{
	{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 },
	{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 },
	{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 },
	{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 },
	{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },
	{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 },
	{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 },
	{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 },
	{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 },
	{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 },
	{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 },
	{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 },
	{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 },
	{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 },
	{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 },
	{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 },
	{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 },
	{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 },
	{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 },
	{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4 }

    };
    
    
    static HashMap<String, Integer> b62 = new HashMap();
    
    public Blosum62() {
        
        for (int i=0; i < 20 ; i++){
            for (int j=0; j < 20 ; j++){
                b62.put(aas[i]+aas[j],blosum[i][j]);
            }
        }
         
    }
    
    static  public int score (String a, String b) {
        
        try {
            Integer score = b62.get(a + b);
            return score;
        }
        catch (NullPointerException ex) {
            //System.out.println("unknown aa: " + a + b);
            return 0;
        }
        
        
    }
    void random_associations() {
	Random rand = new Random();
	RunningStat rs = new RunningStat();
	String aa1,aa2;
	double score;
	for (int i=0; i < 100000000; i++) {
		aa1 = random_aa(rand);
		aa2 = random_aa(rand);
		score = score(aa1,aa2);
		rs.Push(score);
		//System.out.println(aa1+aa2+ " " + String.valueOf(score));
	}
	System.out.println(String.valueOf(rs.Mean()));
	System.out.println(String.valueOf(rs.StandardDeviation()));
    }

    String random_aa(Random rand) {
            double r = rand.nextDouble();
            double p = 0.0;
            int index = 0;
            while (r > p) {
                    p += proportions[index++];
            }
            return aas[index-1];
    }
    
}
