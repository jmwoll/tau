/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tau;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 *
 * @author Jan
 */
public class PACCS {
    
    public Map<String,Double> radii = new HashMap<String,Double>();
    public PACCS(){
        radii.put("H",2.0);
        radii.put("C",2.0);
        radii.put("N",2.0);
        radii.put("O",2.0);
    }
    
    public double ccs(Molecule mol){
        final double numRotamers = 6000; //final int gridSize = 6000;
        final Random randGen = new Random();
        
        double ccsSum = 0.0;
        
        for(int j=0; j < numRotamers; j++){
            final double RandRotX = randGen.nextDouble() * 4 * Math.PI;
            final double RandRotY = randGen.nextDouble() * 4 * Math.PI;
            final double RandRotZ = randGen.nextDouble() * 4 * Math.PI;
            mol.rotate(RandRotX, RandRotY, RandRotZ);
            ccsSum += ccsRotamer(mol);
        }
        return ccsSum / numRotamers;
    }
    
    public double ccsRotamer(Molecule mol){
        final Random randGen = new Random();
        final int numTrials = 8000;
        if(mol.atomXs.size() == 0)System.out.println("could not read atoms");
        final double MinX = Collections.min(mol.atomXs) - 5.0;
        final double MaxX = Collections.max(mol.atomXs) + 5.0;
        final double MinY = Collections.min(mol.atomXs) - 5.0;
        final double MaxY = Collections.max(mol.atomXs) + 5.0;
        final double MaxMinX = MaxX - MinX;
        final double MaxMinY = MaxY - MinY;
        final int NumAtoms = mol.atomSymbols.size();
        final List<String> AtomSymbols = mol.atomSymbols;
        final List<Double> AtomXs = mol.atomXs;
        final List<Double> AtomYs = mol.atomYs;
        int hits = 0;
        for(int outer=0; outer < numTrials; outer++){
            final double RandX = randGen.nextDouble() * MaxMinX + MinX;
            final double RandY = randGen.nextDouble() * MaxMinY + MinY;
            for(int inner=0; inner < NumAtoms; inner++){
                final double DX = RandX - AtomXs.get(inner);
                final double DY = RandY - AtomYs.get(inner);
                final double R = radii.get(AtomSymbols.get(inner));
                if(DX * DX + DY * DY < R * R)
                {
                    hits += 1; // no double / multiple hits
                    break;
                }
            }
        }
        return ( hits / ((double)(numTrials)) ) * (MaxMinX) * (MaxMinY);
    }
}
