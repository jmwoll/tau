/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tau;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Jan
 */
public class Tau {
    
   
        
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        final PACCS p = new PACCS();
        String xyz = "";
        for(int i=0; i < 100; i++){
            xyz += "C "+i+" "+i+" 1\n";
        }
        final Molecule m = Molecule.fromXYZString(xyz);
        System.out.println(p.ccs(m));
    }
    
}
