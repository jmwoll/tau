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
 public class Molecule {
        List<String> atomSymbols = new ArrayList<String>();
        List<Double> atomXs = new ArrayList<Double>();
        List<Double> atomYs = new ArrayList<Double>();
        List<Double> atomZs = new ArrayList<Double>();
        
        public Molecule(List<String> atomSymbols, List<Double> atomXs,
                List<Double> atomYs,
                List<Double> atomZs){
            this.atomSymbols = atomSymbols;
            this.atomXs = atomXs;
            this.atomYs = atomYs;
            this.atomZs = atomZs;
        }
        
        public void rotate(double rotX, double rotY, double rotZ){
            double cx = 0.0; double cy = 0.0; double cz = 0.0;
            final int NumAtoms = atomSymbols.size();
            for(int i=0; i < NumAtoms; i++){
                cx += atomXs.get(i); cy += atomYs.get(i); cz += atomZs.get(i);
            }
            cx = cx / NumAtoms; cy = cy / NumAtoms; cz = cz / NumAtoms;
            for(int i=0; i < NumAtoms; i++){
                double x = atomXs.get(i) - cx; double y = atomYs.get(i) - cy;
                double z = atomZs.get(i) - cz;
                // rotate around x
                final double[] RX = new double[]{y * Math.cos(rotX) - z * Math.sin(rotX),
                    y * Math.sin(rotX) + z * Math.cos(rotX)};
                y = RX[0]; z = RX[1];
                // rotate around y
                final double[] RY = new double[]{x * Math.cos(rotY) + z * Math.sin(rotY),
                    - x * Math.sin(rotY) + z * Math.cos(rotY)};
                x = RY[0]; z = RY[1];
                // rotate around z
                final double[] RZ = new double[]{x * Math.cos(rotZ) - y * Math.sin(rotZ),
                    x * Math.sin(rotZ) + y * Math.cos(rotZ)};
                x = RZ[0]; y = RZ[1];
                // finally submit the rotated coordinate values for the atom
                atomXs.set(i, x + cx); atomYs.set(i, y + cy); atomZs.set(i, z + cz);
            }
        }
        
        public static Molecule fromXYZString(String xyzString){
            final List<String> atomSymbols = new ArrayList<String>();
            final List<Double> atomXs = new ArrayList<Double>();
            final List<Double> atomYs = new ArrayList<Double>();
            final List<Double> atomZs = new ArrayList<Double>();
            xyzString = xyzString.trim().replaceAll("\t", " ").replaceAll("  ", " ");
            while(xyzString.contains("  ")){
                xyzString = xyzString.replaceAll("  ", " ");
            }
            int lneCount = 0;
            for(String lne : xyzString.trim().split("\n")){
                lneCount += 1;
                final String[] atom = lne.trim().split(" ");
                if(atom.length != 4){
                    //if(lneCount > 2)
                        System.out.println("Skipping line: "+lne
                            +"line count= "+atom.length);
                }
                else{
                    atomSymbols.add(atom[0]);
                    atomXs.add(Double.parseDouble(atom[1]));;
                    atomYs.add(Double.parseDouble(atom[2]));;
                    atomZs.add(Double.parseDouble(atom[3]));
                }
            }
            return new Molecule(atomSymbols, atomXs, atomYs, atomZs);
        }
}
