/* TernaryFFT.java
 * Author: Austin Corotan
 * Date: 5/7/2017
 * Description: This program will calculate the fft and the ifft for a
 * polynomial given its respective coefficients.
 */

import java.util.*;
import java.lang.*;
import java.io.*;

public class TernaryFFT {

    // compute the FFT of x[], assuming its length is a power of 3
    public static Complex[] fft(Complex[] x) {
        int n = x.length;

        // base case
        if (n == 1) return new Complex[] { x[0] };

        if (n % 3 != 0) { throw new RuntimeException("n is not a power of 3"); }

        //fft of first terms
        Complex[] first = new Complex[n/3];
        for(int k = 0; k < n/3; k++){
            first[k] = x[3*k];
        }
        Complex[] a = fft(first);

        // fft of second terms
        Complex[] second = new Complex[n/3];
        for (int k = 0; k < n/3; k++) {
            second[k] = x[3*k+1];
        }
        Complex[] b = fft(second);

        // fft of thrid terms
        Complex[] third = new Complex[n/3];
        for (int k = 0; k < n/3; k++) {
            third[k] = x[3*k+2];
        }

        Complex[] c = fft(third);

        // combine
        Complex w = new Complex(1, 0);
        Complex w2 = new Complex(1, 0);
        double kth = 2 * Math.PI / n;
        double kth2 = 4 * Math.PI / n;
        Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
        Complex wk2 = new Complex(Math.cos(kth2), Math.sin(kth2));
        
        double nth = 2 * Math.PI / 3;
        double nth2 = 4 * Math.PI / 3;
        Complex wn = new Complex(Math.cos(nth), Math.sin(nth));
        Complex wn2 = new Complex(Math.cos(nth2), Math.sin(nth2));

        Complex[] combined = new Complex[n];

        for (int k = 0; k < n/3; k++) {
            combined[k]       = (a[k].plus(w.times(b[k]))).plus(w2.times(c[k]));
            combined[k + n/3] = (a[k].plus(w.times(b[k]).times(wn))).plus(w2.times(c[k]).times(wn2));
            combined[k + 2*n/3] = (a[k].plus(w.times(b[k]).times(wn2))).plus(w2.times(c[k]).times(wn));
            w = w.times(wk);
            w2 = w2.times(wk2);
        }
        return combined;
    }

    /* finishes up computation for the recursive ifft */
    public static Complex[] ifft(Complex[] x) {
        int n = x.length;
        Complex[] y = new Complex[n];

        y = Recursiveifft(x);

        for (int i = 0; i < n; i++) {
            y[i] = y[i].times(1.0/n);
        }

        return y;
    }

    /* recursively computes the inverse fast fourier transform */
    public static Complex[] Recursiveifft(Complex[] x) {
        int n = x.length;

        // base case
        if (n == 1) return new Complex[] { x[0] };

        if (n % 3 != 0) { throw new RuntimeException("n is not a power of 3"); }


        //fft of first terms
        Complex[] first = new Complex[n/3];
        for(int k = 0; k < n/3; k++){
            first[k] = x[3*k];
        }
        Complex[] a = Recursiveifft(first);

        // fft of second terms
        Complex[] second = new Complex[n/3];
        for (int k = 0; k < n/3; k++) {
            second[k] = x[3*k+1];
        }
        Complex[] b = Recursiveifft(second);

        // fft of thrid terms
        Complex[] third = new Complex[n/3];
        for (int k = 0; k < n/3; k++) {
            third[k] = x[3*k+2];
        }

        Complex[] c = Recursiveifft(third);

        // combine
        Complex w = new Complex(1, 0);
        Complex w2 = new Complex(1, 0);
        double kth = -2 * Math.PI / n;
        double kth2 = -4 * Math.PI / n;

        Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
        Complex wk2 = new Complex(Math.cos(kth2), Math.sin(kth2));
        
        double nth = 2 * Math.PI / 3;
        double nth2 = 4 * Math.PI / 3;
        Complex wn = new Complex(Math.cos(nth), Math.sin(nth));
        Complex wn2 = new Complex(Math.cos(nth2), Math.sin(nth2));

        Complex wn_inv = wn.reciprocal();
        Complex wn2_inv = wn2.reciprocal();

        Complex[] combined = new Complex[n];

        for (int k = 0; k < n/3; k++) {
            combined[k]       = (a[k].plus(w.times(b[k]))).plus(w2.times(c[k]));
            combined[k + n/3] = (a[k].plus(w.times(b[k]).times(wn_inv))).plus(w2.times(c[k]).times(wn2_inv));
            combined[k + 2*n/3] = (a[k].plus(w.times(b[k]).times(wn2_inv))).plus(w2.times(c[k]).times(wn_inv));
            w = w.times(wk);
            w2 = w2.times(wk2);
        }
        return combined;
    }

    /* reads input text file and stores data in a single array */
    static double[] readInput(File inputFile) throws IOException{
        int size = 0;
        Scanner input = new Scanner (inputFile);
        int i = 0;
        String line;
        while(input.hasNextLine()){
            size++;
            input.nextLine();
        }

        double[] A = new double[size];
        input = new Scanner (inputFile);
        while(input.hasNextLine()){
            line = input.nextLine();
            String[] vals = line.trim().split("\\s+");
            A[i] = Double.parseDouble(vals[0]);
            i++;
        }

        return A;
    }

    /* Accepts the coefficient input text file and passes it to the fft and ifft functions
     * the displays the results of each respective function */
    public static void main(String[] args) { 

        try{
            double[] A;
            File inputFile = new File(args[0]);
            A = readInput(inputFile);

            int n = A.length;
            Complex[] a = new Complex[n];

            // original data
            for (int i = 0; i < n; i++) {
                a[i] = new Complex(A[i], 0);
            }

            // FFT 
            Complex[] b = fft(a);

            // inverse FFT
            Complex[] c = ifft(b);

            //Display Data
            System.out.printf("%7s %12s %10s %10s %10s\n", "a", "Re(b)", "Im(b)", "Re(c)", "Im(c)");
            for (int i = 0; i < a.length; i++) {
                System.out.printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n", a[i].re(), b[i].re(), b[i].im(), c[i].re(), c[i].im());
            }
            System.out.println();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}