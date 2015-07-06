package bielecki;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Scanner;

public class Schrodinger {

	static int N;
	static int n;
	static double L;
	static double dx;
	static double kappa;
	static double tau=0;
	static double dtau;
	static double omega;
	static int Sd;
	static int Sout;
	
	static double Norm;
	static double Xsr;
	static double Energy;
	static double Gsum;
	static double Gp;
	static double Gesum;
	
	
	static double[] arguments;
	static double[] fic;
	static double[] Fir;
	static double[] Fii;
	static double[] Hr;
	static double[] Hi;
	static double[] Gk;
	
	
	public static void main(String[] args) throws IOException {
		
		loadData("dane");
		omega=1*5*omega;
		
		arguments = new double[N+1];
		calculateArguments();
		
		Fii = new double[N+1];
		Fir = new double[N+1];
		initiateFi();
		
		Hi = new double[N+1];
		Hr = new double[N+1];
		calculateHr();
		calculateHi();
		
		makeStep();

		Gk = new double[N+1];
		calculateG();
		calculateNEX();
			
		File fileOuttt = new File("ro2");
	    FileOutputStream foppp = new FileOutputStream(fileOuttt);
	    for (int i=0;i<N;i++)
	    {
	    	foppp.write(String.format("%f\t",Gk[i]).getBytes());
	    }
	    foppp.write(String.format("\n").getBytes());
		
		File fileOut = new File("ro");
        FileOutputStream fop = new FileOutputStream(fileOut);
        saveData(fop);
        
        File fileOutt = new File("energy");
        FileOutputStream fopp = new FileOutputStream(fileOutt);
        
        int m=1;
        double emax=0;
		for(int i=0;i<Sd;i++)
		{		
			makeStep();
			calculateNEX();
			calculateG();
			
			if (emax<Energy)
			{
				emax=Energy;
			}
			
			if (i==m*Sout){
        		saveData(fop);
        		try {
        			    fopp.write(String.format("%f\t%f\t%f\t%f\n",tau,Norm,Xsr,Energy).getBytes());
    			} catch (IOException e) {
    				System.out.println("Blad zapisu fopp");
    			}	
        		m++;
        		
   
        		if (m==600){
        			for (int ii=0;ii<N;ii++)
        			{
        				foppp.write(String.format("%f\t",Gk[ii]).getBytes());
        			}
        			foppp.write(String.format("\n").getBytes());
        		}
        			
        		
        	}
		}
		
		for (int ii=0;ii<N;ii++)
	    {
	    	foppp.write(String.format("%f\t",Gk[ii]).getBytes());
	    }
	    foppp.write(String.format("\n").getBytes());
		
		fop.close();
		fopp.close();
		foppp.close();
		System.out.println(String.format("max<E>: %f", emax));
		System.out.println("Zakonczono wykonywanie programu");
		
		

	}
	
	private static void loadData(String fileName) {
        Scanner scan;
        File fileIn = new File(fileName);
        try {
            scan = new Scanner(fileIn);
            N = scan.nextInt();
            n = scan.nextInt();
            L = scan.nextDouble();
            kappa = scan.nextDouble();
            dtau = scan.nextDouble();
            omega = scan.nextDouble();
            Sd = scan.nextInt();
            Sout = scan.nextInt();
            
        } catch (Exception e1) {
            System.out.println("Invalid input file");
        }

        System.out.println(String.format("Wczytane dane:\nN=%d, n=%d, L=%.3f, kappa=%.3f, dtau=%.4f, omega=%.3f, Sd=%d, Sout=%d", N, n, L, kappa, dtau, omega, Sd, Sout));
    }
	
	private static void saveData(FileOutputStream fop) {
        try {
        	for (int i=0;i<=N;i++)
        	{
                fop.write(String.format("%f\n",Gk[i]).getBytes());
        	}
        	fop.write("\n\n".getBytes());
        } catch (IOException e) {
            System.out.println("Blad zapisu fop");
        }
    }

	
	private static void calculateNEX() {
		
		Norm=dx*Gsum;
		Xsr=dx*Gp;
		Energy=dx*Gesum;
		
	}


	private static void calculateG() {
		Gsum=0;
		Gp=0;
		Gesum=0;
		for (int i=0;i<=N;i++)
		{
			Gk[i]=Fir[i]*Fir[i]+Fii[i]*Fii[i];
			Gsum +=Gk[i];
			Gp+=arguments[i]*Gk[i];
			Gesum+=Fir[i]*Hr[i]+Fii[i]*Hi[i];
		}
		
		
	}

	private static void makeStep() {
		
		for (int i=1;i<N;i++)
		{
			Fir[i]=Fir[i]+Hi[i]*dtau/2;
		}
		calculateHr();
		for (int i=1;i<N;i++)
		{
			Fii[i]=Fii[i]-Hr[i]*dtau;
		}
		calculateHi();
		for (int i=1;i<N;i++)
		{
			Fir[i]=Fir[i]+Hi[i]*dtau/2;
		}
		tau+=dtau;
		
	}


	private static void calculateHr() {
		Hr[0]=Hr[N]=0;
		
		for( int i=1;i<N;i++)
		{
			Hr[i]= -1./2*(Fir[i+1]+Fir[i-1]-2*Fir[i])/(dx*dx)+kappa*(arguments[i]-1./2)*Fir[i]*Math.sin(omega*tau);
			
		}
		
	}
	
	private static void calculateHi() {
		Hi[0]=Hi[N]=0;
		
		for( int i=1;i<N;i++)
		{
			Hi[i]= -1./2*(Fii[i+1]+Fii[i-1]-2*Fii[i])/(dx*dx)+kappa*(arguments[i]-1./2)*Fii[i]*Math.sin(omega*tau);
			
		}
		
	}


	private static void initiateFi() {
		
		Fii[0]=Fii[N]=0;
		Fir[0]=Fir[N]=0;
		
		for(int i=1;i<N;i++)
		{
			Fii[i]=0;
			Fir[i]=Math.sqrt(2)*Math.sin(n*Math.PI*arguments[i]);
		//	System.out.println(Fir[i]);
		}
		
		
	}


	private static void calculateArguments()
	{
		arguments[0]=0;
		dx=1./N;
		for(int i=1;i<=N;i++)
		{
			arguments[i]=arguments[i-1]+dx;
			//System.out.println(arguments[i]);
		}
	}
}
