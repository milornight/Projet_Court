#include "mbed.h"
#include "Grove_LCD_RGB_Backlight.h"
#include <math.h> 
#include <cstdio>
#include <cmath>
#include <complex>
#include <iostream>
#include <valarray>


typedef std::complex<float> Complex;		//definie un type complex qui contient que les nombre float
typedef std::valarray<Complex> CArray;

void spectre(Complex tmp[], int N);

AnalogIn sound(A0);			 //Le son entré du l'entré A0
Grove_LCD_RGB_Backlight i2c(I2C_SDA , I2C_SCL );	// déclaration de l'objet I2C
Serial pc(USBTX, USBRX);	// tx,rx

char str[15];
char str2[15];
char str3[15];
char str4[15];
const float PI = 3.1415926535;

Ticker tick;


float Log2( float num)  
{   
    return log(num)/log(2.0);  
} 

 
void fft(CArray &x)
{
    // TFD
    unsigned int N = x.size(), k = N, n;
    float thetaT = PI / N;
    Complex Wn = Complex(cos(thetaT), -sin(thetaT)), T;  // initier la facteur Wn
	
    while (k > 1)
    {
        n = k;
        k >>= 1;
        Wn *= Wn;
        T = 1.0L;
			
        for (unsigned int l = 0; l < k; l++)
        {
            for (unsigned int a = l; a < N; a += n)
            {
                unsigned int b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= Wn;
        }
    }
		
    //Decimer
    unsigned int m = (unsigned int) Log2(N);
		
    for (unsigned int a = 0; a < N; a++)
    {
        unsigned int b = a;
			
        //Inverser bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
			
        if (b > a)
        {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
}


const int N = 512; //nombre d'échantillonage
const int Fe = 40000; //Frequnce d'echantillonage=2*Bande
Complex tmp[N];
int cpt = 0;
float fs=Fe/N;		// frequence du signal


void liresignal(){
    
    double audio = sound.read();
    tmp[cpt] = audio;
    cpt++;
		//printf("\r\n Re:%f Im:%f",tmp[cpt].real(),tmp[cpt].imag());
	
    if (cpt == N) 
		{
        cpt = 0;
        spectre(tmp,N);
    }
    
}

void spectre(Complex tmp[], int N) 
{
    CArray aux(tmp, N);	// malloc un tableau de taille N
    float tabfreq[N/2]; //vecteur frequence
		float tabPSD[N/2];  //puissance spectrale
		int cpt2 = 0;
    
		fft(aux);
   
	
    for (int i = 1; i <N/2 ; i++)
		{
        aux[i] = (aux[i].imag())*(aux[i].imag())+(aux[i].real())*(aux[i].real()); //norme de vecteur de frequence
				//aux[i+1] = (aux[i+1].imag())*(aux[i+1].imag())+(aux[i+1].real())*(aux[i+1].real());
        tabfreq[i] =((float)i)*(float)Fe/((float)(N)) ;														// la valeur de fréquence pour chaque point
				tabPSD[i]=(aux[i].imag()*aux[i].imag()+aux[i].real()*aux[i].real())/N;		// la puissance pour chaque point
				//tabPSD[i+1]=(aux[i+1].imag()*aux[i+1].imag()+aux[i+1].real()*aux[i+1].real())/N;		
				//printf("\r\n Freq:%.fHz Re:%f Im:%f",tabfreq[i],aux[i].real(),aux[i].imag());
			
        if (aux[i].real() > aux[cpt2].real())
				{
            cpt2 = i;
        }
    }
		
		// Affichage
    wait(0.3);        
		
    int freq = tabfreq[cpt2]; 
		int freq2 =tabfreq[cpt2+1];
		int freq3 =fs*cpt2; //freq3=freq
		int PSD=tabPSD[cpt2];
		int PSD2=tabPSD[cpt2+1];
		
    sprintf(str, "F:%dHz ", freq);
		sprintf(str2, "F2:%dHz ", freq2);
		sprintf(str3, "P:%d", PSD);
		sprintf(str4, "P2:%d", PSD2);
		//printf("%dHz",freq3);
		
		i2c.clear();			// initier l'écran
		i2c.print(str);
		i2c.print(str3);
		i2c.locate(0,1);	// position en 2er ligne 1er colonne
		i2c.print(str2);
		i2c.print(str4);
}



int main() 
{    
    i2c.clear();  
    i2c.setRGB(0,0,255);   // mettre en couleur bleu clair
    tick.attach(&liresignal, 1.0/Fe);
    
    while(1)
		{
			
    }
}

