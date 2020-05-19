/*
Se compila en la maquina virtual con el comandos:
  gcc -o WavReader WavReader.c -lm
  ./complex
  
*/

#include <stdio.h>      /* Standard Library of Input and Output */
#include <complex.h>    /* Standard Library of Complex Numbers */
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int init  = 22; /*Esta variable es importante, lo que hace es evitar que se 
                  guarden las variables de cabeza en el csv generado*/

int init2 = 22; /*Esta variable es importante, lo que hace es similar a la 
                  init evitar que se tome en cuenta las variables de cabeza 
                  en el csv generado*/

#define pi  3.1415926535897

//***********************************************************************
//**********            FUNCIONES PARA LECTURA DE WAV          **********
//***********************************************************************

void to_float(short* buff_short, float* buff_float, int len){
    int i;
    for (i = init2; i < len; i++){
        if(i%16 == 0){
            printf("\n");
        }
        buff_float[i] = (float)(buff_short[i])/32768;
        printf("Num %d: 0x%02x -- %f\n",i,buff_short[i],buff_float[i]);
    }
    init2 = 0; /*Luego de la primera corrida, se ha pasado el encabezado
                 y se baja la bandera*/
}

void save_csv(float* buff_float, int len){
    FILE * file;
    file = fopen ("Data_OUT.csv", "a");

    for (int m=init; m<len; m++){
        fprintf(file, "%f\n", buff_float[m]);
    }
    init = 0; /*Luego de la primera corrida, se ha pasado el encabezado
                y se baja la bandera*/
    fclose(file); 
}

//***********************************************************************
//**********                 FUNCIONES PARA FFT                **********
//***********************************************************************

void bit_rev (double complex *X, int EXP){

    unsigned short N  = 1<<EXP;
    unsigned short N2 = N>>1;
    double complex temp;
    unsigned short j,i,k;
    j = 0;

    for(i=1;i<N-1;i++){
        k = N2;
        while(k <= j){
            j-=k;
            k>>=1;
        }

        j+=k;
        if(i<j){
            temp = X[j];
            X[j] = X[i];
            X[i] = temp;
        }
    }
}


void fft(double complex *X, unsigned short EXP,double complex *W, unsigned short SCALE){
    double temp_re; /* Temporary storage of complex variable */
    double temp_im;

    double U_re; /* Twiddle factor W^k */
    double U_im;

    unsigned short i,j;
    unsigned short id; /* Index for lower point in butterfly */
    unsigned short N=1<<EXP;/* Number of points for FFT */
    unsigned short L; /* FFT stage */
    unsigned short LE; /* Number of points in sub DFT at stage L
    and offset to next DFT in stage */
    unsigned short LE1; /* Number of butterflies in one DFT at
    stage L. Also is offset to lower point
    in butterfly at stage L */
    float scale;
    scale = 0.5;
    if (SCALE == 0){
        scale = 1.0;
    }

    /* FFT butterfly */
    for (L=1; L<=EXP; L++){
        LE=1<<L; /* LE=2^L=points of sub DFT */
        LE1=LE>>1; /* Number of butterflies in sub-DFT */
        U_re = 1.0;
        U_im = 0.;
        for (j=0; j<LE1;j++){
            /* Do the butterflies */
            for(i=j; i<N; i+=LE) {
                id=i+LE1;
                temp_re = (creal(X[id])*U_re - cimag(X[id])*U_im)*scale;
                temp_im = (cimag(X[id])*U_re + creal(X[id])*U_im)*scale;

                X[id] = (creal(X[i])*scale - temp_re) + (cimag(X[i])*scale - temp_im) * I;
                X[i]  = (creal(X[i])*scale + temp_re) + (cimag(X[i])*scale + temp_im) * I;
            }
            /* Recursive compute W^k as U*W^(k-1) */
            temp_re = U_re*creal(W[L-1]) - U_im*cimag(W[L-1]);
            U_im    = U_re*cimag(W[L-1]) + U_im*creal(W[L-1]);
            U_re    = temp_re;
        }
    }
}

void fft_init (double complex *W, unsigned short EXP){
    unsigned short L,LE,LE1;
    for(L=1; L<=EXP; L++){
        LE  = 1<<L;  // LE=2^ points of sub DFT
        LE1 = LE>>1; // Number of butterflies in sub DFT
        W[L-1] = cos(pi/LE1) - sin(pi/LE1) * I;
    }
}


int main(int argc, char *argv[]) {

    //**********************************************************************************
    ///// Prueba numeros complejos
    double complex z1 = 1.0 + 3.0 * I;
    double complex z2 = 1.0 - 4.0 * I;

    printf("Prueba numeros complejos:\n");
    printf("Valores complejos: z1 = %.2f + %.2fi\t z2 = %.2f %+.2fi\n\n\n", creal(z1), cimag(z1), creal(z2), cimag(z2));

    
    //**********************************************************************************
    //// Prueba de funcion bit_rev
    printf("Prueba de inversion de entrada:\n");

    double complex *Set;
    short Num = 16;
    short EXP = 4;
    Set = malloc((int)Num * sizeof(double complex));

    Set[0]  = 0;
    Set[1]  = 1;
    Set[2]  = 2;
    Set[3]  = 3;
    Set[4]  = 4;
    Set[5]  = 5;
    Set[6]  = 6;
    Set[7]  = 7;
    Set[8]  = 8;
    Set[9]  = 9;
    Set[10] = 10;
    Set[11] = 11;
    Set[12] = 12;
    Set[13] = 13;
    Set[14] = 14;
    Set[15] = 15;
    
    bit_rev(Set,EXP);
    
    for (int n=0; n<Num; n++){
        printf("%.2f %+.2fi\n",creal(Set[n]),cimag(Set[n]));
    }
    printf("\n\n");


    //**********************************************************************************
    //// Prueba de funcion fft_init 
    printf("Prueba de generacion de matriz\n");
    EXP = 4;
    double complex *W;
    W = malloc((int)EXP * sizeof(double complex));
    
    fft_init(W, EXP);

    for (int n=0; n<EXP; n++){
        printf("%d: %.2f %+.2fi\n",n,creal(W[n]),cimag(W[n]));
    }
    printf("\n\n");

    
    //**********************************************************************************
    //// Prueba de lectura de archivos wav
    FILE * file;
    file = fopen ("Data_OUT.csv", "w");
    fclose(file);

    int Muestras = 16 + 22;
    float buffer_Float[Muestras];
    short buffer_Short[Muestras];

    FILE *fp = fopen("test.wav", "r");
    
    //Para probar la FFT se comenta el loop y solo se hace pasar las 16 muestras mas los 22 muestras de cabeza
    //if (fp != NULL) {
        size_t byte_read;
    //    do{
            /*Cantidad de muestras leidas = fread(Destino de lectura, 
                                                  tamaño de bit, 
                                                  cantidad de muestras maxima,
                                                  archivo)*/

            byte_read = fread(buffer_Short,sizeof(short),Muestras,fp);
            to_float(buffer_Short, buffer_Float, byte_read);
            save_csv(buffer_Float,byte_read);
    //    }while(byte_read > 0); 
        //Mietras se siga leyendo datos se seguira ejecutando el while

        printf("Se ha escrito el csv!, se llama Data_OUT.csv\n\n\n");
        fclose(fp); //Se cierra el archivo cuando se termina de leer
    //}

    /*
    printf("Num 0 : 0x%02x -- %f\n",buffer_Short[0],buffer_Float[0]);
    printf("Num 5 : 0x%02x -- %f\n",buffer_Short[5],buffer_Float[5]);
    printf("Num 10: 0x%02x -- %f\n",buffer_Short[10],buffer_Float[10]);
    printf("Num 20: 0x%02x -- %f\n",buffer_Short[20],buffer_Float[20]);
    */

    //**********************************************************************************
    //// Prueba de funcion fft 
    printf("Prueba de FFT\n");

    EXP = 4;
    Num = 16;
    double complex *W_FFT;
    double complex *FFT;

    W_FFT = malloc((int)EXP * sizeof(double complex));
    FFT   = malloc((int)Num * sizeof(double complex));

    for (int i = 22; i < (22+Num); i++){
        FFT[i-22] = buffer_Float[i];
    }

    bit_rev(FFT,EXP);
    fft_init(W, EXP);
    fft(FFT, EXP, W, 0);
    
    for (int f=0; f<Num; f++){
        printf("Num %d: %.4f %.4fi\n",f,creal(FFT[f]),cimag(FFT[f]));
    }
    printf("\n\n");

    return 0;
}