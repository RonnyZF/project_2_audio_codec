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
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

char* FFT_TO_WAV  = "Encoder.wav";
char* WAV_TO_CSV  = "Input.csv";
char* FFT_TO_CSV  = "Encoder.csv";
char* IFFT_TO_CSV = "Decoder.csv";
char* OUT = "out.wav";

//Banderas
int FIXED = 8;
//int init3 = 22;
int init2 = 22; /*Esta variable es importante, lo que hace es evitar que se 
                  guarden las variables de cabeza en el csv generado*/

int init  = 22; /*Esta variable es importante, lo que hace es similar a la 
                  init evitar que se tome en cuenta las variables de cabeza 
                  en el csv generado*/

#define pi  3.1415926535897

//***********************************************************************
//********   FUNCIONES PARA CALCULO DE TIEMPO DE PROCESAMIENTO   ********
//***********************************************************************


int frequency_of_primes (int n) {
  int i,j;
  int freq=n-1;
  for (i=2; i<=n; ++i) for (j=sqrt(i);j>1;--j) if (i%j==0) {--freq; break;}
  return freq;
}


//***********************************************************************
//**********         FUNCIONES PARA PASAR A PUNTO FIJO         **********
//***********************************************************************

float fix_to_flo(int x, int e){
    // f es el numero en punto fijo de tipo int
    // e es la cantidad de bits que se le dio la seccion decimal
    double f;
    int sign;
    int c;

    c = abs(x);
    sign = 1;
    
    if (x < 0){
    /* Las siguientes 3 lineas es para devolverlo del complemento a 2 
        si el numero original era negativo */ 

        c = x - 1; 
        c = ~c;
        sign = -1;
    }

    /*Lo que se hace es simplemente multiplicar el numero 
    en punto fijo dividido por 2 elevado a la e*/
    f = (1.0 * c)/(1<<e);
    f = f * sign;
    return f;
}

int flo_to_fix(double f, int e){
    // f es el numero original de tipo double
    // e es la cantidad de bits que se le va a dar a la seccion decimal

    /*Lo que se hace es simplemente multiplicar el numero 
    original multiplicado por 2 elevado a la e*/
    double a = f*(1<<e);

    if (a >= 127){
        a = 127;
    }

    if (a <= -127){
        a = -127;
    }

    int b = (int)(round(a));
    
    if (a < 0){
        /* Las siguientes 3 lineas es para pasarlo a complemento a 2 
        si el numero original es negativo */ 

        b = abs(b);
        b = ~b;
        b = b + 1;
    }
    return b;
}


//***********************************************************************
//******     FUNCIONES PARA LECTURA/ESCRITURA ARCHIVOS WAV/CSV      *****
//***********************************************************************

void to_float(short* buff_short, float* buff_float, int len){
    int i;
    for (i = init; i < len; i++){

        buff_float[i] = fix_to_flo(buff_short[i],15);
        //printf("Num %d: 0x%04x -- %f\n",i,buff_short[i],buff_float[i]);
    }
    init = 0; /*Luego de la primera corrida, se ha pasado el encabezado
                 y se baja la bandera*/
}

void to_complex(char* buff_short, double complex* buff_complex, int len){

    buff_complex[0]   = fix_to_flo(buff_short[0]  ,FIXED) + fix_to_flo(buff_short[1],FIXED)*I;    
    buff_complex[(len-1)/2] = fix_to_flo(buff_short[len-2],FIXED) + fix_to_flo(buff_short[len-1],FIXED)*I;        
    
    for (int i = 2; i < len-2; i = i+2){
        buff_complex[i/2]       = fix_to_flo(buff_short[i],FIXED) + fix_to_flo(buff_short[i + 1],FIXED)*I;
        buff_complex[len-i/2-2] = fix_to_flo(buff_short[i],FIXED) - fix_to_flo(buff_short[i + 1],FIXED)*I;
    }

    /*
    for (int j = 0; j < len - 2; j++){
        printf("Num %d: %f + %fi\n",j,creal(buff_complex[j]), cimag(buff_complex[j]));
    }*/
}

void save_csv(double complex* buff_float, int len){
    FILE * file;
    file = fopen (IFFT_TO_CSV, "a");

    char* progbar="-/|\\";
    int progidx=0;

    for (int m=0; m<len; m++){
        fprintf(file, "%f\n", creal(buff_float[m]));

        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);
    }

    fclose(file); 
}


void save_wav_csv(float* buff_float, int len){
    FILE * file;
    file = fopen (WAV_TO_CSV, "a");

    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len; i++){
        fprintf(file, "%f\n", buff_float[i]);

        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);
    }

    fclose(file); 
}

void save_fft_csv(double complex* buff_complex, int len){
    FILE * file;
    file = fopen (FFT_TO_CSV, "a");
    
    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len/2+1; i++){
        fprintf(file, "%.7f, %.7f\n", creal(buff_complex[i]),cimag(buff_complex[i]));
        
        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);
    }

    fclose(file); 
}

void save_fft(double complex* buff_complex, int len){
    FILE * file;
    file = fopen (FFT_TO_WAV, "a");

    char re;
    char im;
    
    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len/2+1;i++){
        re = (char)flo_to_fix(creal(buff_complex[i]),FIXED);
        im = (char)flo_to_fix(cimag(buff_complex[i]),FIXED);
        fwrite(&re, sizeof(char), 1, file);
        fwrite(&im, sizeof(char), 1, file);

        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);

        //printf("Num %d: 0x%04x = %.5f = %.5f",i,re,fix_to_flo(re,FIXED),creal(buff_complex[i]));
        //printf("  ----  0x%04x = %.5f = %.5f\n",im,fix_to_flo(im,FIXED),cimag(buff_complex[i]));
    }

    //printf("Escritura... ...\n");

    //for (int i=0; i<len; i++){
    //    printf("Num %d: %.5f ",i,creal(buff_complex[i]));
    //    printf("+ %.5fi\n",cimag(buff_complex[i]));
    //}
    fclose(file); 
}


void save_wav(double complex *buff_float, int len){
    FILE * file;
    file = fopen (OUT, "a");
    
    short var;

    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len;i++){
        var = (short)flo_to_fix(creal(buff_float[i]),15);
        fwrite(&var, sizeof(short), 1, file);

        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);
    }

    fclose(file); 
}


//***********************************************************************
//**********                 FUNCIONES PARA FFT                **********
//***********************************************************************

void fft(double complex *X, unsigned short EXP,double complex *W, unsigned short SCALE){
    double temp_re;         /* Temporary storage of complex variable */
    double temp_im;

    double U_re;            /* Twiddle factor W^k */
    double U_im;

    unsigned short i,j;
    unsigned short id;      /* Index for lower point in butterfly */
    unsigned short N=1<<EXP;/* Number of points for FFT */
    unsigned short L;       /* FFT stage */
    
    unsigned short LE;      /* Number of points in sub DFT at stage L and offset 
                               to next DFT in stage */
    
    unsigned short LE1;     /* Number of butterflies in one DFT at stage L. 
                               Also is offset to lower point in butterfly at stage L */
    float scale;
    scale = 0.5;
    if (SCALE == 0){
        scale = 1.0;
    }

    /* FFT butterfly */
    for (L=1; L<=EXP; L++){
        LE=1<<L;            /* LE=2^L=points of sub DFT */
        LE1=LE>>1;          /* Number of butterflies in sub-DFT */
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

void ifft(double complex *X, unsigned short EXP,double complex *W, unsigned short SCALE){

    unsigned short N  = 1<<EXP; /* Number of points for FFT */
    for(int i=0;i<N;i++){
        //printf("Num %d: %f + %f i\n",i,creal(X[i]),cimag(X[i]));
        
        X[i] = conj(X[i])*N;
        
        //printf("Num %d: %f + %f i\n\n",i,creal(X[i]),cimag(X[i]));
    }

    
    fft(X, EXP, W, SCALE);
    /*for(int i=0;i<N;i++){
        printf("Num %d: %f + %f i\n",i,creal(X[i]),cimag(X[i]));
    }*/
}

void fft_init (double complex *W, unsigned short EXP){
    unsigned short L,LE,LE1;
    for(L=1; L<=EXP; L++){
        LE  = 1<<L;  // LE=2^ points of sub DFT
        LE1 = LE>>1; // Number of butterflies in sub DFT
        W[L-1] = cos(pi/LE1) - sin(pi/LE1) * I;
    }
}

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



//***********************************************************************
//**********                  INICIO DEL  MAIN                 **********
//***********************************************************************

int main(int argc, char *argv[]) {
    
    //**********************************************************************************
    // VARIABLES

    int Muestras = 64;
    int MuestrasFFT = Muestras + 2;
    int Num = 64;
    float *buffer_Float;
    short *buffer_Short;
    double complex *W;
    double complex *FFT;
    int EXP = 6;
    char *buffer_ShortFFT;
    double complex *buffer_FloatFFT;

    buffer_ShortFFT = malloc((int)MuestrasFFT   * sizeof(char));
    buffer_FloatFFT = malloc((int)Num * sizeof(double complex));

    buffer_Float    = malloc((int)Muestras * sizeof(float));
    buffer_Short    = malloc((int)Muestras * sizeof(short));
    FFT             = malloc((int)Num * sizeof(double complex));
    W               = malloc((int)EXP * sizeof(double complex));



    //Inicializa contador y clock para ver tiempo de ejecucion.
    clock_t t;
    int f;
    t = clock();
    f = frequency_of_primes (99999);
    printf ("Inicio de reloj...\n");


    //**********************************************************************************
    // ADQUISICION DE ARGUMENTOS

    // structure for the long options. 
    static struct option lopts[] = {
        {"verbose"  ,no_argument,0,'v'},
        {"help"     ,no_argument,0,'h'},
        {0,0,0,0}
    };

    int optionIdx,c;
    char * encoder ="No data";
    while ((c = getopt_long(argc, argv, "cdh",lopts,&optionIdx)) != -1) {
        switch (c) {

        case 'c':
            encoder = "codificar";
        break;
        
        case 'd':
            encoder = "decodifica";
        break;
        
        case 'h':
            printf("Usage: %s [-c] [-d] [-h|--help] archivo.wav\n\n",argv[0]);
            printf("          -c -> Codifica el archivo wav que se entrega como argumento\n");
            printf("          -d -> Decodifica el archivo wav que se entrega como argumento\n");
            printf(" archivo.wav -> Es el documento wav que se desea codificar/decodifica\n\n");
            exit(1);
        
        default:
            printf("Usage: %s [-c] [-d] [-h|--help] archivo.wav\n\n",argv[0]);
            printf("          -c -> Codifica el archivo wav que se entrega como argumento\n");
            printf("          -d -> Decodifica el archivo wav que se entrega como argumento\n");
            printf(" archivo.wav -> Es el documento wav que se desea codificar/decodifica\n\n");
            exit(1);
        }
    }

    if (encoder == "No data")
    {
        printf("Usage: %s [-c] [-d] [-h|--help] archivo.wav\n",argv[0]);
        printf("          -c -> Codifica el archivo wav que se entrega como argumento\n");
        printf("          -d -> Decodifica el archivo wav que se entrega como argumento\n");
        printf(" archivo.wav -> Es el documento wav que se desea codificar/decodifica\n\n");
        exit(1);
    }

    char * Lectura = argv[2];
    printf("Se desea %s el documento -> %s\n\n",encoder,Lectura);


    unsigned short L,LE,LE1;
    for(L=1; L<=EXP; L++){
        LE  = 1<<L;  // LE=2^ points of sub DFT
        LE1 = LE>>1; // Number of butterflies in sub DFT
        W[L-1] = cos(pi/LE1) - sin(pi/LE1) * I;
    }

    //Codificador
    if(encoder == "codificar"){
        
        //**********************************************************************************
        // LIMPIEZA Y CREACION DE ARCHIVOS

        FILE * file;
        file = fopen (WAV_TO_CSV, "w");
        fclose(file);

        FILE * fileFFT;
        fileFFT = fopen (FFT_TO_WAV, "w");
        fclose(fileFFT);

        FILE *fp = fopen(Lectura, "r");
        printf("ENCODER\n");  

        //fft_init(W, EXP);//Se inicializa los coeficientes de la FFT
        //Para probar la FFT se comenta el loop y solo se hace pasar las 16 muestras mas los 22 muestras de cabeza
        int count = 0;
        if (fp != NULL) {
            size_t byte_read;
            printf("Escribiendo:\n");
            do{
                /*Cantidad de muestras leidas = fread(Destino de lectura,
                                                      tamaño de bit,
                                                      cantidad de muestras maxima,
                                                      archivo)*/

                byte_read = fread(buffer_Short,sizeof(short),Muestras,fp);
                to_float(buffer_Short, buffer_Float, byte_read);
                save_wav_csv(buffer_Float,byte_read);

                for (int i = 0; i < Num; i++){
                    FFT[i] = buffer_Float[i];
                }

                bit_rev(FFT,EXP);
                fft(FFT, EXP, W, 1);
                save_fft(FFT, Num);
                //save_fft_csv(FFT, Num);
                count++;
            }while(byte_read > 0);
            printf("Se ha escrito la FFT con el wav codificado!, se llama %s\n\n\n",FFT_TO_WAV);
            //Mietras se siga leyendo datos se seguira ejecutando el while
            
            fclose(fp); //Se cierra el archivo cuando se termina de leer
        }
    }


    //Decodificador

    //**********************************************************************************
    //// Prueba de lectura de archivos fft

    else{
        printf("DECODER\n");

        //**********************************************************************************
        // LIMPIEZA Y CREACION DE ARCHIVOS

        FILE * fileIFFT;
        fileIFFT = fopen (IFFT_TO_CSV, "w");
        fclose(fileIFFT);

        FILE * fileFFT_csv;
        fileFFT_csv = fopen (FFT_TO_CSV, "w");
        fclose(fileFFT_csv);

        FILE * fileOUT;
        fileOUT = fopen (OUT, "w");
        fclose(fileOUT);

        FILE *fpFFT = fopen(Lectura, "r");
        //fft_init(W, EXP);//Se inicializa los coeficientes de la FFT
        //Para probar la FFT se comenta el loop y solo se hace pasar las 16 muestras mas los 22 muestras de cabeza
        if (fpFFT != NULL) {
            size_t byte_readFFT;
            printf("Escribiendo:\n");
            do{
                //Cantidad de muestras leidas = fread(Destino de lectura,
                //                                      tamaño de bit,
                //                                      cantidad de muestras maxima,
                //                                      archivo)

                byte_readFFT = fread(buffer_ShortFFT,sizeof(char),MuestrasFFT,fpFFT);
                to_complex(buffer_ShortFFT, buffer_FloatFFT, byte_readFFT);
                if (byte_readFFT > 0){
                    save_fft_csv(buffer_FloatFFT, Num);
                    bit_rev(buffer_FloatFFT,EXP);
                    ifft(buffer_FloatFFT, EXP, W, 1);
                    save_csv(buffer_FloatFFT,Num);
                    save_wav(buffer_FloatFFT, Num);
                }

            }while(byte_readFFT > 0); 
            //Mietras se siga leyendo datos se seguira ejecutando el while

            printf("Se ha escrito el csv con la FFT decodificada!, se llama %s\n\n\n",IFFT_TO_CSV);
            fclose(fpFFT); //Se cierra el archivo cuando se termina de leer
        }
    }


    //Imprime el tiempo de ejecucion
    t = clock() - t;
    printf ("Le tardo al programa %ld clicks completar el proceso (%f segundos).\n\n\n",t,((float)t)/CLOCKS_PER_SEC);
  

    return 0;
}