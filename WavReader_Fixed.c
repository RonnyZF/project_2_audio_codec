/*
Se compila en la maquina virtual con el comandos:
  gcc -o Fix WavReader_Fixed.c -lm
  ./Fix
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

typedef struct{
    int64_t re;
    int64_t im;
}Complex;

char* FFT_TO_WAV  = "Encoder.wav";
char* WAV_TO_CSV  = "Input.csv";
char* FFT_TO_CSV  = "Encoder.csv";
char* IFFT_TO_CSV = "Decoder.csv";
char* IFFT_TO_WAV = "out.wav";

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

void short_to_complex(short* buff_short, Complex* buff_complex, int len){
    int i;
    for (i = init; i < len; i++){

        buff_complex[i].re = (int)buff_short[i];
        buff_complex[i].im = 0;
        //printf("Num %d: %d + %di\n",i,buff_complex[i].re, buff_complex[i].im);
    }
    init = 0; /*Luego de la primera corrida, se ha pasado el encabezado
                 y se baja la bandera*/
}


void char_to_complex(char* buff_char, Complex* buff_complex, int len){

    double puente_re[len - 2];
    double puente_im[len - 2];

    puente_re[0] = fix_to_flo(buff_char[0],FIXED);
    puente_im[0] = fix_to_flo(buff_char[1],FIXED);

    puente_re[(len-1)/2] = fix_to_flo(buff_char[len-2],FIXED);
    puente_im[(len-1)/2] = fix_to_flo(buff_char[len-1],FIXED);

    for (int i = 2; i < len-2; i = i+2){
        puente_re[i/2] = fix_to_flo(buff_char[i],FIXED);
        puente_im[i/2] = fix_to_flo(buff_char[i + 1],FIXED);

        puente_re[len-i/2-2] =   fix_to_flo(buff_char[i],FIXED);
        puente_im[len-i/2-2] = - fix_to_flo(buff_char[i + 1],FIXED);
    }

    for (int i = 0; i < len-2; i++){
        buff_complex[i].re = flo_to_fix(puente_re[i],15);
        buff_complex[i].im = flo_to_fix(puente_im[i],15);
    }
    /*
    for (int j = 0; j < len - 2; j++){
        printf("Num %d: %f + %fi\n",j,creal(buff_complex[j]), cimag(buff_complex[j]));
    }*/
}

void save_csv(Complex* buff_complex, int len){
    FILE * file;
    file = fopen (IFFT_TO_CSV, "a");

    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len; i++){
        fprintf(file, "%f\n", (double)buff_complex[i].re/(1<<15));
        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);
    }

    fclose(file); 
}


void save_wav(Complex *buff_complex, int len){
    FILE * file;
    file = fopen (IFFT_TO_WAV, "a");
    
    short var;
    double puente;

    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len;i++){
        
        //Por razones que escapan a mi entendimiento, se debe hacer el paso a flotante 
        //y luego tirarlo de vuelta a punto fijo, esto es un misterio, pero funciona
        puente = fix_to_flo(buff_complex[i].re,15);
        var    = (short)flo_to_fix(puente,15);;
        fwrite(&var, sizeof(short), 1, file);
    }

    fclose(file); 
}


void save_fft(Complex* buff_complex, int len){
    FILE * file;
    file = fopen (FFT_TO_WAV, "a");

    char re;
    char im;
    double puente_re;
    double puente_im;
    
    char* progbar="-/|\\";
    int progidx=0;

    for (int i=0; i<len/2+1;i++){

        puente_re = fix_to_flo(buff_complex[i].re,15); 
        puente_im = fix_to_flo(buff_complex[i].im,15);

        if (puente_re > 0.5){
            puente_re = 0.5;
        }
        else if (puente_re < -0.5){
            puente_re = -0.5;
        }

        if (puente_im > 0.5){
            puente_im = 0.5;
        }
        else if (puente_im < -0.5){
            puente_im = -0.5;
        }

        re = (char)flo_to_fix(puente_re,FIXED);
        im = (char)flo_to_fix(puente_im,FIXED);
        fwrite(&re, sizeof(char), 1, file);
        fwrite(&im, sizeof(char), 1, file);

        printf("%c\r",progbar[progidx++&3]);
        fflush(stdout);

        //printf("0x%04x = 0x%08lx -- %.5f = %.5f\n"  ,re,buff_complex[i].re,(double)re/(1<<FIXED),fix_to_flo(buff_complex[i].re,15));
        //printf("0x%04x = 0x%08lx -- %.5f = %.5f\n\n",im,buff_complex[i].im,(double)im/(1<<FIXED),fix_to_flo(buff_complex[i].im,15));
    }

    fclose(file); 
}


//***********************************************************************
//**********                 FUNCIONES PARA FFT                **********
//***********************************************************************

void fft(Complex *X, unsigned short EXP,Complex *W, unsigned short SCALE,int Fix){
    int64_t temp_re;         /* Temporary storage of complex variable */
    int64_t temp_im;

    int64_t U_re;            /* Twiddle factor W^k */
    int64_t U_im;

    unsigned short i,j;
    unsigned short id;      /* Index for lower point in butterfly */
    unsigned short N=1<<EXP;/* Number of points for FFT */
    unsigned short L;       /* FFT stage */
    
    unsigned short LE;      /* Number of points in sub DFT at stage L and offset 
                               to next DFT in stage */
    
    unsigned short LE1;     /* Number of butterflies in one DFT at stage L. 
                               Also is offset to lower point in butterfly at stage L */
    int64_t scale;
    scale = 0.5*(1<<Fix);
    if (SCALE == 0){
        scale = (1<<Fix);
    }

    /* FFT butterfly */
    for (L=1; L<=EXP; L++){
        LE  = 1<<L;            /* LE=2^L=points of sub DFT */
        LE1 = LE>>1;          /* Number of butterflies in sub-DFT */
        
        U_re = 1.0*(1<<Fix);
        U_im = 0;

        for (j=0; j<LE1;j++){
            /* Do the butterflies */
            //IMPORTANTE
            for(i=j; i<N; i+=LE) {
                id=i+LE1;
                temp_re = X[id].re*U_re/(1<<Fix);
                temp_re = temp_re - X[id].im*U_im/(1<<Fix);
                temp_re = temp_re*scale/(1<<Fix);

                temp_im = X[id].im*U_re/(1<<Fix);
                temp_im = temp_im + X[id].re*U_im/(1<<Fix);
                temp_im = temp_im*scale/(1<<Fix);

                X[id].re = X[i].re*scale/(1<<Fix) - temp_re;
                //printf("%d\n",X[id].re);
                //printf("%f\n\n",(double)X[id].re/(1<<Fix));

                X[id].im = X[i].im*scale/(1<<Fix) - temp_im;
                X[i].re  = X[i].re*scale/(1<<Fix) + temp_re;
                X[i].im  = X[i].im*scale/(1<<Fix) + temp_im;
            }

            /* Recursive compute W^k as U*W^(k-1) */
            temp_re = (U_re*W[L-1].re/(1<<Fix) - U_im*W[L-1].im/(1<<Fix));
            U_im    = (U_re*W[L-1].im/(1<<Fix) + U_im*W[L-1].re/(1<<Fix));
            U_re    = temp_re;
        }
    }
}

void ifft(Complex *X, unsigned short EXP,Complex *W, unsigned short SCALE, int Fix){

    unsigned short N  = (1<<EXP); /* Number of points for FFT */
    
    for(int i=0;i<(1<<EXP);i++){        

        X[i].re =  X[i].re*N;
        X[i].im = -X[i].im*N;
    }

    fft(X, EXP, W, SCALE,Fix);
/*
    for(int i=0;i<N;i++){
        //printf("Num %d: %ld + %ld i\n",i,X[i].re,X[i].im);
        printf("Num %d: %f + %f i\n",i,(double)X[i].re/(1<<Fix),(double)X[i].im/(1<<Fix));
    }
*/
}

void fft_init (Complex *W, unsigned short EXP, int Fix){

    unsigned short L,LE,LE1;
    double re;
    double im;

    for(L=1; L<=EXP; L++){
        LE  = 1<<L;         // LE=2^ points of sub DFT
        LE1 = LE>>1;        // Number of butterflies in sub DFT
        re = cos(pi/LE1);
        im = - sin(pi/LE1);
        
        //Se pasan los coeficietes a punto fijo;
        W[L-1].re = re*(1<<Fix);
        W[L-1].im = im*(1<<Fix);
        
        //printf("Coeficiete = %d: %f + %f i\n",L,re,im);
        //printf("Coeficiete = %d: %d + %d i\n",L,W[L-1].re,W[L-1].im);
    }
}

void bit_rev (Complex *X, int EXP){

    unsigned short N  = 1<<EXP;
    unsigned short N2 = N>>1;
    Complex temp;
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
    Complex *buffer_Complex;
    short *buffer_Short;
    Complex *W;
    Complex *FFT;
    int EXP = 6;
    char *buffer_ShortFFT;
    Complex *buffer_ComplexFFT;

    buffer_ShortFFT = malloc((int)MuestrasFFT * sizeof(char));
    buffer_ComplexFFT = malloc((int)Muestras    * sizeof(Complex));

    buffer_Complex    = malloc((int)Muestras * sizeof(Complex));
    buffer_Short    = malloc((int)Muestras * sizeof(short));
    FFT             = malloc((int)Num * sizeof(Complex));
    W               = malloc((int)EXP * sizeof(Complex));


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


    //Codificador
    if(encoder == "codificar"){
        
        //**********************************************************************************
        // LIMPIEZA Y CREACION DE ARCHIVOS

        FILE * fileFFT;
        fileFFT = fopen (FFT_TO_WAV, "w");
        fclose(fileFFT);

        FILE *fp = fopen(Lectura, "r"); //Se abre el wav a codificar
        fft_init(W, EXP, 15);//Se inicializa los coeficientes de la FFT

        printf("ENCODER\n");
        //Lectura del wav original para ser codificado
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
                short_to_complex(buffer_Short, buffer_Complex, byte_read);

                /*
                for (int i = 0; i < Num; i++){
                    FFT[i].re = buffer_Complex[i].re;
                    FFT[i].im = buffer_Complex[i].im;
                }*/

                bit_rev(buffer_Complex,EXP);
                fft(buffer_Complex, EXP, W, 1,15);

                save_fft(buffer_Complex, Num);
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

        FILE * fileOUT;
        fileOUT = fopen (IFFT_TO_WAV, "w");
        fclose(fileOUT);

        FILE *fpFFT = fopen(Lectura, "r");
        fft_init(W, EXP, 15);//Se inicializa los coeficientes de la FFT

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
                char_to_complex(buffer_ShortFFT, buffer_ComplexFFT, byte_readFFT);
                if (byte_readFFT > 0){
                    //save_fft_csv(buffer_ComplexFFT, Num);
                    
                    bit_rev(buffer_ComplexFFT,EXP);
                    ifft(buffer_ComplexFFT, EXP, W, 1,15);
                    save_csv(buffer_ComplexFFT,Num);
                    save_wav(buffer_ComplexFFT, Num);
                }

            }while(byte_readFFT > 0); 
            //Mietras se siga leyendo datos se seguira ejecutando el while

            printf("Se ha escrito el csv con la FFT decodificada!, se llama %s\n\n\n",IFFT_TO_CSV);
            fclose(fpFFT); //Se cierra el archivo cuando se termina de leer
        }
    }

    return 0;
}