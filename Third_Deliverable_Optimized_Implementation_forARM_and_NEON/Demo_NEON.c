/** @mainpage test_NEON - None
 *
 * @author Ronny Zarate <anony@mo.us>
 * @version 1.0.0
**/
#include <stdio.h>
#include <arm_neon.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#ifndef NULL
#define NULL	((void *) 0)
#endif


int32x4_t fix_to_flo_Forced_NEON(int32x4_t x,int e){

	int32x4_t comp = vdupq_n_s32(127);
	int32x4_t c = vabsq_s32(x);
	//printf ("c= %d, %d, %d, %d\n", vgetq_lane_s32(c,0), vgetq_lane_s32(c,1), vgetq_lane_s32(c,2), vgetq_lane_s32(c,3));
	int32x4_t sign = vdupq_n_s32(1);

	uint32x4_t aux = vcgtq_s32(x,comp);
	int32x4_t mask = (int32x4_t)aux;

	//int32x4_t mask = vreinterpret_s32_u32(aux);

	int32x4_t mask_neg = vmvnq_s32(mask);
	//printf ("mask %X, %X, %X, %X\n", vgetq_lane_s32(mask,0), vgetq_lane_s32(mask,1), vgetq_lane_s32(mask,2), vgetq_lane_s32(mask,3));
	//printf ("mask_N %X, %X, %X, %X\n", vgetq_lane_s32(mask_neg,0), vgetq_lane_s32(mask_neg,1), vgetq_lane_s32(mask_neg,2), vgetq_lane_s32(mask_neg,3));

	int32x4_t f_F = vshrq_n_s32(c,e);
	f_F = vmulq_s32(f_F,sign);
	f_F = vandq_s32(f_F,mask_neg);

	comp = vdupq_n_s32(256);
	sign = vdupq_n_s32(-1);
	c = vabdq_s32(vabsq_s32(x),comp);
	//printf ("c_true= %d, %d, %d, %d\n", vgetq_lane_s32(c,0), vgetq_lane_s32(c,1), vgetq_lane_s32(c,2), vgetq_lane_s32(c,3));
	int32x4_t f_T = vshrq_n_s32(c,e);
	f_T = vmulq_s32(f_T,sign);
	f_T = vandq_s32(f_T,mask);

	//printf ("f_F %d, %d, %d, %d\n", vgetq_lane_s32(f_F,0), vgetq_lane_s32(f_F,1), vgetq_lane_s32(f_F,2), vgetq_lane_s32(f_F,3));
	//printf ("f_T %d, %d, %d, %d\n", vgetq_lane_s32(f_T,0), vgetq_lane_s32(f_T,1), vgetq_lane_s32(f_T,2), vgetq_lane_s32(f_T,3));

	int32x4_t f = vaddq_s32(f_F,f_T);
	return f;
	//printf ("f = %d, %d, %d, %d\n", vgetq_lane_s32(f,0), vgetq_lane_s32(f,1), vgetq_lane_s32(f,2), vgetq_lane_s32(f,3));

}
int fix_to_flo_Forced(int32_t x, int e){
    // f es el numero en punto fijo de tipo int
    // e es la cantidad de bits que se le dio la seccion decimal
	int32_t f;
    int sign;
    int c;

    c = abs(x);
    sign = 1;

    if (x > 127){
    /* Las siguientes 3 lineas es para devolverlo del complemento a 2
        si el numero original era negativo */

        c = abs(c - 256);
        sign = -1;
    }

    /*Lo que se hace es simplemente multiplicar el numero
    en punto fijo dividido por 2 elevado a la e*/
    f = (1.0 * c)/(1<<e);
    f = f * sign;
    return f;
}

void print_int32 (int32x4_t data, char* name) {
    int i;
    static int32_t p[4];

    vst1q_s32 (p, data);

    printf ("%s = ", name);
    for (i = 0; i < 4; i++) {
	printf ("%d ", p[i]);
    }
    printf ("\n");
}

int main () {

	srand(time(NULL));

	int N_muestras = 64000;
    int32_t int32_data[N_muestras];
    int32_t int32_out[N_muestras];
    int32_t int32_out_NEON[N_muestras];


    int32_t num;
    int32x4_t v_in;
	int32x4_t v_out;
	int32_t out;

    for(int q=0; q<N_muestras; q++){
    	num = 1 + rand() % (256 - 1);
    	int32_data[q] = num;
    }

    clock_t t;
	int f;
	t = clock();
	printf ("Inicio de reloj para operaciones no vectoriales...\n");

    for(int m=0;m<N_muestras;m++){
		out = fix_to_flo_Forced(int32_data[m],2);
		int32_out[m]=out;
		//printf("%d ",int32_data[m]);
	}

    t = clock() - t;
	printf ("Tiempo de ejecución tomó %ld clicks completar el proceso (%f segundos).\n\n\n",t,((float)t)/CLOCKS_PER_SEC);

    printf("\n");

    // EJECUCION DE MANERA VECTORIAL

    t = clock();
    t = clock() - t;
    printf ("Inicio de reloj para operaciones vectoriales...\n");

    for(int n=0;n<N_muestras;n+=4){
    	v_in = vld1q_s32(int32_data+n);
    	//print_int32 (v_in, "data");
    	v_out = fix_to_flo_Forced_NEON(v_in,2);
    	int32_out_NEON[n]= vgetq_lane_s32(v_out,0);
    	int32_out_NEON[n+1]= vgetq_lane_s32(v_out,1);
    	int32_out_NEON[n+2]= vgetq_lane_s32(v_out,2);
    	int32_out_NEON[n+3]= vgetq_lane_s32(v_out,3);
	}
    printf ("Tiempo de ejecución tomó %ld clicks completar el proceso (%f segundos).\n\n\n",t,((float)t)/CLOCKS_PER_SEC);

    for(int q=0; q<N_muestras; q++){
		//printf("out= %d\n",int32_out[q]);
		//printf("out neon= %d\n",int32_out_NEON[q]);
    	if(int32_out[q]!=int32_out_NEON[q]){
    		printf("Error en la línea %d\n",q);
    	}
	}

    return 0;
}

