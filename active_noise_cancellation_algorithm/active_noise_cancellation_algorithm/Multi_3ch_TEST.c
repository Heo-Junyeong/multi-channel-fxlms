#include "common.h"

#ifdef __3CH__

#define FILTER_LENGTH 64
#define PRI_PATH_LENGTH 256
#define SEC_PATH_LENGTH 256
#define SIGNAL_LENGTH 160000

typedef float dtype;

// Functions
dtype Dot_Product(dtype* a, dtype* b, int size)
{
    dtype result = 0;
    for (int i = 0; i < size; i++)
        result += a[i] * b[i];
    return result; // check memory free
}

// load off-line estimated secondary path
short int p_1[PRI_PATH_LENGTH] = { 0, };
short int p_2[PRI_PATH_LENGTH] = { 0, };
short int p_3[PRI_PATH_LENGTH] = { 0, };

short int s_11[SEC_PATH_LENGTH] = { 0, };
short int s_12[SEC_PATH_LENGTH] = { 0, };
short int s_13[SEC_PATH_LENGTH] = { 0, };
short int s_21[SEC_PATH_LENGTH] = { 0, };
short int s_22[SEC_PATH_LENGTH] = { 0, };
short int s_23[SEC_PATH_LENGTH] = { 0, };
short int s_31[SEC_PATH_LENGTH] = { 0, };
short int s_32[SEC_PATH_LENGTH] = { 0, };
short int s_33[SEC_PATH_LENGTH] = { 0, };

dtype p_1f[SEC_PATH_LENGTH] = { 0, };
dtype p_2f[SEC_PATH_LENGTH] = { 0, };
dtype p_3f[SEC_PATH_LENGTH] = { 0, };

dtype s_11f[SEC_PATH_LENGTH] = { 0, };
dtype s_12f[SEC_PATH_LENGTH] = { 0, };
dtype s_13f[SEC_PATH_LENGTH] = { 0, };
dtype s_21f[SEC_PATH_LENGTH] = { 0, };
dtype s_22f[SEC_PATH_LENGTH] = { 0, };
dtype s_23f[SEC_PATH_LENGTH] = { 0, };
dtype s_31f[SEC_PATH_LENGTH] = { 0, };
dtype s_32f[SEC_PATH_LENGTH] = { 0, };
dtype s_33f[SEC_PATH_LENGTH] = { 0, };


// Signals
short int ref_int[SIGNAL_LENGTH] = { 0, };
dtype ref_float[SIGNAL_LENGTH] = { 0, };

short int primary_noise_int[SIGNAL_LENGTH] = { 0, };
dtype primary_noise_float[SIGNAL_LENGTH] = { 0, };

dtype d1_vec[SEC_PATH_LENGTH] = { 0, };
dtype d2_vec[SEC_PATH_LENGTH] = { 0, };
dtype d3_vec[SEC_PATH_LENGTH] = { 0, };

dtype d1[SIGNAL_LENGTH] = { 0, };
dtype d2[SIGNAL_LENGTH] = { 0, };
dtype d3[SIGNAL_LENGTH] = { 0, };

dtype x_vec_w[FILTER_LENGTH] = { 0, };
dtype x_vec_p[SEC_PATH_LENGTH] = { 0, };

dtype y1_vec[SEC_PATH_LENGTH] = { 0, };
dtype y2_vec[SEC_PATH_LENGTH] = { 0, };
dtype y3_vec[SEC_PATH_LENGTH] = { 0, };

dtype xp11_vec[FILTER_LENGTH] = { 0, };
dtype xp12_vec[FILTER_LENGTH] = { 0, };
dtype xp13_vec[FILTER_LENGTH] = { 0, };
dtype xp21_vec[FILTER_LENGTH] = { 0, };
dtype xp22_vec[FILTER_LENGTH] = { 0, };
dtype xp23_vec[FILTER_LENGTH] = { 0, };
dtype xp31_vec[FILTER_LENGTH] = { 0, };
dtype xp32_vec[FILTER_LENGTH] = { 0, };
dtype xp33_vec[FILTER_LENGTH] = { 0, };

// TEST VECTOR
dtype xp11_vec_t1[FILTER_LENGTH] = { 0, };
dtype xp11_vec_t2[FILTER_LENGTH * 2] = { 0, };
dtype* ptr;

dtype w1[FILTER_LENGTH] = { 0, };
dtype w2[FILTER_LENGTH] = { 0, };
dtype w3[FILTER_LENGTH] = { 0, };

dtype y1, y2, y3;
dtype y1_p, y2_p, y3_p;
dtype e1, e2, e3;
dtype xp_11, xp_12, xp_13,
xp_21, xp_22, xp_23,
xp_31, xp_32, xp_33;

short int e1_debug[SIGNAL_LENGTH] = { 0, };
short int e2_debug[SIGNAL_LENGTH] = { 0, };
short int e3_debug[SIGNAL_LENGTH] = { 0, };

short int e_sum_int[SIGNAL_LENGTH] = { 0, };
dtype e_sum;

int main()
{
    clock_t start1, end1, start2, end2, start3, end3;

    // NLMS Parameters
    dtype alpha = 0.01;
    dtype pwr = 1;
    dtype beta = 1 - alpha;
    dtype tmp;
    dtype mu;
    dtype ep = 1e-6;
    int ptr_idx = FILTER_LENGTH;

    FILE* fp_p1, * fp_p2, * fp_p3;

    FILE* fp_s11, * fp_s12, * fp_s13,
        * fp_s21, * fp_s22, * fp_s23,
        * fp_s31, * fp_s32, * fp_s33;

    FILE* fp_primary_noise;
    FILE* fp_error, * fp_e1, * fp_e2, * fp_e3;

    fp_primary_noise = fopen("sine_250_500_750_au_10s.pcm", "rb");
    fread(primary_noise_int, 2, SIGNAL_LENGTH, fp_primary_noise);
    for (int i = 0; i < SIGNAL_LENGTH; i++)
        primary_noise_float[i] = primary_noise_int[i] / 32768.f;

    fp_error = fopen("error_sum.pcm", "wb"); // sr = 16000 
    fp_e1 = fopen("error_1.pcm", "wb");
    fp_e2 = fopen("error_2.pcm", "wb");
    fp_e3 = fopen("error_3.pcm", "wb");

    /* Load Paths */
    fp_p1 = fopen(".\\paths\\p1_path_board.pcm", "rb");
    fp_p2 = fopen(".\\paths\\p2_path_board.pcm", "rb");
    fp_p3 = fopen(".\\paths\\p3_path_board.pcm", "rb");

    fp_s11 = fopen(".\\paths\\s11_path_board.pcm", "rb");
    fp_s12 = fopen(".\\paths\\s12_path_board.pcm", "rb");
    fp_s13 = fopen(".\\paths\\s13_path_board.pcm", "rb");
    fp_s21 = fopen(".\\paths\\s21_path_board.pcm", "rb");
    fp_s22 = fopen(".\\paths\\s22_path_board.pcm", "rb");
    fp_s23 = fopen(".\\paths\\s23_path_board.pcm", "rb");
    fp_s31 = fopen(".\\paths\\s31_path_board.pcm", "rb");
    fp_s32 = fopen(".\\paths\\s32_path_board.pcm", "rb");
    fp_s33 = fopen(".\\paths\\s33_path_board.pcm", "rb");

    fread(p_1, 2, PRI_PATH_LENGTH, fp_p1);
    fread(p_2, 2, PRI_PATH_LENGTH, fp_p2);
    fread(p_3, 2, PRI_PATH_LENGTH, fp_p3);

    fread(s_11, 2, SEC_PATH_LENGTH, fp_s11);
    fread(s_12, 2, SEC_PATH_LENGTH, fp_s12);
    fread(s_13, 2, SEC_PATH_LENGTH, fp_s13);
    fread(s_21, 2, SEC_PATH_LENGTH, fp_s21);
    fread(s_22, 2, SEC_PATH_LENGTH, fp_s22);
    fread(s_23, 2, SEC_PATH_LENGTH, fp_s23);
    fread(s_31, 2, SEC_PATH_LENGTH, fp_s31);
    fread(s_32, 2, SEC_PATH_LENGTH, fp_s32);
    fread(s_33, 2, SEC_PATH_LENGTH, fp_s33);

    // normalization
    for (int i = 0; i < PRI_PATH_LENGTH; i++) {
        p_1f[i] = p_1[i] / 32768.f;
        p_2f[i] = p_2[i] / 32768.f;
        p_3f[i] = p_3[i] / 32768.f;
    }
    
    for (int i = 0; i < SEC_PATH_LENGTH; i++) {
        s_11f[i] = s_11[i] / 32768.f;
        s_12f[i] = s_12[i] / 32768.f;
        s_13f[i] = s_13[i] / 32768.f;
        s_21f[i] = s_21[i] / 32768.f;
        s_22f[i] = s_22[i] / 32768.f;
        s_23f[i] = s_23[i] / 32768.f;
        s_31f[i] = s_31[i] / 32768.f;
        s_32f[i] = s_32[i] / 32768.f;
        s_33f[i] = s_33[i] / 32768.f;
    }

    ptr = &xp11_vec_t2;

    for (int z = 0; z < 10; z++)
    {
        clock_t s1, s2;
        clock_t t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;

        s1 = clock();

        /* Main Routine */
        for (int i = 0; i < SIGNAL_LENGTH; i++)
        {
            if (i % (SIGNAL_LENGTH / 16) == 0)  t1 = clock();
            
            // Use primary noise as reference signal (for filtering)
            // % ## [section 1]
            ref_float[i] = primary_noise_float[i];
            for (int j = FILTER_LENGTH - 1; j > 0; j--) x_vec_w[j] = x_vec_w[j - 1];
            x_vec_w[0] = ref_float[i];
            

            // Filtering routine to make d1,d2,d3
            for (int j = PRI_PATH_LENGTH - 1; j > 0; j--) d1_vec[j] = d1_vec[j - 1];
            d1_vec[0] = primary_noise_float[i];
            d1[i] = Dot_Product(d1_vec, p_1f, PRI_PATH_LENGTH);

            for (int j = PRI_PATH_LENGTH - 1; j > 0; j--) d2_vec[j] = d2_vec[j - 1];
            d2_vec[0] = primary_noise_float[i];
            d2[i] = Dot_Product(d2_vec, p_2f, PRI_PATH_LENGTH);

            for (int j = PRI_PATH_LENGTH - 1; j > 0; j--) d3_vec[j] = d3_vec[j - 1];
            d3_vec[0] = primary_noise_float[i];
            d3[i] = Dot_Product(d3_vec, p_3f, PRI_PATH_LENGTH);

            // Make filter output y
           
            y1 = Dot_Product(w1, x_vec_w, FILTER_LENGTH);
            y2 = Dot_Product(w2, x_vec_w, FILTER_LENGTH);
            y3 = Dot_Product(w3, x_vec_w, FILTER_LENGTH);

            for (int j = SEC_PATH_LENGTH - 1; j > 0; j--) y1_vec[j] = y1_vec[j - 1];
            y1_vec[0] = y1;
            for (int j = SEC_PATH_LENGTH - 1; j > 0; j--) y2_vec[j] = y2_vec[j - 1];
            y2_vec[0] = y2;
            for (int j = SEC_PATH_LENGTH - 1; j > 0; j--) y3_vec[j] = y3_vec[j - 1];
            y3_vec[0] = y3;

            // Make y prime 
            y1_p = Dot_Product(s_11f, y1_vec, SEC_PATH_LENGTH) + Dot_Product(s_12f, y2_vec, SEC_PATH_LENGTH) + Dot_Product(s_13f, y3_vec, SEC_PATH_LENGTH);
            y2_p = Dot_Product(s_21f, y1_vec, SEC_PATH_LENGTH) + Dot_Product(s_22f, y2_vec, SEC_PATH_LENGTH) + Dot_Product(s_23f, y3_vec, SEC_PATH_LENGTH);
            y3_p = Dot_Product(s_31f, y1_vec, SEC_PATH_LENGTH) + Dot_Product(s_32f, y2_vec, SEC_PATH_LENGTH) + Dot_Product(s_33f, y3_vec, SEC_PATH_LENGTH);

            e1 = d1[i] + y1_p;
            e2 = d2[i] + y2_p;
            e3 = d3[i] + y3_p;

            // Make x prime 
            for (int j = SEC_PATH_LENGTH - 1; j > 0; j--) x_vec_p[j] = x_vec_p[j - 1];
            x_vec_p[0] = ref_float[i];

            xp_11 = Dot_Product(s_11f, x_vec_p, SEC_PATH_LENGTH);
            xp_12 = Dot_Product(s_21f, x_vec_p, SEC_PATH_LENGTH);
            xp_13 = Dot_Product(s_31f, x_vec_p, SEC_PATH_LENGTH);
            xp_21 = Dot_Product(s_12f, x_vec_p, SEC_PATH_LENGTH);
            xp_22 = Dot_Product(s_22f, x_vec_p, SEC_PATH_LENGTH);
            xp_23 = Dot_Product(s_32f, x_vec_p, SEC_PATH_LENGTH);
            xp_31 = Dot_Product(s_13f, x_vec_p, SEC_PATH_LENGTH);
            xp_32 = Dot_Product(s_23f, x_vec_p, SEC_PATH_LENGTH);
            xp_33 = Dot_Product(s_33f, x_vec_p, SEC_PATH_LENGTH);

            //start1 = (double)(clock() / CLOCKS_PER_SEC);
            //memmove(xp11_vec_t1 + 1, xp11_vec_t1, sizeof(dtype) * (FILTER_LENGTH - 1));
            //xp11_vec_t1[0] = xp_11;
            //end1 = (double)(clock() / CLOCKS_PER_SEC);

            //start2 = (double)(clock() / CLOCKS_PER_SEC);
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp11_vec[j] = xp11_vec[j - 1];
            xp11_vec[0] = xp_11;
            //end2 = (double)(clock() / CLOCKS_PER_SEC);

            //start3 = (double)(clock() / CLOCKS_PER_SEC);
            //// Add routine as YC said
            //// xp11_vec_t2[0] = xp_11;
            //*(ptr + ptr_idx - 1) = xp_11;
            //*(ptr + FILTER_LENGTH + ptr_idx - 1) = xp_11;
            ////ptr += FILTER_LENGTH;
            //ptr_idx--;
            //if (ptr_idx == 0) ptr_idx = FILTER_LENGTH;
            //// for (int i = 0; i < size; i++) result += a[i] * b[i];
            //// Pointer Dot Product
            //// (FILTER_LENGTH/2)
            ////xp_11_vec_t2;
            //end3 = (double)(clock() / CLOCKS_PER_SEC);

            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp12_vec[j] = xp12_vec[j - 1];
            xp12_vec[0] = xp_12;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp13_vec[j] = xp13_vec[j - 1];
            xp13_vec[0] = xp_13;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp21_vec[j] = xp21_vec[j - 1];
            xp21_vec[0] = xp_21;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp22_vec[j] = xp22_vec[j - 1];
            xp22_vec[0] = xp_22;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp23_vec[j] = xp23_vec[j - 1];
            xp23_vec[0] = xp_23;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp31_vec[j] = xp31_vec[j - 1];
            xp31_vec[0] = xp_31;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp32_vec[j] = xp32_vec[j - 1];
            xp32_vec[0] = xp_32;
            for (int j = FILTER_LENGTH - 1; j > 0; j--) xp33_vec[j] = xp33_vec[j - 1];
            xp33_vec[0] = xp_33;

            // mu update
            mu = alpha / (FILTER_LENGTH * (pwr + ep));

            // weight update 
            for (int j = 0; j < FILTER_LENGTH; j++) w1[j] = w1[j] - mu * (xp11_vec[j] * e1 + xp12_vec[j] * e2 + xp13_vec[j] * e3);
            for (int j = 0; j < FILTER_LENGTH; j++) w2[j] = w2[j] - mu * (xp21_vec[j] * e1 + xp22_vec[j] * e2 + xp23_vec[j] * e3);
            for (int j = 0; j < FILTER_LENGTH; j++) w3[j] = w3[j] - mu * (xp31_vec[j] * e1 + xp32_vec[j] * e2 + xp33_vec[j] * e3);

            pwr = beta * pwr + (1 - beta) * (xp_11 * xp_11 + xp_12 * xp_12 + xp_13 * xp_13);

            e_sum = (e1 + e2 + e3) / 3;
            e_sum_int[i] = (short int)(e_sum * 32767);

            e1_debug[i] = (short int)(e1 * 32767);
            e2_debug[i] = (short int)(e2 * 32767);
            e3_debug[i] = (short int)(e3 * 32767);

            if (i % (SIGNAL_LENGTH / 16) == 0)
            {
                t2 = clock();
                printf("section 1, clock : %d\n", t2 - t1);
            }

        }

        s2 = clock();

        printf("iter %d, clock : %d\n", z, s2 - s1);
        printf("section 1, clock : %d\n", t2-t1);
        printf("section 2, clock : %d\n", t4-t3);
        printf("section 3, clock : %d\n\n", t6-t5);
    }

    fwrite(e_sum_int, 2, SIGNAL_LENGTH, fp_error);
    fwrite(e1_debug, 2, SIGNAL_LENGTH, fp_e1);
    fwrite(e2_debug, 2, SIGNAL_LENGTH, fp_e2);
    fwrite(e3_debug, 2, SIGNAL_LENGTH, fp_e3);

    //printf("for loop Elapsed time : %lf  \n", (end2 - start2));
    //printf("memmove Elapsed time : %lf  \n", (end1 - start1));

    return 0;
}

#endif