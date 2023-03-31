/*
 * Copyright (C) 2014 David Milligan
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "hist.h"
#include "dualiso.h"
#include "opt_med.h"
#include "wirth.h"
//#include <pthread.h>
#include "../../debayer/debayer.h"
#include <time.h>
#include <omp.h>
#include <assert.h>

#undef STDOUT_SILENT

#define EV_RESOLUTION 65536
#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define COERCE(x,lo,hi) MAX(MIN((x),(hi)),(lo))
#define ABS(a) ((a) > 0 ? (a) : -(a))

//#define LOCK(x) static pthread_mutex_t x = PTHREAD_MUTEX_INITIALIZER; pthread_mutex_lock(&x);
//#define UNLOCK(x) pthread_mutex_unlock(&(x));

#ifdef PERF_INFO
clock_t perf_clock, perf_sub_clock;
#endif

//this is just meant to be fast
int diso_get_preview(uint16_t * image_data,dual_iso_freeze_data_t* iso_data, uint16_t width, uint16_t height, int32_t black, int32_t white, int diso_check)
{
    double a,b;
    uint16_t dark_row_start = -1;

    if (iso_data->freeze < 2){
        //compute the median of the green channel for each multiple of 4 rows
        uint16_t median[4];
        struct histogram * hist[4];
        struct histogram * hist_hi = NULL;
        struct histogram * hist_lo = NULL;
        
        for(int i = 0; i < 4; i++)
            hist[i] = hist_create(white);

        for(uint16_t y = 4; y < height - 4; y += 5)
        {
            hist_add(hist[y % 4], &(image_data[y * width + (y + 1) % 2]), width - (y + 1) % 2, 3);
        }
        
        for(int i = 0; i < 4; i++)
        {
            median[i] = hist_median(hist[i]);
        }
        
        if((median[2] - black) > ((median[0] - black) * 2) &&
        (median[2] - black) > ((median[1] - black) * 2) &&
        (median[3] - black) > ((median[0] - black) * 2) &&
        (median[3] - black) > ((median[1] - black) * 2))
        {
            dark_row_start = 0;
            hist_lo = hist[0];
            hist_hi = hist[2];
        }
        else if((median[0] - black) > ((median[1] - black) * 2) &&
                (median[0] - black) > ((median[2] - black) * 2) &&
                (median[3] - black) > ((median[1] - black) * 2) &&
                (median[3] - black) > ((median[2] - black) * 2))
        {
            dark_row_start = 1;
            hist_lo = hist[1];
            hist_hi = hist[0];
        }
        else if((median[0] - black) > ((median[2] - black) * 2) &&
                (median[0] - black) > ((median[3] - black) * 2) &&
                (median[1] - black) > ((median[2] - black) * 2) &&
                (median[1] - black) > ((median[3] - black) * 2))
        {
            dark_row_start = 2;
            hist_lo = hist[2];
            hist_hi = hist[0];
        }
        else if((median[1] - black) > ((median[0] - black) * 2) &&
                (median[1] - black) > ((median[3] - black) * 2) &&
                (median[2] - black) > ((median[0] - black) * 2) &&
                (median[2] - black) > ((median[3] - black) * 2))
        {
            dark_row_start = 3;
            hist_lo = hist[0];
            hist_hi = hist[2];
        }
        else
        {
    #ifndef STDOUT_SILENT
            err_printf("\nCould not detect dual ISO interlaced lines\n");
    #endif

            for(int i = 0; i < 4; i++)
            {
                hist_destroy(hist[i]);
            }
            return 0;
        }
        
        if(diso_check)
        {
    #ifndef STDOUT_SILENT
            err_printf("\nDetected dual ISO interlaced lines\n");
    #endif

            for(int i = 0; i < 4; i++)
            {
                hist_destroy(hist[i]);
            }
            return 1;
        }

        /* compare the two histograms and plot the curve between the two exposures (dark as a function of bright) */
        const int min_pix = 100;                                /* extract a data point every N image pixels */
        int data_size = (width * height / min_pix + 1);                  /* max number of data points */
        int* data_x = (int *)malloc(data_size * sizeof(data_x[0]));
        int* data_y = (int *)malloc(data_size * sizeof(data_y[0]));
        double* data_w = (double *)malloc(data_size * sizeof(data_w[0]));
        int data_num = 0;
        
        int acc_lo = 0;
        int acc_hi = 0;
        int raw_lo = 0;
        int raw_hi = 0;
        int prev_acc_hi = 0;
        
        int hist_total = hist[0]->count;
        
        for (raw_hi = 0; raw_hi < hist_total; raw_hi++)
        {
            acc_hi += hist_hi->data[raw_hi];
            
            while (acc_lo < acc_hi)
            {
                acc_lo += hist_lo->data[raw_lo];
                raw_lo++;
            }
            
            if (raw_lo >= white)
                break;
            
            if (acc_hi - prev_acc_hi > min_pix)
            {
                if (acc_hi > hist_total * 1 / 100 && acc_hi < hist_total * 99.99 / 100)    /* throw away outliers */
                {
                    data_x[data_num] = raw_hi - black;
                    data_y[data_num] = raw_lo - black;
                    data_w[data_num] = (MAX(0, raw_hi - black + 100));    /* points from higher brightness are cleaner */
                    data_num++;
                    prev_acc_hi = acc_hi;
                }
            }
        }
        
        /**
         * plain least squares
         * y = ax + b
         * a = (mean(xy) - mean(x)mean(y)) / (mean(x^2) - mean(x)^2)
         * b = mean(y) - a mean(x)
         */
        
        double mx = 0, my = 0, mxy = 0, mx2 = 0;
        double weight = 0;
        for (int i = 0; i < data_num; i++)
        {
            mx += data_x[i] * data_w[i];
            my += data_y[i] * data_w[i];
            mxy += (double)data_x[i] * data_y[i] * data_w[i];
            mx2 += (double)data_x[i] * data_x[i] * data_w[i];
            weight += data_w[i];
        }
        mx /= weight;
        my /= weight;
        mxy /= weight;
        mx2 /= weight;
        a = (mxy - mx*my) / (mx2 - mx*mx);
        b = my - a * mx;
        
        free(data_w);
        free(data_y);
        free(data_x);

        for(int i = 0; i < 4; i++)
        {
            hist_destroy(hist[i]);
        }
    } else {
        a = iso_data->a;
        b = iso_data->b;
        dark_row_start = iso_data->dark_row_start;
    }

    if (iso_data->freeze == 1){
        iso_data->freeze = 2;
        iso_data->a = a;
        iso_data->b = b;
        iso_data->dark_row_start = dark_row_start;
    }
    
    //TODO: what's a better way to pick a value for this?
    uint16_t shadow = (uint16_t)(black + 1 / (a * a) + b);
    
    for(int y = 0; y < height; y++)
    {
        int row_start = y * width;
        if (((y - dark_row_start + 4) % 4) >= 2)
        {
            //bright row
            for(int i = row_start; i < row_start + width; i++)
            {
                if(image_data[i] >= white)
                {
                    image_data[i] = y > 2 ? (y < height - 2 ? (image_data[i-width*2] + image_data[i+width*2]) / 2 : image_data[i-width*2]) : image_data[i+width*2];
                }
                else
                {
                    image_data[i] = (uint16_t)(MIN(white,(image_data[i] - black) * a + black + b));
                }
            }
        }
        else
        {
            //dark row
            for(int i = row_start; i < row_start + width; i++)
            {
                if(image_data[i] < shadow)
                {
                    image_data[i] = (uint16_t)(y > 2 ? (y < height - 2 ? (image_data[i-width*2] + MIN(white,(image_data[i+width*2]  - black) * a + black + b)) / 2 : image_data[i-width*2]) : MIN(white,(image_data[i+width*2]  - black) * a + black + b));
                }
                
            }
        }
    }

    return 1;
}


//from cr2hdr 20bit version
//this is not thread safe (yet)

#define BRIGHT_ROW (is_bright[y & 0b11])
#define COUNT(x) ((int)(sizeof(x)/sizeof((x)[0])))

#define raw_get_pixel(x,y) (image_data[(x) + (y) * raw_info.width])
#define raw_get_pixel16(x,y) (image_data[(x) + (y) * raw_info.width])
#define raw_get_pixel20(x,y) ((raw_buffer_32[(x) + (y) * raw_info.width]) & 0xFFFFF)
#define raw_get_pixel_14to20(x,y) ((((uint32_t)image_data[(x) + (y) * raw_info.width]) << 6) & 0xFFFFF)
#define raw_get_pixel32(x,y) (raw_buffer_32[(x) + (y) * raw_info.width])
#define raw_set_pixel32(x,y,value) raw_buffer_32[(x) + (y)*raw_info.width] = value
#define raw_get_pixel_20to16(x,y) ((raw_get_pixel32(x,y) >> 4) & 0xFFFF)
#define raw_set_pixel_20to16_rand(i,value) image_data[i] = COERCE((int)((value) / 16.0 + fast_randn05() + 0.5), 0, 0xFFFF)
#define raw_set_pixel20(x,y,value) raw_buffer_32[(x) + (y) * raw_info.width] = COERCE((value), 0, 0xFFFFF)

static const double fullres_start = 4;
static const double fullres_transition = 4;
static const double fullres_thr = 0.8;

/* trial and error - too high = aliasing, too low = noisy */
static const int ALIAS_MAP_MAX = 15000;

static inline int FC(int row, int col)
{
    if ((row%2) == 0 && (col%2) == 0)
        return 0;  /* red */
    else if ((row%2) == 1 && (col%2) == 1)
        return 2;  /* blue */
    else
        return 1;  /* green */
}
static void white_detect(struct raw_info raw_info, uint16_t * image_data, int* white_dark, int* white_bright, int * is_bright)
{
    /* sometimes the white level is much lower than 15000; this would cause pink highlights */
    /* workaround: consider the white level as a little under the maximum pixel value from the raw file */
    /* caveat: bright and dark exposure may have different white levels, so we'll take the minimum value */
    /* side effect: if the image is not overexposed, it may get brightened a little; shouldn't hurt */
    
    int whites[2]         = {  0,    0};
    int discard_pixels[2] = { 10,   50}; /* discard the brightest N pixels */
    int safety_margins[2] = {100, 1500}; /* use a higher safety margin for the higher ISO */
    /* note: with the high-ISO WL underestimated by 1500, you would lose around 0.15 EV of non-aliased detail */
    
    int* pixels[2];
    int max_pix = raw_info.width * raw_info.height / 2 / 9;
    pixels[0] = malloc(max_pix * sizeof(pixels[0][0]));
    pixels[1] = malloc(max_pix * sizeof(pixels[0][0]));
    memset(pixels[0], 0, sizeof(max_pix * sizeof(pixels[0][0])));
    memset(pixels[1], 0, sizeof(max_pix * sizeof(pixels[0][0])));
    int counts[2] = {0, 0};
    
    /* collect all the pixels and find the k-th max, thus ignoring hot pixels */
    /* change the sign in order to use kth_smallest_int */

#pragma omp parallel for schedule(static) default(none) shared(raw_info, image_data, max_pix, counts, pixels, is_bright) collapse(2)
    for (int y = raw_info.active_area.y1; y < raw_info.active_area.y2; y += 3)
    {
        for (int x = raw_info.active_area.x1; x < raw_info.active_area.x2; x += 3)
        {
            int pix = raw_get_pixel16(x, y);
            
#pragma omp critical
{
            counts[BRIGHT_ROW] = MIN(counts[BRIGHT_ROW], max_pix-1);
            pixels[BRIGHT_ROW][counts[BRIGHT_ROW]] = -pix;
            counts[BRIGHT_ROW]++;
}
        }
    }
    
    whites[0] = -kth_smallest_int(pixels[0], counts[0], discard_pixels[0]) - safety_margins[0];
    whites[1] = -kth_smallest_int(pixels[1], counts[1], discard_pixels[1]) - safety_margins[1];
    
    /* we assume 14-bit input data; out-of-range white levels may cause crash */
    *white_dark = COERCE(whites[0], 5000, 16383);
    *white_bright = COERCE(whites[1], 5000, 16383);
#ifndef STDOUT_SILENT
    printf("White levels    : %d %d\n", *white_dark, *white_bright);
#endif
    free(pixels[0]);
    free(pixels[1]);
}

static void compute_black_noise20(struct raw_info raw_info, uint32_t * raw_buffer_32,int y1, int y2, int dx, int dy, double* out_mean, double* out_stdev)
{
    int x1 = 8;
    int x2 = raw_info.active_area.x1 - 8;

    if (x2 <= 0){
        // We have no access to sensor device "black bars" zone to measure noise (because it's cropped I guess)
        *out_mean = raw_info.black_level;
        *out_stdev = 8+6; /* default to 11 stops of DR */
        return;
    }

    double mean, stdev;

    long long black = 0;
    int num = 0;
    /* compute average level */
    for (int y = y1; y < y2; y += dy)
    {
        for (int x = x1; x < x2; x += dx)
        {
            black += raw_get_pixel20(x, y);
            num++;
        }
    }
    mean = (double) black / num;

    /* compute standard deviation */
    stdev = 0;
    for (int y = y1; y < y2; y += dy)
    {
        for (int x = x1; x < x2; x += dx)
        {
            double dif = raw_get_pixel20(x, y) - mean;
            stdev += dif * dif;
        }
    }
    stdev /= (num-1);
    stdev = sqrt(stdev);
    
    if (num == 0)
    {
        mean = raw_info.black_level;
        /* default to 11 stops (14-log2(8)) of DR */
        stdev = 8 + 6;  // (+6 for 20bits res)
    }
    
    *out_mean = mean;
    *out_stdev = stdev;
}

static int mean2(int a, int b, int white, int* err)
{
    if (a >= white || b >= white)
    {
        if (err) *err = 10000000;
        return white;
    }
    
    int m = (a + b) / 2;
    
    if (err)
        *err = ABS(a - b);
    
    return m;
}

static int mean3(int a, int b, int c, int white, int* err)
{
    int m = (a + b + c) / 3;
    
    if (err)
        *err = MAX(MAX(ABS(a - m), ABS(b - m)), ABS(c - m));
    
    if (a >= white || b >= white || c >= white)
        return MAX(m, white);
    
    return m;
}

/* http://www.developpez.net/forums/d544518/c-cpp/c/equivalent-randn-matlab-c/#post3241904 */

#define TWOPI (6.2831853071795864769252867665590057683943387987502) /* 2 * pi */

/*
 RAND is a macro which returns a pseudo-random numbers from a uniform
 distribution on the interval [0 1]
 */
#define RAND (rand())/((double) RAND_MAX)

/*
 RANDN is a macro which returns a pseudo-random numbers from a normal
 distribution with mean zero and standard deviation one. This macro uses Box
 Muller's algorithm
 */
#define RANDN (sqrt(-2.0*log(RAND))*cos(TWOPI*RAND))

/* anti-posterization noise */
/* before rounding, it's a good idea to add a Gaussian noise of stdev=0.5 */
static float randn05_cache[1024];

void fast_randn_init()
{
    int i;
#pragma omp parallel for schedule(static) default(none) shared(randn05_cache)
    for (i = 0; i < 1024; i++)
    {
        randn05_cache[i] = RANDN / 2;
    }
}

float fast_randn05()
{
    static int k = 0;
    return randn05_cache[(k++) & 1023];
}

static int identify_rggb_or_gbrg(struct raw_info raw_info, uint16_t * image_data)
{
    int w = raw_info.width;
    int h = raw_info.height;
    
    /* build 4 little histograms: one for red, one for blue and two for green */
    /* we don't know yet which channels are which, but that's what we are trying to find out */
    /* the ones with the smallest difference are likely the green channels */
    int hist_size = 16384 * sizeof(int);
    int* hist[4];
    for (int i = 0; i < 4; i++)
    {
        hist[i] = malloc(hist_size);
        memset(hist[i], 0, hist_size);
    }
    
    int y0 = (raw_info.active_area.y1 + 3) & ~3;
    
    /* to simplify things, analyze an identical number of bright and dark lines */
#pragma omp parallel for schedule(static) default(none) shared(raw_info, image_data, y0, h, w,hist) collapse(2)
    for (int y = y0; y < h/4*4; y++)
    {
        for (int x = 0; x < w; x++){
#pragma omp atomic
            hist[(y%2)*2 + (x%2)][raw_get_pixel16(x,y) & 16383]++;
        }
    }
    
    /* compute cdf */
#pragma omp parallel for schedule(static,1) num_threads(4) default(none) shared(hist)
    for (int k = 0; k < 4; k++)
    {
        int acc = 0;
        for (int i = 0; i < 16384; i++)
        {
            acc += hist[k][i];
            hist[k][i] = acc;
        }
    }
    
    /* compare cdf's */
    /* for rggb, greens are at y%2 != x%2, that is, 1 and 2 */
    /* for gbrg, greens are at y%2 == x%2, that is, 0 and 3 */
    double diffs_rggb = 0;
    double diffs_gbrg = 0;
#pragma omp parallel for schedule(static) default(none) shared(hist) reduction(+:diffs_rggb) reduction(+:diffs_gbrg)
    for (int i = 0; i < 16384; i++)
    {
        diffs_rggb += ABS(hist[1][i] - hist[2][i]);
        diffs_gbrg += ABS(hist[0][i] - hist[3][i]);
    }
    
    for (int i = 0; i < 4; i++)
    {
        free(hist[i]); hist[i] = 0;
    }
    
    /* which one is most likely? */
    return diffs_rggb < diffs_gbrg;
}

static int identify_bright_and_dark_fields(struct raw_info raw_info, uint16_t * image_data, int rggb, int * is_bright)
{
    /* first we need to know which lines are dark and which are bright */
    /* the pattern is not always the same, so we need to autodetect it */
    
    /* it may look like this */                       /* or like this */
    /*
     ab cd ef gh  ab cd ef gh               ab cd ef gh  ab cd ef gh
     
     0  RG RG RG RG  RG RG RG RG            0  rg rg rg rg  rg rg rg rg
     1  gb gb gb gb  gb gb gb gb            1  gb gb gb gb  gb gb gb gb
     2  rg rg rg rg  rg rg rg rg            2  RG RG RG RG  RG RG RG RG
     3  GB GB GB GB  GB GB GB GB            3  GB GB GB GB  GB GB GB GB
     4  RG RG RG RG  RG RG RG RG            4  rg rg rg rg  rg rg rg rg
     5  gb gb gb gb  gb gb gb gb            5  gb gb gb gb  gb gb gb gb
     6  rg rg rg rg  rg rg rg rg            6  RG RG RG RG  RG RG RG RG
     7  GB GB GB GB  GB GB GB GB            7  GB GB GB GB  GB GB GB GB
     8  RG RG RG RG  RG RG RG RG            8  rg rg rg rg  rg rg rg rg
     */
    
    /* white level is not yet known, just use a rough guess */
    int white = 10000;
    int black = raw_info.black_level;
    
    int w = raw_info.width;
    int h = raw_info.height;
    
    /* build 4 little histograms */
    int hist_size = 16384 * sizeof(int);
    int* hist[4];
    for (int i = 0; i < 4; i++)
    {
        hist[i] = malloc(hist_size);
        memset(hist[i], 0, hist_size);
    }
    
    int y0 = (raw_info.active_area.y1 + 3) & ~3;
    
    /* to simplify things, analyze an identical number of bright and dark lines */
#pragma omp parallel for schedule(static) default(none) shared(raw_info, image_data, y0,h,w,hist) collapse(2)
    for (int y = y0; y < h/4*4; y++)
    {
        for (int x = 0; x < w; x++)
        {
            if ((x%2) != (y%2))
            {
                /* only check the green pixels */
#pragma omp atomic
                hist[y%4][raw_get_pixel16(x,y) & 16383]++;
            }
        }
    }
    
    int hist_total = 0;
    for (int i = 0; i < 16384; i++)
        hist_total += hist[0][i];
    
    /* choose the highest percentile that is not overexposed */
    /* but not higher than 99.8, to keep a tiny bit of robustness (specular highlights may play dirty tricks) */
    int acc[4] = {0};
    int raw[4] = {0};
    int off[4] = {0};
    int ref;
    int ref_max = hist_total * 0.998;
    int ref_off = hist_total * 0.05;
    for (ref = 0; ref < ref_max; ref++)
    {
        for (int i = 0; i < 4; i++)
        {
            while (acc[i] < ref)
            {
                acc[i] += hist[i][raw[i]];
                raw[i]++;
            }
        }
        
        if (ref < ref_off)
        {
            if (MAX(MAX(raw[0], raw[1]), MAX(raw[2], raw[3])) < black + (white-black) / 4)
            {
                /* try to remove the black offset by estimating it from relatively dark pixels */
                off[0] = raw[0];
                off[1] = raw[1];
                off[2] = raw[2];
                off[3] = raw[3];
            }
        }
        
        if (raw[0] >= white) break;
        if (raw[1] >= white) break;
        if (raw[2] >= white) break;
        if (raw[3] >= white) break;
    }
    
    for (int i = 0; i < 4; i++)
    {
        free(hist[i]); hist[i] = 0;
    }
    
    /* remove black offsets */
    raw[0] -= off[0];
    raw[1] -= off[1];
    raw[2] -= off[2];
    raw[3] -= off[3];
    
    /* very crude way to compute median */
    int sorted_bright[4];
    memcpy(sorted_bright, raw, sizeof(sorted_bright));
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = i+1; j < 4; j++)
            {
                if (sorted_bright[i] > sorted_bright[j])
                {
                    double aux = sorted_bright[i];
                    sorted_bright[i] = sorted_bright[j];
                    sorted_bright[j] = aux;
                }
            }
        }
    }
    double median_bright = (sorted_bright[1] + sorted_bright[2]) / 2;
    
    for (int i = 0; i < 4; i++)
        is_bright[i] = raw[i] > median_bright;
#ifndef STDOUT_SILENT
    printf("ISO pattern     : %c%c%c%c %s\n", is_bright[0] ? 'B' : 'd', is_bright[1] ? 'B' : 'd', is_bright[2] ? 'B' : 'd', is_bright[3] ? 'B' : 'd', rggb ? "RGGB" : "GBRG");
#endif
    if (is_bright[0] + is_bright[1] + is_bright[2] + is_bright[3] != 2)
    {
#ifndef STDOUT_SILENT
        printf("Bright/dark detection error\n");
#endif
        return 0;
    }
    
    if (is_bright[0] == is_bright[2] || is_bright[1] == is_bright[3])
    {
#ifndef STDOUT_SILENT
        printf("Interlacing method not supported\n");
#endif
        return 0;
    }
    return 1;
}

static void apply_correction(double a, double b, int h, int w, struct raw_info raw_info, int black20, int white20, uint32_t* raw_buffer_32, double * corr_ev, int * white_darkened, int * is_bright)
{
    /* apply the correction */
#pragma omp parallel for schedule(static) default(none) shared(black20, is_bright, h, w, b, raw_buffer_32, raw_info, a) collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            int p = raw_get_pixel32(x, y);
            if (p == 0) continue;
            
            if (BRIGHT_ROW)
            {
                /* bright exposure: darken and apply the black offset (fixme: why not half?) */
                p = (p - black20) * a + black20 + b*a;
            }
            else
            {
                p = p - b + b*a;
            }
            
            /* out of range? */
            /* note: this breaks M24-1127 */
            p = COERCE(p, 0, 0xFFFFF);
            
            raw_set_pixel20(x, y, p);
        }
    }
    *white_darkened = (white20 - black20 + b) * a + black20;
    
    double factor = 1/a;
    *corr_ev = log2(factor);

#ifndef STDOUT_SILENT
    printf("ISO difference  : %.2f EV (%d) (a[%f] b[%f])\n", *corr_ev, (int)round(factor*100), a, b);
    printf("Whites  :       : %d(darkened) %d\n", *white_darkened, white20);
    printf("Black delta     : %.2f\n", b/4); /* we want to display black delta for the 14-bit original data, but we have computed it from 16-bit data */
#endif

#ifndef STDOUT_SILENT
    if (factor < 1.2 || !isfinite(factor))
    {
        printf("Doesn't look like interlaced ISO\n");
    }
#endif
}

static int match_exposures(struct raw_info raw_info, uint32_t * raw_buffer_32, double * corr_ev,
                           int * white_darkened, int * is_bright, dual_iso_freeze_data_t* iso_data)
{
    int plot_iso_curve = 0;

    /* guess ISO - find the factor and the offset for matching the bright and dark images */
    int black20 = raw_info.black_level;
    int white20 = MIN(raw_info.white_level, *white_darkened);
    int black = black20/16;
    int white = white20/16;
    int clip0 = white - black;
    int clip  = clip0 * 0.95;    /* there may be nonlinear response in very bright areas */

    int w = raw_info.width;
    int h = raw_info.height;
    int y0 = raw_info.active_area.y1 + 2;

    if (iso_data->freeze == 2){
        apply_correction(iso_data->a, iso_data->b, h, w, raw_info, black20, white20, raw_buffer_32, corr_ev, white_darkened, is_bright);
        return 1;
    }

    /* quick interpolation for matching */
    int* dark   = malloc(w * h * sizeof(dark[0]));
    int* bright = malloc(w * h * sizeof(bright[0]));
    memset(dark, 0, w * h * sizeof(dark[0]));
    memset(bright, 0, w * h * sizeof(bright[0]));
    
    for (int y = y0; y < h-2; y++)
    {
        int* native = BRIGHT_ROW ? bright : dark;
        int* interp = BRIGHT_ROW ? dark : bright;

        for (int x = 0; x < w; x++)
        {
            int pa = raw_get_pixel_20to16(x, y-2) - black;
            int pb = raw_get_pixel_20to16(x, y+2) - black;
            int pn = raw_get_pixel_20to16(x, y) - black;
            int pi = (pa + pb + 1) / 2;
            if (pa >= clip || pb >= clip || pn >= clip)
                pi = pn = clip0;               /* pixel too bright? discard */
            interp[x + y * w] = pi;
            native[x + y * w] = pn;
        }
    }

    /* 
     * Robust line fit (match unclipped data):
     * - divide the image into a MxM grid (M=50)
     * - consider only unclipped data (each ISO has its own white level)
     * - compute low-iso median for each cell (unclipped data only)
     * - compute median difference for each cell (unclipped data only)
     * - line fit: X = low iso medians, Y = low iso medians + median difference
     * - adjust (darken) high ISO to match low ISO
     * - repeat to fine-tune the result
     */

    /* y = a*x + b */
    /* y = bright, x = dark */
    double a = 1;
    double b = 0;

    /* fixme: when to stop? */
    for (int iter = 0; iter < 5; iter++)
    {
        const int M = 50;

        int Darks[3*M*M];
        int Brigs[3*M*M];
        int K = 0;

        for (int ch = 0; ch < 3; ch++)
        {
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    /* cell (i,j) */
                    const int nmax = (1 + w / M) * (1 + h / M);
                    int * darks = malloc(nmax * sizeof(darks[0]));
                    int * diffs = malloc(nmax * sizeof(diffs[0]));
                    int k = 0;
                    int all = 0;

                    for (int y = i * h / M; y < (i + 1) * h / M; y++)
                    {
                        if (y < y0) continue;
                        if (y >= h-2) continue;

                        for (int x = j * w / M; x < (j + 1) * w / M; x++)
                        {
                            /* one data point per channel */
                            /* fixme: scan the image once, not 3 times */
                            if (FC(y,x) != ch)
                            {
                                continue;
                            }

                            int lo = dark[x + y * w];
                            int hi = bright[x + y * w];

                            /* note: clipped pixels were all set to clip0 */
                            if (lo != clip0 && hi != clip0)
                            {
                                assert(k < nmax);
                                assert(lo < clip);
                                assert(hi < clip);
                                darks[k] = lo;
                                diffs[k] = (hi-b)/a - lo;
                                k++;
                            }
                            all++;
                        }
                    }

                    /* more than 25% valid pixels? */
                    if (k > all / 4)
                    {
                        assert(K < COUNT(Darks));
                        Darks[K] = median_int_wirth(darks, k);
                        int dmed = median_int_wirth(diffs, k);
                        Brigs[K] = (dmed + Darks[K]) * a + b;
                        if (0)
                        {
                            /* debug info */
                            /* reuse darks to check the interquartile range of diffs */
                            /* this will be very high at first iteration
                             * but much lower at the next one(s) */
                            int * quart = darks;
                            int q = 0;
                            for (int i = 0; i < k; i++)
                                if (diffs[i] < dmed)
                                    quart[q++] = diffs[i];
                            int p25 = median_int_wirth(quart, q);
                            q = 0;
                            for (int i = 0; i < k; i++)
                                if (diffs[i] > dmed)
                                    quart[q++] = diffs[i];
                            int p75 = median_int_wirth(quart, q);
                            printf("diff range: %d - %d\n", p25, p75);
                        }
                        K++;
                    }

                    free(darks);
                    free(diffs);
                }
            }
        }

        /* linear fit on Diffs vs Darks */

        /*
         * some sort of robust linear fitting
         * median for X, median for Y, median for angle (atan2)
         */
        int * aux = malloc(K * sizeof(aux[0]));
        memcpy(aux, Darks, K * sizeof(aux[0]));
        int mx = median_int_wirth(aux, K);
        memcpy(aux, Brigs, K * sizeof(aux[0]));
        int my = median_int_wirth(aux, K);
        free(aux);

        int * Angles = malloc(K * sizeof(Angles[0]));
        for (int i = 0; i < K; i++)
        {
            double ang = atan2(Brigs[i] - my, Darks[i] - mx);
            while (ang < 0) ang += M_PI;
            Angles[i] = (int)round(ang * 1000000);
        }
        int ma = median_int_wirth(Angles, K);
        free(Angles);

        /* convert to y = ax + b */
        a = tan(ma / 1000000.0);
        b = my - a * mx;

        /* print values normalized to 14-bit */
        /* b is usually negative */
#ifndef STDOUT_SILENT
        printf("Linear fit      : y = %.4f*x - %.2f\n", a, -b/64);
#endif
        if (plot_iso_curve)
        {
            FILE* f = fopen("iso-curve.m", "w");
            fprintf(f, "a = %g\n", a);
            fprintf(f, "b = %g\n", b);
            fprintf(f, "clip = %d\n", clip);
            fprintf(f, "Darks = [ ");
            for (int i = 0; i < K; i++)
                fprintf(f, "%d ", Darks[i]);
            fprintf(f, "];\n");
            fprintf(f, "Brigs = [ ");
            for (int i = 0; i < K; i++)
                fprintf(f, "%d ", Brigs[i]);
            fprintf(f, "];\n");
            fprintf(f, "black_delta = -b / (a + 1);");
            fprintf(f, "graphics_toolkit gnuplot;\n");
            fprintf(f, "plot(Darks - black_delta, Brigs + black_delta, '.r');\n");
            fprintf(f, "axis auto; set(gca,'xscale','log'); set(gca,'yscale','log'); axis tight; ax = axis;\n");
            fprintf(f, "hold on, plot(linspace(1, clip/a, 1000) - black_delta, linspace(1, clip/a, 1000) * a + b + black_delta); axis(ax);\n");
            fprintf(f, "print -dpng iso-curve-%d.png\n", iter);
            fclose(f);
            if(system("octave iso-curve.m"));
        }
    }

#if 0
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            raw_set_pixel16(x, y, dark[x + y*w] * a + b + black);
    save_debug_dng("fit-dark.dng");
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            raw_set_pixel16(x, y, bright[x + y*w] + black);
    save_debug_dng("fit-bright.dng");
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            raw_set_pixel16(x, y, 100 * (bright[x + y*w] - (dark[x + y*w] * a + b)) + white/2);
    save_debug_dng("fit-delta.dng");
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            raw_set_pixel16(x, y, 100 * (bright[x + y*w]/a - b/a - dark[x + y*w]) + white/2);
    save_debug_dng("fit-delta2.dng");
#endif

    free(dark);
    free(bright);

    /* apply the correction */

    /* previous code used different assumption: dark = a * bright + b */
    /* y = ax + b => x = y/a - b/a */
    a = 1/a;
    b *= a; b = -b;

    /*
     * We already know that EXIF black level is the average of the black levels
     * for bright and dark exposures, so their offset is symmetrical around the mean
     * 
     * Usually, the brighter exposure ends up having a lower black level.
     *
     * From fitting: dark-black = a * (bright-black) + b
     * Find black delta so that (dark-black-delta) = a * (bright-black+delta).
     */
    double b20 = b * 16;
    double black_delta20 = b20 / (a + 1);

    if (iso_data->freeze == 1){
        iso_data->freeze = 2;
        iso_data->a = a;
        iso_data->b = b;
    }
    
    apply_correction(a, b20, h, w, raw_info, black20, white20, raw_buffer_32, corr_ev, white_darkened, is_bright);

    return 1;
}

static inline uint32_t * convert_to_20bit(struct raw_info raw_info, uint16_t * image_data)
{
    int w = raw_info.width;
    int h = raw_info.height;
    /* promote from 14 to 20 bits (original raw buffer holds 14-bit values stored as uint16_t) */
    uint32_t * raw_buffer_32 = malloc(w * h * sizeof(raw_buffer_32[0]));
    
#pragma omp parallel for schedule(static) default(none) shared(raw_info,raw_buffer_32, h,w,image_data)
    for (int i = 0; i < w*h; i ++)
        raw_buffer_32[i] = ((uint32_t)image_data[i] << 6) & 0xFFFFF;
    
    return raw_buffer_32;
}

static inline void build_ev2raw_lut(int * raw2ev, int * ev2raw_0, int black, int white)
{
    int* ev2raw = ev2raw_0 + 10*EV_RESOLUTION;
    
#pragma omp parallel for schedule(static) default(none) shared(black, raw2ev)
    for (int i = 0; i < 1<<20; i++)
    {
        double signal = MAX(i/64.0 - black/64.0, -1023);
        if (signal > 0)
            raw2ev[i] = (int)round(log2(1+signal) * EV_RESOLUTION);
        else
            raw2ev[i] = -(int)round(log2(1-signal) * EV_RESOLUTION);
    }
    
#pragma omp parallel for schedule(static) default(none) shared(black, ev2raw)
    for (int i = -10*EV_RESOLUTION; i < 0; i++)
    {
        ev2raw[i] = COERCE(black+64 - round(64*pow(2, ((double)-i/EV_RESOLUTION))), 0, black);
    }
    
#pragma omp parallel for schedule(static) default(none) shared(black, ev2raw, raw2ev, white)
    for (int i = 0; i < 14*EV_RESOLUTION; i++)
    {
        ev2raw[i] = COERCE(black-64 + round(64*pow(2, ((double)i/EV_RESOLUTION))), black, (1<<20)-1);
        
        if (i >= raw2ev[white])
        {
            ev2raw[i] = MAX(ev2raw[i], white);
        }
    }
    
    /* keep "bad" pixels, if any */
    ev2raw[raw2ev[0]] = 0;
    ev2raw[raw2ev[0]] = 0;
    
    /* check raw <--> ev conversion */
    //~ printf("%d %d %d %d %d %d %d *%d* %d %d %d %d %d\n", raw2ev[0],         raw2ev[16000],         raw2ev[32000],         raw2ev[131068],         raw2ev[131069],         raw2ev[131070],         raw2ev[131071],         raw2ev[131072],         raw2ev[131073],         raw2ev[131074],         raw2ev[131075],         raw2ev[131076],         raw2ev[132000]);
    //~ printf("%d %d %d %d %d %d %d *%d* %d %d %d %d %d\n", ev2raw[raw2ev[0]], ev2raw[raw2ev[16000]], ev2raw[raw2ev[32000]], ev2raw[raw2ev[131068]], ev2raw[raw2ev[131069]], ev2raw[raw2ev[131070]], ev2raw[raw2ev[131071]], ev2raw[raw2ev[131072]], ev2raw[raw2ev[131073]], ev2raw[raw2ev[131074]], ev2raw[raw2ev[131075]], ev2raw[raw2ev[131076]], ev2raw[raw2ev[132000]]);
}

static inline double compute_noise20(struct raw_info raw_info, uint32_t * image_data, double * noise_std, double * dark_noise, double * bright_noise, double * dark_noise_ev, double * bright_noise_ev)
{
    // Fixme : We should use black bars to compute noise levels, but they are not present here
    double noise_avg = 0.0;
    for (int y = 0; y < 4; y++)
        compute_black_noise20(raw_info, image_data, raw_info.active_area.y1/4*4 + 20 + y, raw_info.active_area.y2 - 20, 1, 4, &noise_avg, &noise_std[y]);

#ifndef STDOUT_SILENT
    printf("Noise levels    : %.02f %.02f %.02f %.02f (14-bit)\n", noise_std[0], noise_std[1], noise_std[2], noise_std[3]);
#endif

    *dark_noise = MIN(MIN(noise_std[0], noise_std[1]), MIN(noise_std[2], noise_std[3]));
    *bright_noise = MAX(MAX(noise_std[0], noise_std[1]), MAX(noise_std[2], noise_std[3]));
    *dark_noise_ev = log2(*dark_noise);
    *bright_noise_ev = log2(*bright_noise);
    return noise_avg;
}

static inline double * build_fullres_curve(int black)
{
    /* fullres mixing curve */
    static double fullres_curve[1<<20];
    static int previous_black = -1;
    
    if(previous_black == black) return fullres_curve;
    
    previous_black = black;
       
#pragma omp parallel for schedule(static) default(none) shared(black,fullres_curve,fullres_start,fullres_transition)
    for (int i = 0; i < (1<<20); i++)
    {
        double ev2 = log2(MAX(i/64.0 - black/64.0, 1));
        double c2 = -cos(COERCE(ev2 - fullres_start, 0, fullres_transition)*M_PI/fullres_transition);
        double f = (c2+1) / 2;
        fullres_curve[i] = f;
    }
    
    return fullres_curve;
}

/* define edge directions for interpolation */
struct xy { int x; int y; };
const struct
{
    struct xy ack;      /* verification pixel near a */
    struct xy a;        /* interpolation pixel from the nearby line: normally (0,s) but also (1,s) or (-1,s) */
    struct xy b;        /* interpolation pixel from the other line: normally (0,-2s) but also (1,-2s), (-1,-2s), (2,-2s) or (-2,-2s) */
    struct xy bck;      /* verification pixel near b */
}
edge_directions[] = {       /* note: all y coords should be multiplied by s */
    //~ { {-6,2}, {-3,1}, { 6,-2}, { 9,-3} },     /* almost horizontal (little or no improvement) */
    { {-4,2}, {-2,1}, { 4,-2}, { 6,-3} },
    { {-3,2}, {-1,1}, { 3,-2}, { 4,-3} },
    { {-2,2}, {-1,1}, { 2,-2}, { 3,-3} },     /* 45-degree diagonal */
    { {-1,2}, {-1,1}, { 1,-2}, { 2,-3} },
    { {-1,2}, { 0,1}, { 1,-2}, { 1,-3} },
    { { 0,2}, { 0,1}, { 0,-2}, { 0,-3} },     /* vertical, preferred; no extra confirmations needed */
    { { 1,2}, { 0,1}, {-1,-2}, {-1,-3} },
    { { 1,2}, { 1,1}, {-1,-2}, {-2,-3} },
    { { 2,2}, { 1,1}, {-2,-2}, {-3,-3} },     /* 45-degree diagonal */
    { { 3,2}, { 1,1}, {-3,-2}, {-4,-3} },
    { { 4,2}, { 2,1}, {-4,-2}, {-6,-3} },
    //~ { { 6,2}, { 3,1}, {-6,-2}, {-9,-3} },     /* almost horizontal */
};

static inline int edge_interp(float ** plane, int * squeezed, int * raw2ev, int dir, int x, int y, int s)
{
    
    int dxa = edge_directions[dir].a.x;
    int dya = edge_directions[dir].a.y * s;
    int pa = COERCE((int)plane[squeezed[y+dya]][x+dxa], 0, 0xFFFFF);
    int dxb = edge_directions[dir].b.x;
    int dyb = edge_directions[dir].b.y * s;
    int pb = COERCE((int)plane[squeezed[y+dyb]][x+dxb], 0, 0xFFFFF);
    int pi = (raw2ev[pa] * 2 + raw2ev[pb]) / 3;
    
    return pi;
}

static inline void amaze_interpolate(struct raw_info raw_info, uint32_t * raw_buffer_32, uint32_t* dark, uint32_t* bright, int black, int white, int white_darkened, int * is_bright)
{
    int w = raw_info.width;
    int h = raw_info.height;
#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif
    int* squeezed = malloc(h * sizeof(int));
    memset(squeezed, 0, h * sizeof(int));
    
    float** rawData = malloc(h * sizeof(rawData[0]));
    float** red     = malloc(h * sizeof(red[0]));
    float** green   = malloc(h * sizeof(green[0]));
    float** blue    = malloc(h * sizeof(blue[0]));

    int wx = w + 16;
#pragma omp parallel for schedule(static) default(none) shared(h, w, rawData, red, green, blue, wx)
    for (int i = 0; i < h; i++)
    {
        rawData[i] =   malloc(wx * sizeof(rawData[0][0]));
        memset(rawData[i], 0, wx * sizeof(rawData[0][0]));
        red[i]     = malloc(wx * sizeof(red[0][0]));
        green[i]   = malloc(wx * sizeof(green[0][0]));
        blue[i]    = malloc(wx * sizeof(blue[0][0]));
    }
#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t memory allocation took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif

#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif

    /* squeeze the dark image by deleting fields from the bright exposure */
    int yh = -1;
    for (int y = 0; y < h; y ++)
    {
        if (BRIGHT_ROW)
            continue;
        
        if (yh < 0) /* make sure we start at the same parity (RGGB cell) */
            yh = y;
        
#pragma omp parallel for schedule(static) default(none) shared(is_bright, w,yh,y, rawData, raw_buffer_32, raw_info, black)
        for (int x = 0; x < w; x++)
        {
            int p = raw_get_pixel32(x, y);
            
            if (x%2 != y%2) /* divide green channel by 2 to approximate the final WB better */
                p = (p - black) / 2 + black;
            
            rawData[yh][x] = p;
        }
        
        squeezed[y] = yh;
        
        yh++;
    }
    
    /* now the same for the bright exposure */
    yh = -1;
    for (int y = 0; y < h; y ++)
    {
        if (!BRIGHT_ROW)
            continue;
        
        if (yh < 0) /* make sure we start with the same parity (RGGB cell) */
            yh = h/4*2 + y;
#pragma omp parallel for schedule(static) default(none) shared(is_bright, w,yh,y, rawData, raw_buffer_32, raw_info, black)
        for (int x = 0; x < w; x++)
        {
            int p = raw_get_pixel32(x, y);
            
            if (x%2 != y%2) /* divide green channel by 2 to approximate the final WB better */
                p = (p - black) / 2 + black;
            
            rawData[yh][x] = p;
        }
        
        squeezed[y] = yh;
        
        yh++;
        if (yh >= h) break; /* just in case */
    }

#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t squeeze dark/bright images took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#if 0
    void amaze_demosaic_RT(
                           float** rawData,    /* holds preprocessed pixel values, rawData[i][j] corresponds to the ith row and jth column */
                           float** red,        /* the interpolated red plane */
                           float** green,      /* the interpolated green plane */
                           float** blue,       /* the interpolated blue plane */
                           int winx, int winy, /* crop window for demosaicing */
                           int winw, int winh
                           );
#endif
    //IDK if AMaZE is actually thread safe, but I'm just going to assume not, rather than inspecting that huge mess of code
#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif

    int threads;
#pragma omp parallel
{
    #pragma omp master
    {
       threads = omp_get_num_threads();
    }
    //printf("number of threads = %d\n", threads /*omp_get_thread_num()*/);
}
    /* If threads is < 2 just do a normal amaze */
    if (threads < 2)
    {
        /* run the Amaze */
        demosaic(& (amazeinfo_t) {
                     rawData,
                     red,
                     green,
                     blue,
                     0, 0,
                     w, h,
                     0,
                     black /*0*/ }); //Should start with 0, or better with black?
    }
    /* Else do multithreading */
    else
    {
        int startchunk_y[threads];
        int endchunk_y[threads];

        /* How big each thread's chunk is, multiple of 2 - or debayer
         * would start on wrong pixel and magenta stripes appear */
        int chunk_height = h / threads;
        chunk_height -= chunk_height % 2;

        /* To small chunk heights bring AMaZE module to crash */
        while( chunk_height <= 32 )
        {
            if( threads <= 1 ) break;
            threads--;
            chunk_height = h / threads;
            chunk_height -= chunk_height % 2;
        }

        /* Calculate chunks of image for each thread */
        for (int thread = 0; thread < threads; ++thread)
        {
            startchunk_y[thread] = chunk_height * thread;
            endchunk_y[thread] = chunk_height * (thread + 1);
        }

        /* Last chunk must reach end of frame */
        endchunk_y[threads-1] = h;
        amazeinfo_t amaze_arguments[threads];

    #pragma omp parallel
    {
            int thread = omp_get_thread_num();

            /* Amaze arguments */
            amaze_arguments[thread] = (amazeinfo_t) {
                rawData,
                red,
                green,
                blue,
                /* Crop out a part for each thread */
                0, startchunk_y[thread],    /* crop window for demosaicing */
                w, (endchunk_y[thread] - startchunk_y[thread]),
                0,
                0 };
            /* Partial image for this thread */
            demosaic(&amaze_arguments[thread]);
    #pragma omp barrier
    }
}

    //demosaic(& (amazeinfo_t) { rawData, red, green, blue, 0, 0, w, h, 0, 0 });
#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t demosaic took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    
    /* undo green channel scaling and clamp the other channels */
#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif

#pragma omp parallel for schedule(static) default(none) shared(green, red, blue, h, w, black) collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            green[y][x] = COERCE((green[y][x] - black) * 2 + black, 0, 0xFFFFF);
            red[y][x] = COERCE(red[y][x], 0, 0xFFFFF);
            blue[y][x] = COERCE(blue[y][x], 0, 0xFFFFF);
        }
    }
#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t undo green channel scaling and clamp the other channels took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifndef STDOUT_SILENT
    printf("Edge-directed interpolation...\n");
#endif
#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif
    //~ printf("Grayscale...\n");
    /* convert to grayscale and de-squeeze for easier processing */
    uint32_t * gray = malloc(w * h * sizeof(gray[0]));

#pragma omp parallel for schedule(static) default(none) shared(gray, green, red, blue, h, w, squeezed) collapse(2)
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            gray[x + y*w] = green[squeezed[y]][x]/2 + red[squeezed[y]][x]/4 + blue[squeezed[y]][x]/4;
    
    
    uint8_t* edge_direction = malloc(w * h * sizeof(edge_direction[0]));
    int d0 = COUNT(edge_directions)/2;

#pragma omp parallel for schedule(static) default(none) shared(edge_direction, green, red, blue, h, w,d0) collapse(2)
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            edge_direction[x + y*w] = d0;

    double * fullres_curve = build_fullres_curve(black);
    
    //~ printf("Cross-correlation...\n");
    int semi_overexposed = 0;
    int not_overexposed = 0;
    int deep_shadow = 0;
    int not_shadow = 0;
    
    /* for fast EV - raw conversion */
    static int raw2ev[1<<20];   /* EV x EV_RESOLUTION */
    static int ev2raw_0[24*EV_RESOLUTION];
    static int previous_black = -1;
    
    /* handle sub-black values (negative EV) */
    int* ev2raw = ev2raw_0 + 10*EV_RESOLUTION;
    
    if(black != previous_black)
    {
        build_ev2raw_lut(raw2ev, ev2raw_0, black, white);
        previous_black = black;
    }

#pragma omp parallel for schedule(static) default(none) shared(edge_direction, raw2ev,gray,is_bright, white_darkened, h, w, raw_info, raw_buffer_32, black, white, d0, edge_directions, fullres_curve, fullres_thr) reduction(+:not_shadow) reduction(+:deep_shadow) reduction(+:not_overexposed) reduction(+:semi_overexposed)
    for (int y = 5; y < h-5; y ++)
    {
        int s = (is_bright[y%4] == is_bright[(y+1)%4]) ? -1 : 1;    /* points to the closest row having different exposure */
        for (int x = 5; x < w-5; x ++)
        {
            int e_best = INT_MAX;
            int d_best = d0;
            int dmin = 0;
            int dmax = COUNT(edge_directions)-1;
            int search_area = 5;

            /* only use high accuracy on the dark exposure where the bright ISO is overexposed */
            if (!BRIGHT_ROW)
            {
                /* interpolating bright exposure */
                if (fullres_curve[raw_get_pixel32(x, y)] > fullres_thr)
                {
                    /* no high accuracy needed, just interpolate vertically */
                    not_shadow++;
                    dmin = d0;
                    dmax = d0;
                }
                else
                {
                    /* deep shadows, unlikely to use fullres, so we need a good interpolation */
                    deep_shadow++;
                }
            }
            else if (raw_get_pixel32(x, y) < (unsigned int)white_darkened)
            {
                /* interpolating dark exposure, but we also have good data from the bright one */
                not_overexposed++;
                dmin = d0;
                dmax = d0;
            }
            else
            {
                /* interpolating dark exposure, but the bright one is clipped */
                semi_overexposed++;
            }

            if (dmin == dmax)
            {
                d_best = dmin;
            }
            else
            {
                for (int d = dmin; d <= dmax; d++)
                {
                    int e = 0;
                    for (int j = -search_area; j <= search_area; j++)
                    {
                        int dx1 = edge_directions[d].ack.x + j;
                        int dy1 = edge_directions[d].ack.y * s;
                        int p1 = raw2ev[gray[x+dx1 + (y+dy1)*w]];
                        int dx2 = edge_directions[d].a.x + j;
                        int dy2 = edge_directions[d].a.y * s;
                        int p2 = raw2ev[gray[x+dx2 + (y+dy2)*w]];
                        int dx3 = edge_directions[d].b.x + j;
                        int dy3 = edge_directions[d].b.y * s;
                        int p3 = raw2ev[gray[x+dx3 + (y+dy3)*w]];
                        int dx4 = edge_directions[d].bck.x + j;
                        int dy4 = edge_directions[d].bck.y * s;
                        int p4 = raw2ev[gray[x+dx4 + (y+dy4)*w]];
                        e += ABS(p1-p2) + ABS(p2-p3) + ABS(p3-p4);
                    }

                    /* add a small penalty for diagonal directions */
                    /* (the improvement should be significant in order to choose one of these) */
                    e += ABS(d - d0) * EV_RESOLUTION/8;

                    if (e < e_best)
                    {
                        e_best = e;
                        d_best = d;
                    }
                }
            }

            edge_direction[x + y*w] = d_best;
        }
    }

#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t Edge-directed interpolation took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifndef STDOUT_SILENT
    printf("Semi-overexposed: %.02f%%\n", semi_overexposed * 100.0 / (semi_overexposed + not_overexposed));
    printf("Deep shadows    : %.02f%%\n", deep_shadow * 100.0 / (deep_shadow + not_shadow));
#endif
    //~ printf("Actual interpolation...\n");
#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif

#pragma omp parallel for schedule(static) default(none) shared(w, red, green, edge_direction, squeezed, blue, bright,dark, is_bright, h, raw2ev,ev2raw,raw_info,raw_buffer_32)
    for (int y = 2; y < h-2; y ++)
    {
        uint32_t* native = BRIGHT_ROW ? bright : dark;
        uint32_t* interp = BRIGHT_ROW ? dark : bright;
        int is_rg = (y % 2 == 0); /* RG or GB? */
        int s = (is_bright[y%4] == is_bright[(y+1)%4]) ? -1 : 1;    /* points to the closest row having different exposure */

        //~ printf("Interpolating %s line %d from [near] %d (squeezed %d) and [far] %d (squeezed %d)\n", BRIGHT_ROW ? "BRIGHT" : "DARK", y, y+s, yh_near, y-2*s, yh_far);

        for (int x = 2; x < w-2; x += 2)
        {
            for (int k = 0; k < 2; k++, x++)
            {
                float** plane = is_rg ? (x%2 == 0 ? red   : green)
                : (x%2 == 0 ? green : blue );

                int dir = edge_direction[x + y*w];

                /* vary the interpolation direction and average the result (reduces aliasing) */
                int pi0 = edge_interp(plane, squeezed, raw2ev, dir, x, y, s);
                int pip = edge_interp(plane, squeezed, raw2ev, MIN(dir+1, COUNT(edge_directions)-1), x, y, s);
                int pim = edge_interp(plane, squeezed, raw2ev, MAX(dir-1,0), x, y, s);

                interp[x   + y * w] = ev2raw[(2*pi0+pip+pim)/4];
                native[x   + y * w] = raw_get_pixel32(x, y);
            }
            x -= 2;
        }
    }
#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t actual interpolation took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifdef PERF_INFO
    perf_sub_clock = clock();
#endif

#pragma omp parallel for schedule(static) default(none) shared(rawData, red, green, blue, h)
    for (int i = 0; i < h; i++)
    {
        free(rawData[i]);
        free(red[i]);
        free(green[i]);
        free(blue[i]);
    }
    
    free(squeezed); squeezed = 0;
    free(rawData); rawData = 0;
    free(red); red = 0;
    free(green); green = 0;
    free(blue); blue = 0;
    free(gray); gray = 0;
    free(edge_direction);
#ifdef PERF_INFO
    perf_sub_clock = clock()-perf_sub_clock;
    printf("\t free memory took %f seconds\n", ((double) perf_sub_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
}

static inline void mean23_interpolate(struct raw_info raw_info, uint32_t * raw_buffer_32, uint32_t* dark, uint32_t* bright, int black, int white, int white_darkened, int * is_bright)
{
    int w = raw_info.width;
    int h = raw_info.height;
#ifndef STDOUT_SILENT
    printf("Interpolation   : mean23\n");
#endif
    /* for fast EV - raw conversion */
    static int raw2ev[1<<20];   /* EV x EV_RESOLUTION */
    static int ev2raw_0[24*EV_RESOLUTION];
    static int previous_black = -1;
    
    /* handle sub-black values (negative EV) */
    int* ev2raw = ev2raw_0 + 10*EV_RESOLUTION;
    
    if(black != previous_black)
    {
        build_ev2raw_lut(raw2ev, ev2raw_0, black, white);
        previous_black = black;
    }

#pragma omp parallel for schedule(static) default(none) shared(is_bright, dark, bright, white_darkened, h, w, raw_info, raw_buffer_32, black, white, raw2ev, ev2raw)
    for (int y = 2; y < h-2; y ++)
    {
        uint32_t* native = BRIGHT_ROW ? bright : dark;
        uint32_t* interp = BRIGHT_ROW ? dark : bright;
        int is_rg = (y % 2 == 0); /* RG or GB? */
        int white = !BRIGHT_ROW ? white_darkened : raw_info.white_level;

        for (int x = 2; x < w-3; x += 2)
        {

            /* red/blue: interpolate from (x,y+2) and (x,y-2) */
            /* green: interpolate from (x+1,y+1),(x-1,y+1),(x,y-2) or (x+1,y-1),(x-1,y-1),(x,y+2), whichever has the correct brightness */

            int s = (is_bright[y%4] == is_bright[(y+1)%4]) ? -1 : 1;

            if (is_rg)
            {
                int ra = raw_get_pixel32(x, y-2);
                int rb = raw_get_pixel32(x, y+2);
                int ri = mean2(raw2ev[ra], raw2ev[rb], raw2ev[white], 0);

                int ga = raw_get_pixel32(x+1+1, y+s);
                int gb = raw_get_pixel32(x+1-1, y+s);
                int gc = raw_get_pixel32(x+1, y-2*s);
                int gi = mean3(raw2ev[ga], raw2ev[gb], raw2ev[gc], raw2ev[white], 0);

                interp[x   + y * w] = ev2raw[ri];
                interp[x+1 + y * w] = ev2raw[gi];
            }
            else
            {
                int ba = raw_get_pixel32(x+1  , y-2);
                int bb = raw_get_pixel32(x+1  , y+2);
                int bi = mean2(raw2ev[ba], raw2ev[bb], raw2ev[white], 0);

                int ga = raw_get_pixel32(x+1, y+s);
                int gb = raw_get_pixel32(x-1, y+s);
                int gc = raw_get_pixel32(x, y-2*s);
                int gi = mean3(raw2ev[ga], raw2ev[gb], raw2ev[gc], raw2ev[white], 0);

                interp[x   + y * w] = ev2raw[gi];
                interp[x+1 + y * w] = ev2raw[bi];
            }

            native[x   + y * w] = raw_get_pixel32(x, y);
            native[x+1 + y * w] = raw_get_pixel32(x+1, y);
        }
    }

}

static inline void border_interpolate(struct raw_info raw_info, uint32_t * raw_buffer_32, uint32_t* dark, uint32_t* bright, int * is_bright)
{
    int w = raw_info.width;
    int h = raw_info.height;
    
    /* border interpolation */
    for (int y = 0; y < 3; y ++)
    {
        uint32_t* native = BRIGHT_ROW ? bright : dark;
        uint32_t* interp = BRIGHT_ROW ? dark : bright;
        
#pragma omp parallel for schedule(static) default(none) shared(raw_info, raw_buffer_32, interp, native, w, y)
        for (int x = 0; x < w; x ++)
        {
            interp[x + y * w] = raw_get_pixel32(x, y+2);
            native[x + y * w] = raw_get_pixel32(x, y);
        }
    }

    for (int y = h-4; y < h; y ++)
    {
        uint32_t* native = BRIGHT_ROW ? bright : dark;
        uint32_t* interp = BRIGHT_ROW ? dark : bright;
#pragma omp parallel for schedule(static) default(none) shared(raw_info, raw_buffer_32, interp, native, w, y)
        for (int x = 0; x < w; x ++)
        {
            interp[x + y * w] = raw_get_pixel32(x, y-2);
            native[x + y * w] = raw_get_pixel32(x, y);
        }
    }

#pragma omp parallel for schedule(static) default(none) shared(raw_info, raw_buffer_32, w, h, bright, dark, is_bright)
    for (int y = 2; y < h; y ++)
    {
        uint32_t* native = BRIGHT_ROW ? bright : dark;
        uint32_t* interp = BRIGHT_ROW ? dark : bright;
        
        for (int x = 0; x < 2; x ++)
        {
            interp[x + y * w] = raw_get_pixel32(x, y-2);
            native[x + y * w] = raw_get_pixel32(x, y);
        }
        
        for (int x = w-3; x < w; x ++)
        {
            interp[x + y * w] = raw_get_pixel32(x-2, y-2);
            native[x + y * w] = raw_get_pixel32(x-2, y);
        }
    }
}

static inline void fullres_reconstruction(struct raw_info raw_info, uint32_t * fullres, uint32_t* dark, uint32_t* bright, uint32_t white_darkened, int * is_bright, int dark_highlight_threshold)
{
    int w = raw_info.width;
    int h = raw_info.height;
    
    int reversed_dark_highlight_threshold = dark_highlight_threshold > 1048512 ? 0 : 1048512 - dark_highlight_threshold;
    /* reconstruct a full-resolution image (discard interpolated fields whenever possible) */
    /* this has full detail and lowest possible aliasing, but it has high shadow noise and color artifacts when high-iso starts clipping */
#ifndef STDOUT_SILENT
    printf("Full-res reconstruction...\n");
#endif
#pragma omp parallel for schedule(static) default(none) shared(reversed_dark_highlight_threshold, fullres, is_bright, dark, bright, white_darkened, h, w) collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            uint32_t d = dark[x + y*w];
            if (BRIGHT_ROW)
            {
                uint32_t f = bright[x + y*w];
                if (reversed_dark_highlight_threshold && d > f && (int)d > reversed_dark_highlight_threshold){
                    fullres[x + y*w] = d;
                }else{
                    /* if the brighter copy is overexposed, the guessed pixel for sure has higher brightness */
                    fullres[x + y*w] = f < (uint32_t)white_darkened ? f : MAX(f, d);
                }
            }
            else
            {
                fullres[x + y*w] = d;
            }
        }
    }
}

static inline void build_alias_map(struct raw_info raw_info, uint16_t* alias_map, uint32_t* fullres_smooth, uint32_t* halfres_smooth, uint32_t* bright, int dark_noise, int black, int * raw2ev/*, double *fullres_curve*/)
{
    if(!alias_map) return;
    
    int w = raw_info.width;
    int h = raw_info.height;
    
    double * fullres_curve = build_fullres_curve(black);
#ifndef STDOUT_SILENT
    printf("Building alias map...\n");
#endif
    uint16_t* alias_aux = malloc(w * h * sizeof(uint16_t));
    
    /* build the aliasing maps (where it's likely to get aliasing) */
    /* do this by comparing fullres and halfres images */
    /* if the difference is small, we'll prefer halfres for less noise, otherwise fullres for less aliasing */
#pragma omp parallel for schedule(static) default(none) shared(fullres_curve, fullres_smooth, halfres_smooth, h, w, fullres_thr, dark_noise, bright, raw2ev, alias_map) collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            /* do not compute alias map where we'll use fullres detail anyway */
            if (fullres_curve[bright[x + y*w]] > fullres_thr)
                continue;
            
            int f = fullres_smooth[x + y*w];
            int h = halfres_smooth[x + y*w];
            int fe = raw2ev[f];
            int he = raw2ev[h];
            int e_lin = ABS(f - h); /* error in linear space, for shadows (downweights noise) */
            e_lin = MAX(e_lin - dark_noise*3/2, 0);
            int e_log = ABS(fe - he); /* error in EV space, for highlights (highly sensitive to noise) */
            alias_map[x + y*w] = MIN(MIN(e_lin/2, e_log/16), 65530);
        }
    }
    
    memcpy(alias_aux, alias_map, w * h * sizeof(uint16_t));
#ifndef STDOUT_SILENT
    printf("Filtering alias map...\n");
#endif
#pragma omp parallel for schedule(static) default(none) shared(fullres_curve, h, w, fullres_thr, bright, alias_map, alias_aux) collapse(2)
    for (int y = 6; y < h-6; y ++)
    {
        for (int x = 6; x < w-6; x ++)
        {
            /* do not compute alias map where we'll use fullres detail anyway */
            if (fullres_curve[bright[x + y*w]] > fullres_thr)
                continue;
            
            /* use 5th max (out of 37) to filter isolated pixels */
            int neighbours[] = {
                                                                              -alias_map[x-2 + (y-6) * w], -alias_map[x+0 + (y-6) * w], -alias_map[x+2 + (y-6) * w],
                                                 -alias_map[x-4 + (y-4) * w], -alias_map[x-2 + (y-4) * w], -alias_map[x+0 + (y-4) * w], -alias_map[x+2 + (y-4) * w], -alias_map[x+4 + (y-4) * w],
                    -alias_map[x-6 + (y-2) * w], -alias_map[x-4 + (y-2) * w], -alias_map[x-2 + (y-2) * w], -alias_map[x+0 + (y-2) * w], -alias_map[x+2 + (y-2) * w], -alias_map[x+4 + (y-2) * w], -alias_map[x+6 + (y-2) * w], 
                    -alias_map[x-6 + (y+0) * w], -alias_map[x-4 + (y+0) * w], -alias_map[x-2 + (y+0) * w], -alias_map[x+0 + (y+0) * w], -alias_map[x+2 + (y+0) * w], -alias_map[x+4 + (y+0) * w], -alias_map[x+6 + (y+0) * w], 
                    -alias_map[x-6 + (y+2) * w], -alias_map[x-4 + (y+2) * w], -alias_map[x-2 + (y+2) * w], -alias_map[x+0 + (y+2) * w], -alias_map[x+2 + (y+2) * w], -alias_map[x+4 + (y+2) * w], -alias_map[x+6 + (y+2) * w], 
                                                 -alias_map[x-4 + (y+4) * w], -alias_map[x-2 + (y+4) * w], -alias_map[x+0 + (y+4) * w], -alias_map[x+2 + (y+4) * w], -alias_map[x+4 + (y+4) * w],
                                                                              -alias_map[x-2 + (y+6) * w], -alias_map[x+0 + (y+6) * w], -alias_map[x+2 + (y+6) * w],
            };
            alias_aux[x + y * w] = -kth_smallest_int(neighbours, COUNT(neighbours), 5);
        }
    }
#ifndef STDOUT_SILENT
    printf("Smoothing alias map...\n");
#endif
    /* gaussian blur */
#pragma omp parallel for schedule(static) default(none) shared(fullres_curve, h, w, fullres_thr, bright, alias_map, alias_aux) collapse(2)
    for (int y = 6; y < h-6; y ++)
    {
        for (int x = 6; x < w-6; x ++)
        {
            /* do not compute alias map where we'll use fullres detail anyway */
            if (fullres_curve[bright[x + y*w]] > fullres_thr)
                continue;
            
            int c =
            (alias_aux[x+0 + (y+0) * w])+
            (alias_aux[x+0 + (y-2) * w] + alias_aux[x-2 + (y+0) * w] + alias_aux[x+2 + (y+0) * w] + alias_aux[x+0 + (y+2) * w]) * 820 / 1024 +
            (alias_aux[x-2 + (y-2) * w] + alias_aux[x+2 + (y-2) * w] + alias_aux[x-2 + (y+2) * w] + alias_aux[x+2 + (y+2) * w]) * 657 / 1024 +
            (alias_aux[x+0 + (y-2) * w] + alias_aux[x-2 + (y+0) * w] + alias_aux[x+2 + (y+0) * w] + alias_aux[x+0 + (y+2) * w]) * 421 / 1024 +
            (alias_aux[x-2 + (y-2) * w] + alias_aux[x+2 + (y-2) * w] + alias_aux[x-2 + (y-2) * w] + alias_aux[x+2 + (y-2) * w] + alias_aux[x-2 + (y+2) * w] + alias_aux[x+2 + (y+2) * w] + alias_aux[x-2 + (y+2) * w] + alias_aux[x+2 + (y+2) * w]) * 337 / 1024 +
            (alias_aux[x-2 + (y-2) * w] + alias_aux[x+2 + (y-2) * w] + alias_aux[x-2 + (y+2) * w] + alias_aux[x+2 + (y+2) * w]) * 173 / 1024 +
            (alias_aux[x+0 + (y-6) * w] + alias_aux[x-6 + (y+0) * w] + alias_aux[x+6 + (y+0) * w] + alias_aux[x+0 + (y+6) * w]) * 139 / 1024 +
            (alias_aux[x-2 + (y-6) * w] + alias_aux[x+2 + (y-6) * w] + alias_aux[x-6 + (y-2) * w] + alias_aux[x+6 + (y-2) * w] + alias_aux[x-6 + (y+2) * w] + alias_aux[x+6 + (y+2) * w] + alias_aux[x-2 + (y+6) * w] + alias_aux[x+2 + (y+6) * w]) * 111 / 1024 +
            (alias_aux[x-2 + (y-6) * w] + alias_aux[x+2 + (y-6) * w] + alias_aux[x-6 + (y-2) * w] + alias_aux[x+6 + (y-2) * w] + alias_aux[x-6 + (y+2) * w] + alias_aux[x+6 + (y+2) * w] + alias_aux[x-2 + (y+6) * w] + alias_aux[x+2 + (y+6) * w]) * 57 / 1024;
            alias_map[x + y * w] = c;
        }
    }
    
    /* make it grayscale */
#pragma omp parallel for schedule(static) collapse(2)
    for (int y = 2; y < h-2; y += 2)
    {
        for (int x = 2; x < w-2; x += 2)
        {
            int a = alias_map[x   +     y * w];
            int b = alias_map[x+1 +     y * w];
            int c = alias_map[x   + (y+1) * w];
            int d = alias_map[x+1 + (y+1) * w];
            int C = MAX(MAX(a,b), MAX(c,d));
            
            C = MIN(C, ALIAS_MAP_MAX);
            
            alias_map[x   +     y * w] =
            alias_map[x+1 +     y * w] =
            alias_map[x   + (y+1) * w] =
            alias_map[x+1 + (y+1) * w] = C;
        }
    }
    
    free(alias_aux);
}

#define CHROMA_SMOOTH_TYPE uint32_t

#define CHROMA_SMOOTH_2X2
#include "chroma_smooth.c"
#undef CHROMA_SMOOTH_2X2

#define CHROMA_SMOOTH_3X3
#include "chroma_smooth.c"
#undef CHROMA_SMOOTH_3X3

#define CHROMA_SMOOTH_5X5
#include "chroma_smooth.c"
#undef CHROMA_SMOOTH_5X5

static inline void hdr_chroma_smooth(struct raw_info raw_info, uint32_t * input, uint32_t * output, int method, int * raw2ev, int * ev2raw)
{
    int w = raw_info.width;
    int h = raw_info.height;
    int black = raw_info.black_level;
    int white = raw_info.white_level;
    
    switch (method) {
        case 2:
            chroma_smooth_2x2(w, h, input, output, raw2ev, ev2raw, black, white);
            break;
        case 3:
            chroma_smooth_3x3(w, h, input, output, raw2ev, ev2raw, black, white);
            break;
        case 5:
            chroma_smooth_5x5(w, h, input, output, raw2ev, ev2raw, black, white);
            break;
            
        default:
#ifndef STDOUT_SILENT
            err_printf("Unsupported chroma smooth method\n");
#endif
            break;
    }
}

static inline int mix_images(struct raw_info raw_info, uint32_t* fullres, uint32_t* fullres_smooth, uint32_t* halfres, uint32_t* halfres_smooth, uint16_t* alias_map, uint32_t* dark, uint32_t* bright, uint16_t * overexposed, int dark_noise, uint32_t white_darkened, double corr_ev, double lowiso_dr, uint32_t black, uint32_t white, int chroma_smooth_method, int use_fullres)
{
    int w = raw_info.width;
    int h = raw_info.height;
    
    /* mix the two images */
    /* highlights:  keep data from dark image only */
    /* shadows:     keep data from bright image only */
    /* midtones:    mix data from both, to bring back the resolution */
    
    /* estimate ISO overlap */
    /*
     ISO 100:       ###...........  (11 stops)
     ISO 1600:  ####..........      (10 stops)
     Combined:  XX##..............  (14 stops)
     */
    double clipped_ev = corr_ev;
    double overlap = lowiso_dr - clipped_ev;
    
    /* you get better colors, less noise, but a little more jagged edges if we underestimate the overlap amount */
    /* maybe expose a tuning factor? (preference towards resolution or colors) */
    overlap -= MIN(3, overlap - 3);
#ifndef STDOUT_SILENT
    printf("ISO overlap     : %.1f EV (approx)\n", overlap);
#endif
    if (overlap < 0.5)
    {
#ifndef STDOUT_SILENT
        printf("Overlap error\n");
#endif
        return 0;
    }
    else if (overlap < 2)
    {
#ifndef STDOUT_SILENT
        printf("Overlap too small, use a smaller ISO difference for better results.\n");
#endif
    }
#ifndef STDOUT_SILENT
    printf("Half-res blending...\n");
#endif
    /* mixing curve */
    double max_ev = log2(white/64 - black/64);
    double * mix_curve = malloc((1<<20) * sizeof(double));
    
#pragma omp parallel for schedule(static) default(none) shared(corr_ev, max_ev, overlap, black, mix_curve)
    for (int i = 0; i < 1<<20; i++)
    {
        double ev = log2(MAX(i/64.0 - black/64.0, 1)) + corr_ev;
        double c = -cos(MAX(MIN(ev-(max_ev-overlap),overlap),0)*M_PI/overlap);
        double k = (c+1) / 2;
        mix_curve[i] = k;
    }

    /* for fast EV - raw conversion */
    static int raw2ev[1<<20];   /* EV x EV_RESOLUTION */
    static int ev2raw_0[24*EV_RESOLUTION];
    static uint32_t previous_black = -1;

    /* handle sub-black values (negative EV) */
    int* ev2raw = ev2raw_0 + 10*EV_RESOLUTION;

    if(black != previous_black)
    {
        build_ev2raw_lut(raw2ev, ev2raw_0, black, white);
        previous_black = black;
    }

#pragma omp parallel for schedule(static) default(none) shared(bright, dark, h, w, raw2ev, ev2raw, mix_curve, halfres) collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            /* bright and dark source pixels  */
            /* they may be real or interpolated */
            /* they both have the same brightness (they were adjusted before this loop), so we are ready to mix them */
            int b = bright[x + y*w];
            int d = dark[x + y*w];

            /* go from linear to EV space */
            int bev = raw2ev[b];
            int dev = raw2ev[d];

            /* blending factor */
            double k = COERCE(mix_curve[b & 0xFFFFF], 0, 1);

            /* mix bright and dark exposures */
            int mixed = bev * (1-k) + dev * k;
            halfres[x + y*w] = ev2raw[mixed];
        }
    }

    if (chroma_smooth_method)
    {
        if (use_fullres)
        {
            memcpy(fullres_smooth, fullres, w * h * sizeof(uint32_t));
            hdr_chroma_smooth(raw_info, fullres, fullres_smooth, chroma_smooth_method, raw2ev, ev2raw);

        }
        memcpy(halfres_smooth, halfres, w * h * sizeof(uint32_t));
        hdr_chroma_smooth(raw_info, halfres, halfres_smooth, chroma_smooth_method, raw2ev, ev2raw);
    }

    if(alias_map)
    {
        build_alias_map(raw_info, alias_map, fullres_smooth, halfres_smooth, bright, dark_noise, black, raw2ev/*, fullres_curve*/);
    }
    
#pragma omp parallel for schedule(static) default(none) shared(bright, dark, white, overexposed, h, w, white_darkened) collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            overexposed[x + y * w] = bright[x + y * w] >= white_darkened || dark[x + y * w] >= white ? 100 : 0;
        }
    }
    
    /* "blur" the overexposed map */
    uint16_t* over_aux = malloc(w * h * sizeof(uint16_t));
    memcpy(over_aux, overexposed, w * h * sizeof(uint16_t));
    
#pragma omp parallel for schedule(static) default(none) shared(overexposed, h, w, over_aux) collapse(2)
    for (int y = 3; y < h-3; y ++)
    {
        for (int x = 3; x < w-3; x ++)
        {
            overexposed[x + y * w] =
            (over_aux[x+0 + (y+0) * w])+
            (over_aux[x+0 + (y-1) * w] + over_aux[x-1 + (y+0) * w] + over_aux[x+1 + (y+0) * w] + over_aux[x+0 + (y+1) * w]) * 820 / 1024 +
            (over_aux[x-1 + (y-1) * w] + over_aux[x+1 + (y-1) * w] + over_aux[x-1 + (y+1) * w] + over_aux[x+1 + (y+1) * w]) * 657 / 1024 +
            //~ (over_aux[x+0 + (y-2) * w] + over_aux[x-2 + (y+0) * w] + over_aux[x+2 + (y+0) * w] + over_aux[x+0 + (y+2) * w]) * 421 / 1024 +
            //~ (over_aux[x-1 + (y-2) * w] + over_aux[x+1 + (y-2) * w] + over_aux[x-2 + (y-1) * w] + over_aux[x+2 + (y-1) * w] + over_aux[x-2 + (y+1) * w] + over_aux[x+2 + (y+1) * w] + over_aux[x-1 + (y+2) * w] + over_aux[x+1 + (y+2) * w]) * 337 / 1024 +
            //~ (over_aux[x-2 + (y-2) * w] + over_aux[x+2 + (y-2) * w] + over_aux[x-2 + (y+2) * w] + over_aux[x+2 + (y+2) * w]) * 173 / 1024 +
            //~ (over_aux[x+0 + (y-3) * w] + over_aux[x-3 + (y+0) * w] + over_aux[x+3 + (y+0) * w] + over_aux[x+0 + (y+3) * w]) * 139 / 1024 +
            //~ (over_aux[x-1 + (y-3) * w] + over_aux[x+1 + (y-3) * w] + over_aux[x-3 + (y-1) * w] + over_aux[x+3 + (y-1) * w] + over_aux[x-3 + (y+1) * w] + over_aux[x+3 + (y+1) * w] + over_aux[x-1 + (y+3) * w] + over_aux[x+1 + (y+3) * w]) * 111 / 1024 +
            //~ (over_aux[x-2 + (y-3) * w] + over_aux[x+2 + (y-3) * w] + over_aux[x-3 + (y-2) * w] + over_aux[x+3 + (y-2) * w] + over_aux[x-3 + (y+2) * w] + over_aux[x+3 + (y+2) * w] + over_aux[x-2 + (y+3) * w] + over_aux[x+2 + (y+3) * w]) * 57 / 1024;
            0;
        }
    }
    
    free(over_aux); over_aux = 0;
    free(mix_curve);

    return 1;
}

static inline void final_blend(struct raw_info raw_info, uint32_t* raw_buffer_32, uint32_t* fullres, uint32_t* fullres_smooth, uint32_t* halfres_smooth, uint32_t* dark, uint32_t* bright, uint16_t* overexposed, uint16_t* alias_map, int black, int white, int dark_noise, int use_fullres)
{
    /* fullres mixing curve */
    double * fullres_curve = build_fullres_curve(black);
    
    int w = raw_info.width;
    int h = raw_info.height;
    
    /* for fast EV - raw conversion */
    static int raw2ev[1<<20];   /* EV x EV_RESOLUTION */
    static int ev2raw_0[24*EV_RESOLUTION];
    static int previous_black = -1;
    
    /* handle sub-black values (negative EV) */
    int* ev2raw = ev2raw_0 + 10*EV_RESOLUTION;
    
    if(black != previous_black)
    {
        build_ev2raw_lut(raw2ev, ev2raw_0, black, white);
        previous_black = black;
    }
#ifndef STDOUT_SILENT
        printf("Final blending...\n");
#endif

#pragma omp parallel for schedule(static)  collapse(2)
    for (int y = 0; y < h; y ++)
    {
        for (int x = 0; x < w; x ++)
        {
            /* high-iso image (for measuring signal level) */
            int b = bright[x + y*w];

            /* half-res image (interpolated and chroma filtered, best for low-contrast shadows) */
            int hr = halfres_smooth[x + y*w];

            /* full-res image (non-interpolated, except where one ISO is blown out) */
            int fr = fullres[x + y*w];

            /* full res with some smoothing applied to hide aliasing artifacts */
            int frs = fullres_smooth[x + y*w];

            /* go from linear to EV space */
            int hrev = raw2ev[hr];
            int frev = raw2ev[fr];
            int frsev = raw2ev[frs];

            int output = hrev;

            /* blending factor */
            if (use_fullres)
            {
                double f = fullres_curve[b & 0xFFFFF];

                double c = 0;

                if (alias_map)
                {
                    int co = alias_map[x + y*w];
                    c = COERCE(co / (double) ALIAS_MAP_MAX, 0, 1);
                }

                double ovf = COERCE(overexposed[x + y*w] / 200.0, 0, 1);
                c = MAX(c, ovf);

                double noisy_or_overexposed = MAX(ovf, 1-f);

                /* use data from both ISOs in high-detail areas, even if it's noisier (less aliasing) */
                f = MAX(f, c);

                /* use smoothing in noisy near-overexposed areas to hide color artifacts */
                double fev = noisy_or_overexposed * frsev + (1-noisy_or_overexposed) * frev;

                /* limit the use of fullres in dark areas (fixes some black spots, but may increase aliasing) */
                int sig = (dark[x + y*w] + bright[x + y*w]) / 2;
                f = MAX(0, MIN(f, (double)(sig - black) / (4*dark_noise)));

                /* blend "half-res" and "full-res" images smoothly to avoid banding*/
                output = hrev * (1-f) + fev * f;

                /* show full-res map (for debugging) */
                //~ output = f * 14*EV_RESOLUTION;

                /* show alias map (for debugging) */
                //~ output = c * 14*EV_RESOLUTION;

                //~ output = hotpixel[x+y*w] ? 14*EV_RESOLUTION : 0;
                //~ output = raw2ev[dark[x+y*w]];
                /* safeguard */
                output = COERCE(output, -10*EV_RESOLUTION, 14*EV_RESOLUTION-1);
            }

            /* back to linear space and commit */
            raw_set_pixel32(x, y, ev2raw[output]);
        }
    }
}


static inline void convert_20_to_16bit(struct raw_info raw_info, uint16_t * image_data, uint32_t * raw_buffer_32)
{
    int w = raw_info.width;
    int h = raw_info.height;
    /* go back from 20-bit to 16-bit output */
    //raw_info.buffer = raw_buffer_16;
    raw_info.black_level /= 16;
    raw_info.white_level /= 16;
    
#pragma omp parallel for schedule(static) default(none) shared(h,w,image_data, raw_buffer_32, raw_info)
    for (int i = 0; i < w * h; i++)
        raw_set_pixel_20to16_rand(i, raw_buffer_32[i]);
}

/**
 * Fix vertical stripes (banding) from 5D Mark III (and maybe others).
 *
 * These stripes are periodic, they repeat every 8 pixels.
 * It looks like some columns have different luma amplification;
 * correction factors are somewhere around 0.98 - 1.02, maybe camera-specific, maybe depends on
 * certain settings, I have no idea. So, this fix compares luma values within one pixel block,
 * computes the correction factors (using median to reject outliers) and decides
 * whether to apply the correction or not.
 *
 * For speed reasons:
 * - Correction factors are computed from the first frame only.
 * - Only channels with error greater than 0.2% are corrected.
 */
#define FIXP_ONE 65536
#define FIXP_RANGE 65536

static int stripes_coeffs[8] = {0};
static int stripes_correction_needed = 0;

#define PA ((int)(p->a))
#define PB ((int)(p->b_lo | (p->b_hi << 12)))
#define PC ((int)(p->c_lo | (p->c_hi << 10)))
#define PD ((int)(p->d_lo | (p->d_hi << 8)))
#define PE ((int)(p->e_lo | (p->e_hi << 6)))
#define PF ((int)(p->f_lo | (p->f_hi << 4)))
#define PG ((int)(p->g_lo | (p->g_hi << 2)))
#define PH ((int)(p->h))

#define SET_PA(x) { int v = (x); p->a = v; }
#define SET_PB(x) { int v = (x); p->b_lo = v; p->b_hi = v >> 12; }
#define SET_PC(x) { int v = (x); p->c_lo = v; p->c_hi = v >> 10; }
#define SET_PD(x) { int v = (x); p->d_lo = v; p->d_hi = v >> 8; }
#define SET_PE(x) { int v = (x); p->e_lo = v; p->e_hi = v >> 6; }
#define SET_PF(x) { int v = (x); p->f_lo = v; p->f_hi = v >> 4; }
#define SET_PG(x) { int v = (x); p->g_lo = v; p->g_hi = v >> 2; }
#define SET_PH(x) { int v = (x); p->h = v; }

#define RAW_MUL(p, x) ((((int)(p) - raw_info.black_level) * (int)(x) / FIXP_ONE) + raw_info.black_level)
#define F2H(ev) COERCE((int)(FIXP_RANGE/2 + ev * FIXP_RANGE/2), 0, FIXP_RANGE-1)
#define H2F(x) ((double)((x) - FIXP_RANGE/2) / (FIXP_RANGE/2))


static void add_pixel(int hist[8][FIXP_RANGE], int num[8], int offset, int pa, int pb, struct raw_info raw_info)
{
    int a = pa;
    int b = pb;

    if (MIN(a,b) < 32)
        return; /* too noisy */

    if (MAX(a,b) > raw_info.white_level / 1.1)
        return; /* too bright */

    /**
     * compute correction factor for b, that makes it as bright as a
     *
     * first, work around quantization error (which causes huge spikes on histogram)
     * by adding a small random noise component
     * e.g. if raw value is 13, add some uniformly distributed noise,
     * so the value will be between -12.5 and 13.5.
     *
     * this removes spikes on the histogram, thus canceling bias towards "round" values
     */
    double af = a + (rand() % 1024) / 1024.0 - 0.5;
    double bf = b + (rand() % 1024) / 1024.0 - 0.5;
    double factor = af / bf;
    double ev = log2(factor);

    /**
     * add to histogram (for computing the median)
     */
    int weight = log2(a);
    hist[offset][F2H(ev)] += weight;
    num[offset] += weight;
}


static void detect_vertical_stripes_coeffs(struct raw_info raw_info, uint32_t * image_data, int force_correction)
{
    static int hist[8][FIXP_RANGE];
    static int num[8];

    memset(hist, 0, sizeof(hist));
    memset(num, 0, sizeof(num));

    /* compute 7 histograms: b./a, c./a ... h./a */
    /* that is, adjust all columns to make them as bright as a */
    /* process green pixels only, assuming the image is RGGB */

    typedef raw_pixblock_14 raw_pixblock; //TODO: should we change depending on the source bits to raw_pixblock_12 or raw_pixblock_10?
    raw_pixblock * row;
    for (row = (raw_pixblock*)image_data; (void*)row < (void*)image_data + raw_info.pitch * raw_info.height; row += 2 * raw_info.pitch / sizeof(raw_pixblock))
    {
        /* first line is RG */
        raw_pixblock * rg;
        for (rg = row; (void*)rg < (void*)row + raw_info.pitch - sizeof(raw_pixblock); rg++)
        {
            /* next line is GB */
            raw_pixblock * gb = rg + raw_info.pitch / sizeof(raw_pixblock);

            raw_pixblock * p = rg;
            int pb = PB - raw_info.black_level;
            int pd = PD - raw_info.black_level;
            int pf = PF - raw_info.black_level;
            int ph = PH - raw_info.black_level;
            p++;
            int pb2 = PB - raw_info.black_level;
            int pd2 = PD - raw_info.black_level;
            int pf2 = PF - raw_info.black_level;
            int ph2 = PH - raw_info.black_level;
            p = gb;
            //int pa = PA - raw_info.black_level;
            int pc = PC - raw_info.black_level;
            int pe = PE - raw_info.black_level;
            int pg = PG - raw_info.black_level;
            p++;
            int pa2 = PA - raw_info.black_level;
            int pc2 = PC - raw_info.black_level;
            int pe2 = PE - raw_info.black_level;
            int pg2 = PG - raw_info.black_level;

            /**
             * verification: introducing strong banding in one column
             * should not affect the coefficients from the other columns
             **/

            //~ pe = pe * 1.1;
            //~ pe2 = pe2 * 1.1;

            /**
             * Make all columns as bright as a2
             * use linear interpolation, so when processing column b, for example,
             * let bi = (b * 1 + b2 * 7) / (7+1)
             * let ei = (e * 4 + e2 * 4) / (4+4)
             * and so on, to avoid getting tricked by smooth gradients.
             */

            add_pixel(hist, num, 1, pa2, (pb * 1 + pb2 * 7) / 8, raw_info);
            add_pixel(hist, num, 2, pa2, (pc * 2 + pc2 * 6) / 8, raw_info);
            add_pixel(hist, num, 3, pa2, (pd * 3 + pd2 * 5) / 8, raw_info);
            add_pixel(hist, num, 4, pa2, (pe * 4 + pe2 * 4) / 8, raw_info);
            add_pixel(hist, num, 5, pa2, (pf * 5 + pf2 * 3) / 8, raw_info);
            add_pixel(hist, num, 6, pa2, (pg * 6 + pg2 * 2) / 8, raw_info);
            add_pixel(hist, num, 7, pa2, (ph * 7 + ph2 * 1) / 8, raw_info);
        }
    }

    int j,k;

    int max[8] = {0};
    for (j = 0; j < 8; j++)
        for (k = 1; k < FIXP_RANGE-1; k++)
            max[j] = MAX(max[j], hist[j][k]);

    /* compute the median correction factor (this will reject outliers) */
    for (j = 0; j < 8; j++)
    {
        if (num[j] < raw_info.frame_size / 128) continue;
        int t = 0;
        for (k = 0; k < FIXP_RANGE; k++)
        {
            t += hist[j][k];
            if (t >= num[j]/2)
            {
                int c = pow(2, H2F(k)) * FIXP_ONE;
                stripes_coeffs[j] = c;
                break;
            }
        }
    }

#if 0
    /* debug graphs */
    FILE* f = fopen("raw2dng.m", "w");
    fprintf(f, "h = {}; x = {}; c = \"rgbcmy\"; \n");
    for (j = 2; j < 8; j++)
    {
        fprintf(f, "h{end+1} = [");
        for (k = 1; k < FIXP_RANGE-1; k++)
        {
            fprintf(f, "%d ", hist[j][k]);
        }
        fprintf(f, "];\n");

        fprintf(f, "x{end+1} = [");
        for (k = 1; k < FIXP_RANGE-1; k++)
        {
            fprintf(f, "%f ", H2F(k) );
        }
        fprintf(f, "];\n");
        fprintf(f, "plot(log2(%d/%d) + [0 0], [0 %d], ['*-' c(%d)]); hold on;\n", stripes_coeffs[j], FIXP_ONE, max[j], j-1);
    }
    fprintf(f, "for i = 1:6, plot(x{i}, h{i}, c(i)); hold on; end;");
    fprintf(f, "axis([-0.05 0.05])");
    fclose(f);
    system("octave-cli --persist raw2dng.m");
#endif

    stripes_coeffs[0] = FIXP_ONE;

    /* do we really need stripe correction, or it won't be noticeable? or maybe it's just computation error? */
    stripes_correction_needed = 0;
    for (j = 0; j < 8; j++)
    {
        double c = (double)stripes_coeffs[j] / FIXP_ONE;
        if (c < 0.998 || c > 1.002)
            stripes_correction_needed = 1;
    }

    if (stripes_correction_needed || force_correction)
    {
        printf("\n\nVertical stripes correction:\n");
        for (j = 0; j < 8; j++)
        {
            if (stripes_coeffs[j])
                printf("  %.5f", (double)stripes_coeffs[j] / FIXP_ONE);
            else
                printf("    1  ");
        }
        printf("\n");
    }
}

static void apply_vertical_stripes_correction(struct raw_info raw_info, uint32_t * image_data)
{
    /**
     * inexact white level will result in banding in highlights, especially if some channels are clipped
     *
     * so... we'll try to use a better estimation of white level *for this particular purpose*
     * start with a gross under-estimation, then consider white = max(all pixels)
     * just in case the exif one is way off
     * reason:
     *   - if there are no pixels above the true white level, it shouldn't hurt;
     *     worst case, the brightest pixel(s) will be underexposed by 0.1 EV or so
     *   - if there are, we will choose the true white level
     */

    int white = raw_info.white_level * 2 / 3;

    typedef raw_pixblock_14 raw_pixblock; //TODO: should we change depending on the source bits to raw_pixblock_12 or raw_pixblock_10?
    raw_pixblock * row;

    for (row =image_data; (void*)row < (void*)image_data + raw_info.pitch * raw_info.height; row += raw_info.pitch / sizeof(raw_pixblock))
    {
        raw_pixblock * p;
        for (p = row; (void*)p < (void*)row + raw_info.pitch; p++)
        {
            white = MAX(white, PA);
            white = MAX(white, PB);
            white = MAX(white, PC);
            white = MAX(white, PD);
            white = MAX(white, PE);
            white = MAX(white, PF);
            white = MAX(white, PG);
            white = MAX(white, PH);
        }
    }

    int black = raw_info.black_level;
    for (row = image_data; (void*)row < (void*)image_data + raw_info.pitch * raw_info.height; row += raw_info.pitch / sizeof(raw_pixblock))
    {
        raw_pixblock * p;
        for (p = row; (void*)p < (void*)row + raw_info.pitch; p++)
        {
            int pa = PA;
            int pb = PB;
            int pc = PC;
            int pd = PD;
            int pe = PE;
            int pf = PF;
            int pg = PG;
            int ph = PH;

            /**
             * Thou shalt not exceed the white level (the exact one, not the exif one)
             * otherwise you'll be blessed with banding instead of nice and smooth highlight recovery
             *
             * At very dark levels, you will introduce roundoff errors, so don't correct there
             */

            if (stripes_coeffs[0] && pa && pa < white && pa > black + 64) SET_PA(MIN(white, RAW_MUL(pa, stripes_coeffs[0])));
            if (stripes_coeffs[1] && pb && pb < white && pa > black + 64) SET_PB(MIN(white, RAW_MUL(pb, stripes_coeffs[1])));
            if (stripes_coeffs[2] && pc && pc < white && pa > black + 64) SET_PC(MIN(white, RAW_MUL(pc, stripes_coeffs[2])));
            if (stripes_coeffs[3] && pd && pd < white && pa > black + 64) SET_PD(MIN(white, RAW_MUL(pd, stripes_coeffs[3])));
            if (stripes_coeffs[4] && pe && pe < white && pa > black + 64) SET_PE(MIN(white, RAW_MUL(pe, stripes_coeffs[4])));
            if (stripes_coeffs[5] && pf && pf < white && pa > black + 64) SET_PF(MIN(white, RAW_MUL(pf, stripes_coeffs[5])));
            if (stripes_coeffs[6] && pg && pg < white && pa > black + 64) SET_PG(MIN(white, RAW_MUL(pg, stripes_coeffs[6])));
            if (stripes_coeffs[7] && ph && ph < white && pa > black + 64) SET_PH(MIN(white, RAW_MUL(ph, stripes_coeffs[7])));
        }
    }
}
void fix_vertical_stripes_diso(struct raw_info raw_info, uint32_t * image_data, int force_correction)
{
    /* for speed: only detect correction factors from the first frame */
    static int first_time = 1;
    if (first_time)
    {
        detect_vertical_stripes_coeffs(raw_info, image_data, force_correction);
        first_time = 0;
    }

    apply_vertical_stripes_correction(raw_info, image_data);
}

void find_and_fix_cold_pixels(struct raw_info raw_info, uint32_t * raw_buffer_32, int force_analysis)
{
    #define MAX_COLD_PIXELS 200000

    struct xy { int x; int y; };

    static struct xy cold_pixel_list[MAX_COLD_PIXELS];
    static int cold_pixels = -1;

    int w = raw_info.width;
    int h = raw_info.height;

    /* scan for bad pixels in the first frame only, or on request*/
    if (cold_pixels < 0 || force_analysis)
    {
        cold_pixels = 0;

        /* at sane ISOs, noise stdev is well less than 50, so 200 should be enough */
        int cold_thr = MAX(0, raw_info.black_level - 200);

        /* analyse all pixels of the frame */
        for (int y = 0; y < h; y++)
        {
            for (int x = 0; x < w; x++)
            {
                int p = raw_get_pixel32(x, y);
                int is_cold = (p < cold_thr);

                /* create a list containing the cold pixels */
                if (is_cold && cold_pixels < MAX_COLD_PIXELS)
                {
                    cold_pixel_list[cold_pixels].x = x;
                    cold_pixel_list[cold_pixels].y = y;
                    cold_pixels++;
                }
            }
        }
#ifndef STDOUT_SILENT
        printf("\rCold pixels : %d                    \n", (cold_pixels));
#endif
    }

    /* repair the cold pixels */
    for (int p = 0; p < cold_pixels; p++)
    {
        int x = cold_pixel_list[p].x;
        int y = cold_pixel_list[p].y;

        int neighbours[100];
        int k = 0;
        int fc0 = FC(x, y);

        /* examine the neighbours of the cold pixel */
        for (int i = -4; i <= 4; i++)
        {
            for (int j = -4; j <= 4; j++)
            {
                /* exclude the cold pixel itself from the examination */
                if (i == 0 && j == 0)
                {
                    continue;
                }

                /* exclude out-of-range coords */
                if (x+j < 0 || x+j >= w || y+i < 0 || y+i >= h)
                {
                    continue;
                }

                /* examine only the neighbours of the same color */
                if (FC(x+j, y+i) != fc0)
                {
                    continue;
                }

                int p = raw_get_pixel32(x+j, y+i);
                neighbours[k++] = -p;
            }
        }

        /* replace the cold pixel with the median of the neighbours */
        raw_set_pixel32(x, y, -median_int_wirth(neighbours, k));
    }

}

static void find_and_fix_bad_pixels(struct raw_info raw_info, uint32_t * raw_buffer_32, int dark_noise, int bright_noise, int bad_pixels_search_method, int * is_bright, int black, int white)
{
    int w = raw_info.width;
    int h = raw_info.height;

    /* for fast EV - raw conversion */
    static int raw2ev[1<<20];   /* EV x EV_RESOLUTION */
    static int ev2raw_0[24*EV_RESOLUTION];
    static uint32_t previous_black = -1;

    /* handle sub-black values (negative EV) */
    int* ev2raw = ev2raw_0 + 10*EV_RESOLUTION;

    if(black != previous_black)
    {
        build_ev2raw_lut(raw2ev, ev2raw_0, black, white);
        previous_black = black;
    }

#ifndef STDOUT_SILENT
    printf("Looking for hot/cold pixels...\n");
#endif


    /* hot pixel map */
    uint32_t* hotpixel = malloc(w * h * sizeof(uint32_t));
    memset(hotpixel, 0, w * h * sizeof(uint32_t));

    int hot_pixels = 0;
    int cold_pixels = 0;

    /* really dark pixels (way below the black level) are probably noise */
    /* there might be dark pixels not that much below the black level, but they need further checking */
    int cold_thr = MAX(0, black - dark_noise*8);
    int maybe_cold_thr = black + dark_noise*2;

#pragma omp parallel for schedule(static) collapse(2)
    for (int y = 6; y < h-6; y ++)
    {
        for (int x = 6; x < w-6; x ++)
        {
            int p = raw_get_pixel20(x, y);

            int is_hot = 0;
            int is_cold = (p < cold_thr);
            int maybe_cold = (p < maybe_cold_thr);

            /* we don't have no hot pixels on the bright exposure */
            /* but we may have cold pixels */
            if (!BRIGHT_ROW || maybe_cold)
            {
                /* let's look at the neighbours: is this pixel clearly brigher? (isolated) */
                int neighbours[100];
                int k = 0;
                int fc0 = FC(x, y);
                int b0 = is_bright[y%4];
                int max = 0;
                for (int i = -4; i <= 4; i++)
                {
                    /* only look at pixels of the same brightness */
                    if (is_bright[(y+i)%4] != b0)
                        continue;

                    for (int j = -4; j <= 4; j++)
                    {
                        if (i == 0 && j == 0)
                            continue;

                        /* only look at pixels of the same color */
                        if (FC(x+j, y+i) != fc0)
                            continue;

                        int p = raw_get_pixel20(x+j, y+i);
                        neighbours[k++] = -p;
                        max = MAX(max, p);
                    }

                    /* this difference will only get lower, so if it's already too low (see below), stop scanning */
                    /* (don't stop scanning if the pixel is cold, since we'll need this info to interpolate it) */
                    if (raw2ev[p] - raw2ev[max] <= EV_RESOLUTION && !maybe_cold)
                        break;
                }

                is_hot = (raw2ev[p] - raw2ev[max] > EV_RESOLUTION) && (max > black + 8*dark_noise);

                if (maybe_cold)
                {
                    /* there may be cold pixels very close to black level */
                    /* heuristic: if it's much darker than the brightest neighbour, it's a cold pixel */
                    is_cold |= (raw2ev[max] - raw2ev[p] > EV_RESOLUTION * 10);
                }

                if (bad_pixels_search_method == 1)    /* aggressive */
                {
                    int third_max = -kth_smallest_int(neighbours, k, 2);
                    is_hot = ((raw2ev[p] - raw2ev[max] > EV_RESOLUTION/4) && (max > black + 8*dark_noise))
                          || (raw2ev[p] - raw2ev[third_max] > EV_RESOLUTION/2);
                }

                if (is_hot)
                {
                    hot_pixels++;
                    hotpixel[x + y*w] = -kth_smallest_int(neighbours, k, 2);
                }

                if (is_cold)
                {
                    cold_pixels++;
                    hotpixel[x + y*w] = -median_int_wirth(neighbours, k);
                }
            }
        }
    }

    /* apply the correction */
#pragma omp parallel for schedule(static) collapse(2)
    for (int y = 0; y < h; y ++)
        for (int x = 0; x < w; x ++)
            if (hotpixel[x + y*w])
                raw_set_pixel20(x, y, hotpixel[x + y*w]);

#ifndef STDOUT_SILENT
    if (hot_pixels)
        printf("Hot pixels      : %d\n", hot_pixels);

    if (cold_pixels)
        printf("Cold pixels     : %d\n", cold_pixels);
#endif


    free(hotpixel);
}

int diso_get_full20bit(struct raw_info raw_info,dual_iso_freeze_data_t* iso_data, uint16_t * image_data, int interp_method, int use_alias_map, int use_fullres, 
    int chroma_smooth_method, int vertical_stripes_fix, int use_horizontal_stripe_fix, int fix_bad_pixels_dual, int bad_pixels_search_method, 
    int dark_highlight_threshold)
{
    int w = raw_info.width;
    int h = raw_info.height;

    if (w <= 0 || h <= 0) return 0;

#ifdef PERF_INFO
    perf_clock = clock();
#endif
    /* RGGB or GBRG? */
    int rggb = iso_data->rggb;
    if (iso_data->freeze < 2 || iso_data->rggb < 0){
        rggb = identify_rggb_or_gbrg(raw_info, image_data);
        iso_data->rggb = rggb;
    } else {
        printf("Using cache\n");
    }
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("identify_rggb_or_gbrg took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    
    if (!rggb) /* this code assumes RGGB, so we need to skip one line */
    {
        image_data += raw_info.pitch;
        raw_info.active_area.y1++;
        raw_info.active_area.y2--;
        raw_info.height--;
        h--;
    }
    
    int is_bright[4];
#ifdef PERF_INFO
    perf_clock = clock();
#endif

    if (iso_data->freeze < 2 || iso_data->is_bright[0] < 0){
        if(!identify_bright_and_dark_fields(raw_info, image_data, rggb, is_bright)){    
            iso_data->rggb = -1;
            return 0;
        }
    } else {
        is_bright[0] = iso_data->is_bright[0];
        is_bright[1] = iso_data->is_bright[1];
        is_bright[2] = iso_data->is_bright[2];
        is_bright[3] = iso_data->is_bright[3];
    }

    if (iso_data->freeze == 1){
        iso_data->is_bright[0] = is_bright[0];
        iso_data->is_bright[1] = is_bright[1];
        iso_data->is_bright[2] = is_bright[2];
        iso_data->is_bright[3] = is_bright[3];
    }

#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("identify_bright_and_dark_fields %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    int ret = 0;
    
    /* will use 20-bit processing and 16-bit output, instead of 14 */
    raw_info.black_level *= 64;
    raw_info.white_level *= 64;
    
    int black = raw_info.black_level;
    int white = raw_info.white_level;
    
    int white_bright = white;
#ifdef PERF_INFO
    perf_clock = clock();
#endif

    white_detect(raw_info, image_data, &white, &white_bright, is_bright);
    white *= 64;
    white_bright *= 64;
    raw_info.white_level = white;

#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("white_detect took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif

    uint32_t * raw_buffer_32 = convert_to_20bit(raw_info, image_data);
    double noise_std[4];
    double dark_noise, bright_noise, dark_noise_ev, bright_noise_ev;

#ifdef PERF_INFO
    perf_clock = clock();
#endif

    double noise_avg = compute_noise20(raw_info, raw_buffer_32, noise_std, &dark_noise, &bright_noise, &dark_noise_ev, &bright_noise_ev);

#ifdef PERF_INFO

    perf_clock = clock()-perf_clock;
    printf("compute_noise took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifdef PERF_INFO
    perf_clock = clock();
#endif
    /* promote from 14 to 20 bits (original raw buffer holds 14-bit values stored as uint16_t) */

#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("convert_to_20bit took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif

#ifdef PERF_INFO
    perf_clock = clock();
#endif
    /* dark and bright exposures, interpolated */
    uint32_t* dark   = malloc(w * h * sizeof(uint32_t));
    uint32_t* bright = malloc(w * h * sizeof(uint32_t));
    memset(dark, 0, w * h * sizeof(uint32_t));
    memset(bright, 0, w * h * sizeof(uint32_t));
    
    /* fullres image (minimizes aliasing) */
    uint32_t* fullres = malloc(w * h * sizeof(uint32_t));
    memset(fullres, 0, w * h * sizeof(uint32_t));
    uint32_t* fullres_smooth = fullres;
    
    /* halfres image (minimizes noise and banding) */
    uint32_t* halfres = malloc(w * h * sizeof(uint32_t));
    memset(halfres, 0, w * h * sizeof(uint32_t));
    uint32_t* halfres_smooth = halfres;
       
    /* overexposure map */
    uint16_t * overexposed = malloc(w * h * sizeof(uint16_t));
    memset(overexposed, 0, w * h * sizeof(uint16_t));
    
    uint16_t* alias_map = NULL;
    if(use_alias_map)
    {
        alias_map = malloc(w * h * sizeof(uint16_t));
        memset(alias_map, 0, w * h * sizeof(uint16_t));
    }
    
    /* estimate ISO difference between bright and dark exposures */
    double corr_ev = 0;
    int white_darkened = white_bright;
#ifdef PERF_INFO
    perf_clock = clock();
#endif

    int expo_matched = match_exposures(raw_info, raw_buffer_32, &corr_ev, &white_darkened, is_bright, iso_data);

#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("match_exposures took %f seconds\n", ((double) perf_clock) / CLOCKS
       #ifndef STDOUT_SILENT_PER_SEC);
    fflush(stdout);
#endif
    if(expo_matched)
    {
        printf("Exposures matched");
    }
    else
    {
        printf("Exposures not matched");
    }
#endif

    /* estimate dynamic range */
    double lowiso_dr = log2(white - black) - dark_noise_ev;
#ifndef STDOUT_SILENT
    double highiso_dr = log2(white_bright - black) - bright_noise_ev;
    printf("Dynamic range   : %.02f (+) %.02f => %.02f EV (@14 bits, in theory)\n", lowiso_dr-6, highiso_dr-6, highiso_dr + corr_ev-6);
#endif
    /* correction factor for the bright exposure, which was just darkened */
    double corr = pow(2, corr_ev);

    /* update bright noise measurements, so they can be compared after scaling */
    bright_noise /= corr;
    bright_noise_ev -= corr_ev;

    if(vertical_stripes_fix){
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        int force = vertical_stripes_fix > 1;
        fix_vertical_stripes_diso(raw_info, raw_buffer_32, force);

#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("fix_vertical_stripes_diso took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    }

    if (fix_bad_pixels_dual == 1) //Only if auto mode selected. Map mode is handled outside of Dual ISO (and it's the way to go tbh...)
    {
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        /* best done before interpolation */
        find_and_fix_bad_pixels(raw_info, raw_buffer_32, dark_noise, bright_noise, bad_pixels_search_method, is_bright, black, white);

        //find_and_fix_cold_pixels(raw_info, raw_buffer_32, 0); //Simpler alternative, not any better
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("fix bad pixels took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    }

#ifdef PERF_INFO
    perf_clock = clock();
    printf("Interpolation started:\n");
#endif
    if(interp_method == 0)
    {
        amaze_interpolate(raw_info, raw_buffer_32, dark, bright, black, white, white_darkened, is_bright);
    }
    else
    {
        mean23_interpolate(raw_info, raw_buffer_32, dark, bright, black, white, white_darkened, is_bright);
    }
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("interpolation took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifdef PERF_INFO
    perf_clock = clock();
#endif
    border_interpolate(raw_info, raw_buffer_32, dark, bright, is_bright);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("border_interpolate took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    if (use_horizontal_stripe_fix)
    {
#ifdef PERF_INFO
    perf_clock = clock();
#endif
#ifndef STDOUT_SILENT
        printf("Horizontal stripe fix...\n");
#endif

        int* delta = malloc(w * sizeof(delta[0]));

        /* adjust dark lines to match the bright ones */
        for (int y = raw_info.active_area.y1; y < raw_info.active_area.y2; y ++)
        {
            /* apply a constant offset (estimated from unclipped areas) */
            int delta_num = 0;
            for (int x = raw_info.active_area.x1; x < raw_info.active_area.x2; x ++)
            {
                int b = bright[x + y*w];
                int d = dark[x + y*w];
                if (MAX(b,d) < white_darkened)
                {
                    delta[delta_num++] = b - d;
                }
            }

            if (delta_num < 200)
            {
                //~ printf("%d: too few points (%d)\n", y, delta_num);
                continue;
            }

            /* compute median difference */
            int med_delta = median_int_wirth(delta, delta_num);

            if (ABS(med_delta) > 200*16)
            {
#ifndef STDOUT_SILENT
                printf("%d: offset too large (%d)\n", y, med_delta);
#endif
                continue;
            }

            /* shift the dark lines */
            for (int x = 0; x < w; x ++)
            {
                dark[x + y*w] = COERCE(dark[x + y*w] + med_delta, 0, 0xFFFFF);
            }
        }
        free(delta);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("Horizontal stripe fix took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    }

    if (use_fullres){
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        fullres_reconstruction(raw_info, fullres, dark, bright, white_darkened, is_bright, dark_highlight_threshold);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("identify_bright_and_dark_fields took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    }

    //New mem space for the smoothed versions
    if (chroma_smooth_method)
    {
        if (use_fullres)
        {
            fullres_smooth = malloc(w * h * sizeof(uint32_t));
        }
        halfres_smooth = malloc(w * h * sizeof(uint32_t));
    }
#ifdef PERF_INFO
    perf_clock = clock();
#endif
    int mix_res = mix_images(raw_info, fullres, fullres_smooth, halfres, halfres_smooth, alias_map, dark, bright, overexposed, dark_noise, white_darkened, corr_ev, lowiso_dr, black, white, chroma_smooth_method, use_fullres);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("mix_images took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    if(mix_res)
    {
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        /* let's check the ideal noise levels (on the halfres image, which in black areas is identical to the bright one) */
#pragma omp parallel for schedule(static) default(none) shared(raw_info, raw_buffer_32,h,w,bright) collapse(2)
        for (int y = 3; y < h-2; y ++)
            for (int x = 2; x < w-2; x ++)
                raw_set_pixel32(x, y, bright[x + y*w]);

        compute_black_noise20(raw_info, raw_buffer_32, raw_info.active_area.y1 + 20, raw_info.active_area.y2 - 20, 1, 1, &noise_avg, &noise_std[0]);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("compute_black_noise took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifndef STDOUT_SILENT
        double ideal_noise_std = noise_std[0];
#endif
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        final_blend(raw_info, raw_buffer_32, fullres, fullres_smooth, halfres_smooth, dark, bright, overexposed, alias_map, black, white, dark_noise, use_fullres);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("final_blend took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        /* let's see how much dynamic range we actually got */
        compute_black_noise20(raw_info, raw_buffer_32, raw_info.active_area.y1 + 20, raw_info.active_area.y2 - 20, 1, 1, &noise_avg, &noise_std[0]); 
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("compute_black_noise took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
#ifndef STDOUT_SILENT
        printf("Noise level     : %.02f (20-bit), ideally %.02f\n", noise_std[0], ideal_noise_std);
        printf("Dynamic range   : %.02f EV (cooked)\n", log2(white - black) - log2(noise_std[0]));
#endif
#ifdef PERF_INFO
    perf_clock = clock();
#endif
        convert_20_to_16bit(raw_info, image_data, raw_buffer_32);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("convert_20_to_16bit took %f seconds\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
        ret = 1;
    }
#ifdef PERF_INFO
    perf_clock = clock();
#endif
    if (!rggb) /* back to GBRG */
    {
        raw_info.active_area.y1--;
        raw_info.active_area.y2++;
        raw_info.height++;
        h++;
    }
    
    free(dark);
    free(bright);
    free(fullres);
    free(halfres);
    free(alias_map);
    free(overexposed);
    free(raw_buffer_32);
    if (fullres_smooth && fullres_smooth != fullres) free(fullres_smooth);
    if (halfres_smooth && halfres_smooth != halfres) free(halfres_smooth);
#ifdef PERF_INFO
    perf_clock = clock()-perf_clock;
    printf("free memory took %f seconds\n\n", ((double) perf_clock) / CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    return ret;
}

