#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    for (int c = 0; c < im.c; c++) {
        float total = 0;
        for (int x = 0; x < im.w; x++) {
            for (int y = 0; y < im.h; y++) {
                    total = total + get_pixel(im, x, y, c);
            }
        }
        for (int x = 0; x < im.w; x++) {
            for (int y = 0; y < im.h; y++) {
                set_pixel(im, x, y, c, get_pixel(im, x, y, c) / total);
            }
        }
    }
}

image make_box_filter(int w)
{
    // TODO
    image box = make_image(w,w,1);
    l1_normalize(box);
    return box;
}

image convolve_image(image im, image filter, int preserve)
{  
    image convolve;
    if (preserve == 1) {
        convolve = make_image(im.w, im.h, im.c);
        assert(filter.c == 1 || filter.c == im.c);
        for (int x = 0; x < im.w; x++) {
            for (int y = 0; y < im.h; y++) {
                for (int c = 0; c < im.c; c++) {
                    float v = 0;
                    int fx_adj = floor(filter.w / 2);
                    int fy_adj = floor(filter.h / 2);
                    for (int fx = 0; fx < filter.w; fx++) {
                        for (int fy = 0; fy < filter.h; fy++) {
                            int fc = 0;
                            if (filter.c == im.c) {fc = c;}
                            v = v + get_pixel(filter, fx, fy, fc) * get_pixel(im, x + (fx - fx_adj), y + (fy - fy_adj), c);
                        }
                    }
                    set_pixel(convolve, x, y, c, v);
                }
            }
        }
    }

    // TODO: if preserve is 0, should probably return a single channel image, not 3 channel image with same value in each channel
    else {
        convolve = make_image(im.w, im.h, 1);
        assert(filter.c == 1 || filter.c == im.c);
        for (int x = 0; x < im.w; x++) {
            for (int y = 0; y < im.h; y++) {
                float v = 0;
                for (int c = 0; c < im.c; c++) {
                    int fx_adj = floor(filter.w / 2);
                    int fy_adj = floor(filter.h / 2);
                    for (int fx = 0; fx < filter.w; fx++) {
                        for (int fy = 0; fy < filter.h; fy++) {
                            int fc = 0;
                            if (filter.c == im.c) {fc = c;}
                            v = v + get_pixel(filter, fx, fy, fc) * get_pixel(im, x + (fx - fx_adj), y + (fy - fy_adj), c);
                        }
                    }
                }
                set_pixel(convolve, x, y, 0, v);
            }
        }
    }
    return convolve;
}

image make_highpass_filter()
{
    image filter = make_image(3,3,1);
    filter.data[0] = 0;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 4;
    filter.data[5] = -1;
    filter.data[6] = 0;
    filter.data[7] = -1;
    filter.data[8] = 0;
    return filter;
}

image make_sharpen_filter()
{
    image filter = make_image(3,3,1);
    filter.data[0] = 0;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 5;
    filter.data[5] = -1;
    filter.data[6] = 0;
    filter.data[7] = -1;
    filter.data[8] = 0;
    return filter;
}

image make_emboss_filter()
{
    image filter = make_image(3,3,1);
    filter.data[0] = -2;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 1;
    filter.data[5] = 1;
    filter.data[6] = 0;
    filter.data[7] = 1;
    filter.data[8] = 2;
    return filter;
}


// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    int size = ceil(6 * sigma);
    if (size % 2 == 0 ) size = size + 1;
    int center = ceil(size / 2);
    image filter = make_image(size, size, 1);
    for (int x = 0; x < filter.w; x++) {
        for (int y = 0; y < filter.h; y++) {
            int x_center = x - center;
            int y_center = y - center;
            float v = 1/(2*3.14*pow(sigma,2))*exp(-(pow(x_center,2) + pow(y_center,2))/(2*pow(sigma, 2)));
            set_pixel(filter, x, y, 0, v);
        }
    }
    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image make_gx_filter()
{
    image filter = make_image(3,3,1);
    filter.data[0] = -1;
    filter.data[1] = 0;
    filter.data[2] = 1;
    filter.data[3] = -2;
    filter.data[4] = 0;
    filter.data[5] = 2;
    filter.data[6] = -1;
    filter.data[7] = 0;
    filter.data[8] = 1;
    return filter;
}

image make_gy_filter()
{
    image filter = make_image(3,3,1);
    filter.data[0] = -1;
    filter.data[1] = -2;
    filter.data[2] = -1;
    filter.data[3] = 0;
    filter.data[4] = 0;
    filter.data[5] = 0;
    filter.data[6] = 1;
    filter.data[7] = 2;
    filter.data[8] = 1;
    return filter;
}

void feature_normalize(image im)
{
    float maximum = get_pixel(im, 0, 0, 0);
    float minimum = get_pixel(im, 0, 0, 0);
    float range;
    int i,j,k;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                if (maximum < get_pixel(im, i, j, k)) {maximum = get_pixel(im, i, j, k);}
                if (minimum > get_pixel(im, i, j, k)) {minimum = get_pixel(im, i, j, k);}
            }
        }
    }
    range = maximum - minimum;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                if (range == 0) {set_pixel(im, i, j, k, 0);}
                else {set_pixel(im, i, j, k, (get_pixel(im, i, j, k) - minimum) / range);}
            }
        }
    }

}

image *sobel_image(image im)
{
    image *result = calloc(2, sizeof(image));
    image gx = convolve_image(im, make_gx_filter(), 0);
    image gy = convolve_image(im, make_gy_filter(), 0);
    image magnitude = make_image(im.w, im.h, 1);
    image direction = make_image(im.w, im.h, 1);
    int i, j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            float gxv = get_pixel(gx, i, j, 0);
            float gyv = get_pixel(gy, i, j, 0);
            set_pixel(magnitude, i, j, 0, sqrt(pow(gxv,2)+pow(gyv,2)));
            set_pixel(direction, i, j, 0, atan2(gyv, gxv));
        }
    }
    result[0] = magnitude;
    result[1] = direction;
    return result;
}

image colorize_sobel(image im)
{
    image result = make_image(im.w, im.h, 3);
    image *sobel_info = sobel_image(im);
    image gaussian_smooth = make_gaussian_filter(2);
    int i, j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            set_pixel(result, i, j, 0, get_pixel(sobel_info[1], i, j, 0)); // H
            set_pixel(result, i, j, 0, get_pixel(sobel_info[0], i, j, 0)); // S
            set_pixel(result, i, j, 0, get_pixel(sobel_info[0], i, j, 0)); // V
        }
    }
    return convolve_image(result, gaussian_smooth, 1);;
}
