#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c) {
    // TODO Fill this in
    int rx = x >= im.w ? (im.w - 1) : (x < 0 ? 0 : x);
    int ry = y >= im.h ? (im.h - 1) : (y < 0 ? 0 : y);
    int rc = c >= im.c ? (im.c - 1) : (c < 0 ? 0 : c);
    return im.data[rx + (im.w * ry) + (rc * im.w * im.h)];
}

float get_pixel_low_strict(image im, int x, int y, int c) {
    // TODO Fill this in
    // if (x >= im.w || x < 0 || y >= im.h || y < 0 || c <= im.c || c < 0) return 0;
    if (x < 0 || y < 0 || c < 0) return 0;
    int rx = x >= im.w ? (im.w - 1) : x;
    int ry = y >= im.h ? (im.h - 1) : y;
    int rc = c >= im.c ? (im.c - 1) : c;
    return im.data[rx + (im.w * ry) + (rc * im.w * im.h)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    if (x >= im.w || x < 0 || y >= im.h || y < 0 || c >= im.c || c < 0) {return;}
    // float fv = v > 1 ? 1 : (v < 0 ? 0 : v);
    im.data[x + (im.w * y) + (c * im.w * im.h)] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    //std::memcpy(copy.data, im.data, sizeof im.data); 
 
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            for (int c = 0; c < im.c; c++) {
                set_pixel(copy, x, y, c, get_pixel(im, x, y, c));
            }
        }
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    float lumaWeights [3] = {0.299, 0.587, 0.114}; 

    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            float grayScale = 0;
            for (int c = 0; c < im.c; c++) {
                grayScale = grayScale + lumaWeights[c] * get_pixel(im, x, y, c);
            }
            set_pixel(gray, x, y, 0, grayScale);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            set_pixel(im, x, y, c, get_pixel(im, x, y, c) + v);
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    // If a value goes over float value of 1, when it is converted back to a byte from 0 to 255, it will overflow past 255 and start back at 0 at a small value.
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            for (int c = 0; c < im.c; c++) {
                if (get_pixel(im, x, y, c) > 1) {
                    set_pixel(im, x, y, c, 1);
                }
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            float r = get_pixel(im, x, y, 0);
            float g = get_pixel(im, x, y, 1);
            float b = get_pixel(im, x, y, 2);
            float V = three_way_max(r, g, b);
            float _m = three_way_min(r, g, b);
            float c = V - _m;
            float S = V == 0 ? 0 : c / V;

            float _H;
            if (c == 0)      {_H = 0;}
            else if (V == r) {_H = (g-b)/c;}
            else if (V == g) {_H = (b-r)/c + 2;}
            else             {_H = (r-g)/c + 4;} //(V == b)

            float H = _H < 0 ? _H/6 + 1 : _H/6;
            set_pixel(im, x, y, 0, H);
            set_pixel(im, x, y, 1, S);
            set_pixel(im, x, y, 2, V);
        }
    }
}

void hsv_to_rgb(image im)
{
    // for (int x = 0; x < im.w; x++) {
    //     for (int y = 0; y < im.h; y++) {
    //         float h = get_pixel(im, x, y, 0);
    //         float s = get_pixel(im, x, y, 1);
    //         float v = get_pixel(im, x, y, 2);

    //         float _h = 6*h

    //     }
    // }
    // TODO Fill this in
}
