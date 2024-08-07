#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // Nearest Neighbor interpolation
    return get_pixel(im, round(x), round(y), c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image resize = make_image(w,h,im.c);
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            for (int c = 0; c < im.c; c++) {
                float rx = (float) x / w * im.w;
                float ry = (float) y / h * im.h;
                set_pixel(resize, x, y, c, nn_interpolate(im, rx, ry, c));
            }
        }
    }
    return resize;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // what about edge cases? If x is an integer, x - fx = 0, value will always be 0.
    // image coordinate system 0 in top left
    int fx = floor(x);
    int fy = floor(y);
    int cx = ceil(x);
    int cy = ceil(y);
    float top_left =    get_pixel(im, fx, fy, c);
    float top_right =   get_pixel(im, cx, fy, c);
    float bot_left =    get_pixel(im, fx, cy, c);
    float bot_right =   get_pixel(im, cx, cy, c);
    float top_left_d =  x - fx == 0 ? 1 : x - fx;
    float top_right_d = cx - x == 0 ? 0 : cx - x;
    float side_top_d =  y - fy == 0 ? 1 : y - fy;
    float side_bot_d =  cy - y == 0 ? 0 : cy - y;

    float value =
        top_left * top_right_d * side_bot_d + 
        top_right * top_left_d * side_bot_d +
        bot_left * top_right_d * side_top_d +
        bot_right * top_left_d * side_top_d;
    return value;
}

image bilinear_resize(image im, int w, int h)
{
    image resize = make_image(w,h,im.c);
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            for (int c = 0; c < im.c; c++) {
                float rx = (float) x / w * im.w;
                float ry = (float) y / h * im.h;
                set_pixel(resize, x, y, c, bilinear_interpolate(im, rx, ry, c));
            }
        }
    }
    return resize;
}

