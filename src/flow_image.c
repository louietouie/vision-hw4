#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Draws a line on an image with color corresponding to the direction of line
// image im: image to draw line on
// float x, y: starting point of line
// float dx, dy: vector corresponding to line angle and magnitude
void draw_line(image im, float x, float y, float dx, float dy)
{
    assert(im.c == 3);
    float angle = 6*(atan2(dy, dx) / TWOPI + .5);
    int index = floor(angle);
    float f = angle - index;
    float r, g, b;
    if(index == 0){
        r = 1; g = f; b = 0;
    } else if(index == 1){
        r = 1-f; g = 1; b = 0;
    } else if(index == 2){
        r = 0; g = 1; b = f;
    } else if(index == 3){
        r = 0; g = 1-f; b = 1;
    } else if(index == 4){
        r = f; g = 0; b = 1;
    } else {
        r = 1; g = 0; b = 1-f;
    }
    float i;
    float d = sqrt(dx*dx + dy*dy);
    for(i = 0; i < d; i += 1){
        int xi = x + dx*i/d;
        int yi = y + dy*i/d;
        set_pixel(im, xi, yi, 0, r);
        set_pixel(im, xi, yi, 1, g);
        set_pixel(im, xi, yi, 2, b);
    }
}

// Make an integral image or summed area table from an image
// image im: image to process
// returns: image I such that I[x,y] = sum{i<=x, j<=y}(im[i,j])
image make_integral_image(image im)
{
    // TODO: fill in the integral image
    image result = make_image(im.w, im.h, im.c);
    int row, column, channel;
    for (channel = 0; channel < im.c; channel++) {
        for (column = 0; column < im.w; column++) {
            for (row = 0; row < im.h; row++) {
                float self =    get_pixel_low_strict(im, column, row, channel);
                float left =    get_pixel_low_strict(im, column - 1, row, channel);
                float top =     get_pixel_low_strict(im, column, row - 1, channel);
                float diag =    get_pixel_low_strict(im, column - 1, row - 1, channel);
                set_pixel(result, column, row, channel, self + left + top - diag);
            }
        }
    }
    return result;
}

// Apply a box filter to an image using an integral image for speed
// image im: image to smooth
// int s: window size for box filter
// returns: smoothed image
image box_filter_image(image im, int size)
{
    int i,j,k;
    image integ = make_integral_image(im);
    image result = make_image(im.w, im.h, im.c);
    // TODO: fill in S using the integral image.
    int row, column, channel;
    for (channel = 0; channel < im.c; channel++) {
        for (column = 0; column < im.w; column++) {
            for (row = 0; row < im.h; row++) {
                float top_left =    get_pixel_low_strict(integ, column - size, row - size, channel);
                float top_right =   get_pixel_low_strict(integ, column + size, row - size, channel);
                float bot_left =    get_pixel_low_strict(integ, column - size, row + size, channel);
                float bot_right =   get_pixel_low_strict(integ, column + size, row + size, channel);
                set_pixel(result, column, row, channel, bot_right - bot_left - top_left + top_left);
            }
        }
    }
    return result;
}

// Calculate the time-structure matrix of an image pair.
// image im: the input image.
// image prev: the previous image in sequence.
// int s: window size for smoothing.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
image time_structure_matrix(image im, image prev, int size)
{
    int i;
    int converted = 0;
    if(im.c == 3){
        converted = 1;
        im = rgb_to_grayscale(im);
        prev = rgb_to_grayscale(prev);
    }

    // TODO: calculate gradients, structure components, and smooth them
    image Structure = make_image(im.w, im.h, 5);
    image sobel_x = convolve_image(im, make_gx_filter(), 0);
    image sobel_y = convolve_image(im, make_gy_filter(), 0);
    int column, row;
    for (column = 0; column < im.w; column++) {
        for (row = 0; row < im.h; row++) {
            float sxv = get_pixel(sobel_x, column, row, 0);
            float syv = get_pixel(sobel_y, column, row, 0);
            float stv = get_pixel(prev, column, row, 0) - get_pixel(im, column, row, 0); // should this be single pixel or smoothed box?
            set_pixel(Structure, column, row, 0, pow(sxv,2));
            set_pixel(Structure, column, row, 1, pow(syv,2));
            set_pixel(Structure, column, row, 2, sxv*syv);
            set_pixel(Structure, column, row, 3, sxv*stv);
            set_pixel(Structure, column, row, 4, syv*stv);
        }
    }

    if(converted){
        free_image(im); free_image(prev);
    }
    return box_filter_image(Structure, size);
}

// Calculate the velocity given a structure image
// image S: time-structure image
// int stride: only calculate subset of pixels for speed
image velocity_image(image S, int stride)
{
    image result = make_image(S.w/stride, S.h/stride, 3);
    int i, j;
    matrix STS = make_matrix(2,2);
    matrix STT = make_matrix(2,1);
    for(j = (stride-1)/2; j < S.h; j += stride){
        for(i = (stride-1)/2; i < S.w; i += stride){
            // Louis TODO: Isn't this only getting a single pixel? Isn't this underspecified and we need to sum more that 1 to have a fully defined/overdefined equation we can run least-squares on?
            float Ixx = S.data[i + S.w*j + 0*S.w*S.h];
            float Iyy = S.data[i + S.w*j + 1*S.w*S.h];
            float Ixy = S.data[i + S.w*j + 2*S.w*S.h];
            float Ixt = S.data[i + S.w*j + 3*S.w*S.h];
            float Iyt = S.data[i + S.w*j + 4*S.w*S.h];

            // TODO: calculate vx and vy using the flow equation
            // SΔp = T
            // AKA (derivative) * dp = distance
            // AKA (sobel slope of brightness) * (change in x, y position) = (change in brightness)
            // Δp = (S^TS)^-1S^TT

            // check if STS is invertable (determinant != 0)
            // TODO: This determinant might not be correct. Video says to...
            // not check det(STS), but instead check...
            // det(STS - λI).
            // See "The Ancient Secrets of Computer Vision - 08 - Optical Flow" min 49:50
            float vx, vy;
            if (Ixx * Iyy - Ixy * Ixy == 0) {
                vx = 0;
                vy = 0;
            } else {
                STS.data[0][0] = Ixx;
                STS.data[1][0] = Ixy;
                STS.data[0][1] = Ixy;
                STS.data[1][1] = Iyy;
                STT.data[0][0] = Ixt;
                STT.data[1][0] = Iyt;
                matrix velocities = matrix_mult_matrix(matrix_invert(STS), STT);
                vx = velocities.data[0][0];
                vy = velocities.data[1][0];
            }
            set_pixel(result, i/stride, j/stride, 0, vx);
            set_pixel(result, i/stride, j/stride, 1, vy);
        }
    }
    free_matrix(STS);
    free_matrix(STT);
    return result;
}

// Draw lines on an image given the velocity
// image im: image to draw on
// image v: velocity of each pixel
// float scale: scalar to multiply velocity by for drawing
void draw_flow(image im, image v, float scale)
{
    int stride = im.w / v.w;
    int i,j;
    for (j = (stride-1)/2; j < im.h; j += stride) {
        for (i = (stride-1)/2; i < im.w; i += stride) {
            float dx = scale*get_pixel(v, i/stride, j/stride, 0);
            float dy = scale*get_pixel(v, i/stride, j/stride, 1);
            if(fabs(dx) > im.w) dx = 0;
            if(fabs(dy) > im.h) dy = 0;
            draw_line(im, i, j, dx, dy);
        }
    }
}


// Constrain the absolute value of each image pixel
// image im: image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(image im, float v)
{
    int i;
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if (im.data[i] < -v) im.data[i] = -v;
        if (im.data[i] >  v) im.data[i] =  v;
    }
}

// Calculate the optical flow between two images
// image im: current image
// image prev: previous image
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// returns: velocity matrix
image optical_flow_images(image im, image prev, int smooth, int stride)
{
    image S = time_structure_matrix(im, prev, smooth);   
    image v = velocity_image(S, stride);
    constrain_image(v, 6);
    image vs = smooth_image(v, 2);
    free_image(v);
    free_image(S);
    return vs;
}

// Run optical flow demo on webcam
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// int div: downsampling factor for images from webcam
void optical_flow_webcam(int smooth, int stride, int div)
{
#ifdef OPENCV
    CvCapture * cap;
    cap = cvCaptureFromCAM(0);
    image prev = get_image_from_stream(cap);
    image prev_c = nn_resize(prev, prev.w/div, prev.h/div);
    image im = get_image_from_stream(cap);
    image im_c = nn_resize(im, im.w/div, im.h/div);
    while(im.data){
        image copy = copy_image(im);
        image v = optical_flow_images(im_c, prev_c, smooth, stride);
        draw_flow(copy, v, smooth*div);
        int key = show_image(copy, "flow", 5);
        free_image(v);
        free_image(copy);
        free_image(prev);
        free_image(prev_c);
        prev = im;
        prev_c = im_c;
        if(key != -1) {
            key = key % 256;
            printf("%d\n", key);
            if (key == 27) break;
        }
        im = get_image_from_stream(cap);
        im_c = nn_resize(im, im.w/div, im.h/div);
    }
#else
    fprintf(stderr, "Must compile with OpenCV\n");
#endif
}
