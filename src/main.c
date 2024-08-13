#include <math.h>
#include <string.h>
#include "image.h"
#include "test.h"
#include "args.h"
#include "classifier.h"

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    // char *in = find_char_arg(argc, argv, "-i", "data/dog.jpg");
    // char *out = find_char_arg(argc, argv, "-o", "out");
    // //float scale = find_float_arg(argc, argv, "-s", 1);
    // if(argc < 2){
    //     printf("usage: %s [test | grayscale]\n", argv[0]);  
    // } else if (0 == strcmp(argv[1], "test")){
    //     run_tests();
    // } else if (0 == strcmp(argv[1], "grayscale")){
    //     image im = load_image(in);
    //     image g = rgb_to_grayscale(im);
    //     save_image(g, out);
    //     free_image(im);
    //     free_image(g);
    // }
    // return 0;



    // def softmax_model(inputs, outputs):
    //     l = [make_layer(inputs, outputs, SOFTMAX)]
    //     return make_model(l)

    // def neural_net(inputs, outputs):
    //     print(inputs)
    //     l = [   make_layer(inputs, 32, LOGISTIC),
    //             make_layer(32, outputs, SOFTMAX)]
    //     return make_model(l)

    // def make_model(layers):
    //     m = MODEL()
    //     m.n = len(layers)
    //     m.layers = (LAYER*m.n) (*layers)
    //     return m

    data train = load_classification_data("mnist.train", "mnist.labels", 1);
    data test  = load_classification_data("mnist.test", "mnist.labels", 1);

    int batch = 128;
    int iters = 1000;
    double rate = .01;
    double momentum = .9;
    double decay = .0;
    
    //image *result = calloc(2, sizeof(image));
    layer *layers = calloc(1, sizeof(layer));
    layers[0] = make_layer(train.X.cols, train.y.cols, SOFTMAX);
    // layers[1] = make_layer(train.X.cols, train.y.cols, SOFTMAX);
    model m;
    m.layers = layers;
    m.n = 1;

    train_model(m, train, batch, iters, rate, momentum, decay);

    return 0;
}
