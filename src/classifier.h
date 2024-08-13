#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "matrix.h"

layer make_layer(int input, int output, ACTIVATION activation);
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay);