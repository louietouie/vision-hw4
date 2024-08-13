#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "matrix.h"

// Run an activation function on each element in a matrix,
// modifies the matrix in place
// matrix m: Input to activation function
// ACTIVATION a: function to run
void activate_matrix(matrix m, ACTIVATION a)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        double sum = 0;
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            if(a == LOGISTIC){
                // TODO
                // sigmoid function. maps any value to [0,1]
                // prone to vanishing gradient problem
                m.data[i][j] = 1 / (1 + exp(-a));
            } else if (a == RELU){
                // TODO
                m.data[i][j] = max(0, a);
            } else if (a == LRELU){
                // TODO
                m.data[i][j] = max(0.02 * a, a);
            } else if (a == SOFTMAX){
                // TODO
                // often used in the final layer of neural networks
                m.data[i][j] = exp(a);
            }
            sum += m.data[i][j];
        }
        if (a == SOFTMAX) {
            // TODO: have to normalize by sum if we are using SOFTMAX
            // LOUIS: shouldn't this be one level lower? Why only the sum of each row?
            // ANSWER: No. From homework description, "Each row in the matrix is a separate data point so we want to normalize over each data point separately."
            // ANSWER Continued: No. If this matrix m is following what most matrix structures use, every row represents a sample of data, and the columns are the logits for that sample of data.
            // We only want to divide by the sum of the logits for that row.
            // In this case it is an activation function, I think a row represents the outputs from a neural network layer (often final layer) for a single sample of data, and the columns are the nodes of the layer.
            // this is useful if each node in the final layer is a prediction for a class in a multiclass problem (how is this done? error of each node is based on a different class?) Softmax will turn this final layer into probabilities
            // NOTE: tensorflow says don't use softmax as last layer because makes backpropegation hard
            // https://stackoverflow.com/questions/66747564/softmax-function-of-2d-array
            // https://www.singlestore.com/blog/a-guide-to-softmax-activation-function/
            // https://ai.stackexchange.com/questions/20214/why-does-tensorflow-docs-discourage-using-softmax-as-activation-for-the-last-lay#:~:text=softmax%20in%20as%20the%20activation,when%20using%20a%20softmax%20output.
            // https://xeonqq.github.io/machine%20learning/softmax/

            for(j = 0; j < m.cols; ++j){
                m.data[i][j] = m.data[i][j] / sum;
            }
        }
    }
}

// Calculates the gradient of an activation function and multiplies it into the delta for a layer
// matrix m: an activated layer output
// ACTIVATION a: activation function for a layer
// matrix d: delta before activation gradient
void gradient_matrix(matrix m, ACTIVATION a, matrix d)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            // TODO: multiply the correct element of d by the gradient

            if(a == LOGISTIC){
                // f'(x) = f(x) * (1 - f(x))
                d.data[i][j] *= (x * (1 - x));
            } else if (a == RELU){
                // 1 or 0.
                if (x < 0) {
                    d.data[i][j] = 0;
                }
            } else if (a == LRELU){
                // 1 or 0.02
                if (x < 0) {
                    d.data[i][j] *= 0.02;
                }
            } else if (a == SOFTMAX){
                // 1
                continue;
            }

        }
    }
}
// LOUIS NOTES ON ABOVE
// At this point when this method is called, we will know the derivative of the loss function with respect to the output, say dL/dO
// We want the gradient with respect to the nodes input, dL/dI, so we need to know the gradient of the activation function the node runs, dO/dL
// We are able to calculate the gradient of our activation function f(x), f'(x) without x, which is normally necessary (to plug into f'(x)).
// Instead, we only need y. (I think this is because our activation functions are monotonic, AKA constantly increasing, so each y only maps to one possible dx)



// Forward propagate information through a layer
// layer *l: pointer to the layer
// matrix in: input to layer
// returns: matrix that is output of the layer
matrix forward_layer(layer *l, matrix in)
{

    l->in = in;  // Save the input for backpropagation


    // TODO: fix this! multiply input by weights and apply activation function.
    // matrix out = make_matrix(in.rows, l->w.cols); // placeholder
    matrix out = matrix_mult_matrix(l->in, l->w);
    activate_matrix(out, l->activation);


    free_matrix(l->out);// free the old output
    l->out = out;       // Save the current output for gradient calculation
    return out;
}

// Backward propagate derivatives through a layer
// layer *l: pointer to the layer
// matrix delta: partial derivative of loss w.r.t. output of layer
// returns: matrix, partial derivative of loss w.r.t. input to layer
matrix backward_layer(layer *l, matrix delta)
{
    // 1.4.1
    // delta is dL/dy
    // TODO: modify it in place to be dL/d(xw)
    // delta is changed to loss w.r.t input of layer
    gradient_matrix(l->out, l->activation, delta);

    // 1.4.2
    // TODO: then calculate dL/dw and save it in l->dw
    // matrix dw = make_matrix(l->w.rows, l->w.cols); // replace this
    matrix input_x_t = transpose_matrix(l->in);
    matrix dw = matrix_mult_matrix(delta, input_x_t);
    free_matrix(l->dw);
    l->dw = dw;
    free_matrix(input_x_t);

    
    // 1.4.3
    // TODO: finally, calculate dL/dx and return it.
    // matrix dx = make_matrix(l->in.rows, l->in.cols); // replace this
    matrix input_w_t = transpose_matrix(l->w);
    matrix dx = matrix_mult_matrix(delta, input_w_t);
    return dx;
}

// Update the weights at layer l
// layer *l: pointer to the layer
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_layer(layer *l, double rate, double momentum, double decay)
{
    // LOUIS NOTE: this beginner neural network does not have a bias term, only weights. Redmon says for the training we are doing, it isn't necessary.

    // TODO:
    // Calculate Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    // save it to l->v
    scale_matrix(l->dw, rate);
    scale_matrix(l->v, momentum);
    scale_matrix(l->w, decay);
    l->dw = scale_matrix(l->dw, rate) + scale_matrix(l->v, momentum) - scale_matrix(l->w, decay);
    l -> dw = matrix_sub_matrix(l -> dw, l-> w)

    // NEED TO DO dw + v - w

    l->v = l->dw;

    // Update l->w
    l->w = l->w + l->dw;

    // Remember to free any intermediate results to avoid memory leaks

}

// Make a new layer for our model
// int input: number of inputs to the layer
// int output: number of outputs from the layer
// ACTIVATION activation: the activation function to use
layer make_layer(int input, int output, ACTIVATION activation)
{
    layer l;
    l.in  = make_matrix(1,1);
    l.out = make_matrix(1,1);
    l.w   = random_matrix(input, output, sqrt(2./input));
    l.v   = make_matrix(input, output);
    l.dw  = make_matrix(input, output);
    l.activation = activation;
    return l;
}

// Run a model on input X
// model m: model to run
// matrix X: input to model
// returns: result matrix
matrix forward_model(model m, matrix X)
{
    int i;
    for(i = 0; i < m.n; ++i){
        X = forward_layer(m.layers + i, X);
    }
    return X;
}

// Run a model backward given gradient dL
// model m: model to run
// matrix dL: partial derivative of loss w.r.t. model output dL/dy
void backward_model(model m, matrix dL)
{
    matrix d = copy_matrix(dL);
    int i;
    for(i = m.n-1; i >= 0; --i){
        matrix prev = backward_layer(m.layers + i, d);
        free_matrix(d);
        d = prev;
    }
    free_matrix(d);
}

// Update the model weights
// model m: model to update
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_model(model m, double rate, double momentum, double decay)
{
    int i;
    for(i = 0; i < m.n; ++i){
        update_layer(m.layers + i, rate, momentum, decay);
    }
}

// Find the index of the maximum element in an array
// double *a: array
// int n: size of a, |a|
// returns: index of maximum element
int max_index(double *a, int n)
{
    if(n <= 0) return -1;
    int i;
    int max_i = 0;
    double max = a[0];
    for (i = 1; i < n; ++i) {
        if (a[i] > max){
            max = a[i];
            max_i = i;
        }
    }
    return max_i;
}

// Calculate the accuracy of a model on some data d
// model m: model to run
// data d: data to run on
// returns: accuracy, number correct / total
double accuracy_model(model m, data d)
{
    matrix p = forward_model(m, d.X);
    int i;
    int correct = 0;
    for(i = 0; i < d.y.rows; ++i){
        if(max_index(d.y.data[i], d.y.cols) == max_index(p.data[i], p.cols)) ++correct;
    }
    return (double)correct / d.y.rows;
}

// Calculate the cross-entropy loss for a set of predictions
// matrix y: the correct values
// matrix p: the predictions
// returns: average cross-entropy loss over data points, 1/n Σ(-ylog(p))
double cross_entropy_loss(matrix y, matrix p)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < y.rows; ++i){
        for(j = 0; j < y.cols; ++j){
            sum += -y.data[i][j]*log(p.data[i][j]);
        }
    }
    return sum/y.rows;
}


// Train a model on a dataset using SGD
// model m: model to train
// data d: dataset to train on
// int batch: batch size for SGD
// int iters: number of iterations of SGD to run (i.e. how many batches)
// double rate: learning rate
// double momentum: momentum
// double decay: weight decay
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay)
{
    int e;
    for(e = 0; e < iters; ++e){
        data b = random_batch(d, batch);
        matrix p = forward_model(m, b.X);
        fprintf(stderr, "%06d: Loss: %f\n", e, cross_entropy_loss(b.y, p));
        matrix dL = axpy_matrix(-1, p, b.y); // partial derivative of loss dL/dy
        backward_model(m, dL);
        update_model(m, rate/batch, momentum, decay);
        free_matrix(dL);
        free_data(b);
    }
}
