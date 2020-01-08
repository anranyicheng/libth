#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/SELU.c"
#else

void THNN_(SELU_updateOutput)(
          THNNState *state,
          THTensor *input,
          THTensor *output,
          accreal alpha_,
          accreal lambda_,
          bool inplace)
{
  real alpha = TH_CONVERT_ACCREAL_TO_REAL(alpha_);
  real lambda = TH_CONVERT_ACCREAL_TO_REAL(lambda_);
  real la = alpha * lambda;
  if(inplace) {
    TH_TENSOR_APPLY(real, input,
      if(*input_data <= 0) {
        *input_data = (exp(*input_data) - 1) * la;
      } else {
        *input_data *= lambda;
      }
    );
    THTensor_(set)(output, input);
  } else {
    THTensor_(resizeAs)(output, input);
    TH_TENSOR_APPLY2(real, input, real, output,
      *output_data = *input_data <= 0 ? (exp(*input_data)-1)*la : *input_data * lambda;
    );
  }
}

void THNN_(SELU_updateGradInput)(
          THNNState *state,
          THTensor *input,
          THTensor *gradOutput,
          THTensor *gradInput,
          THTensor *output,
          accreal alpha_,
          accreal lambda_,
          bool inplace)
{
  real alpha = TH_CONVERT_ACCREAL_TO_REAL(alpha_);
  real lambda = TH_CONVERT_ACCREAL_TO_REAL(lambda_);
  real la = alpha * lambda;
  THNN_CHECK_NELEMENT(input, gradOutput);
  if(inplace) {
    TH_TENSOR_APPLY2(real, gradOutput, real, output,
      if(*output_data <= 0) {
        *gradOutput_data *= *output_data + la;
      } else {
        *gradOutput_data *= lambda;
      }
    );
    THTensor_(set)(gradInput, gradOutput);
  } else {
    THTensor_(resizeAs)(gradInput, output);
    TH_TENSOR_APPLY3(real, gradInput, real, gradOutput, real, output,
      *gradInput_data = *output_data <= 0 ? *gradOutput_data * (*output_data + la) : *gradOutput_data *lambda;
    );
  }
}

#endif
