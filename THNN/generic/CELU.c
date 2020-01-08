#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/CELU.c"
#else

void THNN_(CELU_updateOutput)(
          THNNState *state,
          THTensor *input,
          THTensor *output,
          accreal alpha_,
          bool inplace)
{
  real alpha = TH_CONVERT_ACCREAL_TO_REAL(alpha_);
  if(inplace) {
    TH_TENSOR_APPLY(real, input,
      if(*input_data <= 0) {
        *input_data = (exp(*input_data / alpha) - 1) * alpha;
      }
    );
    THTensor_(set)(output, input);
  } else {
    THTensor_(resizeAs)(output, input);
    TH_TENSOR_APPLY2(real, input, real, output,
      *output_data = *input_data <= 0 ? (exp(*input_data / alpha)-1)*alpha : *input_data;
    );
  }
}

void THNN_(CELU_updateGradInput)(
          THNNState *state,
          THTensor *input,
          THTensor *gradOutput,
          THTensor *gradInput,
          THTensor *output,
          accreal alpha_,
          bool inplace)
{
  real alpha = TH_CONVERT_ACCREAL_TO_REAL(alpha_);
  THNN_CHECK_NELEMENT(input, gradOutput);
  if(inplace) {
    TH_TENSOR_APPLY3(real, gradOutput, real, output, real, input,
      if(*output_data <= 0) {
        *gradOutput_data *= exp(*input_data / alpha);
      }
    );
    THTensor_(set)(gradInput, gradOutput);
  } else {
    THTensor_(resizeAs)(gradInput, output);
    TH_TENSOR_APPLY3(real, gradInput, real, gradOutput, real, output,
      *gradInput_data = *output_data <= 0 ? *gradOutput_data * (*output_data / alpha + 1) : *gradOutput_data;
    );
  }
}

#endif
