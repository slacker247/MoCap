#include <qam_soft_decision.h>
#include <qam_conv_decoder.h>
using namespace SPUC;
// Create puncture inputs for Viterbi
long* qam_conv_decoder::clear_soft_decision(long rate)
{
  if (rate) {
	for (int i=0;i<2*rate;i++) viterbi_input[i] =  0; 
  } else {
	viterbi_input[0] = 0; // bpsk
  }
  return(viterbi_input);
}
// Data decoder
bool qam_conv_decoder::data_decode(complex<long> data_in)
{
  int j;
  qam_data_demap(rate_index, data_in, soft_decision_level, viterbi_input);

  // check viterbi_ready & viterbi_data for outputs
  for (j=0;j<rx_bits_per_symbol;j++) {
	if (no_conv) {
	  viterbi_ready = 1;
	  viterbi_data = (viterbi_input[j]<0) ? 1:0;
	} else {
	  viterbi_data = viterbi_decoder.depuncture(enc_rate,viterbi_input[j]);
	  viterbi_ready = viterbi_decoder.output_ready;
	}
  }
  bool ready = viterbi_ready;
}
