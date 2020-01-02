#include "rpf.h"

extern const struct rpf librpf_model[];
extern const int librpf_numModels;

const struct rpf *Glibrpf_model;
int Glibrpf_numModels;

void get_librpf_models(int version, int *numModels, const struct rpf **model)
{
  if (version != LIBIFA_RPF_API_VERSION) stop("LIBIFA_RPF binary API version mismatch");
  *numModels = librpf_numModels;
  *model = librpf_model;
}

// [[Rcpp::export]]
void registerCCallable()
{
    Glibrpf_numModels = librpf_numModels;
    Glibrpf_model = librpf_model;
    R_RegisterCCallable("rpf", "get_librpf_model_GPL", (DL_FUNC) get_librpf_models);
}
