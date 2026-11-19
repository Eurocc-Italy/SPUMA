#include "sdv_tracing_local.H"
#include "gamgProfiling.H"

region_t profile_regions[MAX_REGIONS];
int num_profile_regions = 0;
int trace_report_registered = 0;
int trace_rank = -2;
int trace_regions_warned = 0;

bool Foam::gv::enabled = true;
