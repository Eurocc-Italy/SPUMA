#ifndef SDV_TRACING_LOCAL_H
#define SDV_TRACING_LOCAL_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

// Room for the global regions plus the per-level ones (conversion, restrict,
// smooth, smooth_amul, smooth_eig, scale, scale_amul, prolong) of every GAMG
// level; levels are told apart by their number of cells only.
#define MAX_REGIONS 256
#define REGION_NAME_LEN 48

typedef struct {
    char name[REGION_NAME_LEN];
    uint64_t totalCycles;
    uint64_t lastStart;
    uint64_t calls;
    int active; // Open: a nested begin of the same region is ignored
} region_t;

extern region_t profile_regions[MAX_REGIONS];
extern int num_profile_regions;
// Set once the final report is registered: shared by every translation unit
// including this header, so the report runs (CSV + table) only once.
extern int trace_report_registered;
// MPI rank of this process (-1 in a serial run), set from UPstream when the
// first region is registered (gamgTimers.H); -2 while unknown.
extern int trace_rank;
// Set once the MAX_REGIONS warning has been printed.
extern int trace_regions_warned;

// Monotonic time source for the lightweight region profiler: CLOCK_MONOTONIC
// in nanoseconds (the "Cycles" columns of the report hold nanoseconds).
static inline uint64_t read_cycles_raw() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static void trace_report_final() {
    if (num_profile_regions == 0) return;

    // One file per rank in a parallel run, so the ranks do not overwrite
    // each other; only rank 0 prints the summary.
    const int rank = trace_rank;
    char path[64];
    if (rank >= 0) snprintf(path, sizeof path, "profiling_results_rank%d.csv", rank);
    else           snprintf(path, sizeof path, "profiling_results.csv");

    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "Region,TotalCycles,Calls,AvgCycles\n");
        for (int i = 0; i < num_profile_regions; i++) {
            if (profile_regions[i].calls > 0) {
                fprintf(f, "%s,%lu,%lu,%lu\n",
                    profile_regions[i].name, profile_regions[i].totalCycles,
                    profile_regions[i].calls, (profile_regions[i].totalCycles / profile_regions[i].calls));
            }
        }
        fclose(f);
    }

    if (rank > 0) return;

    // The stdout table shows the global regions only; the per-level ones
    // ("<region>@n<nCells>", which overlap them) go to the CSV alone.
    // The regions are nested (gv_vcycle and its parts lie inside gv_solve),
    // so each one is given as a share of gv_solve, not of their sum.
    // gv_construct (solver construction) lies outside gv_solve.
    uint64_t solveTotal = 0;
    for (int i = 0; i < num_profile_regions; i++) {
        if (strcmp(profile_regions[i].name, "gv_solve") == 0)
            solveTotal = profile_regions[i].totalCycles;
    }

    printf("\n==========================================================\n");
    printf("       GAMG REGION SUMMARY (rank 0)                       \n");
    printf("==========================================================\n");
    printf(" %-20s | %-12s | %s\n", "Region", "Avg ns", "% of gv_solve");
    printf("----------------------------------------------------------\n");

    for (int i = 0; i < num_profile_regions; i++) {
        if (profile_regions[i].calls > 0 && !strchr(profile_regions[i].name, '@')) {
            const char* name = profile_regions[i].name;
            const char* mark = (strcmp(name, "gv_construct") == 0) ? " *" : "";
            const int pad = 20 - (int)(strlen(name) + strlen(mark));
            const uint64_t avg = profile_regions[i].totalCycles / profile_regions[i].calls;
            printf(" %s%s%*s | %-12lu | ", name, mark, pad > 0 ? pad : 0, "", avg);
            if (solveTotal > 0)
                printf("%.2f\n", 100.0 * (double)profile_regions[i].totalCycles / (double)solveTotal);
            else
                printf("-\n");
        }
    }
    printf("----------------------------------------------------------\n");
    printf(" %% = region total time / gv_solve total time. Regions are\n");
    printf(" nested (gv_vcycle and its parts lie inside gv_solve), so\n");
    printf(" the column does not add up to 100.\n");
    printf(" * outside gv_solve (solver construction).\n");
    printf(" Per-level regions (<region>@n<cells>): see the CSV files.\n");
    printf("==========================================================\n\n");
}

static inline void register_final_report() {
    if (!trace_report_registered) { atexit(trace_report_final); trace_report_registered = 1; }
}

// Index of the region called `name`, created (with no calls yet) if new; -1
// once MAX_REGIONS are in use, which the begin/end calls below ignore.
static inline int trace_region_id(const char* name) {
    register_final_report();
    for (int i = 0; i < num_profile_regions; i++) {
        if (strcmp(profile_regions[i].name, name) == 0) return i;
    }
    if (num_profile_regions >= MAX_REGIONS) {
        if (!trace_regions_warned) {
            fprintf(stderr, "GAMG profiling: more than %d regions, '%s' and the"
                " next new ones are not recorded\n", MAX_REGIONS, name);
            trace_regions_warned = 1;
        }
        return -1;
    }
    int idx = num_profile_regions++;
    const size_t len = strnlen(name, REGION_NAME_LEN - 1);
    memcpy(profile_regions[idx].name, name, len);
    profile_regions[idx].name[len] = '\0';
    profile_regions[idx].totalCycles = 0;
    profile_regions[idx].calls = 0;
    profile_regions[idx].active = 0;
    return idx;
}

// Begin/end by index (no name lookup): the GV_TIMER macros of gamgTimers.H
// resolve the index once per call site, or once per level.
static inline void trace_begin_id(int id) {
    if (id < 0) return;
    region_t *r = &profile_regions[id];
    // Start only if not already active (GAMG recursion protection)
    if (!r->active) {
        r->lastStart = read_cycles_raw();
        r->active = 1;
        r->calls++;
    }
}

// Close region `id` at time `end` (taken by the caller before any lookup).
static inline void trace_end_at(int id, uint64_t end) {
    if (id < 0) return;
    region_t *r = &profile_regions[id];
    // Close only if it was active
    if (r->active) {
        r->totalCycles += end - r->lastStart;
        r->active = 0;
    }
}

static inline void trace_end_id(int id) {
    trace_end_at(id, read_cycles_raw());
}

#ifdef __cplusplus
}
#endif
#endif
