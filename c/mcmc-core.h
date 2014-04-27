/* mcmc_binner: metagenomic binning utility
Author: Andrey Kislyuk (kislyuk@gatech.edu)
*/

#ifndef MCMC_CORE_H
#define MCMC_CORE_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <fcntl.h>
#include <sys/types.h>
#include <math.h>
#include <float.h>
//#include <limits.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include "mcmc-mats.h"

#define DEFAULT_CHAIN_ORDER 3
#define MAX_CHAIN_ORDER 8
#define DEFAULT_STEPS 20000
#define MAX_STEPS 1e9
#define DEFAULT_RESTARTS 1
#define MAX_RESTARTS 10
#define DEFAULT_SOURCES 2
#define MAX_SOURCES 10000
#define FILENAME_SIZE 1020
#define SEQNAME_SIZE 256
#define MAX_MODEL_DEVIANCE 1e-9
#define MAX_FAILED_PERTURBS 1000
#define DEFAULT_MODEL_VAR 1e-6
#define DEFAULT_FREQ_VAR 5e-5
// Sorensen and Gianola p. 504
#define TARGET_ACCEPT_RATE 0.234

typedef unsigned int uint;

typedef struct {
	char * id;
	uint order;
	uint num_mers;
	double * mer_freqs;
} model;

typedef struct {
	char * id;
	uint order;
	uint num_sources;
	uint num_mers;
	double * source_freqs;
	model ** source_models;
	double logprob;
} ModelSet;

typedef struct {
	ModelSet ** models;
	uint * timestamps;
	uint num_models;
	uint max_models;
	int sp_start;
	int sp_start_timestamp;
} ModelRecord;

typedef struct {
	uint num_seqs;
	uint num_mers;
	char ** source_names;
	uint ** mer_counts;
} SeqStats;

typedef struct {
	uint rows;
	uint cols;
	double * contents;
} basis;

typedef struct {
	uint failed_perturbs;
	uint acceptances;
	uint backward_steps;
	uint double_underflows;
} stats;

typedef struct {
	double model_var;
	double freq_var;
	uint chain_order;
	uint num_sources;
	uint max_failed_perturbs;
	uint steps;
	uint restarts;
	char model_logfile[FILENAME_SIZE+4];
	char model_hintfile[FILENAME_SIZE+4];
	char thread_id;
	stats stats;
} settings;

const double ord_bounds[] = {0, 4, 20, 84, 340, 1364, 5460};
const char nucs[] = {'A', 'T', 'G', 'C'};
const char rc_nucs[] = {'T', 'A', 'C', 'G'};

#endif /* MCMC_CORE_H */
