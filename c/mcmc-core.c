/* mcmc_binner: metagenomic binning utility
Author: Andrey Kislyuk (kislyuk@gatech.edu)
*/

#include "mcmc-core.h"

// NB: lowercase model means a K-mer model representing a single source
// Uppercase ModelSet means an aggregate model of a set of sources. It contains a set of lowercase models

model * newmodel(uint order, char * id) {
	uint i, j;
	model * m = (model *)malloc(sizeof(model)); assert(m);
	m->id = (char *)malloc(sizeof(char) * 128); assert(m->id);
	strncpy(m->id, id, 128); m->id[127] = 0;
	m->order = order;
	m->num_mers = 0; for (i=1; i<=order; i++) { m->num_mers += pow(4, i); }
	m->mer_freqs = (double *)malloc(sizeof(double) * m->num_mers); assert(m->mer_freqs);
	for (i=0; i<order; i++) {
		for (j=ord_bounds[i]; j<ord_bounds[i+1]; j++) {
			m->mer_freqs[j] = 1.0/pow(4, i+1);
		}
	}
	return m;
}

ModelSet * newModelSet(char * id, uint order, uint num_sources) {
	uint i;
	ModelSet * M = (ModelSet *)malloc(sizeof(ModelSet)); assert(M);
	M->id = (char *)malloc(sizeof(char) * 128); assert(M->id);
	strncpy(M->id, id, 128); M->id[127] = 0;
	M->order = order;
	M->num_sources = num_sources;
	M->num_mers = 0; for (i=1; i<=order; i++) { M->num_mers += pow(4, i); }
	M->source_freqs = (double *)malloc(sizeof(double) * num_sources); assert(M->source_freqs);
	M->source_models = (model **)malloc(sizeof(void *) * num_sources); assert(M->source_models);
	for (i=0; i<num_sources; i++) {
		M->source_models[i] = newmodel(order, "");
		M->source_freqs[i] = 1.0/num_sources;
	}
	return M;
}

model * copymodel(model * old_m) {
	uint i, j;
	model * m = (model *)malloc(sizeof(model)); assert(m);
	m->id = (char *)malloc(sizeof(char) * 128); assert(m->id);
	strncpy(m->id, old_m->id, 128); m->id[127] = 0;
	m->order = old_m->order;
	m->num_mers = old_m->num_mers;
	m->mer_freqs = (double *)malloc(sizeof(double) * m->num_mers); assert(m->mer_freqs);
	for (i=0; i<m->order; i++) {
		for (j=ord_bounds[i]; j<ord_bounds[i+1]; j++) {
			m->mer_freqs[j] = old_m->mer_freqs[j];
		}
	}
	return m;
}

ModelSet * copyModelSet(ModelSet * old_ModelSet) {
	uint i;
	ModelSet * M = (ModelSet *)malloc(sizeof(ModelSet)); assert(M);
	M->id = (char *)malloc(sizeof(char) * 128); assert(M->id);
	strncpy(M->id, old_ModelSet->id, 128); M->id[127] = 0;
	M->order = old_ModelSet->order;
	M->num_sources = old_ModelSet->num_sources;
	M->num_mers = old_ModelSet->num_mers;
	M->logprob = old_ModelSet->logprob;
	M->source_freqs = (double *)malloc(sizeof(double) * M->num_sources); assert(M->source_freqs);
	M->source_models = (model **)malloc(sizeof(void *) * M->num_sources); assert(M->source_models);
	for (i=0; i<M->num_sources; i++) {
		M->source_models[i] = copymodel(old_ModelSet->source_models[i]);
		M->source_freqs[i] = old_ModelSet->source_freqs[i];
	}
	return M;
}

// Make a new Model and set all its sources' parameters to those of the overall data mean
ModelSet * newModelSetAtMean(char * id, uint order, uint num_sources, SeqStats * C) {
	uint i, j, k;
	ModelSet * M = newModelSet(id, order, num_sources);
	double mean_mer_freqs[C->num_mers], total_mers[M->order];

	for (i=0; i<C->num_mers; i++) { mean_mer_freqs[i] = 0; }
	for (i=0; i<M->order; i++) { total_mers[i] = 0; }

	for (i=0; i<M->order; i++) {
		for (j=ord_bounds[i]; j<ord_bounds[i+1]; j++) {
			for (k=0; k<C->num_seqs; k++) {
				mean_mer_freqs[j] += C->mer_counts[k][j];
			}
			total_mers[i] += mean_mer_freqs[j];
		}
		for (j=ord_bounds[i]; j<ord_bounds[i+1]; j++) {
			if (total_mers[i] > 0.0) {
				mean_mer_freqs[j] /= total_mers[i];
			}
		}
	}
	for (i=0; i<M->num_sources; i++) {
		for (j=0; j<M->num_mers; j++) {
			M->source_models[i]->mer_freqs[j] = mean_mer_freqs[j];
		}
	}
	return M;
}

ModelRecord * newModelRecord(void) {
	ModelRecord * MR = (ModelRecord *)malloc(sizeof(ModelRecord)); assert(MR);

	MR->max_models = 64;
	MR->models = (void *)malloc(sizeof(void *) * MR->max_models); assert(MR->models);
	MR->timestamps = (uint *)malloc(sizeof(uint) * MR->max_models); assert(MR->timestamps);
	MR->num_models = 0;
	MR->sp_start = MR->sp_start_timestamp = -1;

	return MR;
}

// FIXME
void destroyModelRecord(ModelRecord * MR) {
	int i;
	for (i=0; i<MR->max_models; i++) {
//		destroyModel(MR->models[i]);
	}
//	free(MR->timestamps);
	free(MR);
}

// timestep: step at which model was generated
ModelRecord * addModelSetToRecord(ModelSet * M, int timestep, ModelRecord * MR) {
	if (MR->num_models == MR->max_models) {
		MR->max_models *= 2;
		MR->models = (void *)realloc(MR->models, sizeof(void *) * MR->max_models); assert(MR->models);
		MR->timestamps = (uint *)realloc(MR->timestamps, sizeof(uint) * MR->max_models); assert(MR->timestamps);
	}

	MR->models[MR->num_models] = M;
	MR->timestamps[MR->num_models] = timestep;
	MR->num_models++;
	return MR;
}

/* Performs the following mapping (subject to ordering in nucs[]):
 1 => A, 2 => T, 3 => G, 4 => C,
 5 => AA, 6 => AT, 7 => AG, 8 => AC, 9 => TA, ... 20 => CC,
 21 => AAA, ... */
char * index2nucmer(uint index, char * nucmer) {
	int i, length=0, base_i=0;
	assert(index > 0);
	while (1) {
		int upper_bd = 0;
		for (i=1; i<=length; i++) { upper_bd += pow(4, i); }
		if (index-1 < upper_bd) { break; } else { length++; }
	}
	assert(length <= MAX_CHAIN_ORDER);
	for (i=1; i<length; i++) { base_i += pow(4, i); }
	index -= base_i+1;
	for (i=length; i>0; i--) {
		int stride = pow(4, i-1);
		nucmer[length-i] = nucs[index/stride];
		index %= stride;
	}
	nucmer[length] = '\0';
	return nucmer;
}

// Performs the reverse mapping
uint nucmer2index(char * nucmer) {
	uint index=0, stride=0, i, j;
	for (i=1; i<strlen(nucmer); i++) {
		index += pow(4, i);
	}
	for (i=0; i<strlen(nucmer); i++) {
		for (j=0; j<4; j++) {
			if (nucmer[i] == nucs[j]) {
				stride = pow(4, strlen(nucmer)-i-1);
				index += stride*j;
			}
		}
	}
	return index;
}

// Print a human-readable model 
void printModelSet(ModelSet * M, FILE * stream) {
	uint i, j, n1, n2;
	char nucmer[MAX_CHAIN_ORDER];
	if (stream == NULL) { stream = stderr; }
	for (i=0; i<M->num_sources; i++) {
//		fprintf(stream, "LL=
		fprintf(stream, "Source %d: ID %s, order %d, frequency %f\n", i, M->source_models[i]->id, M->order, M->source_freqs[i]);
		fprintf(stream, "1 2>"); for (j=0; j<4; j++) { fprintf(stream, "\t%c", nucs[j]); } fprintf(stream, "\n");
		for (n1=0; n1<4; n1++) {
			fprintf(stream, "%c\t", nucs[n1]);
			for (n2=0; n2<4; n2++) {
				sprintf(nucmer, "%c%c", nucs[n1], nucs[n2]);
				fprintf(stream, "%.4f\t", M->source_models[i]->mer_freqs[nucmer2index(nucmer)]);
			}
			fprintf(stream, "\n");
		}
		fprintf(stream, "%%G+C = %.4f\n", M->source_models[i]->mer_freqs[nucmer2index("G")] + M->source_models[i]->mer_freqs[nucmer2index("C")]);
	}
}

// Dump a machine readable list of source frequencies and corresponding k-mer frequency vectors
void dumpModelSet(ModelSet * M, FILE * f) {
	uint i, j;
	fprintf(f, "%.16f\n", M->logprob);
	for (i=0; i<M->num_sources; i++) {
		fprintf(f, "%.16f\n", M->source_freqs[i]);
		for (j=0; j<M->source_models[i]->num_mers; j++) {
			fprintf(f, "\t%.16f\n", M->source_models[i]->mer_freqs[j]);
		}
	}
}

// Load a machine readable list of source frequencies and corresponding k-mer frequency vectors
ModelSet * loadModelSet(FILE * f, settings * s) {
	ModelSet * M = newModelSet("im", s->chain_order, s->num_sources);
	uint cur_source, cur_mer;
	float token;

	for (cur_source=0; cur_source < M->num_sources; cur_source++) {
		if (fscanf(f, "%f", &token) == 1) {
			M->source_freqs[cur_source] = (double)token;
//fprintf(stderr, "scanned in %f (in source %d)\n", token, cur_source);
			for (cur_mer=0; cur_mer < M->num_mers; cur_mer++) {
				if (fscanf(f, "%f", &token) == 1) {
//fprintf(stderr, "\tscanned in %f (in source %d, mer %d)\n", token, cur_source, cur_mer);
					M->source_models[cur_source]->mer_freqs[cur_mer] = (double)token;
				} else {
					fprintf(stderr, "[t%c] Error loading model\n", s->thread_id);
					abort();
				}
			}
			if (cur_source > s->num_sources+1) {
				fprintf(stderr, "[t%c] Error loading model: Too many sources in input\n", s->thread_id);
				abort();
			}
		} else {
			fprintf(stderr, "[t%c] Error loading model\n", s->thread_id);
			abort();
		}
	}
	return M;
}


void destroymodel(model * m) {
	free(m->mer_freqs);
	free(m);
}

void destroyModelSet(ModelSet * M) {
	uint i;
	free(M->id);
	free(M->source_freqs);
	for (i=0; i<M->num_sources; i++) {
		destroymodel(M->source_models[i]);
	}
	free(M);
}

// A model data structure with frequencies in log space, for caching purposes
model * logmodel(model * m) {
	uint i;
	model * logm = newmodel(m->order, m->id);
	for (i=0; i<ord_bounds[m->order]; i++) {
		logm->mer_freqs[i] = log(m->mer_freqs[i]);
	}
	return logm;
}

// Get log probability of seeing this sequence given this model
double getSeqLogProb(uint * mer_counts, model * log_m, uint order) {
	uint i;
	double log_prob = 0.0;
	double * mer_freqs = log_m->mer_freqs;
	for (i=ord_bounds[order-1]; i<ord_bounds[order]; i++) {
		log_prob += mer_freqs[i] * mer_counts[i];
	}
	if (order > 1) { // Divide by marginals of the next highest order
		for (i=ord_bounds[order-2]; i<ord_bounds[order-1]; i++) {
			log_prob -= mer_freqs[i] * mer_counts[i];
		}
	}
	return log_prob;
}

// Get aggregate log probability of seeing all given sequences given this Model
// TODO: priors support
double getModelSetLogProb(ModelSet * M, SeqStats * C) {
	uint max_logprob_index=0, i, j;
	double logP = 0, max_logprob, prob_partsum, seq_logprobs[M->num_sources];

	model ** log_models = (model **)malloc(sizeof(void *) * M->num_sources); assert(log_models);
	for (i=0; i<M->num_sources; i++) {
		log_models[i] = logmodel(M->source_models[i]);
	}

	for (i=0; i<C->num_seqs; i++) {
		max_logprob = -1e999; //DBL_MIN;
		for (j=0; j<M->num_sources; j++) {
			seq_logprobs[j] = getSeqLogProb(C->mer_counts[i], log_models[j], M->order);
			if (max_logprob < seq_logprobs[j]) {
				max_logprob = seq_logprobs[j];
				max_logprob_index = j;
			}
		}
		// TODO: FREQUENCIES
		prob_partsum = 0;
		for (j=0; j<M->num_sources; j++) {
			prob_partsum += exp(seq_logprobs[j] - seq_logprobs[max_logprob_index]);
		}
		logP += seq_logprobs[max_logprob_index] + log(prob_partsum);
	}

	for (i=0; i<M->num_sources; i++) {
		destroymodel(log_models[i]);
	}
	free(log_models);
	M->logprob = logP;
	return logP;
}

// Read kmer frequencies from standard input (program input)
SeqStats * readSeqStats(uint order) {
	uint i;
	char buf[1024], seqname[SEQNAME_SIZE];
	char * p;
	FILE * stream = stdin;
	int max_seqs = 1024;
	SeqStats * s = (SeqStats *)malloc(sizeof(SeqStats)); assert(s);

	s->num_seqs = 0;
	s->num_mers = 0; for (i=1; i<=order; i++) { s->num_mers += pow(4, i); }
	s->source_names = malloc(sizeof(void *) * max_seqs); assert(s->source_names);
	s->mer_counts = malloc(sizeof(void *) * max_seqs); assert(s->mer_counts);
	while (fgets(seqname, SEQNAME_SIZE, stream)) {
		p = strrchr(seqname, '\n'); *p = 0; // chomp
		s->source_names[s->num_seqs] = malloc(sizeof(char) * (strlen(seqname)+1)); assert(s->source_names[s->num_seqs]);
		strcpy(s->source_names[s->num_seqs], seqname);
		fgets(buf, 1024, stream);
		p = strrchr(buf, '\n'); *p = 0; // chomp
		s->mer_counts[s->num_seqs] = malloc(sizeof(uint) * s->num_mers); assert(s->mer_counts[s->num_seqs]);

		i = 0;
		p = strtok(buf, " ");
		while (p) {
			assert(i < s->num_mers);
			s->mer_counts[s->num_seqs][i] = atoi(p);
			p = strtok(NULL, " ");
			i++;
		}

		s->num_seqs++;
		if (s->num_seqs >= max_seqs) {
			max_seqs *= 2;
			s->source_names = (void *) realloc(s->source_names, sizeof(void *) * max_seqs); assert(s->source_names);
			s->mer_counts = (void *) realloc(s->mer_counts, sizeof(void *) * max_seqs); assert(s->mer_counts);
		}
	}

	return s;
}

static double leftover_gaussian = 0.0;

// Box-Muller transform (http://www.taygeta.com/random/gaussian.html)
double gaussian(double mean, double variance) {
	double x1, x2, w, y1, y2;
	if (leftover_gaussian) {
		y1 = leftover_gaussian;
		leftover_gaussian = 0;
		return y1;
	}
	if (!variance) { variance = 1; }

	do {
		x1 = 2.0 * ((double)rand() / ((double)(RAND_MAX)+(double)(1))) - 1.0;
		x2 = 2.0 * ((double)rand() / ((double)(RAND_MAX)+(double)(1))) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0);

	w = sqrt((-2.0 * log(w)) / w);
	y1 = (x1 * w) * variance + mean;
	y2 = (x2 * w) * variance + mean;
	leftover_gaussian = y2;
	return y1;
}

// Get an element of the basis matrix A
double getNullMatElt(uint order, uint row, uint col) {
	double * N = Nmats[order-1];
	uint stride = Ndims[order-1][1];
	assert((row * stride) + col < Ndims[order-1][0] * Ndims[order-1][1]);
	return N[(row * stride) + col];
}

basis * newBasis(settings * s) {
	basis * b = malloc(sizeof(basis)); assert(b);
	b->rows = Ndims[s->chain_order-1][0];
	b->cols = Ndims[s->chain_order-1][1];
	b->contents = malloc(sizeof(double) * b->rows * b->cols); assert(b->contents);
	return b;
}

void setBasisElt(basis * b, uint row, uint col, double value) {
	b->contents[(row * b->cols) + col] = value;
}

double getBasisElt(basis * b, uint row, uint col) {
	return b->contents[(row * b->cols) + col];
}

void destroyBasis(basis * b) {
	free(b->contents);
	free(b);
}

model * perturbmodel(model * old_m, model * new_m, settings * s) {
	uint i, j, failed_perturbs=0;
	basis * b = newBasis(s);

	assert(new_m->num_mers == old_m->num_mers);
	for (i=0; i<new_m->num_mers; i++) {
		new_m->mer_freqs[i] = old_m->mer_freqs[i];
	}

	TRY: while (1) {
		if (failed_perturbs > MAX_FAILED_PERTURBS) {
			new_m = old_m; // unable to perturb
			break;
		}
		for (i=0; i<b->rows; i++) {
			double r = gaussian(0.0, s->model_var);
			for (j=0; j<b->cols; j++) { setBasisElt(b, i, j, getNullMatElt(s->chain_order, i, j) * r); }
		}
		assert(b->cols == new_m->num_mers);
		for (i=0; i<b->rows; i++) {
			for (j=0; j<b->cols; j++) {
				assert(getBasisElt(b, i, j) != 0);
//fprintf(stderr, "mer %d, bv %d: %.12f += %.12f\n", j, i, new_m->mer_freqs[j], getBasisElt(b, i, j));
				new_m->mer_freqs[j] += getBasisElt(b, i, j); // (b, i, j);
				if (new_m->mer_freqs[j] <= 0 || new_m->mer_freqs[j] >= 1) {
					failed_perturbs++;
//fprintf(stderr, "Perturbation out of bounds at %d, %d: %f\n", i, j, new_m->mer_freqs[j]);
					goto TRY;
				}
			}
		}
		break;
	}
	destroyBasis(b);
	s->stats.failed_perturbs += failed_perturbs;
	return new_m;
}

// TODO: calibration and basis multiplication
double * perturbfreqs(double * old_freqs, double * new_freqs, settings * s) {
	uint i;

	TRY: while (1) {
		// Unbiased random selection of array index
		uint selected_f = (uint)(((double)rand()/((double)(RAND_MAX)+(double)(1))) * (double)s->num_sources);
		double r = gaussian(0, s->freq_var);
		new_freqs[selected_f] = old_freqs[selected_f] + r;
		for (i=0; i<s->num_sources; i++) {
			if (i == selected_f) { continue; }
			new_freqs[i] = old_freqs[i] - (r / (s->num_sources - 1));
		}
		for (i=0; i<s->num_sources; i++) {
			if (new_freqs[i] <= 0 || new_freqs[i] >= 1) { goto TRY; }
		}
		break;
	}

	//double total=0, tol=1e-100;
	//for (int i=0; i<s.num_sources; i++) { total += new_freqs[i]; assert(total - 1.0 < tol); }
	return new_freqs;
}

ModelSet * perturbModelSet(ModelSet * M, settings * s) {
	uint i;
	ModelSet * newM = newModelSet("", M->order, M->num_sources);
	for (i=0; i<M->num_sources; i++) {
		if (perturbmodel(M->source_models[i], newM->source_models[i], s) == M->source_models[i]) {
			destroyModelSet(newM);
			return M; // unable to perturb
		}
	}
	perturbfreqs(M->source_freqs, newM->source_freqs, s);
	return newM;
}

// Seed the pseudo-random number generator
void seedPRNG(void) {
	int fd, buf;
	if((fd = open("/dev/random", O_RDONLY)) < 0) {
		perror("/dev/random");
		exit(1);
	}
	read(fd, &buf, sizeof(buf));
	close(fd);
	srand(buf);
}

// Parse the command line options
void readOptions(int argc, char ** argv, settings * s) {
	int i;
	char * optspec = "o:n:s:r:m:f:l:i:h:";
	char * opts_with_reqd_args = "onsrmflih";
	extern char *optarg;
	extern int optind, opterr, optopt;
	while ((i = getopt(argc, argv, optspec)) != -1) {
		switch (i) {
			case 'o':
				s->chain_order = atoi(optarg);
				if (s->chain_order < 1 || s->chain_order > MAX_CHAIN_ORDER) {
					fprintf(stderr, "Chain order out of range\n"); abort();
				}
				break;
			case 'n':
				s->num_sources = atoi(optarg);
				if (s->num_sources < 1 || s->num_sources > MAX_SOURCES) {
					fprintf(stderr, "Number of sources out of range\n"); abort();
				}
				break;
			case 's':
				s->steps = atoi(optarg);
				if (s->steps < 1 || s->steps > MAX_STEPS) {
					fprintf(stderr, "Number of steps out of range\n"); abort();
				}
				break;
			case 'r':
				s->restarts = atoi(optarg);
				if (s->restarts > MAX_RESTARTS) {
					fprintf(stderr, "Number of restarts out of range\n"); abort();
				}
				break;
			case 'm':
				s->model_var = atof(optarg);
				break;
			case 'i':
				s->thread_id = optarg[0];
				break;
			case 'f':
				s->freq_var = atof(optarg);
				break;
			case 'l':
				strncpy(s->model_logfile, optarg, FILENAME_SIZE); s->model_logfile[FILENAME_SIZE] = 0;
				break;
			case 'h':
				strncpy(s->model_hintfile, optarg, FILENAME_SIZE); s->model_hintfile[FILENAME_SIZE] = 0;
				break;
			case '?':
				if (strchr(opts_with_reqd_args, optopt)) {
					fprintf(stderr, "Option -%c requires an argument\n", optopt);
				} else {
					fprintf(stderr, "Encountered unknown option\n");
				}
				abort();
			default:
				assert(0);
		}
	}
	if (optind < argc) {
		fprintf(stderr, "Encountered unknown option\n");
		abort();
	}
}

// Compute second-order deviations from identities specified by the k-mer dimension redundancy model
double checkModelSetDeviance(ModelSet * M) {
	double total_deviance=0;
	uint i, j;
	if (M->order < 2) { return 0.0; }
	for (i=0; i<M->num_sources; i++) {
		double * m = M->source_models[i]->mer_freqs;
		double deviances[4];
		deviances[0] = m[nucmer2index("AG")] - (m[nucmer2index("A")] - m[nucmer2index("AA")] - m[nucmer2index("AC")] - m[nucmer2index("AT")]);
		deviances[1] = m[nucmer2index("CC")] - (0.5 - 2*m[nucmer2index("A")] + m[nucmer2index("AA")]
			+ m[nucmer2index("AC")] + m[nucmer2index("AT")] - m[nucmer2index("CA")] - m[nucmer2index("CG")]);
		deviances[2] = m[nucmer2index("TA")] - (2*m[nucmer2index("AC")] + m[nucmer2index("AT")]
			- 2*m[nucmer2index("CA")] - m[nucmer2index("CG")] + m[nucmer2index("GC")]);
		deviances[3] = m[nucmer2index("TC")] - (m[nucmer2index("A")] - m[nucmer2index("AA")] - 2*m[nucmer2index("AC")]
			- m[nucmer2index("AT")] + m[nucmer2index("CA")] + m[nucmer2index("CG")] - m[nucmer2index("GC")]);

		for (j=0; j<4; j++) { total_deviance += fabs(deviances[j]); }
	}
	return total_deviance;
}

// Print a list of models to s->model_logfile in a format directly loadable by Perl
// This is used by the visualization functions only
void printModelRecord(ModelRecord * MR, settings * s) {
	uint i, j, k;
	FILE * record = fopen(s->model_logfile, "w");
	if (record == NULL) { fprintf(stderr, "Warning: Unable to open file \"%s\" for writing\n", s->model_logfile); return; }

	fprintf(record, "$VAR1 = { 'sp_start_timestamp' => %d,\n", MR->sp_start_timestamp);
	for (i=0; i<MR->num_models; i++) {
//		fprintf(record, "\t'%d' => { loglik => %f,\n\t\tmodels => [", model_timestamps[i], Model_record[i]->logprob);
		fprintf(record, "\t'%d' => { loglik => %f,\n\t\tmodels => [", MR->timestamps[i], MR->models[i]->logprob);
		for (j=0; j<MR->models[i]->num_sources; j++) {
			fprintf(record, "\t\t{\n");
			for (k=0; k<MR->models[i]->num_mers; k++) {
				char n[MAX_CHAIN_ORDER+1];
//				fprintf(record, "\t\t\t'%s' => '%.12e',\n", index2nucmer(k+1, n), Model_record[i]->source_models[j]->mer_freqs[k]);
				fprintf(record, "\t\t\t'%s' => '%.12e',\n", index2nucmer(k+1, n), MR->models[i]->source_models[j]->mer_freqs[k]);
			}
			fprintf(record, "\t\t},\n");
		}
		fprintf(record, "\t], freqs => [] },\n");
	}
	fprintf(record, "};\n");
	fclose(record);
	fprintf(stderr, "[t%c] Model record printed to %s\n", s->thread_id, s->model_logfile);
}

// Subtract one model from the other
// Models must be of the same order, source count, and source ordering, so adaptive source jumping will need work
double ModelSetDistance(ModelSet * ModelSet1, ModelSet * ModelSet2) {
	int i, j;
//	Model * M = newModel("diff", Model1->order, Model1->num_sources);
	double sumsqdiffs = 0;
	assert(ModelSet1->num_sources == ModelSet2->num_sources);
	assert(ModelSet1->num_mers == ModelSet2->num_mers);

	for (i=0; i<ModelSet1->num_sources; i++) {
		double m_sumsqdiffs = 0;
//		for (j=0; j<Model1->source_models[i]->num_mers; j++) {
		for (j=ord_bounds[ModelSet1->order - 1]; j<ord_bounds[ModelSet1->order]; j++) {
			double diff = ModelSet2->source_models[i]->mer_freqs[j] - ModelSet1->source_models[i]->mer_freqs[j];
			m_sumsqdiffs += diff * diff;
//fprintf(stderr, "md: %d %d: %e, sum = %e\n", i, j, diff, m_sumsqdiffs);
		}
		sumsqdiffs += m_sumsqdiffs;
	}
	return sumsqdiffs;
}

/* Set and return the first index in Model_record corresponding to a stabilized plateau (stationary phase),
defined as a sequence of models with maximum distance between each not exceeding 5% of the initial.
*/
int findStationaryPhase(ModelRecord * MR, settings * s) {
	int i;
	// TODO: calibrate this number - it should vary with absolute distance between stable freq and init freq
	double max_dist_frac = 0.005;
	double d0 = ModelSetDistance(MR->models[MR->num_models-1], MR->models[0]);

	for (i = MR->num_models/2; i < MR->num_models-1; i += (MR->num_models-i)/2) {
		double cur_dist = ModelSetDistance(MR->models[MR->num_models-1], MR->models[i]);
		if (cur_dist/d0 < max_dist_frac) {
			fprintf(stderr, "Stationary phase clamped at %d of %d (d=%f)\n", i, MR->num_models, cur_dist);
			MR->sp_start = i;
			MR->sp_start_timestamp = MR->timestamps[MR->sp_start];
			return MR->sp_start;
		}
	}
	fprintf(stderr, "[t%c] Warning: unable to clamp stationary phase\n", s->thread_id);
	// failed to clamp
	MR->sp_start = MR->num_models;
	MR->sp_start_timestamp = MR->timestamps[MR->sp_start];
	return MR->sp_start;
}

ModelSet * ModelSetAverage(ModelSet * ModelSet1, ModelSet * ModelSet2) {
	int i, j;
	ModelSet * M = newModelSet("foo", ModelSet1->order, ModelSet1->num_sources);
	assert(ModelSet1->num_sources == ModelSet2->num_sources);

	for (i=0; i<ModelSet1->num_sources; i++) {
		for (j=0; j<ModelSet1->source_models[i]->num_mers; j++) {
			M->source_models[i]->mer_freqs[j] = (ModelSet1->source_models[i]->mer_freqs[j] + ModelSet2->source_models[i]->mer_freqs[j])/2;
		}
		M->source_freqs[i] = (ModelSet1->source_freqs[i] + ModelSet2->source_freqs[i]) / 2;
	}

	// TODO: freqs averaging
	assert(0);

	return M;
}

// Returns a model averaging all models in Model_record corresponding to the stationary phase,
// starting with stationaryStart and up to the end of the simulation
ModelSet * sampleStationaryDistribution(ModelRecord * MR, SeqStats * C, settings * s) {
	int i, j, k, num_samples = 64;
	assert(MR->num_models > 0);
	ModelSet * M = newModelSet("foo", MR->models[0]->order, MR->models[0]->num_sources);
	if (MR->sp_start == -1) {
		findStationaryPhase(MR, s);
	}
	assert(MR->sp_start >= 0 && MR->sp_start < MR->num_models);
	if (num_samples > MR->num_models - MR->sp_start) { num_samples = MR->num_models - MR->sp_start; }

	for (j=0; j<MR->models[0]->num_sources; j++) {
		for (k=0; k<MR->models[0]->source_models[j]->num_mers; k++) {
			M->source_models[j]->mer_freqs[k] = 0;
		}
		M->source_freqs[j] = 0;
	}

	// TODO: sample without replacement, not with
	for (i=0; i<num_samples; i++) {
		int r = MR->sp_start + (rand() % (MR->num_models - MR->sp_start)); // doesn't need much randomness
		for (j=0; j<MR->models[0]->num_sources; j++) {
			for (k=0; k<MR->models[0]->source_models[j]->num_mers; k++) {
				M->source_models[j]->mer_freqs[k] += MR->models[r]->source_models[j]->mer_freqs[k];
			}
			M->source_freqs[j] += MR->models[r]->source_freqs[j];
		}
	}

	for (j=0; j<MR->models[0]->num_sources; j++) {
		for (k=0; k<MR->models[0]->source_models[j]->num_mers; k++) {
			M->source_models[j]->mer_freqs[k] /= num_samples;
		}
		M->source_freqs[j] /= num_samples;
	}

	M->logprob = getModelSetLogProb(M, C);

	return M;
}

ModelSet * getBestModelSet(ModelRecord * MR, SeqStats * C, settings * s) {
//	return sampleStationaryDistribution(MR, C, s); // SP sampling
//	return MR->models[MR->num_models-1]; // last model (not necessarily highest LP)

	double max_modelset_logprob = -1e999;
	int i, best_modelset_i=-1;

	for (i = 0; i < MR->num_models-1; i++) {
		if (max_modelset_logprob < MR->models[i]->logprob) {
			max_modelset_logprob = MR->models[i]->logprob;
			best_modelset_i = i;
		}
	}
	assert(best_modelset_i >= 0);
	return MR->models[best_modelset_i];
}

void calibrateVariance(ModelSet * init_ModelSet, SeqStats * C, settings * s) {
	int i, j;
	settings s_copy = *s;
	s_copy.steps /= 2; // WARNING: HACK
	double best_model_var=DEFAULT_MODEL_VAR, best_acc_delta = 1e999; //DBL_MAX;

	double max_model_var = 1e-3, min_model_var = 1e-8;
	for (s_copy.model_var = max_model_var; s_copy.model_var >= min_model_var; s_copy.model_var /= 4) {
		ModelSet * M = copyModelSet(init_ModelSet);
		double old_prob = getModelSetLogProb(M, C);

		for (i=0; i<s_copy.steps; i++) {
			ModelSet * new_M = perturbModelSet(M, &s_copy);
			if (new_M == M) { // failed to perturb
				destroyModelSet(M);
				M = NULL;
				break;
			}

			double new_prob = getModelSetLogProb(new_M, C);
			double prob_log_ratio = new_prob - old_prob;
			double r;
			// Random double between 0 and 1 with granularity 1/RAND_MAX
			while ((r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)))) == 0.0);
			double rand_cutoff = log(r);

			if (prob_log_ratio > rand_cutoff) {
				destroyModelSet(M);
				M = new_M;
				old_prob = new_prob;
				s_copy.stats.acceptances++;
			} else {
				destroyModelSet(new_M);
			}
		}
		if (!M) {
			fprintf(stderr, "[t%c] var %e too high, skipping\n", s_copy.thread_id, s_copy.model_var);
			continue;
		} // variance too high, continue search at lower var

		double acc_rate = (double)s_copy.stats.acceptances/(double)s_copy.steps;
		double acc_delta = fabs(acc_rate - TARGET_ACCEPT_RATE);
		if (best_acc_delta > acc_delta) {
			fprintf(stderr, "[t%c] Delta %e replaces %e, var %e replaces %e\n", s_copy.thread_id, acc_delta, best_acc_delta, s_copy.model_var, best_model_var);
			best_acc_delta = acc_delta;
			best_model_var = s_copy.model_var;
		} else if (acc_rate > TARGET_ACCEPT_RATE) { // variance lower than optimal, terminate search
			fprintf(stderr, "[t%c] Terminated calibration loop: accept rate %e > %e; optimal var %e, delta %e\n", s_copy.thread_id, acc_rate, TARGET_ACCEPT_RATE, best_model_var, best_acc_delta);
			break;
		}
//		fprintf(stderr, "[t%c] Calibrator: At v=%e: %d / %d = %f\n", s_copy.thread_id, s_copy.model_var, s_copy.stats.acceptances, s_copy.steps, (double)s_copy.stats.acceptances/(double)s_copy.steps);
		s_copy.stats.acceptances = 0; s_copy.stats.failed_perturbs = 0;
	}
//	destroyModel(M);
	fprintf(stderr, "[t%c] Calibrator: Best variance %e, delta %e\n", s_copy.thread_id, best_model_var, best_acc_delta);

	s->model_var = best_model_var;// / 3; // WARNING: HACK. The steps and var estimate divisors depend on external parameters and bias variance estimation.
}

/* Run the MCMC search for s->steps steps.
Return null if the run fails (gets stuck)
*/
ModelRecord * runMCMC(ModelSet * init_ModelSet, ModelRecord * MR, SeqStats * C, settings * s) {
	int i, j;
	addModelSetToRecord(init_ModelSet, 0, MR);
	fprintf(stderr, "[t%c] Initial model:\n", s->thread_id);
	printModelSet(MR->models[0], stderr);

	double old_prob = getModelSetLogProb(MR->models[0], C);

	for (i=0; i<s->steps; i++) {
		ModelSet * new_M = perturbModelSet(MR->models[MR->num_models-1], s);

		if (new_M == MR->models[MR->num_models-1]) { // perturb failed
			return NULL;
		}

		double new_prob = getModelSetLogProb(new_M, C);
		double prob_log_ratio = new_prob - old_prob;
		double r;
		// Random double between 0 and 1 with granularity 1/RAND_MAX
		while ((r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)))) == 0.0);
		double rand_cutoff = log(r);

		if (i % 1000 == 0) { fprintf(stderr, "[t%c] [%d/%d] (LR=%.2f - %.2f = %.4f, cutoff %.2f)\n", s->thread_id, i+1, s->steps, new_prob, old_prob, prob_log_ratio, rand_cutoff); }
		if (prob_log_ratio > rand_cutoff) {
//fprintf(stderr, "[%d/%d] Accepted (LR=%.2f - %.2f = %.4f, cutoff %.2f)\n", i+1, s->steps, new_prob, old_prob, prob_log_ratio, rand_cutoff);
			addModelSetToRecord(new_M, i, MR);
			old_prob = new_prob;
			s->stats.acceptances++;
			if (prob_log_ratio < 0.0) { s->stats.backward_steps++; }
		} else {
			destroyModelSet(new_M);
//fprintf(stderr, "[%d/%d] Rejected (LR=%.2f - %.2f = %.4f, cutoff %.2f)\n", i+1, s->steps, new_prob, old_prob, prob_log_ratio, rand_cutoff);
		}
	}
	return MR;
}

ModelSet * getInitialModelSet(settings * s) {
	ModelSet * IM = NULL;
	if (strcmp(s->model_hintfile, "")) {
		FILE * mhf = fopen(s->model_hintfile, "r");
		if (!mhf) {
			perror("Unable to open model hint file");
			exit(1);
		}
		IM = loadModelSet(mhf, s);
		fclose(mhf);
	} else {
		IM = newModelSet("foo", s->chain_order, s->num_sources);
	}
	assert(IM);
	return IM;
}

/* TODO:
	start condition hinting
	stop parent script on crash in debug mode
	use standard function documentation styles in perl and C
*/
int main(int argc, char ** argv) {
	uint i;
	settings s, init_s_copy;

	s.steps = DEFAULT_STEPS;
	s.chain_order = DEFAULT_CHAIN_ORDER;
	s.num_sources = DEFAULT_SOURCES;
	s.model_var = DEFAULT_MODEL_VAR;
	s.freq_var = DEFAULT_FREQ_VAR;
	s.restarts = DEFAULT_RESTARTS;
	s.thread_id = '0';
	strcpy(s.model_logfile, "");
	strcpy(s.model_hintfile, "");
	readOptions(argc, argv, &s);

	if (strcmp(s.model_logfile, "")) { sprintf(strrchr(s.model_logfile, '\0'), ".%c", s.thread_id); }

	s.max_failed_perturbs = 1000;
	s.stats.failed_perturbs = 0;
	s.stats.acceptances = 0;
	s.stats.backward_steps = 0;

	SeqStats * C = readSeqStats(s.chain_order);

	seedPRNG();

	ModelSet * start_M = getInitialModelSet(&s);
	
	double deviance = checkModelSetDeviance(start_M);
	if (deviance > MAX_MODEL_DEVIANCE) {
		fprintf(stderr, "[t%c] Model deviance %f exceeds maximum allowed deviance %f\n", s.thread_id, deviance, MAX_MODEL_DEVIANCE);
	}

	fprintf(stderr, "[t%c] Calibrating variance...\n", s.thread_id);

	calibrateVariance(start_M, C, &s);

	ModelRecord * best_MR = NULL;
	ModelSet * best_ModelSet = NULL;
	float best_ModelSet_prob = -1e999; //DBL_MIN;
	// restart N times, perturbing from mean with different variances
	// select the model with highest log prob, saving its model record
	for (i=1; i<=s.restarts; i++) {
		ModelRecord * cur_MR = newModelRecord();
		if (!runMCMC(start_M, cur_MR, C, &s)) {
			fprintf(stderr, "[t%c] Failed to run perturbations. Lowering variance...\n", s.thread_id);
			s.model_var /= 2;
			s.freq_var /= 2;
			// TODO: this is hacky
			if (MAX_RESTARTS > s.restarts) { s.restarts++; }
			continue;
		}
		ModelSet * cur_BM = getBestModelSet(cur_MR, C, &s);

		if (best_ModelSet_prob < cur_BM->logprob) {
			fprintf(stderr, "[t%c] Restart %d: LL=%e replacing LL=%e\n", s.thread_id, i, cur_BM->logprob, best_ModelSet_prob);
			if (best_MR) { destroyModelRecord(best_MR); }
			best_MR = cur_MR;
			best_ModelSet_prob = cur_BM->logprob;
			best_ModelSet = cur_BM;
		} else {
			destroyModelRecord(cur_MR);
		}
		fprintf(stderr, "[t%c] Restart %d done\n", s.thread_id, i);
	}

	fprintf(stderr, "[t%c] Restarts done\n", s.thread_id);

	assert(best_MR);
	assert(best_ModelSet);

//	Model * best_Model = getBestModel(best_MR, &s);
	deviance = checkModelSetDeviance(best_ModelSet);
	if (deviance > MAX_MODEL_DEVIANCE) {
		fprintf(stderr, "[t%c] Model deviance %f exceeds maximum allowed deviance %f\n", s.thread_id, deviance, MAX_MODEL_DEVIANCE);
//		abort();
	}

	fprintf(stderr, "[t%c] Algorithm statistics:\n", s.thread_id);
	fprintf(stderr, "[t%c] \tAcceptance rate: %d / %d = %f\n", s.thread_id, s.stats.acceptances, s.steps, (double)s.stats.acceptances/(double)s.steps);
	if (strcmp(s.model_logfile, "")) { printModelRecord(best_MR, &s); }

	dumpModelSet(best_ModelSet, stdout);

	return EXIT_SUCCESS;
}
