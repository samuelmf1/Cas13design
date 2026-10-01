/* SPDX-License-Identifier: GPL-2.0-or-later
 *
 * RNAhybrid core copyright (C) 2004 Marc Rehmsmeier, Peter Steffen,
 * Matthias Hoechsmann. This prototype links and adapts RNAhybrid 2.1.2 code.
 *
 * Batched Cas13design RNAhybrid runner.
 *
 * Read target,query pairs (comma- or tab-delimited), one pair per line, and
 * evaluate them in one persistent RNAhybrid 2.1.2 process.  Output deliberately
 * uses RNAhybrid's existing compact format so the production R parser can be
 * tested unchanged.
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "energy.h"
#include "globals.h"
#include "hybrid_core.h"
#include "input.h"
#include "minmax.h"

void tableAlloc(int target_length, int query_length);

static float maximal_duplex_energy(const char *seq_arg)
{
    const int seq_len = (int)strlen(seq_arg);
    char *seq = calloc((size_t)seq_len + 1, sizeof(char));
    char *rc = calloc((size_t)seq_len + 1, sizeof(char));
    float result = 0.0f;
    int i;

    if (seq == NULL || rc == NULL) {
        fprintf(stderr, "allocation failure\n");
        exit(2);
    }

    for (i = 0; i < seq_len; i++) {
        const char c = seq_arg[i];
        if (c == 'a' || c == 'A') seq[i] = A;
        else if (c == 'c' || c == 'C') seq[i] = C;
        else if (c == 'g' || c == 'G') seq[i] = G;
        else if (c == 'u' || c == 'U' || c == 't' || c == 'T') seq[i] = U;
        else seq[i] = N;
    }

    for (i = 0; i < seq_len; i++) {
        const char c = seq[seq_len - i - 1];
        if (c == A) rc[i] = U;
        else if (c == C) rc[i] = G;
        else if (c == G) rc[i] = C;
        else if (c == U) rc[i] = A;
        else rc[i] = N;
    }
    rc[seq_len] = '\0';

    for (i = 0; i < seq_len - 1; i++) {
        result += (float)stack_dg_ar
            [(int)rc[seq_len - i - 1]]
            [(int)rc[seq_len - (i + 1) - 1]]
            [(int)seq[i + 1]]
            [(int)seq[i]];
    }

    free(seq);
    free(rc);
    return result;
}

static void strip_space(char *seq)
{
    char *src = seq;
    char *dst = seq;
    while (*src != '\0') {
        if (!isspace((unsigned char)*src)) *dst++ = *src;
        src++;
    }
    *dst = '\0';
}

int main(void)
{
    enum { MAX_PAIR_LEN = 30, LINE_LEN = 256 };
    char line[LINE_LEN];
    char target[MAX_PAIR_LEN + 1];
    char query[MAX_PAIR_LEN + 1];
    char *separator;
    char *raw_target;
    char *raw_query;

    iloop_upper_limit = ILOOPUPPERLIMITDEFAULT;
    bloop_upper_limit = BLOOPUPPERLIMITDEFAULT;
    helix_start = 0;
    helix_end = 0;

    init_constants();
    init_energies();
    tableAlloc(MAX_PAIR_LEN, MAX_PAIR_LEN);

    r1 = t1 = calloc(2 * MAX_PAIR_LEN, sizeof(char));
    r2 = t2 = calloc(2 * MAX_PAIR_LEN, sizeof(char));
    r3 = t3 = calloc(2 * MAX_PAIR_LEN, sizeof(char));
    r4 = t4 = calloc(2 * MAX_PAIR_LEN, sizeof(char));
    x = calloc(MAX_PAIR_LEN + 2, sizeof(char));
    y = calloc(MAX_PAIR_LEN + 2, sizeof(char));
    if (r1 == NULL || r2 == NULL || r3 == NULL || r4 == NULL ||
        x == NULL || y == NULL) {
        fprintf(stderr, "allocation failure\n");
        return 2;
    }

    while (fgets(line, sizeof(line), stdin) != NULL) {
        float mde;
        float xi;
        float theta;

        separator = strchr(line, ',');
        if (separator == NULL) separator = strchr(line, '\t');
        if (separator == NULL) {
            fprintf(stderr, "expected target,query pair: %s", line);
            return 2;
        }
        *separator = '\0';
        raw_target = line;
        raw_query = separator + 1;
        strip_space(raw_target);
        strip_space(raw_query);
        if (strlen(raw_target) > MAX_PAIR_LEN || strlen(raw_query) > MAX_PAIR_LEN) {
            fprintf(stderr, "sequence length exceeds %d\n", MAX_PAIR_LEN);
            return 2;
        }
        strcpy(target, raw_target);
        strcpy(query, raw_query);

        m = (int)strlen(target);
        n = (int)strlen(query);
        if (m < 1 || n < 1 || m > MAX_PAIR_LEN || n > MAX_PAIR_LEN) {
            fprintf(stderr, "sequence length outside 1..%d: %d,%d\n",
                    MAX_PAIR_LEN, m, n);
            return 2;
        }

        strcpy(x, " ");
        strcat(x, target);
        convert_x();

        y[0] = ' ';
        for (int i = 0; i < n; i++) y[i + 1] = query[n - i - 1];
        y[n + 1] = '\0';
        convert_y();

        mde = maximal_duplex_energy(query);
        xi = XI_SLOPE_3UTR_HUMAN * mde + XI_INTERCEPT_3UTR_HUMAN;
        theta = THETA_SLOPE_3UTR_HUMAN * mde + THETA_INTERCEPT_3UTR_HUMAN;

        mainloop(0, 1, 0, 0.0f, 0, 1.0f, 1,
                 "command_line", target, "command_line", query,
                 0, 0, xi, theta);
    }

    return ferror(stdin) ? 2 : 0;
}
