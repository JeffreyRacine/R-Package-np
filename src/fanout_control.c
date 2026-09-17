/* Private cooperative fan-out mailbox. No R vector is borrowed by MPI.
 * A pool owns one duplicated control context, reused across transactions.
 * The registry deliberately survives R unwinds: only complete requests and
 * the R owner's CLOSED-receipt proof permit retirement. No GC finalizer calls
 * MPI or frees pending memory. General Rmpi nonblocking APIs stay disabled. */
#include "fanout_control.h"
#include <R.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

enum { NP_CANCEL = 1, NP_RELEASE = 2, NP_CONTROL_TAG = 0 };
typedef struct {
    int data[3];
    MPI_Request request;
    int posted;
} np_control_slot;

typedef struct np_control_pool {
    MPI_Comm base, channel;
    int size, rank, generation, active, poisoned, cancelled, released;
    np_control_slot *out;
    np_control_slot in;
    SEXP owner;
    struct np_control_pool *next;
} np_control_pool;

static np_control_pool *np_control_pools = NULL;

static void np_control_check(np_control_pool *p, int status)
{
    if (status != MPI_SUCCESS) {
        char message[MPI_MAX_ERROR_STRING];
        int length = 0;
        p->poisoned = 1;
        MPI_Error_string(status, message, &length);
        Rf_error("MPI fan-out control transport failed; pool quarantined: %.*s",
                 length, message);
    }
}

static np_control_pool *np_control_find(MPI_Comm base)
{
    np_control_pool *p;
    for (p = np_control_pools; p; p = p->next)
        if (p->base == base) return p;
    return NULL;
}

static np_control_pool *np_control_active(MPI_Comm base)
{
    np_control_pool *p = np_control_find(base);
    if (!p || !p->active)
        Rf_error("MPI fan-out control has no active transaction");
    if (p->poisoned)
        Rf_error("MPI fan-out control transport is quarantined after a native failure");
    return p;
}

static void np_control_post_receive(np_control_pool *p)
{
    np_control_check(p, MPI_Irecv(p->in.data, 3, MPI_INT, 0,
                                 NP_CONTROL_TAG, p->channel, &p->in.request));
    p->in.posted = 1;
}

SEXP np_fanout_control_begin(MPI_Comm base, SEXP owner)
{
    np_control_pool *p = np_control_find(base);
    if (base == MPI_COMM_NULL) Rf_error("MPI fan-out requires an active pool");
    if (!p) {
        int size = 0, rank = 0, status;
        status = MPI_Comm_size(base, &size);
        if (status != MPI_SUCCESS || size < 2)
            Rf_error("MPI fan-out requires a healthy worker pool");
        status = MPI_Comm_rank(base, &rank);
        if (status != MPI_SUCCESS)
            Rf_error("MPI fan-out cannot determine pool rank");
        if ((size_t) size > SIZE_MAX / (2 * sizeof(np_control_slot)))
            Rf_error("MPI fan-out control allocation exceeds addressable storage");
        p = calloc(1, sizeof(*p));
        if (!p) Rf_error("cannot allocate MPI fan-out control owner");
        p->out = rank == 0 ? calloc((size_t) size * 2, sizeof(*p->out)) : NULL;
        if (rank == 0 && !p->out) {
            free(p);
            Rf_error("cannot allocate MPI fan-out control slots");
        }
        p->base = base;
        p->channel = MPI_COMM_NULL;
        p->size = size;
        p->rank = rank;
        p->owner = R_NilValue;
        p->next = np_control_pools;
        np_control_pools = p;
        /* Partial activation is retained, never silently retried. */
        p->poisoned = 1;
        np_control_check(p, MPI_Comm_dup(base, &p->channel));
        np_control_check(p, MPI_Comm_set_errhandler(p->channel, MPI_ERRORS_RETURN));
        p->poisoned = 0;
    }
    if (p->active || p->poisoned)
        Rf_error("an earlier MPI fan-out control transaction is not retired");
    if (p->generation == INT_MAX)
        Rf_error("MPI fan-out control generation exhausted; close this pool");
    if (p->rank == 0 && TYPEOF(owner) != ENVSXP)
        Rf_error("MPI fan-out master requires a metadata owner");
    if (p->rank != 0 && owner != R_NilValue)
        Rf_error("MPI fan-out worker must not retain a master owner");
    /* Preserve before starting requests; R allocation failure cannot strand
       an unrooted R owner behind a successful nonblocking MPI call. */
    if (owner != R_NilValue) R_PreserveObject(owner);
    p->owner = owner;
    p->active = 1;
    p->generation++;
    p->cancelled = p->released = 0;
    p->in.request = MPI_REQUEST_NULL;
    p->in.posted = 0;
    if (p->rank == 0) {
        for (size_t i = 0; i < (size_t) p->size * 2; ++i) {
            p->out[i].request = MPI_REQUEST_NULL;
            p->out[i].posted = 0;
        }
    } else np_control_post_receive(p);
    return R_NilValue;
}

SEXP np_fanout_control_send(MPI_Comm base, int rank, int command)
{
    np_control_pool *p = np_control_active(base);
    np_control_slot *slot;
    if (p->rank != 0 || rank < 1 || rank >= p->size ||
        (command != NP_CANCEL && command != NP_RELEASE))
        Rf_error("invalid private MPI fan-out control command");
    slot = &p->out[(size_t) rank * 2 + command - 1];
    if (slot->posted) return R_NilValue;
    if (command == NP_CANCEL && p->out[(size_t) rank * 2 + 1].posted)
        Rf_error("MPI fan-out cannot cancel a released worker");
    slot->data[0] = p->generation;
    slot->data[1] = command;
    slot->data[2] = rank;
    /* Mark initiated before MPI: a failed send must never be retried or have
       its storage reused. A compile-only adversarial build uses Issend. */
    slot->posted = 1;
#ifdef NP_FANOUT_SYNCHRONOUS_CONTROL
    np_control_check(p, MPI_Issend(slot->data, 3, MPI_INT, rank,
                                  NP_CONTROL_TAG, p->channel, &slot->request));
#else
    np_control_check(p, MPI_Isend(slot->data, 3, MPI_INT, rank,
                                 NP_CONTROL_TAG, p->channel, &slot->request));
#endif
    return R_NilValue;
}

SEXP np_fanout_control_poll(MPI_Comm base)
{
    np_control_pool *p = np_control_active(base);
    if (p->rank == 0) {
        int complete = 1;
        for (size_t i = 2; i < (size_t) p->size * 2; ++i) {
            np_control_slot *slot = &p->out[i];
            if (slot->request != MPI_REQUEST_NULL) {
                int done = 0;
                np_control_check(p, MPI_Test(&slot->request, &done, MPI_STATUS_IGNORE));
                if (!done) complete = 0;
            }
        }
        return Rf_ScalarInteger(complete);
    }
    if (p->in.request != MPI_REQUEST_NULL) {
        int done = 0;
        np_control_check(p, MPI_Test(&p->in.request, &done, MPI_STATUS_IGNORE));
        if (done) {
            int command = p->in.data[1];
            if (p->in.data[0] != p->generation || p->in.data[2] != p->rank ||
                (command != NP_CANCEL && command != NP_RELEASE) ||
                (command == NP_CANCEL && p->cancelled)) {
                p->poisoned = 1;
                Rf_error("MPI fan-out control received an invalid generation/command");
            }
            if (command == NP_CANCEL) {
                p->cancelled = 1;
                np_control_post_receive(p);
            } else p->released = 1;
        }
    }
    return Rf_ScalarInteger(p->cancelled + 2 * p->released);
}

SEXP np_fanout_control_finish(MPI_Comm base)
{
    np_control_pool *p = np_control_active(base);
    if (p->rank == 0) {
        for (int rank = 1; rank < p->size; ++rank) {
            if (!p->out[(size_t) rank * 2 + 1].posted ||
                p->out[(size_t) rank * 2].request != MPI_REQUEST_NULL ||
                p->out[(size_t) rank * 2 + 1].request != MPI_REQUEST_NULL)
                Rf_error("MPI fan-out has pending control retirement");
        }
    } else if (!p->released || p->in.request != MPI_REQUEST_NULL) {
        Rf_error("MPI fan-out worker has not consumed its release");
    }
    p->active = 0;
    if (p->owner != R_NilValue) R_ReleaseObject(p->owner);
    p->owner = R_NilValue;
    return R_NilValue;
}

SEXP np_fanout_control_owner(MPI_Comm base)
{
    np_control_pool *p = np_control_find(base);
    return p ? p->owner : R_NilValue;
}

void np_fanout_control_close(MPI_Comm base)
{
    np_control_pool **link = &np_control_pools, *p;
    while (*link && (*link)->base != base) link = &(*link)->next;
    p = *link;
    if (!p) return;
    if (p->active || p->poisoned)
        Rf_error("MPI pool has pending/quarantined fan-out transport; close was not attempted");
    np_control_check(p, MPI_Comm_free(&p->channel));
    *link = p->next;
    free(p->out);
    free(p);
}

void np_fanout_control_finalize(void)
{
    np_control_pool *p;
    for (p = np_control_pools; p; p = p->next)
        if (p->active || p->poisoned)
            Rf_error("MPI has pending/quarantined fan-out transport; finalization was not attempted");
    while (np_control_pools) np_fanout_control_close(np_control_pools->base);
}
