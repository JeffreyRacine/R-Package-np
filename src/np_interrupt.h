#ifndef NP_INTERRUPT_H
#define NP_INTERRUPT_H

/* Failure-only adapter; never returns to a released native search. */
void np_raise_real_interrupt(void);

#endif
